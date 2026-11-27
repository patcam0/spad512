function [I, R, nFramesTotal, info] = spad_roi_correlations_fast( ...
    current_path, subfolders, im_size, hot_pixels_lin, sparseMode, readerOpts, roi, varargin)
%SPAD_ROI_CORRELATIONS_FAST Accumulate intensity and same-frame coincidences.
%
% roi = [rowMin colMin nRows nCols], in MATLAB (one-based) coordinates.
% I is nRows-by-nCols. R is P-by-P, with P = nRows*nCols and pixels ordered
% as I(:). R includes its diagonal and is an unnormalized coincidence sum.
% Empty frames contribute to nFramesTotal even if they contain no ROI hits.
% The direct sparse reader omits empty ROI columns from the multiplication.
%
% Example:
%   [I,R,nFramesTotal,info] = spad_roi_correlations_fast( ...
%       current_path, subfolders(1:n_folders), im_size, hot_pixels_lin, ...
%       spad_cam_obj.sparse, opts, [rowMin colMin nRows nCols], ...
%       'UseGPU','auto', 'BatchFrames',8192, 'MatrixFormat','auto');
%
% Name/value options:
%   UseGPU           'auto' (default), true, or false. 'auto' uses a GPU
%                    when supported and a conservative memory estimate fits.
%   BatchFrames      Maximum frames per multiplication (default 8192).
%                    Capped at 2^24 to protect binary single-precision sums.
%   MatrixFormat     'auto', 'dense', or 'sparse' (default 'auto'). The auto
%                    choice is a density heuristic; benchmark both formats.
%   SparseThreshold  Auto uses sparse matrices at density <= this (0.02).
%   ComputePrecision 'auto', 'single', or 'double' (default 'auto'). Auto
%                    uses single for logical frames, double otherwise.
%                    I and R ALWAYS accumulate in double precision.
%   ReturnOnGPU      Keep I and R on GPU if one is used (default false).
%   Verbose          Print one message per folder (default false).
%   ReadSparseROI    Read sparse events directly into the ROI matrix instead
%                    of a full frame cube (default true). Reuses the format
%                    rules and helpers of read_sparse_frames_fast.
%   FrameReader      Optional replacement reader handle, with signature
%                    [data,nFrames] = f(filename,im_size,hot_pixels_lin,opts).
%                    Supplying it bypasses the direct sparse reader.
%
% Requires the project's existing listFilenames_fast and reader dependencies.
% The packed path calls the existing read_packed_frames_fast. GPU execution
% requires Parallel Computing Toolbox; single sparse matrices need R2025a+.
%
% GPU batching bounds numeric ROI working memory, not the size of R. R alone
% uses 8*P^2 bytes. The direct sparse path still reads one file's encoded
% values at a time. The packed/fallback path still reads one full frame cube.
% Counts from logical data remain exact up to the double integer limit
% (2^53): each binary batch has at most 2^24 frames and is promoted to double
% before it is added. Explicit single on nonbinary data can introduce error.

    p = inputParser;
    p.FunctionName = mfilename;
    addParameter(p, 'UseGPU', 'auto');
    addParameter(p, 'BatchFrames', 8192, @(x) isnumeric(x) && isscalar(x) ...
        && isreal(x) && isfinite(x) && x >= 1 && x == fix(x) && x <= 2^24);
    addParameter(p, 'MatrixFormat', 'auto');
    addParameter(p, 'SparseThreshold', 0.02, @(x) isnumeric(x) && isscalar(x) ...
        && isreal(x) && isfinite(x) && x >= 0 && x <= 1);
    addParameter(p, 'ComputePrecision', 'auto');
    addParameter(p, 'ReturnOnGPU', false, @is_scalar_flag);
    addParameter(p, 'Verbose', false, @is_scalar_flag);
    addParameter(p, 'ReadSparseROI', true, @is_scalar_flag);
    addParameter(p, 'FrameReader', [], @(x) isempty(x) || isa(x, 'function_handle'));
    parse(p, varargin{:});
    cfg = p.Results;
    cfg.MatrixFormat = validatestring(cfg.MatrixFormat, {'auto','dense','sparse'});
    cfg.ComputePrecision = validatestring(cfg.ComputePrecision, {'auto','single','double'});

    validateattributes(im_size, {'numeric'}, ...
        {'real','finite','integer','positive','numel',2}, mfilename, 'im_size');
    validateattributes(roi, {'numeric'}, ...
        {'real','finite','integer','positive','numel',4}, mfilename, 'roi');
    validateattributes(sparseMode, {'numeric','logical'}, ...
        {'real','scalar','binary'}, mfilename, 'sparseMode');
    im_size = double(im_size(:).');
    roi = double(roi(:).');
    nRows = roi(3);
    nCols = roi(4);
    if roi(1) + nRows - 1 > im_size(1) || roi(2) + nCols - 1 > im_size(2)
        error('spad_roi_correlations_fast:ROI', 'The ROI extends outside im_size.');
    end
    rows = roi(1) + (0:nRows-1);
    cols = roi(2) + (0:nCols-1);
    nPix = nRows*nCols;
    if iscell(subfolders)
        subfolders = cellfun(@char, subfolders, 'UniformOutput', false);
    else
        subfolders = cellstr(subfolders);
    end
    subfolders = subfolders(:);
    basePath = fullfile(char(current_path), 'data', 'intensity_images');

    % A lookup table maps stream pixel IDs straight to ROI indices, including
    % the hot-pixel mask. This is prepared once, instead of once per file.
    directSparse = logical(sparseMode) && logical(cfg.ReadSparseROI) ...
        && isempty(cfg.FrameReader);
    if directSparse
        [sparseOpts, roiMap, bytesPerValue, thresholdValue] = ...
            prepare_sparse_roi(im_size, hot_pixels_lin, readerOpts, rows, cols);
    elseif isempty(cfg.FrameReader)
        if sparseMode
            cfg.FrameReader = @read_sparse_frames_fast;
        else
            cfg.FrameReader = @read_packed_frames_fast;
        end
    end

    [useGPU, batchFrames, device, deviceReason] = select_device(cfg, nPix);
    if useGPU
        R = gpuArray.zeros(nPix, nPix, 'double');
        Ivec = gpuArray.zeros(nPix, 1, 'double');
    else
        R = zeros(nPix, nPix, 'double');
        Ivec = zeros(nPix, 1, 'double');
    end
    nFramesTotal = 0;
    info = struct('UsedGPU', useGPU, 'Device', device, ...
        'DeviceReason', deviceReason, 'BatchFrames', batchFrames, ...
        'FilesRead', 0, 'SparseFiles', 0, 'DenseFiles', 0, 'EmptyFiles', 0, ...
        'DirectSparseROI', directSparse);

    for folderIdx = 1:numel(subfolders)
        folderPath = fullfile(basePath, subfolders{folderIdx});
        filenames = listFilenames_fast(folderPath);
        if cfg.Verbose
            fprintf('Folder %d/%d (%d files): %s\n', ...
                folderIdx, numel(subfolders), numel(filenames), subfolders{folderIdx});
        end

        for fileIdx = 1:numel(filenames)
            filename = fullfile(folderPath, filenames{fileIdx});
            if directSparse
                [Xraw, nFrames] = read_sparse_roi( ...
                    filename, sparseOpts, roiMap, bytesPerValue, thresholdValue, nPix);
            else
                [data, nFrames] = cfg.FrameReader( ...
                    filename, im_size, hot_pixels_lin, readerOpts);
                validateattributes(nFrames, {'numeric'}, ...
                    {'real','finite','scalar','integer','nonnegative'});
                if size(data,1) ~= im_size(1) || size(data,2) ~= im_size(2) ...
                        || size(data,3) ~= nFrames || ndims(data) > 3
                    error('spad_roi_correlations_fast:ReaderShape', ...
                        'Reader must return [im_size(1),im_size(2),nFrames]: %s', filename);
                end
                % Each column is one frame. No frame/pixel permutation.
                Xraw = reshape(data(rows, cols, :), nPix, nFrames);
                clear data
            end
            info.FilesRead = info.FilesRead + 1;
            nFramesTotal = nFramesTotal + double(nFrames);
            nonzerosInFile = nnz(Xraw);
            if nFrames == 0 || nonzerosInFile == 0
                info.EmptyFiles = info.EmptyFiles + 1;
                clear Xraw
                continue;
            end
            nMatrixFrames = size(Xraw, 2);

            if strcmp(cfg.MatrixFormat, 'auto')
                useSparse = nonzerosInFile <= cfg.SparseThreshold * double(nPix) * double(nMatrixFrames);
            else
                useSparse = strcmp(cfg.MatrixFormat, 'sparse');
            end
            if useSparse
                info.SparseFiles = info.SparseFiles + 1;
            else
                info.DenseFiles = info.DenseFiles + 1;
            end
            precision = cfg.ComputePrecision;
            if strcmp(precision, 'auto')
                if islogical(Xraw)
                    precision = 'single';
                else
                    precision = 'double';
                end
            end

            for firstFrame = 1:batchFrames:nMatrixFrames
                lastFrame = min(firstFrame + batchFrames - 1, nMatrixFrames);
                raw = Xraw(:, firstFrame:lastFrame);

                if useSparse
                    % Logical sparse arrays are not supported on GPU. Cast
                    % while sparse on the CPU, then transfer only nonzeros.
                    X = cast(sparse(raw), precision);
                    if useGPU
                        X = gpuArray(X);
                    end
                else
                    raw = full(raw);
                    if useGPU
                        % For logical frames this uploads one byte/pixel;
                        % the numeric conversion happens on the GPU.
                        raw = gpuArray(raw);
                    end
                    X = cast(raw, precision);
                end
                clear raw

                Ivec = Ivec + double(full(sum(X, 2)));
                gram = X * transpose(X);
                R = R + double(full(gram));
                clear X gram
            end
            clear Xraw
        end
    end

    I = reshape(Ivec, nRows, nCols);
    if useGPU && ~cfg.ReturnOnGPU
        I = gather(I);
        R = gather(R);
    end
end


function tf = is_scalar_flag(x)
    tf = (islogical(x) || isnumeric(x)) && isreal(x) && isscalar(x) ...
        && (x == 0 || x == 1);
end


function [useGPU, batchFrames, deviceName, reason] = select_device(cfg, nPix)
    autoGPU = (ischar(cfg.UseGPU) || (isstring(cfg.UseGPU) && isscalar(cfg.UseGPU))) ...
        && strcmpi(cfg.UseGPU, 'auto');
    if ~autoGPU && ~is_scalar_flag(cfg.UseGPU)
        error('spad_roi_correlations_fast:UseGPU', 'UseGPU must be ''auto'', true, or false.');
    end
    batchFrames = double(cfg.BatchFrames);
    deviceName = 'CPU';
    useGPU = false;
    reason = 'CPU requested.';
    if ~autoGPU && ~logical(cfg.UseGPU)
        return;
    end
    if ~canUseGPU()
        reason = 'No supported GPU available.';
        if ~autoGPU
            error('spad_roi_correlations_fast:GPU', '%s', reason);
        end
        return;
    end

    d = gpuDevice;
    % Include R, the temporary product (possibly sparse), conversion and
    % addition buffers. Leave 25% free for library workspaces and other use.
    % This is a conservative estimate, not a device allocation guarantee.
    budget = 0.75 * double(d.AvailableMemory);
    reserve = 6 * 8 * double(nPix)^2 + 32 * double(nPix);
    possibleBatch = floor((budget - reserve) / (32 * double(nPix)));
    fits = possibleBatch >= 1 && double(nPix)^2 <= double(intmax('int32'));
    if ~fits
        reason = sprintf(['GPU estimate needs at least %.3g GiB for R and temporaries. ' ...
            'Use a smaller ROI or UseGPU=false.'], reserve/2^30);
        if ~autoGPU
            error('spad_roi_correlations_fast:GPUMemory', '%s', reason);
        end
        return;
    end
    batchFrames = min([batchFrames, possibleBatch, floor(double(intmax('int32'))/nPix)]);
    useGPU = true;
    deviceName = d.Name;
    reason = 'GPU available and memory estimate fits.';
end


function [opts, roiMap, bytesPerValue, thresholdValue] = prepare_sparse_roi( ...
    im_size, hot_pixels_lin, readerOpts, rows, cols)
    % Retain the original reader's option normalization and stream rules.
    opts = parse_reader_options_fast(readerOpts);
    opts.sparse = true;
    opts.format = 'sparse';
    opts.frameCountLocation = normalize_frame_count_location_fast(opts.frameCountLocation);
    nPixFull = prod(im_size);
    threshold = sparse_separator_threshold_fast(opts, nPixFull);
    bytesPerValue = bytes_per_value_fast(opts.indexClass);
    thresholdValue = cast(threshold, opts.indexClass);
    keepVec = make_keep_mask_fast(hot_pixels_lin, im_size, opts.pixelOrder);

    [rr, cc] = ndgrid(rows, cols);
    switch normalize_pixel_order_fast(opts.pixelOrder)
        case 'row-major'
            streamIdx = (rr - 1)*im_size(2) + cc;
        case 'matlab'
            streamIdx = rr + (cc - 1)*im_size(1);
    end
    streamIdx = streamIdx(:);
    localIdx = (1:numel(streamIdx)).';
    keep = keepVec(streamIdx);
    roiMap = zeros(nPixFull, 1);
    roiMap(streamIdx(keep)) = localIdx(keep);
end


function [X, nFrames] = read_sparse_roi( ...
    filename, opts, roiMap, bytesPerValue, thresholdValue, nPixROI)
    fileInfo = dir(filename);
    if isempty(fileInfo)
        error('spad_roi_correlations_fast:File', 'File not found: %s', filename);
    end
    if mod(fileInfo.bytes, bytesPerValue) ~= 0
        error('spad_roi_correlations_fast:FileSize', ...
            'File size is not a whole number of %s values: %s', opts.indexClass, filename);
    end
    fid = fopen(filename, 'rb', opts.endian);
    if fid < 0
        error('spad_roi_correlations_fast:Open', 'Could not open file: %s', filename);
    end
    closeFile = onCleanup(@() fclose(fid)); %#ok<NASGU>
    vals = fread(fid, Inf, ['*' opts.indexClass]);
    if isempty(vals)
        nFrames = 0;
        X = sparse([], [], false(0,1), nPixROI, 0);
        return;
    end
    vals = vals(:);
    isMarker = vals >= thresholdValue;

    % These three branches preserve read_sparse_frames_fast's interpretation
    % of separators, header counts, footer counts and trailing frames.
    switch opts.frameCountLocation
        case 'none'
            isHit = ~isMarker;
            frameNumsAll = cumsum(double(isMarker)) + 1;
            frameNums = frameNumsAll(isHit);
            nFrames = nnz(isMarker) + double(~isMarker(end));
        case 'start'
            [isHeader, ~] = mark_sparse_count_headers_fast(isMarker, true);
            isHit = ~isMarker & ~isHeader;
            frameNumsAll = cumsum(double(isHeader));
            frameNums = frameNumsAll(isHit);
            nFrames = nnz(isHeader);
            if any(frameNums == 0)
                error('spad_roi_correlations_fast:Header', ...
                    'Hits before first frame count header: %s', filename);
            end
        case 'end'
            isFooter = false(size(isMarker));
            if numel(vals) > 1
                isFooter(1:end-1) = ~isMarker(1:end-1) & isMarker(2:end);
            end
            if ~isMarker(end)
                isFooter(end) = true;
            end
            isHit = ~isMarker & ~isFooter;
            frameNumsAll = cumsum(double(isMarker)) + 1;
            frameNums = frameNumsAll(isHit);
            nFrames = nnz(isFooter);
    end
    if nFrames == 0
        X = sparse([], [], false(0,1), nPixROI, 0);
        return;
    end
    idx = double(vals(isHit));
    if opts.pixelBase == 0
        idx = idx + 1;
    end
    bad = idx < 1 | idx > numel(roiMap) | idx ~= floor(idx);
    if any(bad)
        error('spad_roi_correlations_fast:PixelID', ...
            'Pixel ID outside the image range: %s', filename);
    end
    if any(frameNums < 1 | frameNums > nFrames)
        error('spad_roi_correlations_fast:FrameID', ...
            'Inconsistent frame separators/counts: %s', filename);
    end
    localIdx = roiMap(idx);
    keep = localIdx ~= 0;
    localIdx = localIdx(keep);
    frameNums = frameNums(keep);
    % Empty ROI frames have no effect on I or R. Renumber the remaining
    % columns consecutively, while returning the original nFrames above.
    if isempty(frameNums)
        X = sparse([], [], false(0,1), nPixROI, 0);
        return;
    end
    compactFrames = cumsum([1; double(diff(frameNums) ~= 0)]);
    % Logical values combine repeated (pixel,frame) entries with OR, matching
    % the original reader's assignment into logical frames, rather than
    % counting duplicate photon IDs more than once.
    X = sparse(localIdx, compactFrames, true(size(localIdx)), nPixROI, compactFrames(end));
end
