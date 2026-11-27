function [I, R, nFrames, info] = read_sparse_intensity_correlation( ...
    filename, im_size, hot_pixels_lin, readerOptions, processing)
%READ_SPARSE_INTENSITY_CORRELATION Accumulate sparse frames without images.
%
% [I,R,nFrames,info] = read_sparse_intensity_correlation( ...
%     filename, im_size, hot_pixels_lin, readerOptions, processing)
%
% Inputs
%   im_size         : [numberOfRows numberOfColumns] before cropping/binning.
%   hot_pixels_lin  : Same hot-pixel input as read_sparse_frames_fast; [] skips.
%   readerOptions   : The same options struct passed to your existing reader.
%   processing      : Optional struct with these fields:
%       ROI                 [] = full image, or [rowFirst rowLast colFirst colLast].
%                           Bounds are inclusive, 1-based sensor coordinates.
%       BinSize             [1 1], or [binRows binColumns], or a scalar.
%       BinMode             'binary' (default) or 'sum'.
%       FramesPerBlock      100000 (limits each sparse frame matrix).
%       ComputeCorrelation  true; false returns R=[] without allocating it.
%       SparseCorrelation   false; true returns a sparse double R.
%
% Outputs
%   I       : Summed, cropped, binned intensity image (double).
%   R       : Unnormalized same-frame correlations, including the diagonal.
%             For each binned frame x: I(:) += x; R += x*x.'.
%   nFrames : All decoded frames, including empty or completely excluded ones.
%   info    : ROI, BinSize, ImageSize, BinMode and NumberOfBins.
%
% BinMode='binary': x(b)=1 if bin b contains any hit in that frame. This
%   implements R(b,b)=R(b,b)+1 for the unique occupied bins of each frame.
% BinMode='sum': x(b) is the number of hit records in bin b in that frame.
%   Several hits in one bin are retained; correlations use their products.
%   As with a binary SPAD frame, each native pixel should occur at most once
%   per frame. Repeated native-pixel records are counted separately in 'sum'.
%
% Bins start at the ROI's top-left corner; partial edge bins are retained.
% R uses MATLAB linear indexing into I, regardless of the input pixel order.
% The file is read into memory, as in the original reader; only the sparse
% frame matrices are blocked. A full double R alone uses 8*numel(I)^2 bytes.
%
% Requires the SAME format/option helper files as read_sparse_frames_fast.
% No full pixel-by-frame or rows-by-columns-by-frames array is constructed.

    if nargin < 3, hot_pixels_lin = []; end
    if nargin < 4, readerOptions = []; end
    if nargin < 5, processing = struct(); end

    validateattributes(im_size, {'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',2}, ...
        mfilename, 'im_size');
    im_size = double(im_size(:).');

    % Reuse the original parser so the file-format defaults stay identical.
    if isempty(readerOptions)
        opts = parse_reader_options_fast();
    else
        opts = parse_reader_options_fast(readerOptions);
    end
    opts.sparse = true;
    opts.format = 'sparse';
    opts.frameCountLocation = ...
        normalize_frame_count_location_fast(opts.frameCountLocation);
    processing = processing_options(processing, im_size);

    % One lookup per hit handles cropping, binning AND hot-pixel rejection.
    [pixelToBin, imageSize] = make_bin_lookup(im_size, opts, processing);
    nBins = prod(imageSize);

    values = read_values(filename, opts);
    [pixelIdx, frameNums, nFrames] = decode_hits(values, opts, prod(im_size));
    clear values

    I = zeros(imageSize);
    if ~processing.ComputeCorrelation
        R = [];
    elseif processing.SparseCorrelation
        R = sparse(nBins, nBins);
    else
        R = zeros(nBins, nBins);
    end
    info = struct('ROI', processing.ROI, 'BinSize', processing.BinSize, ...
        'ImageSize', imageSize, 'BinMode', processing.BinMode, ...
        'NumberOfBins', nBins);

    if isempty(pixelIdx), return; end

    % The helper returns a keep-mask in the same order as the stored IDs.
    if ~isempty(hot_pixels_lin)
        keepPixels = make_keep_mask_fast(hot_pixels_lin, im_size, opts.pixelOrder);
        pixelToBin(~keepPixels(:)) = 0;
    end
    binIdx = pixelToBin(pixelIdx);
    keep = binIdx > 0;
    binIdx = binIdx(keep);
    frameNums = frameNums(keep);
    clear pixelIdx pixelToBin keep
    if isempty(binIdx), return; end

    % Hit/frame vectors are already in frame order. Locate block boundaries
    % once, rather than scanning every hit again for every block.
    blockNumber = floor((frameNums - 1) / processing.FramesPerBlock);
    blockStart = [1; find(diff(blockNumber)) + 1];
    blockStop = [blockStart(2:end) - 1; numel(binIdx)];
    frameOffset = blockNumber(blockStart) * processing.FramesPerBlock;
    clear blockNumber

    % Logical sparse values merge hits in the same bin/frame using OR.
    % Double values instead add them, giving photon-count binning.
    if strcmp(processing.BinMode, 'binary')
        hitValue = true;
    else
        hitValue = 1;
    end

    intensity = zeros(nBins, 1);
    for block = 1:numel(blockStart)
        hits = blockStart(block):blockStop(block);
        localFrame = frameNums(hits) - frameOffset(block);
        nColumns = min(processing.FramesPerBlock, nFrames - frameOffset(block));

        X = double(sparse(binIdx(hits), localFrame, hitValue, nBins, nColumns));
        intensity = intensity + full(sum(X, 2));

        if processing.ComputeCorrelation
            correlations = X * X.';
            if processing.SparseCorrelation
                R = R + correlations;
            else
                R = R + full(correlations);
            end
        end
    end
    I = reshape(intensity, imageSize);
end


function settings = processing_options(settings, im_size)
% Keep processing options separate from the existing file-format options.
    validateattributes(settings, {'struct'}, {'scalar'}, mfilename, 'processing');
    parser = inputParser;
    parser.FunctionName = 'read_sparse_intensity_correlation';
    parser.PartialMatching = false;
    addParameter(parser, 'ROI', []);
    addParameter(parser, 'BinSize', [1 1]);
    addParameter(parser, 'BinMode', 'binary');
    addParameter(parser, 'FramesPerBlock', 100000);
    addParameter(parser, 'ComputeCorrelation', true);
    addParameter(parser, 'SparseCorrelation', false);
    parse(parser, settings);
    settings = parser.Results;

    if isempty(settings.ROI)
        settings.ROI = [1 im_size(1) 1 im_size(2)];
    end
    validateattributes(settings.ROI, {'numeric'}, ...
        {'real','finite','integer','positive','vector','numel',4}, ...
        mfilename, 'processing.ROI');
    settings.ROI = double(settings.ROI(:).');
    roi = settings.ROI;
    if roi(1) > roi(2) || roi(3) > roi(4) || ...
            roi(2) > im_size(1) || roi(4) > im_size(2)
        error('read_sparse_intensity_correlation:InvalidROI', ...
            'ROI must be [rowFirst rowLast colFirst colLast] inside im_size.');
    end

    validateattributes(settings.BinSize, {'numeric'}, ...
        {'real','finite','integer','positive','vector'}, ...
        mfilename, 'processing.BinSize');
    if isscalar(settings.BinSize)
        settings.BinSize = [settings.BinSize settings.BinSize];
    end
    if numel(settings.BinSize) ~= 2
        error('read_sparse_intensity_correlation:InvalidBinSize', ...
            'BinSize must be a scalar or [binRows binColumns].');
    end
    settings.BinSize = double(settings.BinSize(:).');
    settings.BinMode = validatestring(settings.BinMode, {'binary','sum'});
    validateattributes(settings.FramesPerBlock, {'numeric'}, ...
        {'real','finite','integer','positive','scalar'}, ...
        mfilename, 'processing.FramesPerBlock');
    settings.FramesPerBlock = double(settings.FramesPerBlock);
    validateattributes(settings.ComputeCorrelation, {'logical'}, {'scalar'});
    validateattributes(settings.SparseCorrelation, {'logical'}, {'scalar'});
end


function [pixelToBin, imageSize] = make_bin_lookup(im_size, opts, settings)
    roi = settings.ROI;
    rows = roi(1):roi(2);
    cols = roi(3):roi(4);
    imageSize = ceil([numel(rows) numel(cols)] ./ settings.BinSize);

    % Implicit expansion forms the ROI lookup in MATLAB image orientation.
    binRow = floor((0:numel(rows)-1).' / settings.BinSize(1));
    binCol = floor((0:numel(cols)-1) / settings.BinSize(2));
    pixelToBin = zeros(im_size);
    pixelToBin(rows, cols) = 1 + binRow + imageSize(1) * binCol;

    % Convert the lookup to the file's pixel order, not the output order.
    if strcmp(normalize_pixel_order_fast(opts.pixelOrder), 'row-major')
        pixelToBin = pixelToBin.';
    end
    pixelToBin = pixelToBin(:);
end


function values = read_values(filename, opts)
    filename = char(filename);
    fileInfo = dir(filename);
    if isempty(fileInfo)
        error('read_sparse_intensity_correlation:FileNotFound', ...
            'File not found: %s', filename);
    end
    bytesPerValue = bytes_per_value_fast(opts.indexClass);
    if mod(fileInfo.bytes, bytesPerValue) ~= 0
        error('read_sparse_intensity_correlation:IncompleteValue', ...
            'File size is not a whole number of %s values: %s', ...
            opts.indexClass, filename);
    end
    fid = fopen(filename, 'rb', opts.endian);
    if fid < 0
        error('read_sparse_intensity_correlation:OpenFailed', ...
            'Could not open file: %s', filename);
    end
    closeFile = onCleanup(@() fclose(fid)); %#ok<NASGU>
    values = fread(fid, Inf, ['*' opts.indexClass]);
    message = ferror(fid);
    if ~isempty(message)
        error('read_sparse_intensity_correlation:ReadFailed', ...
            'Could not read the complete file: %s\n%s', filename, message);
    end
end


function [idx, frameNums, nFrames] = decode_hits(values, opts, nPix)
% Same frame decoding as read_sparse_frames_fast, without forming images.
    if isempty(values)
        idx = zeros(0, 1);
        frameNums = zeros(0, 1);
        nFrames = 0;
        return;
    end
    threshold = sparse_separator_threshold_fast(opts, nPix);
    isMarker = values >= cast(threshold, opts.indexClass);

    switch opts.frameCountLocation
        case 'none'
            isHit = ~isMarker;
            allFrames = cumsum(double(isMarker)) + 1;
            nFrames = nnz(isMarker) + double(~isMarker(end));

        case 'start'
            [isHeader, ~] = mark_sparse_count_headers_fast(isMarker, true);
            isHit = ~isMarker & ~isHeader;
            allFrames = cumsum(double(isHeader));
            nFrames = nnz(isHeader);

        case 'end'
            % The last non-marker before a separator is the footer count.
            isFooter = false(size(isMarker));
            if numel(values) > 1
                isFooter(1:end-1) = ~isMarker(1:end-1) & isMarker(2:end);
            end
            if ~isMarker(end)
                isFooter(end) = true;
            end
            isHit = ~isMarker & ~isFooter;
            allFrames = cumsum(double(isMarker)) + 1;
            nFrames = nnz(isFooter);
    end

    frameNums = allFrames(isHit);
    idx = double(values(isHit));
    if opts.pixelBase == 0
        idx = idx + 1;
    end
    if any(idx < 1 | idx > nPix | idx ~= floor(idx))
        error('read_sparse_intensity_correlation:InvalidPixelID', ...
            'Sparse pixel IDs are outside the valid image range.');
    end
    if any(frameNums < 1 | frameNums > nFrames)
        error('read_sparse_intensity_correlation:InvalidFrameNumber', ...
            'Hit records do not belong to a valid frame. Check header/footer options.');
    end
end
