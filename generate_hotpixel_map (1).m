function [hot_pixels,info] = generate_hotpixel_map(host,port,nFrames,exposureUs,thresholds,statsMode)
%GENERATE_HOTPIXEL_MAP Dark calibration: Poisson variance/mean threshold scan.
%   [hot_pixels,info] = generate_hotpixel_map('127.0.0.1',9999,1e6,4)
%   Cover the sensor completely before calling. Uses the existing spad_cam.
%   hot_pixels is 512x512 logical: 1=usable, 0=hot, before ROI/binning.
%   thresholds: optional count thresholds; [] tests every distinct count level.
%   statsMode: 'retained' (default) omits rejected pixels from statistics;
%              'zeroed' includes rejected pixels as artificial zeros.
%   Both modes retain genuine zero-count pixels. No Statistics Toolbox needed.
%   Save with: writematrix(uint8(hot_pixels),'spad_hot_pixels_512.txt','Delimiter','tab')
    if nargin < 1 || isempty(host), host = '127.0.0.1'; end
    if nargin < 2 || isempty(port), port = 9999; end
    if nargin < 3 || isempty(nFrames), nFrames = 1e6; end
    if nargin < 4 || isempty(exposureUs), exposureUs = 4; end
    if nargin < 5, thresholds = []; end
    if nargin < 6 || isempty(statsMode), statsMode = 'retained'; end
    validateattributes(nFrames,{'numeric'},{'scalar','integer','positive','finite'});
    validateattributes(exposureUs,{'numeric'},{'scalar','real','positive','finite'});
    if ~isempty(thresholds)
        validateattributes(thresholds,{'numeric'},{'vector','real','nonnegative','finite'});
    end
    statsMode = validatestring(statsMode,{'retained','zeroed'});
    minKeepFraction = 0.5; % Reject fits that discard most of the sensor; adjustable.

    cam = spad_cam(host,port);
    cleanup = onCleanup(@() delete(cam));
    cam.stream = 1; cam.bits = 1; cam.sparse = 0;
    cam.measurement_time = exposureUs; cam.overlap = 1;
    cam.frames = nFrames; cam.frame_cap = 100000; cam.batch_size = 64;
    cam.image_width = 512; cam.roi = []; cam.binning = [1 1];
    cam.hot_pixels_lin = [];       % Calibration must observe previously hot pixels too.
    cam.correlation_mask = []; cam.mask_intensity = false;
    cam.correlate = false; cam.live_plot = false; cam.store_frames = false;
    cam.use_gpu = false;           % Can be set true if a supported GPU is available.
    cam.precision = 'double'; cam.verbose = true;
    fprintf('Dark calibration: %d frames, %.6g us exposure, full 512x512 sensor.\n', ...
        nFrames,exposureUs);
    [~,measurement] = cam.imageAcquisition();
    if ~measurement.complete || measurement.nFrames ~= nFrames
        error('generate_hotpixel_map:Incomplete','Dark acquisition did not complete.');
    end
    darkSum = double(measurement.I);

    % Sorting and prefix sums evaluate all thresholds without rescanning the image.
    counts = sort(darkSum(:));
    pixelCount = numel(counts);
    if isempty(thresholds)
        [thresholds,nKept] = unique(counts,'last');
    else
        thresholds = unique(floor(double(thresholds(:))));
        bins = histcounts(counts,[-Inf; thresholds+0.5; Inf]);
        nKept = cumsum(bins(1:end-1)).';
    end
    sumCounts = [0; cumsum(counts)];
    sumSquares = [0; cumsum(counts.^2)];
    s1 = sumCounts(nKept+1); s2 = sumSquares(nKept+1);
    if strcmp(statsMode,'retained')
        n = nKept;
    else
        n = pixelCount*ones(size(nKept));
    end
    mu = s1./max(n,1);
    variance = max(0,(s2-s1.^2./max(n,1))./max(n-1,1)); % Sample variance: var(x,0).
    fano = variance./mu;           % Poisson target: variance/mean = 1.
    score = abs(fano-1);
    valid = nKept >= max(2,ceil(minKeepFraction*pixelCount)) & n >= 2 & ...
        mu > 0 & isfinite(score);
    score(~valid) = Inf;           % Excludes all-zero and nearly-all-masked solutions.
    bestScore = min(score);
    if ~isfinite(bestScore)
        error('generate_hotpixel_map:NoValidThreshold', ...
            'No valid positive-mean fit. Increase exposure/frames or widen the threshold range.');
    end
    best = find(score == bestScore,1,'last'); % Prefer fewer exclusions for tied fits.
    threshold = thresholds(best);
    hot_pixels = darkSum <= threshold;
    info = struct('darkSum',darkSum,'threshold',threshold,'mean',mu(best), ...
        'variance',variance(best),'fano',fano(best),'relativeMismatch',bestScore, ...
        'hotPixelCount',nnz(~hot_pixels),'nFrames',measurement.nFrames, ...
        'exposureUs',exposureUs,'statsMode',statsMode,'minKeepFraction',minKeepFraction);
    info.scan = table(thresholds,mu,variance,fano,nKept/pixelCount,score, ...
        'VariableNames',{'Threshold','Mean','Variance','Fano','RetainedFraction','Score'});
    fprintf('Threshold: %g counts; mean %.6g; variance %.6g; Fano %.6g.\n', ...
        threshold,info.mean,info.variance,info.fano);
    fprintf('Hot pixels: %d / %d (%.3f%%). Statistics mode: %s.\n', ...
        info.hotPixelCount,pixelCount,100*info.hotPixelCount/pixelCount,statsMode);

    figure('Name','Dark hot-pixel calibration','NumberTitle','off');
    layout = tiledlayout(1,2);
    ax = nexttile(layout);
    plotValid = valid & thresholds > 0 & fano > 0;
    loglog(ax,thresholds(plotValid),fano(plotValid),'LineWidth',1.2);
    yline(ax,1,'k--','Poisson'); xline(ax,threshold,'r--','Selected');
    xlabel(ax,'Threshold (summed counts)'); ylabel(ax,'Variance / mean'); grid(ax,'on');
    title(ax,sprintf('Best Fano = %.4g',info.fano));
    ax = nexttile(layout);
    imagesc(ax,hot_pixels); axis(ax,'image'); colormap(ax,'gray'); ax.CLim = [0 1];
    colorbar(ax); title(ax,sprintf('1 = usable; 0 = hot (%d pixels)',info.hotPixelCount));
end
