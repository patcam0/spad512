hot_pixels_lin = reshape(hot_pixels,512*512,1);
im_size = [512 512];

spad_cam_obj.measurement_time = 0.02;
spad_cam_obj.frames = 1000000;
spad_cam_obj.image_width=512;

% Take an intensity image:
spad_cam_obj.imageAcquisition()

% Full-resolution preview without allocating R:
preview = struct('ComputeCorrelation', false);
Ipreview = read_sparse_intensity_correlation( ...
    filename, im_size);

% Draw/adjust the rectangle, then double-click inside to accept:
roi = select_intensity_roi(Ipreview);
if isempty(roi), return; end

% Apply ROI and 2-by-2 binning:
processing = struct('ROI', roi, ...
                    'BinSize', [2 2], ...
                    'BinMode', 'binary');

[I, R, nFrames, info] = read_sparse_intensity_correlation( ...
    filename, im_size, hot_pixels_lin, opts, processing);