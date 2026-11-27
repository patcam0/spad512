%
%   The idea is to create a virtual disk in the PC memory and test the
%   write and read time from the camera. It seems that direct streaming to
%   matlab is slow, and limited to dense frame encoding, so perhaps it is
%   more efficient to do it this way.
%


%%
% clear all
addpath(genpath('C:\Users\Lab\Dropbox\Patrick PhD\Code\MATLAB\Functions'))
addpath(genpath('D:\Patrick\Code\Matlab\equipment'))
addpath(genpath('C:\Users\Lab\Dropbox\Postdoc Naples\Code\Matlab\Superkick\superkick functions'))
hot_pixels = readmatrix('\\NAS-Cardano-NEW\Instruments Control\SPAD512\spad_hot_pixels_512.txt');

%%
host = '127.0.0.1';  % More reliable than 'localhost' if the server uses IPv4
port = 9999;
spad_cam_obj = spad_cam(host,port);
%%
spad_cam_obj.verbose = true;
% Define acquisition parameters
spad_cam_obj.bits = 1;
spad_cam_obj.measurement_time = 0.2;
spad_cam_obj.frames = 100000;
spad_cam_obj.frame_trigger = 0;
spad_cam_obj.overlap = 1;
spad_cam_obj.sparse = 1;
spad_cam_obj.stream = 0;

ram_disk_path = ['F:\'];
spad_cam_obj.setSaveFolder(ram_disk_path);

%% 1. Separate read and write speeds, no processing.

width_range = 2.^(1:9);
eff_write_times = zeros(size(width_range));
eff_read_times = zeros(size(width_range));
try
    rmdir([ram_disk_path 'data'],'s')
end

counter=0;
for width = width_range
    counter = counter + 1;
    % Set image width
    spad_cam_obj.image_width=width;

    % Capture
    tic
    spad_cam_obj.imageAcquisition()
    time_acq_and_write = toc;
    eff_write_times(counter) = time_acq_and_write;
    pause(0.2)

    % Read
    filenames = listFilenames_fast([spad_cam_obj.save_path 'data\intensity_images\acq00000']);
    tic
    for fileIdx = 1:numel(filenames)
        filename = fullfile([spad_cam_obj.save_path 'data\intensity_images\acq00000'], filenames{fileIdx});
        [frames, nFrames] = read_sparse_frames_fast(filename, [512 spad_cam_obj.image_width]);
        % frames = frames.^2;
        size(frames)
        clear frames
    end
    eff_read_times(counter) = toc;
    % if counter~=length(width_range)
    pause(0.1)
    rmdir([ram_disk_path 'data'],'s')
    % end
    disp(['Time to acquire and write ' num2str(spad_cam_obj.frames) ' frames: ' num2str(time_acq_and_write,'%.2f') ' s.'])
    disp(['Effective frame rate: ' num2str(spad_cam_obj.frames/time_acq_and_write,'%.2f') ' fps.'])
    disp(['Time to read ' num2str(spad_cam_obj.frames) ' frames: ' num2str(1000*eff_read_times(counter),'%.2f') ' ms.'])

end
%%
time_max = max(max(eff_write_times(:)),max(eff_read_times(:)));
time_min = min(min(eff_write_times(:)),min(eff_read_times(:)));
figure

subplot(2,1,1)
semilogx(width_range,eff_write_times./spad_cam_obj.frames,'linewidth',2,'DisplayName','Write')
hold on
semilogx(width_range,eff_read_times./spad_cam_obj.frames,'linewidth',2,'DisplayName','Read')
ylabel('Time per frame (s)')
xlabel('Image width')
title('Time per frame')

ylim([0 time_max./spad_cam_obj.frames])
legend

subplot(2,1,2)
loglog(width_range,spad_cam_obj.frames./eff_write_times,'linewidth',2,'DisplayName','Write fps')
hold on
loglog(width_range,spad_cam_obj.frames./eff_read_times,'linewidth',2,'DisplayName','Read fps')
xlabel('Image width')
ylabel('Effective frame rate (Hz)')
ylim([0 spad_cam_obj.frames./time_min])
legend
title('Frame rate')
sgtitle('Virtual ram disk')

%% 2. Processing overhead time
%% 2a. prepare acquisition: hot pixels map
spad_cam_obj.measurement_time = 100;
spad_cam_obj.frames = 10000;
spad_cam_obj.image_width=512;

try
    rmdir([ram_disk_path 'data'],'s')
end
% Take an intensity image:
spad_cam_obj.imageAcquisition()
filenames = listFilenames_fast([spad_cam_obj.save_path 'data\intensity_images\acq00000']);
intensity_short = zeros([512,spad_cam_obj.image_width]).';
for fileIdx = 1:numel(filenames)
    filename = fullfile([spad_cam_obj.save_path 'data\intensity_images\acq00000'], filenames{fileIdx});
    [frames, nFrames] = read_sparse_frames_fast(filename, [512 spad_cam_obj.image_width],'row-major');
    intensity_short=intensity_short+squeeze(sum(single(frames),3));
end
intensity_short(1:2,:) = 0;
%%
mask = zeros(size(intensity_short));
mask(intensity_short>22) = 1;
% intensity_short(intensity_short>0.02*max(intensity_short(:))) = 0;
% [roi, rectangleHandle] = select_intensity_roi(intensity_short);
figure
subplot(2,1,1)
imagesc(intensity_short)
axis image
colormap hot

subplot(2,1,2)
imagesc(mask)
axis image
colormap hot
sum(mask,'all')/512^2
%%
writematrix(1-mask,'\\NAS-Cardano-NEW\Instruments Control\SPAD512\spad_hot_pixels_512.txt')

%%


width_range = 2.^(1:9);
eff_write_times = zeros(size(width_range));
eff_read_times = zeros(size(width_range));
try
    rmdir([ram_disk_path 'data'],'s')
end

counter=0;
for width = width_range
    counter = counter + 1;
    % Set image width
    spad_cam_obj.image_width=width;

    % Capture
    tic
    spad_cam_obj.imageAcquisition()
    time_acq_and_write = toc;
    eff_write_times(counter) = time_acq_and_write;
    pause(0.2)

    % Read
    filenames = listFilenames_fast(spad_cam_obj.save_path);
    tic
    for fileIdx = 1:numel(filenames)
        filename = fullfile(spad_cam_obj.save_path, filenames{fileIdx});
        [frames, nFrames] = read_sparse_frames_fast(filename, spad_cam_obj.image_width);
        frames = frames.^2;
        clear frames
    end
    eff_read_times(counter) = toc;
    % if counter~=length(width_range)
    pause(0.1)
    rmdir([ram_disk_path 'data'],'s')
    % end
    disp(['Time to acquire and write ' num2str(spad_cam_obj.frames) ' frames: ' num2str(time_acq_and_write,'%.2f') ' s.'])
    disp(['Effective frame rate: ' num2str(spad_cam_obj.frames/time_acq_and_write,'%.2f') ' fps.'])
    disp(['Time to read ' num2str(spad_cam_obj.frames) ' frames: ' num2str(1000*eff_read_times(counter),'%.2f') ' ms.'])

end
%%
filecap = 0;
im_size = [128 512];
hot_pixels = readmatrix('C:\Users\Lab\Dropbox\Postdoc Naples\Code\Matlab\Superkick\Superkick data and analysis\spad_hot_pixels_256.txt');
hot_pixels = hot_pixels(:,64+1:end-64);
%%
hot_pixels_lin = reshape(hot_pixels,im_size(1)*im_size(2),1);
% Process data:
subfolders = listSubfolders_fast([ram_disk_path 'data\intensity_images\']);
n_folders = length(subfolders);
im_stack = zeros([n_folders im_size]);
for folder = 1:10
    disp(['______ Loading folder ' num2str(folder) ' of ' num2str(n_folders) ' ______'])
    if spad_cam_obj.sparse==1
        im_stack(folder,:,:) = read_multiple_sets_bin_fast([char(ram_disk_path) 'data\intensity_images\\' char(subfolders(folder))],im_size,hot_pixels_lin,filecap,'sparse');
    else
        im_stack(folder,:,:) = read_multiple_sets_bin_fast([char(ram_disk_path) 'data\intensity_images\\' char(subfolders(folder))],im_size,hot_pixels_lin,filecap);
    end
end

% Remove first row as it contains strange artefact:
rowsCut = 5;
if spad_cam_obj.sparse==1
    im_stack(:,1:rowsCut,:) = 0;
end
I = squeeze(sum(im_stack,1));
figure
imagesc(I)
axis image
colormap hot

%% Select ROI:
I = squeeze(sum(im_stack,1));
figure
imagesc(squeeze(sum(im_stack,1)))
axis image
colormap hot

% Desired size in pixels
nRows = 128;
nCols = 128;

% Initial position; only translation is allowed
roi = drawrectangle(gca, ...
    'Position', [1 1 nCols nRows], ...
    'InteractionsAllowed', 'translate');

wait(roi);  % Drag to position, then double-click

% Snap position to integer indices and keep within image I
colMin = min(round(roi.Position(1)));
rowMin = min(round(roi.Position(2)));

vecimage_full = 1:prod(im_size);

close gcf

figure
imagesc(I(rowMin:rowMin+nRows-1,colMin:colMin+nCols-1));
axis image
colormap hot


%% Try coincidence extraction:
opts = parse_reader_options_fast;

filecap = 0;

I = zeros(nRows,nCols);
R = zeros(nRows*nCols,nRows*nCols);
tic
frames_per_folder = 100000;
n_frames_tot = 0;
% folderPath = current_path;
for folder = 1:2
    disp(['______ Loading folder ' num2str(folder) ' of ' num2str(n_folders) ' ______'])

    folderPath = [char(ram_disk_path) 'data\intensity_images\\' char(subfolders(folder))];
    filenames = listFilenames_fast(folderPath);

    nFiles = numel(filenames);
    for fileIdx = 1:nFiles
        fprintf('Reading file %d of %d: %s\n', fileIdx, nFiles, filenames{fileIdx});

        filename = fullfile(folderPath, filenames{fileIdx});

        if spad_cam_obj.sparse
            [data, nFrames] = read_sparse_frames_fast(filename, im_size, hot_pixels_lin, opts);
        else
            [data, nFrames] = read_packed_frames_fast(filename, im_size, hot_pixels_lin, opts);
        end
        n_frames_tot = n_frames_tot + nFrames;
        I = I + squeeze(sum(double(data(rowMin:rowMin+nRows-1,colMin:colMin+nCols-1,:)),3));
        data_lin = double(reshape(data(rowMin:rowMin+nRows-1,colMin:colMin+nCols-1,:),[nFrames,nRows*nCols]));
        R = R + data_lin.'*data_lin;
    end
end
read_time = toc;
disp(['Read framerate = ' num2str(n_frames_tot/read_time,'.2f') ' fps'])
% R = R/spad_cam_obj.frames;
% I = I/spad_cam_obj.frames;

%%
tic
[I, R, nFramesTotal, perf] = spad_roi_correlations_fast( ...
    ram_disk_path, subfolders(1:2), im_size, hot_pixels_lin, ...
    spad_cam_obj.sparse, opts, [rowMin colMin nRows nCols], ...
    'UseGPU', 'auto', ...
    'BatchFrames', 4096, ...
    'MatrixFormat', 'auto', ...
    'ReturnOnGPU', false);
t = toc;
time_per_frame = t/nFramesTotal;
disp(['Average processing time for one frame = ' num2str(1e6*time_per_frame,'%.2f') ' us.'])
disp(['Effective frame rate = ' num2str(1/time_per_frame,'%.2f') ' fps.'])

%% Remove bright spots

threshold = 0.12;
mask = ones(size(I));
mask(I<threshold*max(I(:))) = 0;
vec_mask = find(reshape(mask,[1,nRows*nCols])==1);
figure
subplot(1,2,1)
imagesc(I)
axis image
colormap hot
subplot(1,2,2)
imagesc(mask)
axis image
%% Show correlation projections

vec_mask = find(reshape(mask,[1,nRows*nCols])==1);
% vec_mask = ones(1,nRows*nCols);
% I_lin = reshape(I,[1,nRows*nCols]);
R_corr = R;
R_corr(diag(ones(1,nRows*nCols))==1) = 0;
R_corr(vec_mask,:)=0;
R_corr(:,vec_mask) =0;
G2 = R_corr-I_lin(vec_mask)*I_lin(vec_mask).';
G2 = G2-mean(G2,'all');
R_corr = R_corr - mean(R_corr,'all');

% R = R - mean(R,'all');

sum_proj =  Convolution_reader(G2,I);
% sum_proj = imresize(sum_proj,1/2);
% minus_proj = Correlation_reader2(R_corr,I);
% minus_proj(nRows,nCols) = 0;
figure
subplot(1,2,1)
imagesc(I.*(1-mask))
axis image
colormap hot
title('Intensity')

subplot(1,2,2)
imagesc(sum_proj);
axis image
colormap hot
title('Sum-coordinate projection')


%%


R_corr_binned = Macropix4D_2(R_corr,I,nCols,nRows,nCols/4,nRows/4);
% [Imarg,Irow,Icol,Xs,Ys]=Coincidence_and_1dJPD_reader(R_corr,I.*(1-mask),1:nRows*nCols);
[Imarg,Xs,Ys]=conditionalReader(R_corr_binned,imresize(I.*(1-mask),1/4),1:nRows*nCols/16);