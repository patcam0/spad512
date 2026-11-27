function [S,R,Ntotal,totalTime] = spad_correlation_example(hermesCam,hotPixelsCol,nexposure,nframes,Nset,Nsave,savepath,vecimage1,vecimage2)
% CORRELATIONSPADMPD_HERMES_SNAP JPD measurement with the MPD Hermes SPAD in
% snap mode
%
%   hermesCam: Camera object
%   hotPixelsCol: Location of hot pixels in column form
%   nexposure: exposure time of each frame = nexposure*10ns
%   nframes: number of frames per set
%   Nset: number of sets. With slow cable 1 set ~ 3seconds. With fast
%         cable, 1 ~ 1second
%   Nsave: Save data every Nsave sets. If Nsave = 0, do not save
%   savepath: location to save data
%   vecimage1: Region of photon 1 on the camera
%   vecimage2: Region of photon 2 on the camera 


nintegframes = 1;   %number of integrated frames
ncounters = 1;  
halfar = false;     %this is also possbile / halfar = 'Disabled';
signdata = false;
forc8bit = false;
mode = 'None';

% Change settings and apply them to camera
hermesCam.SetCameraPar(nexposure,nframes,nintegframes, ncounters, forc8bit, halfar, signdata );
hermesCam.ApplySettings();

xdim = 64;
ydim = 32;

S=zeros(length(vecimage1),length(vecimage2));
R=zeros(ydim*xdim,1);

% Its faster to do not use a matlab function here  ~65000 frames/s 
t1=tic;
for kk=1:Nset

% Get snap image
hermesCam.SnapPrepare();
hermesCam.SnapAcquire();

data_ext = hermesCam.SnapGetImageBuffer(); % Get image from buffer

data=double(reshape(data_ext,ydim*xdim,nframes)); %reshape data to get a matrix with Nframes column

data(hotPixelsCol==1,:) = 0; % Set hot pixels to 0
%process data 
data=sparse(data);  %As much pixel are 0 its faster to use the saprse matrix version to do the matrix multiplication
S=S+data(vecimage1,:)*data(vecimage2,:)';
R=R+sum(data,2);
kk

if Nsave
    if mod(kk,Nsave)==0
        save(savepath,'S','R','kk','nexposure','nframes','Nset','xdim','ydim','vecimage1','vecimage2')
    end
end
end
totalTime=toc(t1);
Ntotal=nframes*Nset;
framerate=Ntotal/totalTime;
totalTimeDisp=seconds(totalTime);
totalTimeDisp.Format='hh:mm:ss';
disp(['Total acquisition time (hh:mm:ss): '])
disp(totalTimeDisp)
disp(['Framerate: ' num2str(framerate/1000) 'kpfs'])
R = R/Ntotal;
S = S/Ntotal;
end

