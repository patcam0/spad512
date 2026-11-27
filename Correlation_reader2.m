%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CrossCorrelation_reader.m                 
%   Program to read Imarg from I4D directly    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D4D]=Correlation_reader2(I4D,Idis)

% Inputs
% - I4D is the 4D matrix, represented under a form of a matrix
% - Idis is just an direct image, under the form of a matrix (It just use
% to get the dimensions of the ROI)


%% Settings
[DY,DX]=size(Idis);


%% Calculate the cross correlation image from the different matrices and contribution

KK=0;
D4D=zeros(2*DY-1,2*DX-1);

for K1=1:1:DX
    K11=0;
    KK=KK+1;
    for K2=DX:-1:(DX-K1+1)
        v1=((K1+K11-1)*DY+1):1:(DY*(K1+K11));
        v2=((K2-1)*DY+1):1:(DY*(K2));        
    D4D(:,KK)=D4D(:,KK)+Sumdiag(I4D(v1,v2))';
        K11=K11-1;
    end
end
for K1=2:1:DX
    K11=0;
    KK=KK+1;
    for K2=1:1:(DX-K1+1)
        v1=((K1+K11-1)*DY+1):1:(DY*(K1+K11));
        v2=((K2-1)*DY+1):1:(DY*(K2));        
    D4D(:,KK)=D4D(:,KK)+Sumdiag(I4D(v1,v2))';
        K11=K11+1;
    end
end

 
  end

 