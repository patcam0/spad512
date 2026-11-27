%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convolution_reader.m                   Created by: Hugo Defienne    4/15/15
%   v. 0.1                                  Latest revision:    4/15/15
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D4D]=Convolution_reader(I4D,Idis)

%% Settings
[DY,DX]=size(Idis);
I4D=flipud(I4D);

%% GUI

%% Calculate the cross convolution image from the different matrices and contribution

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

 