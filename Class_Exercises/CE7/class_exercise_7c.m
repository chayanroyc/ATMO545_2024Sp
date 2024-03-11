%-----------------------------------------------------
%   Class Exercise 7C
%   (On OI for temperature)
%-----------------------------------------------------

clear all;

% stations: ELP, BGS, OKC

% get observations (by interpolation)
yo=[-2 -4.8 -5.2]'; 

% get forecast (by interpolation)
xb=[-0.5 -1.0 -3.6]';

% since we already did the interpolation, construct H as
H=[1 0 0; 0 1 0; 0 0 1];

% calculate distances between stations
%d12=6378*acos(cos(31.8*pi/180)*cos(32.39*pi/180)*cos(-106.4*pi/180+101.48*pi/180)+sin(31.8*pi/180)*sin(32.39*pi/180));
d12=468; % ELP to BGS
d13=926.3; % ELP to OKC
d23=520.9; % BGS to OKC

d12=468; % ELP to BGS
d13=908.5; % ELP to OKC
d23=490.6; % BGS to OKC

% forecast error variance
sigFs=1^2;
% observation error variance
sigRs=0.5^2;

% length scale of background error covariance function
L=556;

% construct B
B=zeros(3,3);
for i=1:3
    for j=1:3
        if (i==j)
            B(i,j)=sigFs;
        elseif (i==1 & j==2)
            B(i,j)=sigFs*exp(-0.5*(d12/L)^2);
        elseif (i==1 & j==3)
            B(i,j)=sigFs*exp(-0.5*(d13/L)^2);
        elseif (i==2 & j==1)
            B(i,j)=sigFs*exp(-0.5*(d12/L)^2);
        elseif (i==2 & j==3)
            B(i,j)=sigFs*exp(-0.5*(d23/L)^2);
        elseif (i==3 & j==1)
            B(i,j)=sigFs*exp(-0.5*(d13/L)^2);
        elseif (i==3 & j==2)
            B(i,j)=sigFs*exp(-0.5*(d23/L)^2);
        end
    end
end

% construct observation error covariance (assume diagonal)
R=sigRs*eye(3);

% calculate the weight matrix
W=B*H'*inv(H*B*H'+R);

% calculate your analysis
xa=xb+B*H'*inv(H*B*H'+R)*(yo-H*xb);
disp(['xa ','xb ','yo ']);
disp([xa,xb,yo]);

% calculate your analysis error
Pa=(eye(3)-W*H)*B;
disp('Pa');
disp(Pa)
disp(['Analysis and Background Errors']);
disp([diag(Pa).^0.5,diag(B).^0.5])
