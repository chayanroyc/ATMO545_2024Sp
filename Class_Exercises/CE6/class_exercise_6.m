
%-----------------------------------------------------
%   Class Exercise 6
%   (On constructing matrices and vectors in OI)
%-----------------------------------------------------

clear all;

% let x be true temperature states
xt = [294 292 295]';

% background temperature states
xb=[298.322,294.569,293.819]';
eb = abs(xb-xt);

% construct background error covariance
sigmab=5;
rhob=0.5;

sigmab_e=sigmab;
sigmab_f=sigmab;
sigmab_g=sigmab;

rho_e_f=rhob;
rho_e_g=0;
rho_f_g=rhob;

B  = [sigmab_e^2 sigmab_e*sigmab_f*rho_e_f sigmab_e*sigmab_g*rho_e_g; ...
      sigmab_e*sigmab_f*rho_e_f sigmab_f^2 sigmab_f*sigmab_g*rho_f_g; ...
      sigmab_e*sigmab_g*rho_e_g sigmab_f*sigmab_g*rho_f_g sigmab_g^2];

% positions of states in km
% ix=[1 2 3];
i_e=1; i_f=2; i_g=3;

% positions of obs in km
% iy=[1.3, 2.4];
i_1 = 1.3; i_2=2.4;

% construct H
H=[(i_f-i_1)/(i_f-i_e),(i_1-i_e)/(i_f-i_e), 0 ; ...
    0, (i_g-i_2)/(i_g-i_f), (i_2-i_f)/(i_g-i_f)];

yt=H*xt; %true state in obs space
yb=H*xb; %background state in obs space

% construct observation vector
yo=[294.273,292.762]';

% construct observation error covariance
sigma_1=1;
sigma_2=1;
R=[sigma_1^2 0; 0 sigma_2^2];

% calculate weights
BHT=B*H';
HBHT=H*B*H';

% OI equations
W=BHT*inv(HBHT+R);
xh=xb+W*(yo-yb);
Ph=(eye(3)-W*H)*B;

% some diagnostics (chi-square)
Jh = (1/6)*( (yo-H*xh)'*inv(R)*(yo-H*xh) + (xh-xb)'*inv(Ph)*(xh-xb));
Jbo = (1/6)*( (yo-H*xb)'*inv(R)*(yo-H*xb) );
Jbb = (1/6)*( (xt-xb)'*inv(B)*(xt-xb) );
Jb = Jbo+Jbb;

Jo=(yo-H*xt)'*inv(R)*(yo-H*xt);
Jho=(yo-H*xh)'*inv(R)*(yo-H*xh);

disp([' Jbo',' Jbb',' Jb']);
disp([Jbo,Jbb,Jb]);

disp([' Jh',' Jb',' Jo',' Jho']);
disp([Jh,Jb,Jo,Jho]);

eh=abs(xh-xt);

disp('estimate, truth, prior ');
disp([xh, xt, xb]);
disp('post error prior error');
disp([eh,eb]);
disp('estimation error');
disp(diag(Ph).^0.5);

