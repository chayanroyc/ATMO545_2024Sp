%-----------------------------------------------------
%   Class Exercise 7B
%   (On Error Covariances)
%-----------------------------------------------------
clear all;

% wind speed obs in m/s
yo=[9,14]';     
% wind speed guess in m/s
xb=[8,12]';     

% case a
stdB=1.5;
stdR=1;
R=stdR^2*eye(2);
B=stdB^2*eye(2);
H=[1 0; 0 1];
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb,yo]);


% case b
rhoB=[1 0.5;0.5 1];
D=stdB*eye(2);
B=D*rhoB*D;
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb,yo]);

% case c
B=stdB^2*eye(2);
rhoR=rhoB;
D=stdR*eye(2);
R=D*rhoR*D;
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb,yo]);

% case d
rhoB=[1 0.9;0.9 1];
R=stdR^2*eye(2);
D=stdB*eye(2);
B=D*rhoB*D;
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb,yo]);

% case e
rhoR=[1 0.9;0.9 1];
B=stdB^2*eye(2);
D=stdR*eye(2);
R=D*rhoR*D;
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb,yo]);

% case f
yo=9;
H=[1 0];
B=stdB^2*eye(2);
R=stdR^2;
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb]);

% case g
rhoB=[1 0.5;0.5 1];
D=stdB*eye(2);
B=D*rhoB*D;
W=B*H'*inv(H*B*H'+R);
xa=xb+W*(yo-H*xb);
Pa=(eye(2)-W*H)*B;
disp([xa,sqrt(diag(Pa)),xb]);