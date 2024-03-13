clear all;

% Class Exercise 5 (Lect 016)
% Numerical Example 02/22/2022 

% Section 2.3.1 (Kriging) of Wikle & Berliner, 2007
% Problem 2 (Case 1)

H=[0 1 0; ...
   0 0 1];

x_mean = [18 18 18]';
y = [16 23]';

P=[ 1 0.61 0.22; ...
    0.61 1 0.37; ...
    0.22 0.37 1];

I2=diag(ones(2,1));
I3=diag(ones(3,1));
R=(0.5)*I2;

K=P*H'*inv(R+H*P*H');

x_est_mean_1 = x_mean + K*(y-H*x_mean);
x_est_var_1  = (I3-K*H)*P;

disp('P');
disp(P);

disp('K');
disp(K)

disp('x_est_1');
disp(x_est_mean_1);
disp('x_est_var_1');
disp(x_est_var_1);

x_est_mean_2= inv(H'*inv(R)*H+inv(P))*(H'*inv(R)*y+inv(P)*x_mean);
x_est_var_2=inv(H'*inv(R)*H+inv(P));
disp('x_est_2');
disp(x_est_mean_2);
disp('x_est_var_2');
disp(x_est_var_2);

% Problem 2 Section 2.3.1 Kriging (case 2 P=I)
H=[0 1 0; ...
   0 0 1];

x_mean = [18 18 18]';
y = [16 23]';

%P=[ 1 0.61 0.22; ...
%    0.61 1 0.37; ...
%    0.22 0.37 1];

P=diag(ones(3,1));

I2=diag(ones(2,1));
I3=diag(ones(3,1));
R=(0.5)*I2;

K=P*H'*inv(R+H*P*H');

x_est_mean_1 = x_mean + K*(y-H*x_mean);
x_est_var_1  = (I3-K*H)*P;

disp('P');
disp(P);

disp('K');
disp(K)

disp('x_est_1');
disp(x_est_mean_1);
disp('x_est_var_1');
disp(x_est_var_1);

% let's do this using a different (equivalent) expression
x_est_mean_2= inv(H'*inv(R)*H+inv(P))*(H'*inv(R)*y+inv(P)*x_mean);
x_est_var_2=inv(H'*inv(R)*H+inv(P));
disp('x_est_2');
disp(x_est_mean_2);
disp('x_est_var_2');
disp(x_est_var_2);