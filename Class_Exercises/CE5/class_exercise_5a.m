clear all;

% Class Exercise 5 (Lect 009/010)
% Numerical Example 02/7/2024 

% Section 2.1.1 of Wikle & Berliner, 2007
% Problem 1A: 
%  Try to calculate the posterior mean and variance using two
%  approaches (1):  least-squares (2): Bayessian formulation
% Problem 1B:
%  Try to reproduce Figure 1 & 2

% case 1
x_mean = 20;
x_var  = 3;
x_err  = sqrt(x_var);

yobs   = [19, 23]';
n      =length(yobs);

y_var  = 1;
y_err  = sqrt(y_var);

y_mean = mean(yobs);

% from least squares
wy = x_var / ( n*x_var + y_var );
wx = y_var/ (n*x_var + y_var);
x_est_mean_ls   = wy*yobs(1) + wy*yobs(2) + wx*x_mean;
x_est_var_ls    = 1 / (n/y_var + 1/x_var); 

% same as taking the mean of y
%w1=(1/y_var)/( (1/y_var)+(1/y_var) );
%w2=1-w1;
%y_est_mean_ls=w1*yobs(1)+w2*yobs(2); (average)
%y_var_mean = ( (1/y_var)+(1/y_var) )^-1;
%w=(1/y_var_mean)/( (1/x_var)+(1/y_var_mean) );
%x_est_mean_ls   = x_mean + w*(y_mean-x_mean);
%x_est_var_ls = ((1/x_var)+(1/y_var_mean))^-1;

% same as
%
% K_1=n*x_var/(y_var+n*x_var);
% x_est_mean_ls=x_mean + K_1*(y_est_mean_ls-x_mean);
% x_est_var_ls    = (1-K_1)*x_var;

% from bayesian formulation
nsamples=1e3;
x_x=linspace(x_mean-3*x_err, x_mean+3*x_err,nsamples);
p_x=normpdf(x_x,x_mean,x_err);
p_y_x_1=normpdf(x_x,yobs(1),y_err);
p_y_x_2=normpdf(x_x,yobs(2),y_err);
p_y_x = p_y_x_1.*p_y_x_2;
p_x_y_bayes = p_x.*p_y_x;
p_x_y_ls = normpdf(x_x,x_est_mean_ls,(x_est_var_ls));
x_est_mean_bayes = x_x(p_x_y_bayes==max(p_x_y_bayes));

disp('x_est_mean_ls, x_est_mean_bayes');
disp([x_est_mean_ls, x_est_mean_bayes]);

sigma_est_bayes = sqrt( sum( ( (x_x-x_est_mean_bayes).^2).*  (p_x_y_bayes/sum(p_x_y_bayes)) ) );
p_x_y_bayes1=normpdf(x_x, x_est_mean_bayes, sigma_est_bayes^2);

figure(1)
clf(1)
nx=length(x_x);
plot(x_x,p_x./sum(p_x),'k:','LineWidth',2);
hold on;
plot(x_x,p_y_x./sum(p_y_x),'k-.','LineWidth',2);
hold on;
plot(x_x,p_x_y_bayes1./sum(p_x_y_bayes1),'k-','LineWidth',2);
hold on;
plot(x_x,p_x_y_ls./sum(p_x_y_ls),'k-','LineWidth',1);
xlim([14 26]);
ylim([0 0.01]);
set(gca,'Fontsize',14);
xlabel('x','Fontsize',16);
ylabel('pdf','Fontsize',16);
legend('Prior','Likelihood','Posterior','LS');
set(gca,'Fontsize',14);


%--------------------------------------------------------------------------
% case 2
x_mean = 20;
x_var  = 3;
x_err  = sqrt(x_var);

yobs   = [19, 23]';
n      =length(yobs);

y_var  = 10;
y_err  = sqrt(y_var);

y_mean = mean(yobs);

% from least squares
wy = x_var / ( n*x_var + y_var );
wx = y_var/ (n*x_var + y_var);
x_est_mean_ls   = wy*yobs(1) + wy*yobs(2) + wx*x_mean;
x_est_var_ls    = 1 / (n/y_var + 1/x_var); 

% same as taking the mean of y
%w1=(1/y_var)/( (1/y_var)+(1/y_var) );
%w2=1-w1;
%y_est_mean_ls=w1*yobs(1)+w2*yobs(2); (average)
%y_var_mean = ( (1/y_var)+(1/y_var) )^-1;
%w=(1/y_var_mean)/( (1/x_var)+(1/y_var_mean) );
%x_est_mean_ls   = x_mean + w*(y_mean-x_mean);
%x_est_var_ls = ((1/x_var)+(1/y_var_mean))^-1;

% same as
%
% K_2=n*x_var/(y_var+n*x_var);
% x_est_mean_ls=x_mean + K_2*(y_est_mean_ls-x_mean);
% x_est_var_ls    = (1-K_2)*x_var;

% from bayesian formulation
nsamples=1e3;
x_x=linspace(x_mean-3*x_err, x_mean+3*x_err,nsamples);
p_x=normpdf(x_x,x_mean,x_err);
p_y_x_1=normpdf(x_x,yobs(1),y_err);
p_y_x_2=normpdf(x_x,yobs(2),y_err);
p_y_x = p_y_x_1.*p_y_x_2;
p_x_y_bayes = p_x.*p_y_x;
p_x_y_ls = normpdf(x_x,x_est_mean_ls,(x_est_var_ls));
x_est_mean_bayes = x_x(p_x_y_bayes==max(p_x_y_bayes));
sigma_est_bayes = sqrt( sum( ( (x_x-x_est_mean_bayes).^2).*  (p_x_y_bayes/sum(p_x_y_bayes)) ) );
p_x_y_bayes1=normpdf(x_x, x_est_mean_bayes, sigma_est_bayes^2);


disp('x_est_mean_ls, x_est_mean_bayes');
disp([x_est_mean_ls, x_est_mean_bayes]);


figure(2)
clf(2)
nx=length(x_x);
plot(x_x,p_x./sum(p_x),'k:','LineWidth',2);
hold on;
plot(x_x,p_y_x./sum(p_y_x),'k-.','LineWidth',2);
hold on;
plot(x_x,p_x_y_bayes1./sum(p_x_y_bayes1),'k-','LineWidth',2);
hold on;
plot(x_x,p_x_y_ls./sum(p_x_y_ls),'k-','LineWidth',1);
xlim([14 26]);
ylim([0 2.5e-3]);
set(gca,'Fontsize',14);
xlabel('x','Fontsize',16);
ylabel('pdf','Fontsize',16);
legend('Prior','Likelihood','Posterior','LS');
set(gca,'Fontsize',14);