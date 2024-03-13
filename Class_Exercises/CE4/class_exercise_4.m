clear all;

% class exercise 4
% least squares and bayes theorem

% our truth of state x
x_t = 0.5;

% our first guess of state x
x_b = 3;    
% our estimate of the error of our first guess of x
sigma_b = 3; 

% our observation
y_o = 0;
% our estimate of its error
sigma_o = 1;

%==========================================================================
% find x_h and sigma_h using least squares with 1 obs
%==========================================================================
var_b = sigma_b^2;
var_o = sigma_o^2;

a_b = (1/var_b)/( (1/var_b)+(1/var_o) );
a_o = 1-a_b;

x_h     = a_b*x_b + a_o*y_o;
invvar  = 1/var_b + 1/var_o;
var_h   = inv(invvar);
sigma_h = sqrt(var_h);

disp('LS state estimates');
disp([x_t, x_b, y_o, x_h]);
disp('LS error estimates');
disp([sigma_b, sigma_o, sigma_h]);

%==========================================================================
% find x_h and sigma_h using bayes theorem with 1 obs
%==========================================================================
% how many samples?
nsamples=1000;

% draw your state space
x = linspace(-10, 10, nsamples);

% find p(x) --our initial (prior) belief
p_x = normpdf(x,x_b,sigma_b);
% find p(y|x) obs likelihood  -- our observations
p_y_x = normpdf(x, y_o, sigma_o);
% find p(x|y) -- product of our prior and likelihood (posterior)
p_x_y = p_y_x.*p_x;

% E(x|y) is the peak of the distribution (normal)
x_h = x(p_x_y == max(p_x_y));
% or the summation of the product of x and p(x|y)
x_h = sum(x.*(p_x_y/sum(p_x_y)));
% the standard deviation is the summation of the product of square
% deviation and p(x|y)
sigma_h = sqrt( sum( ( (x-x_h).^2).*  (p_x_y/sum(p_x_y)) ) );

disp('Bayes state estimates');
disp([x_t, x_b, y_o, x_h]);
disp('Bayes error estimates');
disp([sigma_b, sigma_o, sigma_h]);

% plot
figure(1)
clf(1)
C=[0,136,55; 202,0,32; 56,108,176]/255; 
plot(x,p_x/sum(p_x),'Color',C(1,:),'Linewidth',2);
hold on;
plot(x,p_y_x/sum(p_y_x),'r','Linewidth',2);
hold on;
plot(x,p_x_y/sum(p_x_y),'b','Linewidth',2);
grid on;
set(gca,'Fontsize',20);
xlabel('x','Fontsize',22);
ylabel('pdf','Fontsize',22);
title('Bayes','Fontsize',22);
legend('Prior','Obs Likelihood','Posterior')
xlim([-10 10]);

%---------------------------------------------
% now what if we have nobs observations of x
% initialize seed (for reproducibity purposes)
j=1;
randn('state',j);

% if first, generate samples (obs) of x_t (y_o)
nobs = 4;
gen_obs =0;
do_f=0;
do_order=0;

if (gen_obs == 1)
    e_o = 0+sigma_o*randn(nobs,1);
    y_o = x_t+e_o;

    if (nobs==4)
        save y_o_4.dat y_o -ASCII;
    elseif (nobs==1000 & sigma_o ==1)
        save y_o_nobs.dat y_o -ASCII;
    elseif (nobs==1000 & sigma_o ==0.5)
        save y_o_nobs_05.dat y_o -ASCII; 
    end

else
    if (nobs==4)
        load y_o_4.dat;
        y_o=y_o_4;
    elseif (nobs==1000 & sigma_o ==1)
        load y_o_nobs.dat;
        y_o=y_o_nobs;
    elseif (nobs==1000 & sigma_o ==0.5)
        load y_o_nobs_05.dat;
        y_o=y_o_nobs_05;        
    end
end

%==========================================================================
% find x_h using the least squares and nobs
%==========================================================================
x_b_old = x_b(1);
sigma_b_old = sigma_b(1);
var_o = sigma_o^2;

x_b = zeros(nobs,1);
sigma_b = zeros(nobs,1);
x_h = zeros(nobs,1);
sigma_h = zeros(nobs,1);

if (do_order==1)
    y_o=y_o( randperm(nobs) );
end

% loop through all obs
for i=1:nobs
    if do_f==1
        switch i
            case{1,2}
                var_o = 1^2;
            case{3,4}
                var_o = 0.5^2;
        end
    end
    
    var_b = sigma_b_old^2;   
    a_b = (1/var_b)/( (1/var_b)+(1/var_o) );
    a_o = 1-a_b;

    x_h(i)     = a_b*x_b_old + a_o*y_o(i);
    invvar  = 1/var_b + 1/var_o;
    var_h   = inv(invvar);
    sigma_h(i) = sqrt(var_h);
    
    % update first guess
    x_b(i) = x_b_old;
    sigma_b(i) = sigma_b_old;
    
    x_b_old = x_h(i);
    sigma_b_old = sigma_h(i);
end
    
disp('LS state estimates');
disp([x_t, x_b(nobs), y_o(nobs), x_h(nobs)]);
disp('LS error estimates');
disp([sigma_b(1), sigma_o, sigma_h(nobs)]);


%==========================================================================
% find x_h using the bayes with nobs
%==========================================================================
x_b_old = x_b(1);
sigma_b_old = sigma_b(1);
var_o = sigma_o^2;

x_b = zeros(nobs,1);
sigma_b = zeros(nobs,1);
x_h = zeros(nobs,1);
sigma_h = zeros(nobs,1);

x = linspace(-10, 10, nsamples);

% loop through all obs
figure(2)
clf(2)
do_plot=1;
if nobs==4
    do_plot=1;
end
    
for i=1:nobs
    if do_f==1
        switch i
            case{1,2}
                sigma_o = 1;
            case{3,4}
                sigma_o = 0.5;
        end
    end
    p_x = normpdf(x,x_b_old,sigma_b_old);
    p_y_x = normpdf(x, y_o(i), sigma_o);
    p_x_y = p_y_x.*p_x;
 
    x_h(i) = sum(x.*(p_x_y/sum(p_x_y)));
    sigma_h(i) = sqrt( sum( ( (x-x_h(i)).^2).*  (p_x_y/sum(p_x_y)) ) );

    % update first guess
    x_b(i) = x_b_old;
    sigma_b(i) = sigma_b_old;
    
    x_b_old = x_h(i);
    sigma_b_old = sigma_h(i);
    
    
    switch i
        case (1)
            ii=1;
            do_plot=1;
        case (2)
            ii=2;
            do_plot=1;
        case (nobs-1)
            ii=3;
            do_plot=1;
        case (nobs)
            ii=4;
            do_plot=1;
        otherwise
            do_plot=0;
    end
        
    if (do_plot==1)
        
            subplot(2,2,ii)
            plot(x,p_x/sum(p_x),'Color',C(1,:),'Linewidth',2);
            hold on;
            plot(x,p_y_x/sum(p_y_x),'r','Linewidth',2);
            hold on;
            plot(x,p_x_y/sum(p_x_y),'b','Linewidth',2);
            grid on;
            set(gca,'Fontsize',20);
            xlabel('x','Fontsize',22);
            ylabel('pdf','Fontsize',22);
            title(['Obs ',num2str(i)],'Fontsize',22);
            legend('Prior','Obs Likelihood','Posterior','Location','NorthWest');
            xlim([-10 10]);
            ylim([0 0.03]);
    end
        

end

disp('Bayes state estimates');
disp([x_t, x_b(1), y_o(nobs), x_h(nobs)]);
disp('Bayes error estimates');
disp([sigma_b(1), sigma_o, sigma_h(nobs)]);
