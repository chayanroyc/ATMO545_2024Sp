% Class Exercise 1a
% Numerical Representation of Toy Example (same as CE1 but see line 125)
% 
% See Problem Set in Lecture 003

% Tt is the true value of temperature in the room
% T1 is the measured temperature using thermometer 1 (taiwo's)
% T2 is the measured temperature using thermometer 2 (sam's)

%--------------------------------------------------------------------------
% editable variables
%--------------------------------------------------------------------------

% Set Tt
Tt = 297; % degrees Kelvin (my kind of room temp :) )

% Set e1 and e2 (note: e1 and e2 are random, unbiased and normally
% distributed)

mu1=0;  % mean of e1
mu2=0;  % mean of e2

sigma1 = 1; %2.5; % standard deviation of e1
sigma2 = 2; % standard deviation of e2

nsamples = 1e3; % number of samples (realization)
%--------------------------------------------------------------------------

% generate nsamples of e1 and e2 given their characteristics

% initialize seed (for reproducibility purposes)
j=1;
randn('state',j);

e1 = mu1 + sigma1*randn(nsamples,1);
e2 = mu2 + sigma2*randn(nsamples,1);


% generate nsamples of T1 and T2 (our observations)
T1 = Tt + e1;
T2 = Tt + e2;


% Diagnostics
%--------------------------------------------------------------------------
% check histogram of T1 & T2
figure(1)
clf(1)
subplot(1,2,1);
[n1,x1]=hist(T1,50);
h=stairs(x1,n1./sum(n1));
set(h,'Linewidth',2);
hold on;
plot(Tt,0,'k.','Markersize',50);
xlim([290 305]);
ylim([0 0.1]);
set(gca,'Fontsize',16);
xlabel('Temperature','Fontsize',20);
ylabel('Probability','Fontsize',20);
title('Histogram of T1','Fontsize',25);

subplot(1,2,2);
[n2,x2]=hist(T2,50);
h=stairs(x2,n2./sum(n2));
set(h,'Linewidth',2);
hold on;
plot(Tt,0,'k.','Markersize',50);
xlim([290 305]);
ylim([0 0.1]);
set(gca,'Fontsize',16);
xlabel('Temperature','Fontsize',20);
ylabel('Probability','Fontsize',20);
title('Histogram of T2','Fontsize',25);

figure(2)
clf(2)
plot(T1,T2,'b.','Markersize',25);
hold on;
plot(Tt,Tt,'k.','Markersize',50);
xlim([290 305]);
ylim([290 305]);
grid on;
set(gca,'Fontsize',16);
xlabel('T_1 (deg Kelvin)','Fontsize',20);
ylabel('T_2 (deg Kelvin)','Fontsize',20);
title('T_1 vs T_2','Fontsize',25);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
set(gca,'TIckLength',[0.05 0.05]);
% check statistics
disp('Mean and Std Deviation of e1');
disp([mean(e1), std(e1)]);
disp('Mean and Std Deviation of e2');
disp([mean(e2), std(e2)]);
disp('Correlation between T1 and T2');
disp(corr(T1,T2));


%--------------------------------------------------------------------------

% Estimate the true temperature as a linear combination of 2 pieces of
% information, i.e.

% Ta = a1*T1 + a2*T2 or Ta = T2 + a1*(T1-T2)
% Since we assume e1 and e2 to be unbiased, it follows that 
% E(T1)=E(T2)=E(Ta)=Tt
% E(Ta)=a1*E(T1)+a2*E(T2)
% Tt = a1*Tt + a2*Tt
% 1 = a1 + a2

% Ta will be the best estimate of Tt if coefficients minimize mean square
% error:
% sigmaa^2 = E[(Ta-Tt)^2]
% 
% E[(Ta-Tt)^2]= E[ ( a1*(T1-Tt) + (1-a1)*(T2-Tt) )^2]
% Solve a1, a2, Ta and sigmaa 
% a1 = sigma2^2/(sigma1^2+sigma2^2);
% a2 = sigma1^2/(sigma1^2+sigma2^2);
% Ta = a1*T1+a2*T2;
% sigmaa = sqrt( (sigma1^2*sigma2^2)/(sigma1^2+sigma2^2);

%==========================================================================

% pretend we dont know sigma2, now assume a different sigma2
inf=2; % change o higher/lower than 1 to reflect over/under estimation
sigma2a=sigma2*inf;

a1 = sigma2a^2/(sigma1^2+sigma2a^2);
a2 = sigma1^2/(sigma1^2+sigma2a^2);

Ta = a1*T1+a2*T2;
sigmaa = sqrt( (sigma1^2*sigma2^2)/(sigma1^2+sigma2^2) );

disp('a1, a2 and sigmaa')
disp([a1, a2, sigmaa]);

% Diagnostics
figure(3)
clf(3)
C=[0,136,55; 202,0,32; 56,108,176]/255; 
plot(T1,'.','Color',C(1,:)','Markersize',20);
hold on;
plot(T2,'r.','Markersize',20);
hold on;
plot(Ta,'b.','Markersize',20);
hold on;
plot([1:length(Ta)],Tt*ones(length(Ta),1),'k-');
ylim([290 305]);
grid on;
set(gca,'Fontsize',16);
ylabel('Temperature (deg Kelvin) ','Fontsize',20);
xlabel('Realization','Fontsize',20);
hl=legend('T_1','T_2','T_a','T_t','Location','BestOutside');
set(hl,'Fontsize',16);
set(hl,'Box','off');
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
title('DA Performance','Fontsize',25);

% Check if our estimate is more accurate
disp(['T1, ',' T2,',' Ta,',' Tt']);
disp([mean(T1),mean(T2),mean(Ta), Tt]);
disp('root mean square errors for T1, T2 and Ta');
disp([ sqrt(mean((T1-Tt).^2)),sqrt(mean((T2-Tt).^2)), sqrt(mean((Ta-Tt).^2)) ]);

%--------------------------------------------------------------------------
% we can also find the coefficients by
X=[T1,T2];
Y=Ta;
num1=cov(X(:,1),Y);
b1 = num1(2,1)/var(X(:,1));
num2=cov(X(:,2),Ta);
b2 = num2(2,1)/var(X(:,2));

% or
betas=X\Y;

% disp('betas')
% disp([b1,b2,betas']);
%--------------------------------------------------------------------------