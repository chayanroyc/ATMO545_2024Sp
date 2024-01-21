% Class Exercise 2
% Numerical Representation of Toy Example (influence of observations over
% space)
% 

% Tt_1 is the true value of temperature in the room (note: this is grid 0)
% Tt_2 is the true value of temperature in the other room (note: this is grid 1)
% To_2 is the measured temperature in the other room
% Tb_1 is the guess temperature in the room
% Tb_2 is the guess temperature in the other room

% Set Tt_1 & Tt_2
Tt_1 = 297; % degrees Kelvin (my kind of room temp :) )
Tt_2 = 0.98*Tt_1;

% Set eb_1 and eb_2 (note: eb_1 and eb_2 are random, unbiased and normally
% distributed -- but now they are correlated)
% Also set eo_1

sigmab = 2.5; % standard deviation of eb_1 = standard deviation of eb_2
sigmao = 1.0; % standard deviation of eo_1
rhob_12 = 0.9; % correlation of eb_1 and eb_2

nsamples = 1000; % number of samples (realization)

% initialize seed (for reproducibity purposes)
j=1;
randn('state',j);

% generate nsamples of eb_1 and eb_2 and er_2 given their characteristics
eo_2 = 0 + sigmao*randn(nsamples,1);

mu=[0 0];
Sigma=[sigmab^2 rhob_12*sigmab^2; rhob_12*sigmab^2 sigmab^2];
R=chol(Sigma);

z = repmat(mu,nsamples,1) + randn(nsamples,2)*R;
eb_1=z(:,1);
eb_2=z(:,2);

% generate nsamples of To_2, Tb_1 and Tb_2
To_2 = Tt_2 + eo_2;
Tb_1 = Tt_1 + eb_1;
Tb_2 = Tt_2 + eb_2;
%--------------------------------------------------------------------------

% Estimate the true temperature in the room Ta_1 as a linear combination
% of 2 pieces of information, see solution in lecture/handout i.e.

alpha= sigmao.^2/sigmab.^2;
w1 = rhob_12/(1+alpha);
Ta_1 = Tb_1 + w1*(To_2-Tb_2);
sigmaa_1_square = sigmab.^2 * ( 1 - (rhob_12^2)/(1+alpha) );
sigmaa_1 = sqrt(sigmaa_1_square);


% now assume that our knowledge of rhob_12 is wrong (infa=1 is correct)
infa=1; % 0.5
rhob_12a=rhob_12*infa;

% now assume that our knowledge of sigmab is wrong (infb=1 is correct)
infb=1;
sigmaba=sqrt(infb)*sigmab;

alpha= sigmao.^2/sigmaba.^2;
w1 = rhob_12a/(1+alpha);
Ta_1 = Tb_1 + w1*(To_2-Tb_2);
sigmaa_1_square = sigmaba.^2 * ( 1 - (rhob_12a^2)/(1+alpha) );
sigmaa_1 = sqrt(sigmaa_1_square);

disp('==================================================================');

figure(1)
clf(1)
plot(Tb_1,Tb_2,'k.','Markersize',10);
ylabel('Guess Temperature at Grid 1, (deg K)','Fontsize',20);
xlabel('Guess Temperature at Grid 2, (deg K)','Fontsize',20);
grid on;
set(gca,'Fontsize',16);
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
title('Sanity Check','Fontsize',25);

% Diagnostics
figure(2)
clf(2)
C=[0,136,55; 202,0,32; 56,108,176]/255; 
plot(Tb_2,'m.','Markersize',10);
hold on;
plot(To_2,'r.','Markersize',10);
hold on;
plot(Tb_1,'.','Color',C(1,:)','Markersize',10);
hold on;
plot(Ta_1,'b.','Markersize',10);
hold on;
plot([1:length(Ta_1)],Tt_1*ones(length(Ta_1),1),'k-','Linewidth',2);
hold on;
plot([1:length(Ta_1)],Tt_2*ones(length(Ta_1),1),'k:','Linewidth',2);
ylim([280 305]);
grid on;
set(gca,'Fontsize',14);
ylabel('Temperature (deg Kelvin) ','Fontsize',20);
xlabel('Realization','Fontsize',20);
hl=legend('Tb_2','To_2','Tb_1','Ta_1','Tt_1','Tt_2','Location','BestOutside');
set(hl,'Fontsize',16);
set(hl,'Box','off');
set(gca,'YMinorTick','on');
set(gca,'XMinorTick','on');
title('DA Performance','Fontsize',25);

% Check if our estimate is more accurate
% Check if our estimate is more accurate
disp('root mean square errors for Tb_1 and Ta_1');
disp([sqrt( mean((Tb_1-Tt_1).^2)), sqrt(mean((Ta_1-Tt_1).^2))]);
disp('sample Tb_1, sample Ta_1, Tt_1');
disp([Tb_1(1,1), Ta_1(1,1), Tt_1]);
disp('std deviation of error estimates, sigmab, sigmaba, sigmao, sigmaa_1, rhob_12, rhob_12a');
disp([sigmab, sigmaba, sigmao, sigmaa_1, rhob_12, rhob_12a]);