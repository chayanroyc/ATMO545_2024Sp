% Class Exercise 3
% Numerical Representation of Toy Example 
% (influence of observations over space)
% 
% See Problem 1.6 in Chapter 1 Handout

% Tt_0 is the true value of temperature in room 0
% Tt_1 is the true value of temperature in room 1
% Tt_2 is the true value of temperature in room 2

% To_0 is the measured temperature in room 0
% To_1 is the measured temperature in room 1
% To_2 is the measured temperature in room 2

% Tb_0 is the guess temperature in room 0
% Tb_1 is the guess temperature in room 1
% Tb_2 is the guess temperature in room 2

% Set Tt_1 & Tt_1 & Tt_2
Tt_0 = 295.5; % degrees Kelvin (my kind of room temp :) )
Tt_1 = 0.98*Tt_0;
Tt_2 = 1.02*Tt_0;


generate_files = 0;
if (generate_files == 1)

% Set eb_0, eb_1 and eb_2 (note: these are random, unbiased and normally
% distributed -- but now they are correlated)
% Set eo_0, eo_1, and eo_2 as well

sigmab = 3.5; % standard deviation of background
sigmao = 1.0; % standard deviation of observation

rhob_01 = 0.6; % correlation of eb_1 and eb_2
rhob_02 = 0.8; % correlation of eb_0 and eb_2
rhob_12 = 0.4; % correlation of eb_1 and eb_2

nsamples = 10; % number of samples (realization)

% initialize seed (for reproducibity purposes)
j=1;
randn('state',j);

% generate nsamples of eb and eo given their characteristics


eo_0 = 0 + sigmao*randn(nsamples,1);
eo_1 = 0 + sigmao*randn(nsamples,1);
eo_2 = 0 + sigmao*randn(nsamples,1);
eo_11 = 0 + sigmao*randn(nsamples,1);


mu=[0 0 0];
Sigma=[sigmab^2 rhob_01*sigmab^2 rhob_02*sigmab^2; ...
       rhob_01*sigmab^2 sigmab^2 rhob_12*sigmab^2; ...
       rhob_02*sigmab^2 rhob_12*sigmab^2 sigmab^2];
   
R=chol(Sigma);
z = repmat(mu,nsamples,1) + randn(nsamples,3)*R;

eb_0=z(:,1);
eb_1=z(:,2);
eb_2=z(:,3);

% generate nsamples of To_0, To_1, To_2, Tb_0, Tb_1, and Tb_2
To_0 = Tt_0 + eo_0;
To_1 = Tt_1 + eo_1;
To_2 = Tt_2 + eo_2;
To_11 = Tt_1 + eo_11;

Tb_0 = Tt_0 + eb_0;
Tb_1 = Tt_1 + eb_1;
Tb_2 = Tt_2 + eb_2;

else
    
load To_0.dat;
load To_1.dat;
load To_11.dat;
load To_2.dat;

load Tb_0.dat;
load Tb_1.dat;
load Tb_2.dat;

sigmab=std(Tb_0); 
sigmao=std(To_0);

rhob_01 = corr(Tb_0, Tb_1);
rhob_02 = corr(Tb_0, Tb_2);
rhob_12 = corr(Tb_1, Tb_2);

end
    
%--------------------------------------------------------------------------

% Estimate the true temperature in the room Ta_0
alpha = sigmao^2/sigmab^2;

% prob 1 (2 obs 3 model guesses)
% given: To_1, To_2, Tb_0, Tb_1, Tb_2;
w1_1 = ( rhob_01*(1+alpha)-rhob_02*rhob_12) / ( (1+alpha)^2 - rhob_12^2 );
w2_1 = ( rhob_02*(1+alpha)-rhob_01*rhob_12) / ( (1+alpha)^2 - rhob_12^2 );
Ta_0_1 = Tb_0 + w1_1*(To_1-Tb_1) + w2_1*(To_2-Tb_2);

k =  ((1+alpha)*(rhob_01^2 + rhob_02^2)-2*rhob_01*rhob_02*rhob_12);
sigmaa_0_1 = sqrt ( sigmab^2* (1 - k /((1+alpha)^2-rhob_12^2) ) ); 

% prob 2 (1 obs 2 model guesses)
% given To_1, Tb_0, Tb_1
w1_2 = rhob_01/(1+alpha);
Ta_0_2 = Tb_0 + w1_2*(To_1-Tb_1);
sigmaa_0_2 = sqrt( sigmab^2*(1- rhob_01^2/(1+alpha)) );

% prob 0 (1 obs collocated with 1 model guess)
% given To_1, Tb_0

w1_3 = 1/(1+alpha);
Ta_0_3 = Tb_0 + w1_3*(To_0-Tb_0);
sigmaa_0_3 = sqrt( sigmab^2*(1-w1_3) );

% prob 3 (2 collocated obs with 2 model guesses)
% given To_1, To_11, Tb_0, Tb_1
w1_4 = rhob_01/(2+alpha);
w2_4 = w1_4;
Ta_0_4 = Tb_0 + w1_4*(To_1-Tb_1) + w2_4*(To_11-Tb_1);
sigmaa_0_4 = sqrt( sigmab^2*(1- (2*rhob_01^2)/(2+alpha)) );

% prob 4 (2 'wrongly assumed isolated' obs + 2 guesses;
% given To_1, To_2, Tb_0, Tb_1, Tb_2 (assumed rhob_12=0 rhob_01=rhob_02)
w1_5 = rhob_01/(1+alpha);
w2_5 = w1_5;
Ta_0_5 = Tb_0 + w1_5*(To_1-Tb_1) + w2_5*(To_2-Tb_2);
sigmaa_0_5 = sqrt ( sigmab^2*( 1 - (2*rhob_01^2)/(1+alpha)) );



% Diagnostics
figure(1)
clf(1)
plot(Tb_0,'g-','LineWidth',2);
hold on;
plot(To_0,'r-','LineWidth',2);
hold on;
plot(Ta_0_1,'b-','LineWidth',2);
hold on;
plot(Ta_0_2,'b:','LineWidth',2);
hold on;
plot(Ta_0_3,'m-','LineWidth',2);
hold on;
plot(Ta_0_4,'m:','LineWidth',2);
hold on;
plot(Ta_0_5,'c-','LineWidth',2);
hold on;
plot([1:length(Ta_0_1)],Tt_0*ones(length(Ta_0_1),1),'k-','LineWidth',2);
ylim([280 305]);
grid on;
set(gca,'Fontsize',14);
ylabel('Temperature (deg Kelvin) ','Fontsize',14);
xlabel('Realization','Fontsize',14);
legend('Tb_0','To_0','Ta_0_1','Ta_0_2','Ta_0_3','Ta_0_4','Ta_0_5','Tt_0','Location','BestOutside');
title('Temperature Estimates');

% Check if our estimate is more accurate
disp('==============================================');
disp('mean square errors for Tb_0 and Ta_0_1 relative to Tt');
disp([sqrt(mean((Tb_0-Tt_0).^2)),(mean((Ta_0_1-Tt_0).^2))]);
disp(' mean square errors for Tb_0 and Ta_0_1 relative to To_0');
disp([(mean((Tb_0-To_0).^2)),(mean((Ta_0_1-To_0).^2))]);
disp('variance of Tb_0 and Ta_0_1');
disp([sigmab^2, sigmaa_0_1^2]);
disp('mean Tb_0, mean Ta_0_1');
disp([mean(Tb_0), mean(Ta_0_1)]);

disp('==============================================');
disp('mean square errors for Tb_0 and Ta_0_2');
disp([mean((Tb_0-Tt_0).^2),mean((Ta_0_2-Tt_0).^2)]);
disp('mean square errors for Tb_0 and Ta_0_2 relative to To_0');
disp([mean((Tb_0-To_0).^2),mean((Ta_0_2-To_0).^2)]);
disp('variance of Tb_0 and Ta_0_2');
disp([sigmab^2, sigmaa_0_2^2]);
disp('mean Tb_0, mean Ta_0_2');
disp([mean(Tb_0), mean(Ta_0_2)]);

disp('==============================================');
disp('mean square errors for Tb_0 and Ta_0_3');
disp([mean((Tb_0-Tt_0).^2),mean((Ta_0_3-Tt_0).^2)]);
disp('mean square errors for Tb_0 and Ta_0_3 relative to To_0');
disp([mean((Tb_0-To_0).^2),mean((Ta_0_3-To_0).^2)]);
disp('variance of Tb_0 and Ta_0_3');
disp([sigmab^2, sigmaa_0_3^2]);
disp('mean Tb_0, mean Ta_0_3');
disp([mean(Tb_0), mean(Ta_0_3)]);

disp('==============================================');
disp('mean square errors for Tb_0 and Ta_0_4');
disp([mean((Tb_0-Tt_0).^2),mean((Ta_0_4-Tt_0).^2)]);
disp('mean square errors for Tb_0 and Ta_0_4 relative to To_0');
disp([mean((Tb_0-To_0).^2),mean((Ta_0_4-To_0).^2)]);
disp('variance of Tb_0 and Ta_0_4');
disp([sigmab^2, sigmaa_0_4^2]);
disp('mean Tb_0, mean Ta_0_4');
disp([mean(Tb_0), mean(Ta_0_4)]);

disp('==============================================');
disp('mean square errors for Tb_0 and Ta_0_5');
disp([mean((Tb_0-Tt_0).^2),mean((Ta_0_5-Tt_0).^2)]);
disp('mean square errors for Tb_0 and Ta_0_5 relative to To_0');
disp([mean((Tb_0-To_0).^2),mean((Ta_0_5-To_0).^2)]);
disp('variance of Tb_0 and Ta_0_5');
disp([sigmab^2, sigmaa_0_5^2]);
disp('mean Tb_0, mean Ta_0_5');
disp([mean(Tb_0), mean(Ta_0_5)]);


if (generate_files == 1)

save Tb_0.dat Tb_0 -ASCII;
save Tb_1.dat Tb_1 -ASCII;
save Tb_2.dat Tb_2 -ASCII;

save To_0.dat To_0 -ASCII;
save To_1.dat To_1 -ASCII;
save To_2.dat To_2 -ASCII;
save To_11.dat To_11 -ASCII;

end


rmseb1a = sqrt(mean((Tb_0-Tt_0).^2));
rmsea1a = sqrt(mean((Ta_0_1-Tt_0).^2));
rmseb1b = sqrt(mean((Tb_0-To_0).^2));
rmsea1b = sqrt(mean((Ta_0_1-To_0).^2));

rmseb2a = sqrt(mean((Tb_0-Tt_0).^2));
rmsea2a = sqrt(mean((Ta_0_2-Tt_0).^2));
rmseb2b = sqrt(mean((Tb_0-To_0).^2));
rmsea2b = sqrt(mean((Ta_0_2-To_0).^2));

rmseb3a = sqrt(mean((Tb_0-Tt_0).^2));
rmsea3a = sqrt(mean((Ta_0_3-Tt_0).^2));
rmseb3b = sqrt(mean((Tb_0-To_0).^2));
rmsea3b = sqrt(mean((Ta_0_3-To_0).^2));

rmseb4a = sqrt(mean((Tb_0-Tt_0).^2));
rmsea4a = sqrt(mean((Ta_0_4-Tt_0).^2));
rmseb4b = sqrt(mean((Tb_0-To_0).^2));
rmsea4b = sqrt(mean((Ta_0_4-To_0).^2));

rmseb5a = sqrt(mean((Tb_0-Tt_0).^2));
rmsea5a = sqrt(mean((Ta_0_5-Tt_0).^2));
rmseb5b = sqrt(mean((Tb_0-To_0).^2));
rmsea5b = sqrt(mean((Ta_0_5-To_0).^2));

meanTb = mean(Tb_0);
meanTa1 = mean(Ta_0_1);
meanTa2 = mean(Ta_0_2);
meanTa3 = mean(Ta_0_3);
meanTa4 = mean(Ta_0_4);
meanTa5 = mean(Ta_0_5);

mbias1a = meanTa1-Tt_0;
mbias2a = meanTa2-Tt_0;
mbias3a = meanTa3-Tt_0;
mbias4a = meanTa4-Tt_0;
mbias5a = meanTa5-Tt_0;

disp('Statistics for Case 0 Problem 1b, 2b, 3b, 4b');
disp([meanTb, meanTa3,mbias3a, rmsea3a,rmsea3b,sigmaa_0_3,sigmaa_0_3/rmsea3a, sigmaa_0_3/rmsea3b]);
disp('Statistics for Case 2 Problem 1a');
disp([meanTb, meanTa1,mbias1a, rmsea1a,rmsea1b,sigmaa_0_1,sigmaa_0_1/rmsea1a, sigmaa_0_1/rmsea1b]);
disp('Statistics for Case 1 Problem 2a');
disp([meanTb, meanTa2,mbias2a, rmsea2a,rmsea2b,sigmaa_0_2,sigmaa_0_2/rmsea2a, sigmaa_0_2/rmsea2b]);
disp('Statistics for Case 3 Problem 3a');
disp([meanTb, meanTa4,mbias4a, rmsea4a,rmsea4b,sigmaa_0_4,sigmaa_0_4/rmsea4a, sigmaa_0_4/rmsea4b]);
disp('Statistics for Case 4 Problem 4a');
disp([meanTb, meanTa5,mbias5a, rmsea5a,rmsea5b,sigmaa_0_5,sigmaa_0_5/rmsea5a, sigmaa_0_5/rmsea5b]);