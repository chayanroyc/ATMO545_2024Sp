%-----------------------------------------------------
%   Class Exercise 7A
%   (On Weight Matrix)
%   Based on Daley's Fig. 4.7
%-----------------------------------------------------
clear all;
clf;

%   Define parameters
x1=-2.0;				% position of obs 1
x0=0.0;                 % position of analysis point
x=-4.0:0.1:4.0;         % x-grid
x2=x;					% position of obs 2
alfa=0.25;              % alfa = var(obs err)/var(bkgd err)

%   Preliminary calculations
alfp1=1.0+alfa;
alfp1sq=alfp1^2;
d10=abs(x1-x0);         % normalized distance between points 1 and 0
d12=abs(x1-x2);         % normalized distance between points 1 and 2
d20=abs(x2-x0);         % normalized distance between points 2 and 0


%
%   Calculate error correlations of background at obs 1 and obs 2 locations
%
rho10=(1+d10)*exp(-d10);
rho20=(1+d20).*exp(-d20);
rho12=(1+d12).*exp(-d12);

%
%   Calculate weights for obs 1 and obs 2
%
w1a=rho10*alfp1 - rho12.*rho20;
w2a=rho20*alfp1 - rho12*rho10;
denom=alfp1sq-rho12.*rho12;
w1=w1a./denom;
w2=w2a./denom;

%
%   Now calculate the analysis error variance
%   normalized by the background error variance
%
tmp=(alfp1*(rho10^2+(rho20.*rho20))-2*rho10*rho20.*rho12)./denom;
epssq=1-tmp;

%
%   Plot results
%
plot(x,w1,'b-','Linewidth',2);
hold on;
plot(x,w2,'r-','Linewidth',2);
hold on;
plot(x,epssq,'k-','Linewidth',2)
xlabel('x/L','Fontsize',14);
ylabel('normalized error','Fontsize',14);
set(gca, 'Fontsize',14);

hold on;
n=size(x,2);
zer=zeros(n,1);
plot(x,zer,'k:','Linewidth',1);
yzero=-0.2:0.01:1.2;
xzero=yzero*0;
plot(xzero,yzero,'k--','Linewidth',1);
xzero=xzero+x1;
plot(xzero,yzero,'k--','Linewidth',1);
legend('w1','w2','analy err var','zero line','x1 position','x0 position ');
grid on;

figure(2)
clf(2)
plot(x,rho20,'k-','Linewidth',2);
hold on;
plot(x,rho12,'r-','Linewidth',2);
hold on;
plot(x1,rho10,'b.','Markersize',30);
xlabel('x/L','Fontsize',14);
ylabel('correlation','Fontsize',14);
set(gca, 'Fontsize',14);
legend('rho_{20}','rho_{12}','rho_{10}');
grid on;

