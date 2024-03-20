% advect -- program to solve the advection equation
clear all;

% spatial domain (1D)
La=-180;        % lower bound
Lb=180;         % upper bound
dx=1.0;         % delta x
x=La:dx:Lb-dx;  % x axis
N=length(x);    % number of gridpoints

% advection parameters
dt=1.0;             % time step (1 s)
u=1.0;              % advection speed in m/s
C=u*dt/(2.*dx);     % coefficient
nStep=300;          % length of integration ( s)

% construct system noise covariance (assume no noise for now)
var_true=0;                       % true system variance
Ld=0.25*abs(La);                  % correlation length
[Rt,De]=gcorr('gauss',Lb-La,Ld,N,1);   % construct true error correlation
Qt=var_true*Rt;                   % true error covariance

%==========================================================================
% generate truth xt~N(mu,Qt), xf~N(xt,Q)
j=1;
randn('state',j);

% mean of initial true state is a square wave
mu = sqrwv ( x, La, Lb );
xt=mvnrnd(mu,Qt)';                % initial true state
%==========================================================================
% assume noise of the forecast model initial states
varQ=0.5;                         % model error variance  
Q=varQ*Rt;                        % model error covariance  
xf=mvnrnd(xt,Q)';                 % initial forecast state 
Pf=Q;                             % initial forecast error covariance
%==========================================================================

% initialize matrices to be saved
xt_save=zeros(N,nStep);           % true state
xf_save=zeros(N,nStep);           % forecast state

% dynamic model matrix (this is a banded matrix based on different methods)
% Forward in Time Central in Space (finite difference)
% a(i,k+1)=a(i,k) + (-u*dt/2*dx)*( a(i+1,k)-a(i-1,k) )
% Lax Method
% a(i,k+1)=1/2*[a(i+1,k)+a(i-1,k)] + [-u*dt/(2*dx)]*[a(i+1,k)-a(i-1,k)];
% Lax-Wendroff (let C=(u*dt)/(2*dx)
% a(i,k+1)=(C+2*C^2)*a(i-1,k) + (1-4*C^2)*a(i,k) + (2*C^2-C)*a(i+1,k)
%==========================================================================
% initialize
M=zeros(N,N);

% first grid point (periodic domain)
M(1,N)=C+2*C^2;
M(1,1)=1-4*C^2;
M(1,2)=2*C^2-C;
% central grid points
for i=2:N-1
    M(i,i-1)=C+2*C^2;
    M(i,i)=1-4*C^2;
    M(i,i+1)=2*C^2-C;
end
% last grid point
M(N,N-1)=C+2*C^2;
M(N,N)=1-4*C^2;
M(N,1)=2*C^2-C;

% plot M
figure(4);
clf(4);
pcolor(M);
shading flat;
set(gca,'Fontsize',14);
title('Transition Matrix, M');
colorbar;
xlabel('x-axis');
ylabel('y-axis');
set(gca,'TickLength',[0.05 0.05]);
grid on;

%==========================================================================
% MAIN LOOP (time loop)
xa=xf;
for istep=1:nStep
    disp(istep);

    % true state time update
    xt=mvnrnd(M*xt,Qt)';
    xt_save(:,istep)=xt;
    
    % forecast state time update
    xf=M*xa;
    xf_save(:,istep)=xf;

    % without observations, xf is equal to xa
    xa=xf;
    
end
%==========================================================================
% make a movie of the true and forecast state
figure(12)
clf(12)
for istep=1:nStep
    subplot(2,1,1)
    plot(xt_save(:,istep),'k-','Linewidth',2);
    ylabel('State','Fontsize',14);
   
    
    ylim([-2 2]);
    xlim([0 361]);
    set(gca,'Xtick',[0:60:360]);
    set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
    set(gca, 'Fontsize',14);
    set(gca,'TickLength',[0.05 0.05]);
    grid on;
    title(['True State for Time # ',num2str(istep)]);
    
    subplot(2,1,2)
    plot(xf_save(:,istep),'Color',[0 0.5 0],'Linewidth',2);
    ylabel('State','Fontsize',14);
     xlabel('Position','Fontsize',14);
    ylim([-2 2]);
    xlim([0 361]);
    set(gca,'Xtick',[0:60:360]);
    set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
    set(gca, 'Fontsize',14);
    set(gca,'TickLength',[0.05 0.05]);
    grid on;
    title(['Forecast State for Time # ',num2str(istep)]);
    
    if (istep==1)
        pause;
    end
    Movme(istep)=getframe(gcf);
end

% plot results for different times
figure(11)
clf(11)

% for time=1
time=1;
subplot(4,1,1);
plot(xt_save(:,time),'k-', 'Linewidth',2);
hold on;
plot(xf_save(:,time),'Color',[0 0.5 0], 'Linewidth',2);
hold on;
%xlabel('x-domain','Fontsize',14);
ylabel('State','Fontsize',14);
ylim([-1.5 1.5]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;
title(['Time # ',num2str(time)]);

% for time=5
time=5;
subplot(4,1,2);
plot(xt_save(:,time),'k-', 'Linewidth',2);
hold on;
plot(xf_save(:,time),'Color',[0 0.5 0], 'Linewidth',2);
hold on;
%xlabel('x-axis','Fontsize',14);
ylabel('State','Fontsize',14);
ylim([-2 2]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;
title(['Time # ',num2str(time)]);

% for time=100
time=100;
subplot(4,1,3);
plot(xt_save(:,time),'k-', 'Linewidth',2);
hold on;
plot(xf_save(:,time),'Color',[0 0.5 0], 'Linewidth',2);
hold on;
ylabel('State','Fontsize',14);
xlim([0 361]);
ylim([-2 2]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;
title(['Time # ',num2str(time)]);


% for time=240
time=240;
subplot(4,1,4);
plot(xt_save(:,time),'k-', 'Linewidth',2);
hold on;
plot(xf_save(:,time),'Color',[0 0.5 0], 'Linewidth',2);
hold on;
xlabel('x-domain','Fontsize',14);
ylabel('State','Fontsize',14);
ylim([-2 2]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;
title(['Time # ',num2str(time)]);
legend('True','Forecast','Location','Best');




