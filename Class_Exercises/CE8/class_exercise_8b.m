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
Rt=gcorr('gauss',Lb-La,Ld,N,0);   % construct true error correlation
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
% add model error (play around with model error)
add_model_error=0;                  
if add_model_error == 0
    Q=0;                          % no model error
else
    Q=varQ*Rt;                    % with model error
end

% set up observations
nobs=20;                          % number of obs in the spatial domain
dobs=5;                           % get obs every dobs time step
var_y=0.01;                       % obs error variance
y_count=0;                        % initial count for obs
y_save=zeros(nobs, nStep/dobs);   % obs taken every 5 time step 
R=var_y*eye(nobs);                % obs error covariance

% do analysis?
analysis=1;

% initialize matrices to be saved
xt_save=zeros(N,nStep);           % true state
xf_save=zeros(N,nStep);           % forecast state
Pf_save=zeros(N,N,nStep);         % forecast error covariance
xa_save=zeros(N,nStep);           % analysis state
Pa_save=zeros(N,N,nStep);         % analysis error covariance

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
%==========================================================================    
% set initial analysis to be the forecast states
xa=xf;
Pa=Pf;

%==========================================================================
% load files instead
%load xo_f.dat;
%xf=xo_f;

%load xo_t.dat;
%xt=xo_t;

%load Q.dat;
%load M.dat;


% construct observation operator
H=zeros(nobs,N);
col=[4:N/nobs:N];
for i=1:nobs
    H(i,col(i))=1;
end
%==========================================================================
% construct identity matrix (NxN)
I=eye(N);

% d=0;
% MAIN LOOP (time loop)
for istep=1:nStep
    disp(istep);

    % true state time update
    xt=mvnrnd(M*xt,Qt)';
    xt_save(:,istep)=xt;
    
    % forecast state time update
    xf=M*xa;
    xf_save(:,istep)=xf;
    
    % forecast error covariance time update 
    Pf=M*Pa*M' + Q;
    Pf_save(:,:,istep)=Pf;
%==========================================================================    
    % is there obs?
    if (mod(istep,dobs)==0)
        
        % get obs
        yt=H*xt;
        % add noise
        y=yt+sqrt(var_y)*randn(nobs,1);
        % obs counter
        y_count=y_count+1;
        y_save(:,y_count)=y;
%========================================================================== 
        if analysis==1
            % compute Kalman Gain
            K=Pf*H'*inv(H*Pf*H'+R);
            % analysis state update
            xa=xf+K*(y-H*xf);
            % analysis error covariance update
            Pa=(I-K*H)*Pf*(I-K*H)' + K*R*K';
        else
            xa=xf;
            Pa=Pf;
        end
%==========================================================================
    else
        xa=xf;
        Pa=Pf;
    end
    xa_save(:,istep)=xa;
    Pa_save(:,:,istep)=Pa;
    Var_save(:,istep)=diag(Pa);
    
    
        % alright, let's make a movie out of this
    figure(10)
    clf(10)

    plot(xt_save(:,istep),'k-','Linewidth',2);
    hold on;
    
    if (analysis==1)
        plot(xf_save(:,istep),'Color',[0 0.5 0],'Linewidth',3);
        hold on;
        plot(xa_save(:,istep),'b','Linewidth',2);
        hold on;
    end
    if (mod(istep,dobs)==0)
        plot(col,y_save(:,y_count),'ro');
        hold on;
    end
    %plot(Var_save(:,istep).^0.5,'b:');
    xlabel('x-domain','Fontsize',14);
    ylabel('solution','Fontsize',14);
    ylim([-2 2]);
    xlim([0 361]);
    set(gca,'Xtick',[0:60:360]);
    set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
    set(gca, 'Fontsize',14);
    set(gca,'TickLength',[0.05 0.05]);
    grid on;
    title(['Iteration # ',num2str(istep)]);
    

    if (mod(istep,dobs)==0)
        if (analysis==1)
            legend('True','Forecast','Analysis','Obs','Location','NorthEastOutside');
        else
            legend('True','Obs','Location','NorthEastOutside');
        end
    else
        if (analysis==1)
            legend('True','Forecast','Analysis','Location','NorthEastOutside');
        else
            legend('True','Location','NorthEastOutside');
        end
    
    end
    
    % pause;
    %Movme(istep)=getframe(gcf);
    
    
    
    
end
%==========================================================================
% plot results
figure(1)
clf(1)

if (varQ==0)
    ymin=-1.5;
    ymax=1.5;
else
    ymin=-1.5;
    ymax=1.5;
end

% for time=5
subplot(3,1,1);
plot(xt_save(:,5),'k-','Linewidth',2);
hold on;
plot(xf_save(:,5),'Color',[0 0.5 0],'Linewidth',2);
hold on;
plot(xa_save(:,5),'b','Linewidth',2);
hold on;
if (nobs>0)
    plot(col,y_save(:,1),'ro');
    hold on;
end
%plot(Var_save(:,5).^0.5,'b:');
xlabel('x-domain','Fontsize',14);
ylabel('solution','Fontsize',14);
ylim([ymin ymax]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;
title(['Iteration # ',num2str(5)]);
legend('True','Forecast','Analysis','Obs','Location','NorthEastOutside');

% for time=100
subplot(3,1,2);
plot(xt_save(:,nStep/2),'k-','Linewidth',2);
hold on;
plot(xf_save(:,nStep/2),'Color',[0 0.5 0],'Linewidth',2);
hold on;
plot(xa_save(:,nStep/2),'b','Linewidth',2);
hold on;
if (nobs>0)
plot(col,y_save(:,nStep/dobs/2),'ro');
hold on;
end

xlabel('x-axis','Fontsize',14);
ylabel('solution','Fontsize',14);
ylim([ymin ymax]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;

title(['Iteration # ',num2str(nStep/2)]);
legend('True','Forecast','Analysis','Obs','Location','NorthEastOutside');

% for time=240
subplot(3,1,3);
plot(xt_save(:,nStep),'k-','Linewidth',2);
hold on;
plot(xf_save(:,nStep),'Color',[0 0.5 0],'Linewidth',2);
hold on;
plot(xa_save(:,nStep),'b','Linewidth',2);
hold on;
if (nobs>0)
    plot(col,y_save(:,nStep/dobs),'ro');
    hold on;
end
xlabel('x-axis','Fontsize',14);
ylabel('solution','Fontsize',14);
ylim([ymin ymax]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
legend('True','Forecast','Analysis','Obs','Location','NorthEastOutside');
set(gca, 'Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;

title(['Iteration # ',num2str(nStep)]);


figure(2)
clf(2)
subplot(3,1,1)
plot(Var_save(:,5).^0.5,'k-','Linewidth',2);
set(gca, 'Fontsize',14);
xlim([0 361]);
ylim([ymin ymax]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
title(['Analysis Std for Iteration # ',num2str(5)]);
xlabel('x-domain','Fontsize',14);
ylabel('solution','Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;


subplot(3,1,2)
plot(Var_save(:,nStep/2).^0.5,'b-','Linewidth',2);
set(gca, 'Fontsize',14);
title(['Analysis Std for  Iteration # ',num2str(nStep/2)]);
xlim([0 361]);
ylim([ymin ymax]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
ylim([ymin ymax]);
set(gca, 'Fontsize',14);
xlabel('x-domain','Fontsize',14);
ylabel('solution','Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;


subplot(3,1,3)
plot(Var_save(:,nStep).^0.5,'r-','Linewidth',2);
set(gca, 'Fontsize',14);
title(['Analysis Std for Iteration # ',num2str(nStep)]);
ylim([ymin ymax]);
xlim([0 361]);
set(gca,'Xtick',[0:60:360]);
set(gca,'XTickLabel',{'-180','-120',' -60','  0 ','  60','  120','  180'}); 
set(gca, 'Fontsize',14);
xlabel('x-domain','Fontsize',14);
ylabel('solution','Fontsize',14);
set(gca,'TickLength',[0.05 0.05]);
grid on;
