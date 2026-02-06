% Spacecraft Guidance And Navigation
% Academic Year 2023/2024
% Assignment #1 exercise #1
% Mutti Marcello 220252 10698636

%% EX1
clear all; close all; clc;

% Set the default font size for axes labels and tick labels
set(0, 'DefaultAxesFontSize', 20);

% Set the default font size for legend
set(0, 'DefaultLegendFontSize', 20);

% Set the default linewidth
set(0, 'DefaultLineLineWidth', 3);

% Circular Restricted Three Body Problem Mass Ratio
mu=0.012150;

% Libration Points L1 & L2 are located along the x-axis
xL1=lagrange(mu,1,1);
xL2=lagrange(mu,2,0);

rL1=[xL1 0 0]';
rL2=[xL2 0 0]';

fprintf('Libration point L1 location [-]: [%.10f, %.2f, %.2f] \n', rL1);

%% EX2
close all; clc;

% Initialization of the periodic halo orbit individuation process
x0=[1.08892819445324 0 0.0591799623455459 0 0.257888699435051 0];   % Initial state guess
Phi0=reshape(eye(6),[1,36]);                                        % STM initial conditions
y0=[x0 Phi0];

% CR3BP symmetry property
S=diag([1 -1 1]);

% Modified options for ODE integrators are set
abs_tol=1e-9;   % Absolute tolerance
rel_tol=1e-12;  % Relative tolerance
options = odeset('reltol', rel_tol, 'abstol', abs_tol,'Events',@(t,y) plane_intersection(t,y));

% Numerical Half-Orbit Propagation
[tt, yy]=ode113(@(t,x) CR3BP(t,x,mu),[0 4],y0,options);

XX=yy(:,1:6);                           % Collections of integrated state
xx=[XX(:,1:3); XX(end-1:-1:1,1:3)*S];   % Orbit completed exploiting symmetry S

figure
hold on
grid on
view(30,45)
plot3(xx(:,1),xx(:,2),xx(:,3),'b')
plot3(xx(size(XX,1),1),xx(size(XX,1),2),xx(size(XX,1),3),'ob')
plot3(1-mu,0,0,'ko',xL1,0,0,'r*',xL2,0,0,'r*');
text(1-mu+0.01,0,0.01,'Moon','FontSize',15)
text(xL1+0.01,0,0,'L1','FontSize',15)
text(xL2+0.01,0,0,'L2','FontSize',15)
hold off
legend('Initial Guess Orbit Propagation','State Discontinuity','Location','best')
xlabel('x [km/LU]')
ylabel('y [km/LU]')
zlabel('z [km/LU]')

% An iterative process is set up to correct the initial state x0,
% according to the problem constraints

err=1;              % Error initialization
n_max=200;          % Maximum number of iterations           
Y=[y0(1) y0(5)]';   % Initial state components to be corrected, y0 previously defined
LU=384000;          % Earth-Moon Semimajor Axis [km]
TU=27.3*24*3600;    % Moon Revolution Period [s]
tol=TU/LU*1e-5;     % Tolerance

it=1;
while err>tol && it<=n_max
    y0(1)=Y(1);
    y0(5)=Y(2);
    
    % Numerical Half-Orbit Propagation
    [tt, yy]=ode113(@(t,x) CR3BP(t,x,mu),[0 4],y0,options);
    
    % Final STM at time of event te: STM(te), and reduced STM
    Phi_f=reshape(yy(end,7:end),[6 6])
    Phi_r=Phi_f([4 6],[1 5]);
    
    % Initial conditions update
    Y=Y-Phi_r\[yy(end,4) yy(end,6)]';

    % Stopping criterion set as final perpendicolarity constraint
    err=abs(yy(end,4))+abs(yy(end,6));

    it=it+1;
end

fprintf('Applied differential STM correction in %.0f iterations with final error %.8f [m/s]\n',it-1,1000*err*LU/TU)

% Corrected initial conditions
y00=yy(1,:);
fprintf('Corrected state initial conditions: [%.16f %f %.16f %f %.16f %f]\n',y00(1),y00(2),y00(3),y00(4),y00(5),y00(6))

XX=yy(:,1:6);                           % Collections of integrated state
xx=[XX(:,1:3); XX(end-1:-1:1,1:3)*S];   % Orbit completed exploiting symmetry S

% Number of represented Halo Orbits
n_h=10;                             
% Collection of propagated positions
halos=cell(n_h,1);
halos{1}=xx;

figure
hold on
grid on
view(30,45)
plot3(xx(:,1),xx(:,2),xx(:,3),'b')
plot3(1-mu,0,0,'ko',xL1,0,0,'r*',xL2,0,0,'r*');
text(1-mu+0.01,0,0.01,'Moon','FontSize',15)
text(xL1+0.01,0,0,'L1','FontSize',15)
text(xL2+0.01,0,0,'L2','FontSize',15)
hold off
legend('Periodic Halo Orbit','Location','best')
xlabel('x [km/LU]')
ylabel('y [km/LU]')
zlabel('z [km/LU]')
% 
% % Collision
% d=zeros(length(xx),1);
% rm=cspice_bodvrd('Moon','RADII',3);
% for i=1:length(xx)
%     d(i)=norm(xx(i,:)-[1-mu,0,0])-mean(rm)/LU;
% end
% figure
% plot(d)

%% EX3
close all; clc;

% Initial Condition updated as last iteration initial condition for
% numerical continuation
y0=y00;

% Again an iterative process is set up to correct the initial states of a
% family of evenly spaced halo orbits along z
z0new=y00(3);
z0_h=linspace(z0new,0.034,n_h);

% Collection of corrected initial conditions
Y0=zeros(n_h,6);
Y0(1,:)=y00(1:6);

figure
hold on
grid on
view(30,45)
plot3(xx(:,1),xx(:,2),xx(:,3),'b')
plot3(1-mu,0,0,'ko',xL1,0,0,'r*',xL2,0,0,'r*');
text(1-mu+0.01,0,0.01,'Moon','FontSize',15)
text(xL1+0.01,0,0,'L1','FontSize',15)
text(xL2+0.01,0,0,'L2','FontSize',15)

for k=2:n_h
    y0(3)=z0_h(k);      % Initial z coordinate constraint
    
    err=1;              % Error initialization
    n_max=200;          % Maximum number of iterations           
    Y=[y0(1) y0(5)]';   % Initial state components to be corrected, y0 previously defined

    it=1;
    while err>tol && it<=n_max
        y0(1)=Y(1);
        y0(5)=Y(2);
    
        % Numerical Half-Orbit Propagation
        [tt, yy]=ode113(@(t,x) CR3BP(t,x,mu),[0 40],y0,options);
        
        % Final STM at time of event te: STM(te), and reduced STM
        Phi_f=reshape(yy(end,7:end),[6 6]);
        Phi_r=Phi_f([4 6],[1 5]);
        
        % Initial conditions update
        Y=Y-Phi_r\[yy(end,4) yy(end,6)]';

        % Stopping Criterion
        err=abs(yy(end,4))+abs(yy(end,6)+abs(yy(end,2)));
        it=it+1;
    end
    fprintf('Applied differential STM correction in %.0f iterations with final error %.8f [m/s]\n',it-1,1000*err*LU/TU)
    
    % Initial Condition updated as last iteration initial condition for
    % numerical continuation
    y0=yy(1,:);
    Y0(k,:)=y0(1:6);

    XX=yy(:,1:6);                           % Collections of integrated state
    xx=[XX(:,1:3); XX(end-1:-1:1,1:3)*S];   % Orbit completed exploiting symmetry S

    plot3(xx(:,1),xx(:,2),xx(:,3),'b')

    halos{k}=xx;    % Halo Orbits Family
end
hold off
legend('Halo Orbit Family','Location','best')
xlabel('x [km/LU]')
ylabel('y [km/LU]')
zlabel('z [km/LU]')

figure
hold on
grid on
view(30,45)
plot3(xL2,0,0,'r*');
text(xL2,-0.01,-0.01,'L2','FontSize',15)
for i=1:n_h
    plot3(halos{i}(:,1),halos{i}(:,2),halos{i}(:,3),'b')       % 3D mapping
    plot3(halos{i}(:,1),halos{i}(:,2),-0.1*ones(size(halos{i}(:,3))),'Color',[.7 .7 .7],'LineWidth',1)    % XY projection
    plot3(halos{i}(:,1),0.15*ones(size(halos{i}(:,2))),halos{i}(:,3),'Color',[.7 .7 .7],'LineWidth',1)    % XZ projection
    plot3(1.08*ones(size(halos{i}(:,1))),halos{i}(:,2),halos{i}(:,3),'Color',[.7 .7 .7],'LineWidth',1)    % YZ projection
end
hold off
ylim([-0.15 0.15])
zlim([-0.1 0.07])
xlabel('x [km/LU]')
ylabel('y [km/LU]')
zlabel('z [km/LU]')

%% FUNCTIONS

function [xLi] = lagrange(mu, i, d)
%     Evaluates x coordinate of Libration Points Li for i=1,2
%     Example: [xLi] = lagrange(mu, i, d)
%     INPUTS:
%         mu  [1x1] CR3BP mass ratio constant
%         i   [1x1] Libration point identifier: i=1 or i=2
%         d   [1x1] Solver display: d=1 display on, d=0 display off
%     OUTPUTS:
%         xLi [1x1] adimentional x coordinate of Libration Points Li
   
    % Potential derivative along x
    dUdx=@(rx) rx-(1-mu).*(rx+mu)./(abs(rx+mu).^3)-mu.*(rx+mu-1)./(abs(rx+mu-1).^3);

    % Li guess region definition
    if i==2
        x0=[1-mu 10]*(1+1e-6);
    else
        x0=[-mu 1-mu]*(1-1e-6);
    end
    
    % Solver Display Settings
    if d==1
        opt=optimset('Display','iter-detailed','TolX',1e-11);
    else
        opt=optimset('Display','off','TolX',1e-11);
    end

    % x coordinate found as zero of potential derivative along x
    [xLi,~,ex_flag]=fzero(dUdx,x0,opt);

    if ex_flag~=1   % failure safety
        error('ERROR lagrange')
    end
end

function [rhs] = CR3BP(~, x, mu)
%     Computes CR3BP equations of motion RHS together with relative STM
%     propagation RHS
%     Example: [rhs] = CR3BP(t, x, mu)
%     INPUTS:
%         t  [1x1]  reference epoch (omissible)
%         x  [1x42] state and state transition matrix
%         mu [1x1]  CR3BP mass ratio constant
%     OUTPUTS:
%         rhs [42x1] equations of motion and STM propagation RHS

    % State and STM verticality check
    [m,n]=size(x);
    if n>m
        x=x';
    end

    % State Variables Extraction
    rx=x(1);
    ry=x(2);
    rz=x(3);
    vx=x(4);
    vy=x(5);
    vz=x(6);
    
    % STM Extraction and Reshaping 
    Phi = reshape(x(7:end),6,6);

    % Two points relative distance
    r1 = sqrt((rx + mu)^2 + ry^2 + rz^2);
    r2 = sqrt((rx + mu - 1)^2 + ry^2 + rz^2);

    % Potential derivative along cartesian axis
    dUdx = rx - (1-mu)/r1^3*(mu+rx) + mu/r2^3*(1-mu-rx);
    dUdy = ry - (1-mu)/r1^3*ry - mu/r2^3*ry;
    dUdz = -(1-mu)/r1^3*rz-mu/r2^3*rz;

    % STM integration matrix A(t)=df/dx assembly
    dfdx =[0, 0, 0, 1, 0, 0;
           0, 0, 0, 0, 1, 0;
           0, 0, 0, 0, 0, 1;
           (mu - 1)/((mu + rx)^2 + ry^2 + rz^2)^(3/2) - mu/((mu + rx - 1)^2 + ry^2 + rz^2)^(3/2) + (3*mu*(2*mu + 2*rx - 2)*(mu + rx - 1))/(2*((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2)) - (3*(2*mu + 2*rx)*(mu + rx)*(mu - 1))/(2*((mu + rx)^2 + ry^2 + rz^2)^(5/2)) + 1, (3*mu*ry*(mu + rx - 1))/((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*ry*(mu + rx)*(mu - 1))/((mu + rx)^2 + ry^2 + rz^2)^(5/2), (3*mu*rz*(mu + rx - 1))/((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*rz*(mu + rx)*(mu - 1))/((mu + rx)^2 + ry^2 + rz^2)^(5/2), 0, 2, 0;
           (3*mu*ry*(2*mu + 2*rx - 2))/(2*((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2)) - (3*ry*(2*mu + 2*rx)*(mu - 1))/(2*((mu + rx)^2 + ry^2 + rz^2)^(5/2)), (mu - 1)/((mu + rx)^2 + ry^2 + rz^2)^(3/2) - mu/((mu + rx - 1)^2 + ry^2 + rz^2)^(3/2) - (3*ry^2*(mu - 1))/((mu + rx)^2 + ry^2 + rz^2)^(5/2) + (3*mu*ry^2)/((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2) + 1, (3*mu*ry*rz)/((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*ry*rz*(mu - 1))/((mu + rx)^2 + ry^2 + rz^2)^(5/2), -2, 0, 0;
           (3*mu*rz*(2*mu + 2*rx - 2))/(2*((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2)) - (3*rz*(2*mu + 2*rx)*(mu - 1))/(2*((mu + rx)^2 + ry^2 + rz^2)^(5/2)), (3*mu*ry*rz)/((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2) - (3*ry*rz*(mu - 1))/((mu + rx)^2 + ry^2 + rz^2)^(5/2), (mu - 1)/((mu + rx)^2 + ry^2 + rz^2)^(3/2) - mu/((mu + rx - 1)^2 + ry^2 + rz^2)^(3/2) - (3*rz^2*(mu - 1))/((mu + rx)^2 + ry^2 + rz^2)^(5/2) + (3*mu*rz^2)/((mu + rx - 1)^2 + ry^2 + rz^2)^(5/2), 0, 0, 0];
    
    % STM propagation RHS
    Phidot = dfdx*Phi;

    % RHS assembly
    rhs = zeros(42,1);
    rhs(1:3) = x(4:6);
    rhs(4) = 2*vy+dUdx;
    rhs(5) = -2*vx+dUdy;
    rhs(6) = dUdz;
    rhs(7:end) = Phidot(:);
    
end

function [value, isterminal, direction] = plane_intersection(~, y)
%     Verifies intersection of propagated orbit with y=0 plane
%     Example: [value,isterminal,direction] = plane_intersection(~, y)
%     INPUTS:
%         t   [1x1]  reference epoch (omissible)
%         y   [1x42] state and state transition matrix
%     OUTPUTS:
%         value      [1x1] value to be brought to zero
%         isterminal [1x1] if isterminal=1 the propagation is stopped
%                          if isterminal=0 the propagation continues
%         direction  [1x1] direction of approach of value=0
%                          if direction=-1 from above
%                          if direction=+1 from below
%                          if direction=0 either
                             
    % y(2) identifies the y coordinate of the state
    value = y(2);
    isterminal = 1;
    direction = -1;
end
