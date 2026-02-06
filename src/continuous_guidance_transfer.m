% Spacecraft Guidance And Navigation
% Academic Year 2023/2024
% Assignment #1 exercise #3
% Mutti Marcello 220252 10698636

%% EX2
clear all; close all; clc;

% Set the default font size for axes labels and tick labels
set(0, 'DefaultAxesFontSize', 20);

% Set the default font size for legend
set(0, 'DefaultLegendFontSize', 20);

% Set the default linewidth
set(0, 'DefaultLineLineWidth', 2);

% All necessary kernel are loaded
cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\gm_de432.tpc');
cspice_furnsh('kernels\de432s.bsp');
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

% Modified options for ODE integrators are set
abs_tol=1e-9;   % Absolute tolerance
rel_tol=1e-12;  % Relative tolerance

% Launch epoch
t0=cspice_str2et('2023-May-28 14:13:09.0000 UTC');

% Spacecraft Properties
m0=1000;        % Mass [kg]
Tmax=800e-6;    % Thrust [kg*km/s^2]
Isp=3120;       % Specific Impulse [s]

% Earth Standard Gravity
g0=9.80665e-3;  % [km/s^2]

% Adimensional Units
LU=cspice_convrt(1,'AU','KM');              % Length Unit [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % Time Unit [s]
MU=1000;                                    % Mass Unit [kg]
AD=[LU, MU, TU];

% Adimensionalized Units
TIG=[Tmax*TU^2/(MU*LU) Isp/TU g0*TU^2/LU];
m0=m0/MU;
fprintf('Adimensionalized Tmax %.4e, Isp %.4e, g0 %.4e\n',TIG(1),TIG(2),TIG(3))

% Earth state at launch
xxe0=cspice_spkezr('Earth',t0,'ECLIPJ2000','NONE','Sun');

% Adimensional position and velocity
rre0=xxe0(1:3)/LU;
vve0=xxe0(4:6)/(LU/TU);
fprintf('Adimensionalized r0: [%.4e %.4e %.4e]\n',rre0(1),rre0(2),rre0(3))
fprintf('Adimensionalized v0: [%.4e %.4e %.4e]\n',vve0(1),vve0(2),vve0(3))

%% EX3
close all; clc;

% Zero-finding problem set up
ex_flag=-1;
TIG(1)=Tmax*TU^2/(MU*LU);

% State initial conditions
x0=[rre0; vve0; m0];

% Modified Options for fsolve
ConTol=1e-4;
fsopt=optimoptions('fsolve','Display','iter-detailed','OptimalityTolerance',ConTol);

% Modified ODE options
odeopt=odeset('RelTol',rel_tol,'AbsTol',abs_tol);

% Correction attempted iteratively for different initial conditions
while ex_flag<=0 || U800(8)<0

%     % Random Costate initial conditions guess [-20 20]
%     ll0g=-20+40*rand(7,1);
%     ll0g(7)=abs(ll0g(7));
% 
%     % Random adimensional final time guess [0 2pi]
%     tfg=2*pi*rand(1,1);
%     
%     % Variables initial guess
%     u0=[ll0g; tfg];

    u0=[-9.26244714410867e+000
        -9.68615319549581e+000
        -6.73339045029483e+000
        -13.9106394854821e+000
        -6.07969361135546e+000
        -15.1336618276910e+000
         15.3661223099864e+000
         592.368595235588e-003];
    
    % Optimal Solution
    [U800,~,ex_flag]=fsolve(@(u) timeopt_conv(u,x0,t0,AD,TIG),u0,fsopt);
    
end

% Optimal solution
fprintf('[Tmax=800mN]\n')
fprintf('Costate initial conditions: [%.4e %.4e %.4e %.4e %.4e %.4e %.4e]\n',U800(1:7))
fprintf('Time of Flight: %.4f [days]\n',U800(8)*AD(3)/(24*3600))
fprintf('Date of Arrival: %s UTC\n',cspice_et2utc(U800(8)*AD(3)+t0,'C',3))

% Final Defects given optimal initial conditions
[df800]=timeopt_conv(U800,x0,t0,AD,TIG);
fprintf('Final Position Defect:  %.4e [km]\n',norm(df800(1:3))*AD(1))
fprintf('Final Velocity Defect:  %.4e [m/s]\n',1e3*norm(df800(4:6))*AD(1)/AD(3))

% Transfer Trajectory Plot
[tt,yy]=ode113(@(t,y) LE(t,y,TIG),[0 U800(8)],[x0; U800(1:7)],odeopt);
ttp=linspace(0,tt(end),100);
rr_e=cspice_spkpos('Earth',ttp*AD(3)+t0,'ECLIPJ2000','NONE','Sun'); % Earth position
rr_v=cspice_spkpos('Venus',ttp*AD(3)+t0,'ECLIPJ2000','NONE','Sun'); % Venus position

figure
hold on
grid on
view(45,20)
plot3(yy(:,1)*AD(1),yy(:,2)*AD(1),yy(:,3)*AD(1),'b')
plot3(rr_e(1,:),rr_e(2,:),rr_e(3,:),'g')
plot3(rr_v(1,:),rr_v(2,:),rr_v(3,:),'r')
plot3(yy(1,1)*AD(1),yy(1,2)*AD(1),yy(1,3)*AD(1),'sb')
plot3(yy(end,1)*AD(1),yy(end,2)*AD(1),yy(end,3)*AD(1),'db')
legend('SC','Earth','Venus')
title('Earth-Venus Transfer (@Sun ECLIPJ2000)')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% Optimal Hamiltonian, Switching function
H800=zeros(size(tt));
St800=zeros(size(tt));
for i=1:length(tt)
    ff=LE(tt(i),yy(i,:),TIG);
    H800(i)=1+dot(yy(i,8:14),ff(1:7));
    St800(i)=-(norm(yy(i,11:13))*TIG(2)*TIG(3))/(yy(i,7))-yy(i,14);
end
fprintf('Hamiltonian maximum relative variation: %.4e%% \n',100*abs((max(H800)-min(H800))/min(H800)))

% Hamiltonian plot
figure
plot(tt,H800)
grid on
xlabel('t [-]')
ylabel('H(t) [-]')

% Switching function plot
figure
plot(tt,St800)
grid on
xlabel('t [-]')
ylabel('S_t(t) [-]')

%% EX4
close all; clc;

% Initial Condition updated as last iteration initial condition for
% numerical continuation
u0=U800;

% New Spacecraft Thrust
Tmaxn=500-6;    % [kg*km/s^2]

% Family of Solutions
n_cont=13;
Tmax_v=linspace(Tmax,Tmaxn,n_cont);

figure
hold on
grid on
view(45,20)
for i=1:n_cont
    
    % Updated Thrust parameter
    TIG(1)=Tmax_v(i)*TU^2/(MU*LU);

    ex_flag=0;
    while ex_flag<=0 || U500(8)<0
        
        % Optimal solution
        [U500,~,ex_flag]=fsolve(@(u) timeopt_conv(u,x0,t0,AD,TIG),u0,fsopt);
        [tt,yy]=ode113(@(t,y) LE(t,y,TIG),[0 U500(8)],[x0; U500(1:7)],odeopt);

    end

    % Initial Condition updated as last iteration initial condition for
    % numerical continuation
    u0=U500;
    ttp=linspace(0,tt(end),100);
    
    % Transfer Trajectory Plot
    if i==1 || i==n_cont
        if i==1
            pt='b--';
        else
            pt='b';
        end
        plot3(yy(1,1)*AD(1),yy(1,2)*AD(1),yy(1,3)*AD(1),'sb')
        plot3(yy(end,1)*AD(1),yy(end,2)*AD(1),yy(end,3)*AD(1),'db')
        plot3(yy(:,1)*AD(1),yy(:,2)*AD(1),yy(:,3)*AD(1),pt)
    end
    

end

rr_e=cspice_spkpos('Earth',ttp*AD(3)+t0,'ECLIPJ2000','NONE','Sun'); % Earth position
rr_v=cspice_spkpos('Venus',ttp*AD(3)+t0,'ECLIPJ2000','NONE','Sun'); % Venus position
plot3(rr_e(1,:),rr_e(2,:),rr_e(3,:),'g')
plot3(rr_v(1,:),rr_v(2,:),rr_v(3,:),'r')
legend('','','SC (800mN)','','','SC (500mN)','Earth','Venus')
title('Earth-Venus Transfer (@Sun ECLIPJ2000)')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

% Optimal solution
fprintf('[Tmax=500mN]\n')
fprintf('Costate initial conditions: [%.4e %.4e %.4e %.4e %.4e %.4e %.4e]\n',U500(1:7))
fprintf('Time of Flight: %.4f [days]\n',U500(8)*AD(3)/(24*3600))
fprintf('Date of Arrival: %s UTC\n',cspice_et2utc(U500(8)*AD(3)+t0,'C',3))

% Final Defects given optimal initial conditions
[df500]=timeopt_conv(U500,x0,t0,AD,TIG);
fprintf('Final Position Defect:  %.4e [km]\n',norm(df500(1:3))*AD(1))
fprintf('Final Velocity Defect:  %.4e [m/s]\n',1e3*norm(df500(4:6))*AD(1)/AD(3))

% Optimal Hamiltonian, Switching function
H500=zeros(size(tt));
St500=zeros(size(tt));
for i=1:length(tt)
    ff=LE(tt(i),yy(i,:),TIG);
    H500(i)=1+dot(yy(i,8:14),ff(1:7));
    St500(i)=-(norm(yy(i,11:13))*TIG(2)*TIG(3))/(yy(i,7))-yy(i,14);
end
fprintf('Hamiltonian maximum relative variation: %.4e%% \n',100*abs((max(H500)-min(H500))/min(H500)))

% Hamiltonian plot
figure
plot(tt,H500)
grid on
xlabel('t [-]')
ylabel('H(t) [-]')

% Switching function plot
figure
plot(tt,St500)
grid on
xlabel('t [-]')
ylabel('S_t(t) [-]')

%% FUNCTIONS

function [rhs] = LE(~, y, TIG)
%     Two Body Proble adimentional Lagrange-Euler equations RHS for Optimal
%     Control Problem
%     Example: [rhs] = LE(~, y, TIG)
%     INPUTS:
%         t   [1x1]  reference epoch (omissible)
%         y   [14x1] state and costate
%         TIG [1x3]  adimensional parameters TIG=[Tmax Isp g0]
%     OUTPUTS:
%         rhs [14x1] state and costate dynamics RHS
        
    % Verticality check
    [m,n]=size(y);
    if n>m
        y=y';
    end
    
    % Variables extraction
    rr=y(1:3);
    vv=y(4:6);
    m=y(7);
    lrr=y(8:10);
    lvv=y(11:13);
    r=norm(rr);
    lv=norm(lvv);

    % RHS initialization and assembly
    rhs=zeros(14,1);
    rhs(1:3)=vv;
    rhs(4:6)=-rr/r^3-TIG(1)/m*lvv/lv;
    rhs(7)=-TIG(1)/(TIG(2)*TIG(3));
    rhs(8:10)=-3*dot(rr,lvv)*rr/r^5+lvv/r^3;
    rhs(11:13)=-lrr;
    rhs(14)=-lv*TIG(1)/m^2;

end

function [df] = timeopt_conv(u, x0, t0, AD, TIG)
%     Computes the adimensional constraints residuals of the time optimal 
%     control problem
%     Example: [df] = timeopt_conv(u, x0, t0, AD, TIG)
%     INPUTS:
%         u   [8x1] problem adimensional variables u=[costate, tf]
%         x0  [6x1] state initial adimesional conditions
%         t0  [1x1] initial dimensional reference epoch
%         AD  [3x1] adimensional units AD=[LU MU TU]
%         TIG [1x3] adimensional parameters TIG=[Tmax Isp g0]
%     OUTPUTS:
%         df  [8x1] adimesional constraints redisuals
    
    % Variables verticality check
    [m,n]=size(u);
    if n>m
        u=u';
    end
    % Initial state verticality check
    [m,n]=size(x0);
    if n>m
        x0=x0';
    end
    
    % Variables extraction
    ll0=u(1:7);         % costate
    tf=u(8);            % final time
    tfd=tf*AD(3)+t0;    % dimensional final time

    % State and costate initial conditions
    y0=[x0; ll0];
    
    % Modified ODE options
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);

    % State and costate propagation
    [~,yy]=ode113(@(t,y) LE(t,y,TIG),[0 tf],y0,options);
    
    % Hamiltonian computation at tf
    ff=LE(tf,yy(end,:),TIG);
    ff=ff(1:7);
    ll=yy(end,8:14);
    Hf=1+dot(ll,ff);
    
    % Venus state extraction at tf
    xxv=cspice_spkezr('Venus',tfd,'ECLIPJ2000','NONE','Sun');
    rrv=xxv(1:3)/AD(1);         % adimensional position
    vvv=xxv(4:6)/(AD(1)/AD(3)); % adimensional velocity
    aav=-rrv/norm(rrv)^3;       % adimensional instantaneous acceleration

    % State final constraint derivative
    psi_dot=[vvv; aav; 0];

    % Constraints residuals initialization and assembly
    df=zeros(8,1);
    df(1:3)=yy(end,1:3)'-rrv;   % position defect
    df(4:6)=yy(end,4:6)'-vvv;   % velocity defect
    df(7)=yy(end,14);           % lambda_m final constraint
    df(8)=Hf-dot(ll,psi_dot);   % transversality condition
end