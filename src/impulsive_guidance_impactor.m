% Spacecraft Guidance And Navigation
% Academic Year 2023/2024
% Assignment #1 exercise #2
% Mutti Marcello 220252 10698636

%% EX1
clear all; close all; clc;

% Set the default font size for axes labels and tick labels
set(0, 'DefaultAxesFontSize', 20);

% Set the default font size for legend
set(0, 'DefaultLegendFontSize', 20);

% Set the default linewidth
set(0, 'DefaultLineLineWidth', 2);

% All necessary kernel are loaded
cspice_furnsh('kernels\pck00010.tpc');
cspice_furnsh('kernels\naif0012.tls');
cspice_furnsh('kernels\gm_de432.tpc');
cspice_furnsh('kernels\de432s.bsp');
cspice_furnsh('kernels\20099942_Apophis.bsp');
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

% Close encounter observability window
et0_ce=cspice_str2et('2029-Jan-1 00:00:00.0000 UTC');
etf_ce=cspice_str2et('2029-Jul-31 00:00:00.0000 UTC');

m=1e4;
ett_ce=linspace(et0_ce,etf_ce,m);
d_ap_sun=zeros(1,m);
d_ap_earth=zeros(1,m);
d_ap_moon=zeros(1,m);
beta=zeros(1,m);
lon_gt=zeros(1,m);
lat_gt=zeros(1,m);

for i=1:m
    % Aphophis-Sun relative position and distance
    rr_ap_sun=cspice_spkpos('20099942',ett_ce(i),'ECLIPJ2000','NONE','Sun');
    d_ap_sun(i)=norm(rr_ap_sun);

    % Aphophis-Earth relative position and distance
    rr_ap_earth=cspice_spkpos('20099942',ett_ce(i),'ECLIPJ2000','NONE','Earth');
    d_ap_earth(i)=norm(rr_ap_earth);

    % Earth-Aphophis-Sun Angle
    cb=dot(-rr_ap_sun,-rr_ap_earth)/(d_ap_sun(i)*d_ap_earth(i));
    sb=norm(cross(-rr_ap_sun,-rr_ap_earth))/(d_ap_sun(i)*d_ap_earth(i));
    beta(i)=rad2deg(atan2(sb,cb));
    
    % Aphophis-Moon relative distance
    d_ap_moon(i)=norm(cspice_spkpos('20099942',ett_ce(i),'ECLIPJ2000','NONE','Moon'));
end

% TCA computed as epoch of constrained minimum of relative Apophis-Earth
% distance on the close encounter window
eTCA_guess=ett_ce(d_ap_earth==min(d_ap_earth));

options=optimoptions(@fmincon,'Display','iter-detailed');
eTCA=fmincon(@(t) norm(cspice_spkpos('20099942',t,'ECLIPJ2000','NONE','Earth')),eTCA_guess,[],[],[],[],et0_ce,etf_ce,[],options);

% TCA and DCA prior SC impact
eTCA_UTC=cspice_et2utc(eTCA,'C',3);
fprintf('Time of Closest Approach (prior impact): %s \n',eTCA_UTC)
fprintf('Distance of Closest Approach (prior impact):  %.4f [km]\n\n',norm(cspice_spkpos('20099942',eTCA,'ECLIPJ2000','NONE','Earth')))

% Earth Flatness Factor Computation
earth_radii=cspice_bodvrd('Earth','RADII',3);
r_eq=earth_radii(1);
r_pol=earth_radii(3);
earth_f=(r_eq-r_pol)/r_eq;

% TCA centered 12h window ground track
ett_TCA=linspace(eTCA-6*3600,eTCA+6*3600,m+1);
for i=1:m
    rr_ap_earth_gt=cspice_spkpos('20099942',ett_TCA(i),'IAU_earth','NONE','Earth');
    [lo,la,~]=cspice_recgeo(rr_ap_earth_gt,r_eq,earth_f);
    lon_gt(i)=rad2deg(lo);
    lat_gt(i)=rad2deg(la);
end

figure
subplot(2,1,1)
plot(ett_ce,cspice_convrt(d_ap_sun,'KM','AU'))
grid on
xlabel('Reference epoch [s]')
ylabel('Distance [AU]')
legend('r_{Aph-Sun}','Location','best')
subplot(2,1,2)
plot(ett_ce,cspice_convrt(d_ap_moon,'KM','AU'),ett_ce,cspice_convrt(d_ap_earth,'KM','AU'))
grid on
xlabel('Reference epoch [s]')
ylabel('Distance [AU]')
legend('r_{Aph-Moon}','r_{Aph-Earth}','Location','best')

figure
plot(ett_ce,cspice_convrt(d_ap_moon-d_ap_earth,'KM','AU'))
grid on
xlabel('Reference epoch [s]')
ylabel('Distance [AU]')
legend('r_{Aph-Moon}-r_{Aph-Earth}','Location','best')

figure
subplot(1,2,1)
plot(ett_ce,beta)
grid on
title('Sun-Aphophis-Earth Angle')
xlabel('Ephemerides epoch [s]')
ylabel('\beta_{Sun-Aphophis-Earth} [deg]')
subplot(1,2,2)
plot(ett_ce(m/2-m/25:m/2+m/40),beta(m/2-m/25:m/2+m/40))
grid on
title('Zoom')
xlabel('Ephemerides epoch [s]')
ylabel('\beta_{Sun-Aphophis-Earth} [deg]')

figure
hold on
grid on
plot(lon_gt,lat_gt,'r.')
plot(lon_gt(1),lat_gt(1),'r')
plot(lon_gt(m/2+1),lat_gt(m/2+1),'ko')
hold off
legend('','Ground Track','Closest Approach')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
xlim([-180 180])
ylim([-90 90])

%% EX2
close all; clc;

% Minimization problem set up
% Variables Bounds
t0_LB=cspice_str2et('2024-10-01 00:00:00.0000 UTC');    % Launch
t0_UB=cspice_str2et('2025-02-01 00:00:00.0000 UTC');
t_dsm_LB=cspice_str2et('2025-04-01 00:00:00.0000 UTC'); % DSM
t_dsm_UB=cspice_str2et('2026-08-01 00:00:00.0000 UTC');
t_imp_LB=cspice_str2et('2028-08-01 00:00:00.0000 UTC'); % Impact
t_imp_UB=cspice_str2et('2029-02-28 00:00:00.0000 UTC');
UB=1e12*ones(21,1);
LB=-1e12*ones(21,1);
UB(19:21)=[t0_UB t_dsm_UB t_imp_UB]';
LB(19:21)=[t0_LB t_dsm_LB t_imp_LB]';

% Linear time constraints
A=zeros(2,21);
A(1,19)=1;
A(1,20)=-1;
A(2,20)=1;
A(2,21)=-1;
B=zeros(2,1);

% Modified options for ODE integrators are set
abs_tol=1e-9;   % Absolute Tolerance
rel_tol=1e-12;  % Relative Tolerance
options = odeset('reltol', rel_tol, 'abstol', abs_tol);

% % Initial Guess Generation Process
% % For the purpose of the excercise y0 is fixed 
% y0=zeros(21,1);
% 
% y0(19)=t0_LB+(t0_UB-t0_LB)*rand(1);
% y0(20)=t_dsm_LB+(t_dsm_UB-t_dsm_LB)*rand(1);
% y0(21)=t_imp_LB+(t_imp_UB-t_imp_LB)*rand(1);
% y0(1:6)=cspice_spkezr('Earth',y0(19),'ECLIPJ2000','NONE','Sun'); 
% 
% [~,y0dsm]=ode113(@(t,x) SC_rhs(t,x),[y0(19) y0(20)],y0(1:6),options);
% y0(7:12)=y0dsm(end,:)';   
% 
% [~,y0imp]=ode113(@(t,x) SC_rhs(t,x),[y0(20) y0(21)],y0(7:12),options);
% y0(13:18)=y0imp(end,:)';

% Initial feasible guess
y0=[49.2169459712239e+006
    139.029206142530e+006
   -7.70738675524294e+003
   -24.9356712816935e+000
    12.7828246193314e+000
    1.73489673747209e+000
   -90.4387624985575e+006
    123.927291118605e+006
    7.90561119119632e+006
   -20.5817574004446e+000
   -17.0864426689424e+000
    857.124084480135e-003
    59.3329568400330e+006
    133.257338829496e+006
   -723.718029360025e+003
   -23.9915964332561e+000
    15.1494536013268e+000
    1.72850446845409e+000
    786.416823182466e+006
    817.879910298779e+006
    915.607432187826e+006];

%% EX3
close all; clc;

% Maximum Initial Constraint Violation
[C0,Ceq0]=non_lcon(y0);
c0=norm([C0,Ceq0'],inf);

% Multiple Shooting Minimization Problem Solution
ConTol=1e-10;
fmin_opt=optimoptions('fmincon','Display','iter-detailed','ConstraintTolerance',ConTol);
Y=fmincon(@(y) avoidance(y),y0,[],[],[],[],LB,UB,@(y) non_lcon(y),fmin_opt);  % omitted linear time constraints
% Y=fmincon(@(y) avoidance(y),y0,A,B,[],[],LB,UB,@(y) non_lcon(y),fmin_opt);  % with linear time constraints

% Final Constraint Violation
[Cf,Ceqf]=non_lcon(Y);

% Earth State at Launch
x_ie=cspice_spkezr('Earth',Y(19),'ECLIPJ2000','NONE','Sun');

% Flow from x1 to x2
[~,xx12]=ode113(@(t,x) SC_rhs(t,x),[Y(19) Y(20)],Y(1:6),options);

% Flow from x2 to x3
[~,xx23]=ode113(@(t,x) SC_rhs(t,x),[Y(20) Y(21)],Y(7:12),options);

% Launch
fprintf('Time of Launch: %s UTC\n',cspice_et2utc(Y(19),'C',3))
fprintf('Launch DV: %.4f, [%.4f %.4f %.4f] [km/s]\n\n', norm(x_ie(4:6)-Y(4:6)),x_ie(4)-Y(4),x_ie(5)-Y(5),x_ie(6)-Y(6))

% Deep Space Maneuver
fprintf('Time of DSM: %s UTC\n',cspice_et2utc(Y(20),'C',3))
fprintf('DSM DV: %.4f, [%.4f %.4f %.4f] [km/s]\n\n', norm(xx12(end,4:6)'-Y(10:12)),xx12(end,4)-Y(10),xx12(end,5)-Y(11),xx12(end,6)-Y(12))

% Impact
fprintf('Time of Impact: %s UTC\n\n',cspice_et2utc(Y(21),'C',3))

% DCA and TCA post SC impact
[DCA,TCA] = avoidance(Y);
fprintf('Time of Closest Approach: %s \n',cspice_et2utc(TCA,'C',3))
fprintf('Distance of Closest Approach: %.4f [km], %.4f Re\n',-DCA,-DCA/mean(earth_radii))

% Transfer Trajectory Plot
Earth_eph=cspice_spkpos('Earth',linspace(t0_LB,TCA+2.5e6,m/10),'ECLIPJ2000','NONE','Sun');
Earth_TCA=cspice_spkpos('Earth',TCA,'ECLIPJ2000','NONE','Sun');
Aph_eph=cspice_spkpos('20099942',linspace(t_imp_LB,Y(21),m/10),'ECLIPJ2000','NONE','Sun');
% Initial post impact Aphophis state
v_imp=Y(16:18);
x0_ap=cspice_spkezr('20099942',Y(21),'ECLIPJ2000','NONE','Sun');
x0_ap(4:6)=x0_ap(4:6)+5*1e-5*v_imp;
% Aphophis post impact state propagation
[~, Aph_nb]=ode113(@(t,x) nbody_rhs(t,x),[Y(21) TCA+2.5e6],x0_ap,options);

figure
hold on 
grid on
view(30,45)
plot3(Earth_eph(1,:),Earth_eph(2,:),Earth_eph(3,:),'k')
plot3(xx12(:,1),xx12(:,2),xx12(:,3),'r')
plot3(xx23(:,1),xx23(:,2),xx23(:,3),'r')
plot3(Aph_eph(1,:),Aph_eph(2,:),Aph_eph(3,:),'b')
plot3(Aph_nb(:,1),Aph_nb(:,2),Aph_nb(:,3),'b--')

plot3(xx12(1,1),xx12(1,2),xx12(1,3),'dr')
plot3(xx23(1,1),xx23(1,2),xx23(1,3),'sr')
plot3(xx23(end,1),xx23(end,2),xx23(end,3),'*r')
plot3(Earth_TCA(1),Earth_TCA(2),Earth_TCA(3),'ko')
title('Planetary Protection (@Sun ECLIPJ2000)')
legend('Earth','SC','','Aphophis^{-}','Aphophis^{+}','Launch','DSM','Impact','Closest Approach','Location','best')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')

%% FUNCTIONS

function [d_TCA, t_TCA_n] = avoidance(y)
%     Objective function of the multiple shooting minimization problem
%     Example: [d_TCA, t_TCA_n] = avoidance(y)
%     INPUTS:
%         y [21x1] problem variables y=[x1 x2 x3 t1 t2 t3]'
%     OUTPUTS:
%         d_TCA   [1x1] objective function, DCA
%         t_TCA_n [1x1] TCA

    % Variables verticality check
    [m,n]=size(y);
    if n>m
        y=y';
    end
    
    % SC impact velocity
    v_imp=y(16:18);
    
    % Aphophis State prior impact
    x0_ap=cspice_spkezr('20099942',y(21),'ECLIPJ2000','NONE','Sun');

    % Aphophis State post impact
    x0_ap(4:6)=x0_ap(4:6)+5*1e-5*v_imp;

    % Modified options for ODE integrator
    options = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(t,x) close_encounter(t,x));

    % Aphophis post impact state propagation
    [~, ~, t_TCA_n, xx_TCA]=ode113(@(t,x) nbody_rhs(t,x),[y(21) 10*y(21)],x0_ap,options);

    % Earth position at TCA
    r_earth=cspice_spkpos('Earth',t_TCA_n,'ECLIPJ2000','NONE','Sun');

    % DCA
    d_TCA=-norm(r_earth-xx_TCA(1:3)');
end

function [C, Ceq] = non_lcon(y)
%     Non-linear constraints of the bjective function
%     Example: [C, Ceq]= non_lcon(y)
%     INPUTS:
%         y [21x1] problem variables y=[x1 x2 x3 t1 t2 t3]'
%     OUTPUTS:
%         C   [1x1]  Non-linear inequality constraint
%         Ceq [15x1] Non-linear equality constraint

    % Variables verticality check    
    [m,n]=size(y);
    if n>m
        y=y';
    end

    % Nodes Extraction
    x1=y(1:6);
    x2=y(7:12);
    x3=y(13:18);
    
    % Earth State at Launch
    x_ie=cspice_spkezr('Earth',y(19),'ECLIPJ2000','NONE','Sun');

    % Aphophis Position at Impact
    r_faph=cspice_spkpos('20099942',y(21),'ECLIPJ2000','NONE','Sun');

    % Modified options for ODE integrator
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [~, xx12]=ode113(@(t,x) SC_rhs(t,x),[y(19) y(20)],x1,options); % flow from x1 to x2
    [~, xx23]=ode113(@(t,x) SC_rhs(t,x),[y(20) y(21)],x2,options); % flow from x2 to x3

    C=norm(x_ie(4:6)-x1(4:6))+norm(xx12(end,4:6)'-x2(4:6))-5;  % dv limit    [

    Ceq=[(x_ie(1:3)-x1(1:3));         % initial node position continuity            [km. km/s]
         (xx12(end,1:3)'-x2(1:3));    % middle node position flow continuity
         (xx23(end,:)'-x3);           % final node state flow continuity
         (r_faph-x3(1:3))];           % final node position continuity
end

function [value, isterminal, direction] = close_encounter(t, x)
%     Verifies minimun of Earth-Aphophis distance
%     Example: [value,isterminal,direction] = close_encounter(t, x)
%     INPUTS:
%         t   [1x1] reference epoch
%         x   [6x1] Aphophis state
%     OUTPUTS:
%         value      [1x1] value to be brought to zero
%         isterminal [1x1] if isterminal=1 the propagation is stopped
%                          if isterminal=0 the propagation continues
%         direction  [1x1] direction of approach of value=0
%                          if direction=-1 from above
%                          if direction=+1 from below
%                          if direction=0 either

    % Earth state at reference epoch t
    x_earth=cspice_spkezr('Earth',t,'ECLIPJ2000','NONE','Sun');

    % State Defect
    dx=x-x_earth;
    rho=dx(1:3);        % Position defect
    rho_dot=dx(4:6);    % Velocity defect
    
    value=2*dot(rho,rho_dot);
    isterminal=1;
    direction=1;
end

function [rhs] = SC_rhs(~, x)
%     Propagates SC state only considering the attraction of the Sun
%     Example: [rhs] = SC_rhs(t, x)
%     INPUTS:
%         t   [1x1] reference epoch
%         x   [6x1] SC state
%     OUTPUTS:
%         rhs [6x1] equations of motion RHS

    % State verticality check
    [m,n]=size(x);
    if n>m
        x=x';
    end

    % Sun Gravitational Parameter
    GM=cspice_bodvrd('Sun','GM',1);

    % RHS initialization and assembly
    rr=x(1:3);
    vv=x(4:6);
    rhs=zeros(6,1);
    rhs(1:3)=vv;
    rhs(4:6)=-GM*rr/(norm(rr)^3);
end

function [rhs] = nbody_rhs(t, x)
%     Propagates Aphophis state considering Sun centerednon-inertial n-body
%     model
%     Example: [rhs] = nbody_rhs(t, x)
%     INPUTS:
%         t   [1x1] reference epoch
%         x   [6x1] Aphophis state
%     OUTPUTS:
%         rhs [6x1] equations of motion RHS

    % n-body problem initialization
    labels={'Sun'; 'Mercury'; 'Venus'; 'Earth'; 'Moon'; 'Mars Barycenter';
                  'Jupiter Barycenter'; 'Saturn Barycenter'; 'Uranus Barycenter';
                  'Neptune Barycenter'; 'Pluto Barycenter'};
    bodies=cell(size(labels));
    for i=1:length(labels)
        bodies{i}.name=labels{i};
        bodies{i}.GM=cspice_bodvrd(labels{i},'GM',1);
    end
    
    % RHS initialization
    rhs=zeros(6,1);
    rhs(1:3)=x(4:6);
    
    % Aphophis position with respect to the Sun and norm
    rr=x(1:3);
    r=sqrt(dot(rr,rr));

    % Gravitational constant of the Sun
    GM0=bodies{1}.GM;

    % Contribution of central body
    rhs(4:6)=-GM0*rr/r^3;

    % The contributions of remaining bodies are iteratively added
    for i=1:length(bodies)
        if i==1
            continue;
        end

        % i-th body position with respect to the Sun
        rr_s = cspice_spkpos(bodies{i}.name, t, 'ECLIPJ2000', 'NONE', 'Sun');

        % Aphophis position with respect to i-th body and norm
        dd=rr-rr_s;
        d=sqrt(dot(dd,dd));

        % Non-inertial coefficients computation
        q=dot(rr,rr-2*rr_s)/dot(rr_s,rr_s);
        f=q*(3+3*q+q^2)/(1+(1+q)^(3/2));

        % Non-inertial acceleration correction
        aa=-bodies{i}.GM/d^3*(rr+rr_s*f);
        rhs(4:6)=rhs(4:6)+aa;
    end
end
