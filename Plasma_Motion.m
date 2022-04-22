function Plasma_Motion
clc

%GLOBALS 
global Re;
global m_e;
global m_p;
global q_e;
global q_p;
global M;
global E;

%CONSTANTS
Re=6378;%km Equatorial, polar is 6357
m_e=9.1093897*10^-31;%kg mass of electron 
m_p=1.6726231*10^-27;%kg mass of proton 
q_e=-1.602*10^-19;%C(As) Charge of electron
q_p=1.602*10^-19;%C(As) Charge of proton
%TEMPORARY CONSTAINTS
global m;
m0=1; %kg
global q;
q=1; %C

%VARIABLES
%Inner Van Allen Belt range from Earth's surface 0.2-2 Re
%Outer Van Allen Belt range from Earth's surface 3-10 Re
%(greatest intensity 4-5Re)
%Magnetic Field Strength at Surface of Earth 25-65 x10^-6 T
%Equatorial=30x10^-6 T, 60 at poles, 25 at South Anomoly 
%Assignment 21 of PH545 relates equatorial field strength to L shell(B0 eq)
%use L=2 (from center of earth) for inner or L=5.5 for outer
%B_in = 3.71*10^-6 B_out = 0.178*10^-6
%B=[0,0,3.71*10^-6];%T(kg/s^2A) Magnetic Field Intensity 
E=0; %V/m or N/C (mkg/s^3A) Electric Field, Assumed neglible
M=1000; %Tm^3 from paper 31000.*10^-9
c=3*10^8; %m/s

%INITIAL CONDITIONS 
%2D Magnetic Equitorial Plane 
r0=[7, 0, 0];%m 
for i=1 %i=1 
v0=[0, i*3 0];%m/s %v0=[0, 3, 0];%m/s 
v_abs(i,:)=sqrt(v0(1)^2+v0(2)^2+v0(3)^2); %used for ode45 guiding center 
Tmin=0;Tmax=53;%74 for e

%Energy MeV
m=m0/sqrt(1-(v0(2)^2/c^2)); %check if relativitc 
energy=((m-m0)*c^2)/(1.602E-19*1.0E6); %relatictic 

A=q*M/m;

r_st=r0(1)/(1-r0(1)^2*v0(2)/A);
v_st=A/r_st^2;
%v_st*0.204
%v_st/4
r0_mag=sqrt(sum(r0.^2,2));
v0_mag=sqrt(sum(v0.^2,2));
r_max=2*r_st/(1+sqrt(1-4*v0_mag/v_st));
r_min=2*r_st/(1+sqrt(1+4*v0_mag/v_st));


%Call ODE45 Functions
options=odeset('Events',@events,'RelTol',1.0e-10,'AbsTol',1.0e-14);
[t,y]=ode45(@Lorentz,[Tmin,Tmax],[r0,v0],options);
%to disable options change isterminal to 0


%RESULTS
%Particle Motion Lorentz 
r=[y(:,1), y(:,2), y(:,3)];%m 
v=[y(:,4), y(:,5), y(:,6)];%m/s
r_mag=sqrt(sum(r.^2,2));%m 
v_mag=sqrt(sum(v.^2,2));%m/s

%Call Magntic Field Function 
[Bx, By, Bz] = mag_vec(r);
B=[Bx, By, Bz];
B_mag=sqrt(sum(B.^2,2));%Bz is B_mag 

%Radius of gyration/Larmor/Cyclotron Equation
%r_g=m*v_perp/(|q|*B)
%use vectorized form to get correct direction
%r_gr=r+m*v_perp/(q*B)[B/B] 
%r_gr=r+m/(q*B^2)[v_perp x B] 
r_gr=r+(m./(q.*B_mag.^2)).*cross(v,B);
r_grx=r_gr(:,1);% don't need but keeping for mental vizualization
r_gry=r_gr(:,2);
r_gr_mag=sqrt(sum(r_gr.^2,2));

%Iterative gyrocenter approximation 
r_gc=r_gr;%feed intial condition 
for iterations=1:24 
[B_gcx, B_gcy, B_gcz] = mag_vec(r_gc); %extract magnetic field at new iterarion
r_gc=r+(m./(q.*B_gcz)).*cross(v,B./B_mag); %corssed values at particle postion
end

%Ode45 for Guiding Center 
Tmax2=t(end);
options2=odeset('RelTol',1.0e-10,'AbsTol',1.0e-14);
[t1,z]=ode45(@GuidingCenter,[Tmin,Tmax2],[r(1,:)],options2,v_abs(i,:));
[t2,z1]=ode45(@GuidingCenter,[Tmin,Tmax2],[r_gr(1,:)],options2,v_abs(i,:)); %can pass uniterated input this way
[t3,z2]=ode45(@GuidingCenter,[Tmin,Tmax2],[r_gc(1,:)],options2,v_abs(i,:));
%Particle Guiding Center
r_trace=[z(:,1), z(:,2), z(:,3)];%m
r_gr_trace=[z1(:,1), z1(:,2), z1(:,3)];%m 
r_gc_trace=[z2(:,1), z2(:,2), z2(:,3)];%m

%Guiding Center ode45 Tracing Error 
%Tracing error Arc
arc=acos((dot(r0,r(end,:),2))./(r0_mag*r_mag(end,:)));%Rad
z_mag=sqrt(sum(z.^2,2));
z1_mag=sqrt(sum(z1.^2,2));
z2_mag=sqrt(sum(z2.^2,2));
arcT=acos((dot(z(1,:),z(end,:),2))./(z_mag(1)*z_mag(end)));%Rad
arcT1=acos((dot(z1(1,:),z1(end,:),2))./(z1_mag(1)*z1_mag(end)));%Rad
arcT2=acos((dot(z2(1,:),z2(end,:),2))./(z2_mag(1)*z2_mag(end)));%Rad
distance_arcTrace=2*pi*r(1,1)*(arcT/(2*pi));
distance_arcTrace1=2*pi*r_gr(1,1)*(arcT1/(2*pi));
distance_arcTrace2=2*pi*r_gc(1,1)*(arcT2/(2*pi));
distance_arcActual=2*pi*r(1,1)*(arc/(2*pi));
distance_arcActual1=2*pi*r_gr(1,1)*(arc/(2*pi));
distance_arcActual2=2*pi*r_gc(1,1)*(arc/(2*pi));
distance_pos(i,:)=abs(distance_arcActual-distance_arcTrace);
distance_gr(i,:)=abs(distance_arcActual1-distance_arcTrace1);
distance_gc(i,:)=abs(distance_arcActual2-distance_arcTrace2);
%}

%Calculate magnetic moment mu
%M=m*v_perpendicular^2/(2*B);% Am^2 magnetic moment
mu=m.*(v_mag.^2./(Bz));
%Magnetic moment of gyrocenter/radius
%create new Magntic Field Function using r_gr
[B_grx, B_gry, B_grz] = mag_vec(r_gr);
B_gr=[B_grx, B_gry, B_grz];
B_gr_mag=sqrt(sum(B_gr.^2,2));
mu_gr=m.*(v_mag.^2./(B_grz));
%Magnetic moment for guding center 
mu_gc=m.*(v_mag.^2./(B_gcz));
%Second Order Approximation for Adiabatic Invarient %eq A1 from 1.5054594
%n=3 for equatorial magnetic dipole
%v_phi below is phi in cylindrical coorinates(rho, phi, z)
%tet=atan2(y(:,2),y(:,1));
%vtet=-y(:,4).*sin(tet)+y(:,5).*cos(tet);
v_phi=(cross(r,v))./(r_mag); 
I=(m)*(v_mag.^2./B_mag).*(1+((3*v_phi(:,3)*m)./(q*r_mag.*B_mag))+((m/(q*r_mag.*B_mag)).^2)*(6*v_phi(:,3).^2+3.75*v_mag.^2)); %using Lorentz trajecory position for B
I2=(m)*(v_mag.^2./B_gcz).*(1+((3*v_phi(:,3)*m)./(q*r_mag.*B_gcz))+((m/(q*r_mag.*B_gcz)).^2)*(6*v_phi(:,3).^2+3.75*v_mag.^2));%using gc position for B
I_mag=sqrt(sum(I.^2,2));

%Calculate Canonical Angular Momentum P
P=m.*r_mag.*v_phi(:,3)-(q.*M./r_mag);%From cjp-2017-0415

%Exact anayltical value of mu
k=(4*q*M/m).*(v_mag)./((P/m).^2); %unitless
[E_type1, E_type2] = ellipke(2.*k./(1+k));
mu_exact=abs((q/pi)*P/m).*(pi-2.*sqrt(1+k).*E_type2);

%Magnetic Moment Ecaxt Value - Mean for each velocity 
MUexact_pos(i,:)=max(abs(mu_exact(1)-max(mu)),abs(mu_exact(1)-min(mu)));
MUexact_gr(i,:)=max(abs(mu_exact(1)-max(mu_gr)),abs(mu_exact(1)-min(mu_gr)));
MUexact_gc(i,:)=max(abs(mu_exact(1)-max(mu_gc)),abs(mu_exact(1)-min(mu_gc)));
MUexact_I(i,:)=max(abs(mu_exact(1)-max(I)),abs(mu_exact(1)-min(I)));
MUexact_I2(i,:)=max(abs(mu_exact(1)-max(I2)),abs(mu_exact(1)-min(I2)));

%Where to start gyrocenter approcimation? compare:
%@ what position is local Magnetic Moment same as exact Value
%anaylically use mu equation, isolate for r_mag, feed mu exact 
%v_mag is contant so can use any point, value
%-min(r_mag))/(max(r_mag)-min(r_mag)) is normalization
r_mag_mu(i,:)=(((mu_exact(1)*M)/(m*v_mag(1)^2)).^(1/3)-min(r_mag))/(max(r_mag)-min(r_mag)); %analyitally
%@ what position is local gyroperiod same as exact Value
%inverse of gyrofrequency which is freq=abs(q)*B/m, dont forget 2pi factor!
r_mag_gyroperiod(i,:)=((t(end)*abs(q)*M/(m*2*pi)).^(1/3)-min(r_mag))/(max(r_mag)-min(r_mag)); %analyitally
%@ what position is local drift period same as exact Value
%inverse of driftfrequency which is freq=(3/2)*(v^2*r)/A where A=abs(q)*M/m, dont forget 2pi factor!
stop_angle=abs(atan2(y(end,2),y(end,1)));
r_mag_driftperiod(i,:)=((2*2*pi/3)*(abs(q)*M/(m*v_mag(1)^2*(t(end)*2*pi/stop_angle)))-min(r_mag))/(max(r_mag)-min(r_mag));%analyitally
%}
end

%figure 17
MagFieldLines

%Plot 1 drift 
figure; %figure 20 %Set isterminal to 0 in events 
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
plot(y(:,1),y(:,2),'b')
plot(r_grx,r_gry,'r')
tet2=-0.1:0.1:2*pi;
plot(0+0.3*cos(tet2),5+0.3*sin(tet2),'k' )
plot(0+0.00*cos(tet2),5+0.00*sin(tet2),'k.' )
text(0.3, 5,'{\bf B}')
text(-0.5, 0,'{\nabla B}')
text(-0.1, -5,'{v_D}')

%Plot 1 gyroperiod
figure; %figure 21
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
plot(y(:,1),y(:,2),'b')
plot(r_grx,r_gry,'r')

%Gyrocenter iterations 
figure; %figure 22 23 FOR 23 SET V0(2) LINE 43 TO 2.5 AND 3.5 And iterations line 93 to iterations=1:10 
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
plot(y(:,1),y(:,2),'b')
plot(r_grx,r_gry,'r')
plot(r_gc(:,1),r_gc(:,2),'k')
Guiding_Center_Iterations
legend('Trajectory','1st Iteration','Converged Gyrocenter','2nd Iteration','3rd Iteration')
%Gyrocenter markers 
%figure 22
Plot_w_Markers

%Plot GC approximation 
figure; %figure 24
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
plot(y(:,1),y(:,2),'b')
plot(r_grx,r_gry,'r')
plot(r_gc(:,1),r_gc(:,2),'k')
plot(r_trace(:,1),r_trace(:,2),'b','LineWidth',2)
plot(r_gr_trace(:,1),r_gr_trace(:,2),'r','LineWidth',2)
plot(r_gc_trace(:,1),r_gc_trace(:,2),'k','LineWidth',2)

%Magntic Moments
figure; %figure 27 FOR right panel set v0(2)-(LINE 43) to 0.3
xlabel('Time/Gyroperiod');
ylabel('\mu/\mu_{Exact}');
xlim([0 1])
grid on;
hold on;
plot(t/t(end),mu./mu_exact,'b');
plot(t/t(end),mu_gr./mu_exact,'r');
plot(t/t(end),mu_gc./mu_exact,'k');
plot(t/t(end),mu_exact./mu_exact, 'c');
legend('Trajectory','1st Gyrocenter','Converged Gyrocenter','Exact')

%Magntic Moments
figure; %figure 28 
xlabel('Time/Gyroperiod');
ylabel('\mu/\mu_{Exact}');
xlim([0 1])
grid on;
hold on;
plot(t/t(end),mu./mu_exact,'b');
plot(t/t(end),mu_gr./mu_exact,'r');
plot(t/t(end),mu_gc./mu_exact,'k');
plot(t/t(end),I./mu_exact, 'm');
plot(t/t(end),I2./mu_exact, 'g');
plot(t/t(end),mu_exact./mu_exact, 'c');
legend('Trajectory','1st Gyrocenter','Converged Gyrocenter','2nd Order Taylor','2nd Order Taylor @ gc','Exact');

%Overlay Exact Postions 
%figure 31
Overlay_exact_positions
end

function [Bx, By, Bz] = mag_vec(r)
global M;
%Tsyganeko model 
SPS=0;
CPS=1;
P=r(:,1).*r(:,1); %m^2
T=r(:,2).*r(:,2); %m^2
U=r(:,3).*r(:,3); %m^2
V=3*r(:,3).*r(:,1); %m^2
Q=M./sqrt(P+T+U).^5; %Tm^2

Bx=Q.*((T+U-2.*P).*SPS-V.*CPS); %T
By=-3.*r(:,2).*Q.*(r(:,1).*SPS+r(:,3).*CPS); %T
Bz=Q.*((P+T-2.*U).*CPS-V.*SPS); %T
%{
r_mag=sqrt(sum(r.^2,2));
Bz=M./(r_mag.^3);%T
n=size(Bz);
Bx=zeros(n);
By=zeros(n);
%E=0; %V/m or N/C (mkg/s^3A) Electric Field, Assumed neglible
%}
end

function [delBx, delBy, delBz] = position_vec(r)
global M;
%Partial Derivative of Tarsagov?(fix spelling)model
%need to do parital derivative of B not Bx By Bz
%simlified by changing SPS=0; CPS=1;
%wolfram alpha 
P=r(:,1).*r(:,1); %m^2
T=r(:,2).*r(:,2); %m^2
U=r(:,3).*r(:,3); %m^2
Q=-3*M./sqrt(P+T+U).^7; %Tm^2

delBx=Q*((r(:,1)^3)-(4*P*r(:,3))+(r(:,1)*(T-5*r(:,2)*r(:,3)-4*U))+(r(:,3)*(T+U))); %T
delBy=Q*((P*(r(:,2)+r(:,3)))-(5*r(:,1)*r(:,2)*r(:,3))+(r(:,2)^3)-(4*T*r(:,3))-(4*r(:,2)*U)+(r(:,3)^3)); %T
delBz=Q*((r(:,1)^3)+(P*(r(:,2)+3*r(:,3)))+(r(:,1)*(T-4*U))+(r(:,2)^3)+(3*T*r(:,3))-(4*r(:,2)*U)-(2*r(:,3)^3)); %T
end

function dxdt = Lorentz(t,x)
%global Re;
%global m_e;
%global m_p;
%global q_e;
%global q_p;
global m;
global q;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %m
v=[x(4), x(5), x(6)]; %m
r_mag=sqrt(sum(r.^2,2));
v_mag=sqrt(sum(v.^2,2));

%Call Magnetic field Function
[Bx,By,Bz]=mag_vec(r);
B=[Bx,By,Bz];

%DIFFERENTIAL EQUATION 
%Lorentz Equation
%F=q(E + v x B)
%a=(q/m)(E + v x B)
a_L=(q/m)*(E + cross(v,B));
a_Lx=a_L(1);
a_Ly=a_L(2);
a_Lz=a_L(3);

%Magnetic Gradient Drift, guided center approach
%v_GradDrift= (m*v_perp^2/(2q*B^3))*(B x grad_B)

% e cross b-negligible as e of 1 particle vs 
%b of earth magnetic field, curvature-not yet assuming vertical field line, 
%no change in magnetic field or electric-so no polarization drift,

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; 
end

function dxdt = GuidingCenter(t,x,v_abs)
global m;
global q;

%global v_mag; %alternative way to get v_abs into ode45

%variables: x=x(1); y=x(2); z=x(3);  
r=[x(1), x(2), x(3)]; %m
r_mag=sqrt(sum(r.^2,2));

%Call Magnetic field Function
[Bx,By,Bz]=mag_vec(r);
B=[Bx,By,Bz];
B_mag=sqrt(sum(B.^2,2));
%Call paartial derivative of Magnetic field Function
[delBx,delBy,delBz]=position_vec(r);
delB=[delBx,delBy,delBz];

%DIFFERENTIAL EQUATION
%Vgc same direction as Vphi, break up into Vx and Vx 
Vgc=(m/q)*(v_abs^2/(2*B_mag^3))*cross(B,delB);
%return 3 column array make magntiude 
%vtet=-y(:,4).*sin(tet)+y(:,5).*cos(tet); rearrange this 

dxdt=[Vgc(1);Vgc(2);0.0]; %no x(4) for velocity so just set to 0 manualy
end

function [value,isterminal,direction] = events(t,y)
global q;
% Locate the time when height passes through zero in a decreasing direction
% and stop integration. 
tet=atan2(y(2),y(1));
%vtet=-y(4)*sin(tet)+y(5)*cos(tet);
%value = vtet*sign(-10); %stop when theta velocity zero again 
vrad=y(4)*cos(tet)+y(5)*sin(tet);
value = vrad*sign(-10); %stop when radial velocity zero again 
isterminal = 1;% 1 to stop the integration
if q>0 %direction -1 for ion, 1 for electon 
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
else
direction = 1;
end
end




