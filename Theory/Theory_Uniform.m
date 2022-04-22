function Theory_Uniform
clc

%GLOBALS 
global Re;
global mp;
global me;
global qp;
global qe;
global E;
global Beta;
global B0;

%PREPERATION
Re=6378;%km Earths Equatorial Radius
me0=9.109*10^-31;%kg mass of electron 
mp0=1.673*10^-27;%kg mass of proton;% 1850 roughly mass ratio of proton to electron 
qp=1.602*10^-19; %C
qe=-1.602*10^-19; %C
E=0; %V/m  Electric Field
Beta=1*10^-2; %T/m Gradient Magnetic Field Strength 
B0=1*10^-1; %T Baseline Magnetic Field Strength 
c=3*10^8; %m/s speed of light 

%IN TERMS OF ENERGY  
energy_p=1*1.0E6; %Mev %CHANGE TO 4*1.0E6 for FIG 3 
%apply %m=m0/sqrt(1-(v0_mag^2/c^2))
mp=(energy_p*1.602E-19)/c^2+mp0; %kg relativistic mass of proton 
vp_mag=sqrt((1-(mp0/mp)^2)*c^2); %m/s velocity of proton 
energy_e=10*1.0E6;%Mev   %CHANGE TO 40*1.0E6 for FIG 3 
me=(energy_e*1.602E-19)/c^2+me0; %kg %relativistic mass of electron 
ve_mag=sqrt((1-(me0/me)^2)*c^2);

Tp=2*pi*mp/(B0*abs(qp)); %gyroperiod of proton 
Te=2*pi*me/(B0*abs(qe)); %gyroperiod of electron 
rgp=mp*vp_mag/(B0*abs(qp)); %gyroradius of proton
rge=me*ve_mag/(B0*abs(qe)); %gyroradius of electron 

%INITIAL CONDITIONS PROTON 
r0p=[10, 0, 0];%m  
v0P=[0, vp_mag, 0];%m/s 
Tmin=0;Tmax=0.000002;
%CALL ODE45 FUNCTION PROTON 
[t,y]=ode45(@Lorentz,[Tmin,Tmax],[r0p,v0P],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14));
%RESULTS 
rp=[y(:,1), y(:,2), y(:,3)];%m 
vp=[y(:,4), y(:,5), y(:,6)];%m
%CALL MAGNETIC FIELD PROTON
[Bx, By, Bz] = mag_vec(rp); %T
B=[Bx, By, Bz]; %T
B_mag=sqrt(sum(B.^2,2));%T

%INITIAL CONDITIONS ELECTRON 
r0e=[15, 0, 0];%m  
v0e=[0, ve_mag, 0];%m/s 
Tmax1=0.0000001;
%Call ODE45 FUNCTION ELECTRON 
[t1,y1]=ode45(@Lorentz1,[Tmin,Tmax1],[r0e,v0e],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14));
%RESULTS 
re=[y1(:,1), y1(:,2), y1(:,3)];%m 
ve=[y1(:,4), y1(:,5), y1(:,6)];%m 
%CALL MAGNETIC FIELD ELECTRON 
[Bx1, By1, Bz1] = mag_vec(re); %T
B1=[Bx1, By1, Bz1];%T
B_mag1=sqrt(sum(B1.^2,2));%T

%INITIAL CONDITIONS PROTON 3D
r2=[10, 0, 0];%m  
v2=[0, vp_mag*0.9, vp_mag*0.1];%m/s 
Tmin2=0;Tmax2=0.000003;
%CALL ODE45 FUNCTION PROTON 
[t2,y2]=ode45(@Lorentz,[Tmin2,Tmax2],[r2,v2],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14));


%FIGURES
%Plot Trajectory in Unform Field 
figure; %figure 1
xlabel('X (m)');
ylabel('Y (m)');
grid on;
axis equal;
hold on;
plot(y(:,1),y(:,2),'b')
plot(y1(:,1),y1(:,2),'r')
legend('Proton', 'Electron','AutoUpdate','off')
tet2=-0.1:0.1:2*pi;
plot(14+0.1*cos(tet2),1+0.1*sin(tet2),'k' )
plot(14+0.00*cos(tet2),1+0.00*sin(tet2),'k.' )
text(14.1, 1,'{\bf B}')

%Plot Trajectory
figure; %figure 2
plot3(y2(:,1),y2(:,2),y2(:,3),'b')
hold on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;
axis equal;
end

function [Bx, By, Bz] = mag_vec(r)
global Beta;
global B0;
Bz=B0;%+Beta*r(:,1); %T
n=size(Bz);
Bx=zeros(n); %T
By=zeros(n); %T
end

function dxdt = Lorentz(t,x)
global mp;
global qp;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see below 
r=[x(1), x(2), x(3)]; %m
v=[x(4), x(5), x(6)]; %m

%Call magentic field function 
[Bx,By,Bz]=mag_vec(r);
B=[Bx,By,Bz];

%DIFFERENTIAL EQUATION 
%Lorentz Equation
%F=q(E + v x B)
%a=(q/m)(E + v x B)
a_L=(qp/mp)*(E + cross(v,B));
a_Lx=a_L(1);
a_Ly=a_L(2);
a_Lz=a_L(3);

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; 
end

function dxdt = Lorentz1(t,x)
global me;
global qe;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %m
v=[x(4), x(5), x(6)]; %m

%Call Magnetic field Function
[Bx,By,Bz]=mag_vec(r);
B=[Bx,By,Bz];

%DIFFERENTIAL EQUATION 
%Lorentz Equation
%F=q(E + v x B)
%a=(q/m)(E + v x B)
a_L=(qe/me)*(E + cross(v,B));
a_Lx=a_L(1);
a_Ly=a_L(2);
a_Lz=a_L(3);

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; 
end






