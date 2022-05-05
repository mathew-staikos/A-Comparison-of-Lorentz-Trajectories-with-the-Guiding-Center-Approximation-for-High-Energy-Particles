function Sections_6dot2_to_6dot5
clear
global m q M Re E;
Re=6378;%km Equatorial Radius of Earth
m0=1.6726231*10^-27;%kg rest mass of proton 
q=1.602*10^-19; %C charge of proton 
M=31000.*10^-9; %TRe^3 Magnetic Field Strength of Earth 
E=0; %V/m or N/C (mkg/s^3A) Electric Field

%Initial Conditions 
Pangle=pi/4;  %Rad Pitch Angle 
phase=0; %Rad
v0_mag=1; %Re/s velocity magnitude 
%position components 
x0=7; %Re
y0=0; %Re
z0=0; %Re
r0=sqrt(x0^2+y0^2+z0^2); %Re
%velocity comp0nents 
vx0=v0_mag*sin(phase); %Re/s
vy0=v0_mag*sin(Pangle)*cos(phase); %Re/s
vz0=v0_mag*cos(Pangle); %Re/s
v0=sqrt(vx0^2+vy0^2+vz0^2); %Re/s
%relativistic mass 
m=m0*(1/sqrt(1-((vy0*Re)^2/(3*10^5)^2))); %kg 

%Call ODE 45  
Tmax=400; %Sufficiently high to reach stopping conditon
options=odeset('Events',@events,'RelTol',1.0e-10,'AbsTol',1.0e-14); %stopping condition
[t,y]=ode45(@Lorentz,[0 Tmax],[x0,y0,z0,vx0,vy0,vz0],options);
rvec=[y(:,1), y(:,2), y(:,3)]; %Re matrix of position vectors 
vvec=[y(:,4), y(:,5), y(:,6)]; %Re/s matrix of velocity vectors 
r=sqrt(y(:,1).^2+y(:,2).^2+y(:,3).^2); %Re matrix of position magnitudes 
v=sqrt(y(:,4).^2+y(:,5).^2+y(:,6).^2); %Re/s matrix of velocity magnitudes 
rho=sqrt(y(:,1).^2+y(:,2).^2); %Re matrix of position magnitudes in XY plane

%Magnetic field at particle positions
[Bx, By, Bz] = mag_vec(rvec); %T
B=[Bx, By, Bz]; %T
B_mag=sqrt(sum(B.^2,2)); %T 

%Additonal 
alpha = acos((dot(vvec,B,2))./(v.*B_mag));%Rad Pitch angle at every position
v_para = v.*cos(alpha); %Re/s Parallel velocity at every position
v_perp = v.*sin(alpha); %Re/s Perpendicular velocity at every position

%Magnetic Moment 
mu=m.*((v_perp*Re*1.0E3).^2)./B_mag; % Am^2  

%Canonical Angular Momentum 
v_phi_t=(y(:,1).*y(:,5)-y(:,2).*y(:,4))./rho; %Re/s velocity in phi direction at every position
tet=acos(y(:,3)./r); %Rad Collatitude at every position
p=(m*r.*sin(tet).*v_phi_t-q*M*(sin(tet)).^2./r);%kgRe2/s  Canonical Angular Momentum 
%Radial Bounds

rr_min=( p(1)+sqrt(p(1)^2+4*q*M*m*v0.*(sin(tet)).^3))./(2*m*v0*sin(tet)); %Re Minimum Radial Bounds 
rr_max=(-p(1)-sqrt(p(1)^2-4*q*M*m*v0.*(sin(tet)).^3))./(2*m*v0*sin(tet)); %Re Maximum Radial Bounds

%Plots
figure(1) %figure 61 left panel 
plot(t/t(end),(p-p(1))/p(1))
xlabel('Time/Half Bounce Period')
ylabel('(P-P_{0})/P_{0}')
grid on 

figure(2) %figure 61 right panel 
plot(t/t(end),(v-v(1))/v(1))
xlabel('Time/Half Bounce Period')
ylabel('(v-v_{0})/v_{0}')
grid on 

figure(3); %Figure 62
plot3(y(:,1),y(:,2),y(:,3),'b');
hold on;
xlabel('X (Re)');
ylabel('Y (Re)');
zlabel('Z (Re)');
grid on; 
axis equal;

figure(4); %Figure 63
plot(t/t(end),mu,'b')
xlabel('Time/Half Bounce Period')
ylabel('\mu (Am^2)');
grid on; 

figure(5) %figure 60 + 64(see lineS 107 + 108)
plot(rho,y(:,3),'b',rr_min.*sin(tet),rr_min.*cos(tet),'k',rr_max.*sin(tet),rr_max.*cos(tet),'k')
hold on
grid on
xlabel('\rho (Re)');
ylabel('Z (Re)');


%Phase Distribution Attempt 
%New phase, Same pitch angle
phase2=pi; %Rad
%Inital Conditions
y1=0; %Re
z1=0; %Re
vx1=v0_mag*sin(phase2); %Re/s
vy1=v0_mag*sin(Pangle)*cos(phase2); %Re/s
vz1=v0_mag*cos(Pangle); %Re/s
x1=(p(1)+sqrt(p(1)^2+4*q*M*m*vy1))./(2*m*vy1); %Re
%Call ODE 45
[t,y]=ode45(@Lorentz,[0 Tmax],[x1,y1,z1,vx1,vy1,vz1],options);
%Remove '%' for figure 64
%plot(sqrt(y(:,1).^2+y(:,2).^2), y(:,3),'r')
end

function dxdt = Lorentz(t,x)
global m;
global q;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %Re
v=[x(4), x(5), x(6)]; %Re

%Call Magnetic field Function
[Bx,By,Bz]=mag_vec(r); %T
B=[Bx,By,Bz]; %T

%DIFFERENTIAL EQUATION 
%Lorentz Equation
a_L=(q/m)*(E + cross(v,B)); %Re/s^2
a_Lx=a_L(1); %Re/s^2
a_Ly=a_L(2); %Re/s^2
a_Lz=a_L(3); %Re/s^2

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; %Output new positon and velcoity
end

function [Bx, By, Bz] = mag_vec(r)
global M;
SPS=0; %No tilt of Earth relative to field lines 
CPS=1; %No tilt of Earth relative to field lines 
P=r(:,1).*r(:,1); %Re^2
T=r(:,2).*r(:,2); %Re^2
U=r(:,3).*r(:,3); %Re^2
V=3*r(:,3).*r(:,1); %Re^2
Q=M./sqrt(P+T+U).^5; %T(Re^-2)

Bx=Q.*((T+U-2.*P).*SPS-V.*CPS); %T
By=-3.*r(:,2).*Q.*(r(:,1).*SPS+r(:,3).*CPS); %T
Bz=Q.*((P+T-2.*U).*CPS-V.*SPS); %T
end

function [value,isterminal,direction] = events(t,y)
value = y(3); %stop when hit equitorial plane  
isterminal = 1;%stop the integration
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end
