function VaryVelocity_RhoZ_Plots
clc
%GLOBALS 
global Re;
global m;
global q;
global M;
global E;

%CONSTANTS
Re=6378;%km Equatorial Radius of Earth
m0=1.6726231*10^-27;%kg rest mass of proton 
q=1.602*10^-19; %C charge of proton 
M=31000.*10^-9; %TRe^3 Magnetic Field Strength of Earth 
E=0; %V/m or N/C (mkg/s^3A) Electric Field
c=3.0E8; %speed of light m/s
 
%INITIAL CONDITIONS
%position
r0=[7, 0, 0];%Re
r0_mag=sqrt(sum(r0.^2,2));%Re
rho0_mag=sqrt(r0(1)^2 + r0(2)^2);%Re
%velocity
Pangle=pi/4; %Rad Pitch Angle 
v0_mag=9.52481; %Re/s
v0=[0, v0_mag*sin(Pangle), v0_mag*cos(Pangle)];%Re/s
%Other
m=m0*(1/sqrt(1-((v0_mag*Re)^2/(3*10^5)^2)));%relativistic mass
Energy=(m-m0)*c*c/(1.602E-19)/1.0E6; % Kinetic energy in MeV

%Call ODE45 Functions
Tmax=1300/v0_mag; 
if v0_mag == 9.52481
Tmax=1000; %Sufficiently high for stopping condition
end
options=odeset('Events',@events,'RelTol',1.0e-10,'AbsTol',1.0e-14); %stops when particle hits Earth's atmohsphere
[t1,y1]=ode45(@Lorentz,[0,Tmax],[r0,v0],options);
t1(end)
%RESULTS
r=[y1(:,1), y1(:,2), y1(:,3)];%Re
size=[(length(r)),1]; %size of matrix
rho=[y1(:,1), y1(:,2), zeros(size)];%Re
rho_mag = sqrt(sum(rho.^2,2));%Re

%RHO-Z PLOTS FOR VELCOITIES OF INTEREST
figure;%FIGURE 70 SET LINE 25 TO v0_mag=9.52481; %Re/s
       %FIGURE 71 LEFT PANEL SET LINE 25 TO v0_mag=4.33731; %Re/s
       %FIGURE 71 RIGHT PANEL SET LINE 25 TO v0_mag=9; %Re/s
plot(rho_mag,y1(:,3),'b');
hold on;
xlabel('\rho (Re)');
ylabel('Z (Re)');
grid on; 
if v0_mag == 9.52481
t_circle=0:0.01:pi;
rho_circleA=1.015679*sin(t_circle); 
z_circleA=1.015679*cos(t_circle); 
plot(rho_circleA,z_circleA,'k');
end

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
%Stop Integration when particle reachs 100km above Earths atmosphere 
value = 1.015679 -sqrt(y(1)^2+y(2)^2+y(3)^2); %In terms of Re
isterminal = 1;% 1 to stop the integration
direction = 0;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end




