function phase_shift %FIGURE 45 
clear
global m;
global q;
global M;

Re=6378;%km Equatorial Radius of Earth
m0=1.6726231*10^-27; %kg rest mass of proton 
q=1.602*10^-19; %C q_p=1.602*10^-19;%C(As) Charge of proton
M=31000.*10^-9; %TRe^3 Magnetic Field Strength of Earth 
E=0; %V/m or N/C (mkg/s^3A) Electric Field
c=3.0E8; %speed of light m/s

%Initial Position 
x0=7; %Re
y0=0; %Re
%Initial Velocity 
vx0=0; %Re/s
vy0=6; %Re/s
%relativistic mass 
m=m0*(1/sqrt(1-((vy0*Re)^2/(3*10^5)^2)));%kg relativistic mass
A=q*M/m; %simplyfy notation  

%Call ode 45
tmax=1000; %sufficiently high to allow events to stop integration
optionsGP=odeset('Events',@eventsGP,'RelTol',1.0e-10,'AbsTol',1.0e-14);%stop after gyroperiod
[t1,y1]=ode45(@lorentz, [0 tmax], [x0 y0 vx0 vy0],optionsGP);

%Phase Distribution preperations
r0=sqrt(x0^2+y0^2); 
v0=sqrt(vx0^2+vy0^2);
r_st=r0/(1-r0^2*vy0/A);
v_st=A/r_st^2;
r_max=2*r_st/(1+sqrt(1-4*v0/v_st));
r_min=2*r_st/(1+sqrt(1+4*v0/v_st));
p=(m*r0*vy0)-(q*M/r0);

%change in phi over gyroperiod
ph_range=atan2(y1(:,2),y1(:,1));
max(ph_range);
min(ph_range);
ph=linspace(min(ph_range),max(ph_range));

start = linspace(r_min,r_max);

%Plot particle trajectory, rmin and rmax bounds over phi range, and
%starting line for all phases between rmin and rmax along positive x axis
plot(y1(:,1),y1(:,2),'b',r_min*cos(ph),r_min*sin(ph),'k--',r_max*cos(ph),r_max*sin(ph),'k--',start,zeros(size(start)),'k')
grid on 
axis('equal')
xlabel('X (Re)')
ylabel('Y (Re)')
hold on

%Plot 6 local velocity vectors between phases 0 and 2pi, along with their
%phases shefted to align with positve x-axis starting line and new velocity
%vectos 
s=size(y1);
iterations=6;
for i=1:iterations
v_tangents=linspace(1,s(1),iterations);
kp=v_tangents(i);
kp_integer=round(kp); %remove decimals for indexing
if kp_integer==203 %manual shifting to increase visual appeal of plot
kp_integer=240;
end    
if kp_integer==608 %manual shifting to increase visual appeal of plot
kp_integer=660;
end    
quiver(y1(kp_integer,1),y1(kp_integer,2), y1(kp_integer,3),y1(kp_integer,4), 0.05,'k');%all tangential velocity
quiver(y1(1,1),y1(1,2), y1(1,3),y1(1,4), 0.05,'r'); %phase shifted velocity for phase = 0
hold on
if i==1
text(6.6,0.0,'{\delta} = 0.0'); 
end

if i==2
phi_loc=atan2(y1(kp_integer,2),y1(kp_integer,1)); %calculate change in phi at local positon 
hold on
r_mag=sqrt(y1(kp_integer,1)^2+y1(kp_integer,2)^2); %calculate radius  at local positon 
phi_vec=linspace(0,phi_loc); %turn change in phi to vector to grapgh 
plot(r_mag*cos(phi_vec),r_mag*sin(phi_vec),'r--'); %draw phase shift to x-axis
hold on
delta=(kp_integer/s(1))*2*pi; %use to detrmine valuye for text below
text(7.6,0.63,'{\delta} = 1.5'); %look and change (real_phi)
vx=v0*sin(delta); %new velocity vector 
vy=v0*cos(delta); %new velocity vector 
v_phi=(cross([y1(kp_integer,1),y1(kp_integer,2), 0],[y1(kp_integer,3), y1(kp_integer,4), 0]))./(r_mag); %v_phi 
x=(p+sqrt(p^2+4*m*v_phi(3)*q*M))/(2*m*v_phi(3)); %use canaical momentum to get new x positon
quiver(x, 0, vx, vy, 0.05,'r'); %phase shifted velocity
end

if i==3
phi_loc=atan2(y1(kp_integer,2),y1(kp_integer,1));   
text(8.15,0.42,'{\delta} = 2.5'); 
hold on
r_mag=sqrt(y1(kp_integer,1)^2+y1(kp_integer,2)^2);
phi_vec=linspace(0,phi_loc);
plot(r_mag*cos(phi_vec),r_mag*sin(phi_vec),'r--');
hold on
delta=(kp_integer/s(1))*2*pi; 
vx=v0*sin(delta);
vy=v0*cos(delta);
v_phi=(cross([y1(kp_integer,1),y1(kp_integer,2), 0],[y1(kp_integer,3), y1(kp_integer,4), 0]))./(r_mag); 
x=(p+sqrt(p^2+4*m*v_phi(3)*q*M))/(2*m*v_phi(3)); 
quiver(x, 0, vx, vy, 0.05,'r');
end

if i==4
phi_loc=atan2(y1(kp_integer,2),y1(kp_integer,1));
text(8.25,-1.95,'{\delta} = 4.1'); 
hold on
r_mag=sqrt(y1(kp_integer,1)^2+y1(kp_integer,2)^2);
phi_vec=linspace(0,phi_loc);
plot(r_mag*cos(phi_vec),r_mag*sin(phi_vec),'r--');
hold on
delta=(kp_integer/s(1))*2*pi; 
vx=v0*sin(delta);
vy=v0*cos(delta);
v_phi=(cross([y1(kp_integer,1),y1(kp_integer,2), 0],[y1(kp_integer,3), y1(kp_integer,4), 0]))./(r_mag); 
x=(p+sqrt(p^2+4*m*v_phi(3)*q*M))/(2*m*v_phi(3)); 
quiver(x, 0, vx, vy, 0.05,'r');
end

if i==5
phi_loc=atan2(y1(kp_integer,2),y1(kp_integer,1)); 
text(7.33,-1.83,'{\delta} = 5.0'); 
hold on
r_mag=sqrt(y1(kp_integer,1)^2+y1(kp_integer,2)^2);
phi_vec=linspace(0,phi_loc);
plot(r_mag*cos(phi_vec),r_mag*sin(phi_vec),'r--');
hold on
delta=(kp_integer/s(1))*2*pi;
vx=v0*sin(delta);
vy=v0*cos(delta);
v_phi=(cross([y1(kp_integer,1),y1(kp_integer,2), 0],[y1(kp_integer,3), y1(kp_integer,4), 0]))./(r_mag);
x=(p+sqrt(p^2+4*m*v_phi(3)*q*M))/(2*m*v_phi(3)); 
quiver(x, 0, vx, vy, 0.05,'r');
end
end
end

function dydt = lorentz(t,y)
global m q;
r=sqrt(y(1)^2+y(2)^2); %radial position of particle 
Bz=magneticfield(r);%T magnetic field at positon r
dydt=[y(3);y(4); (q/m)*(y(4)*Bz); (q/m)*(-y(3)*Bz)];%Outputs new positon and velocity vectors
end

function Bz=magneticfield(r)
global M;
Bz=M/r^3; %T
end

function [value,isterminal,direction] = eventsGP(t,y)
% stop intgration after gyroperiod
tet=atan2(y(2),y(1));
vrad=y(3)*cos(tet)+y(4)*sin(tet);
value = vrad*sign(-10); %stop when radial velocity zero again 
isterminal = 1;% 1 to stop the integration
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both

end


