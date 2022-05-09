function Phases
clc
%GLOBALS 
global Re;
global E;
global m;
global q;
global M;

Re=6378;%km Equatorial Radius of Earth
m0=1.6726231*10^-27; %kg rest mass of proton 
q=1.602*10^-19; %C q_p=1.602*10^-19;%C(As) Charge of proton
M=31000.*10^-9; %TRe^3 Magnetic Field Strength of Earth 
E=0; %V/m or N/C (mkg/s^3A) Electric Field
c=3.0E8; %speed of light m/s

%INITIAL CONDITIONS 
x0=7.0; %Re
y0=0.0; %Re
vx0=0.0; %Re
vy0=3.5; %Re/s %CHNAGE TO 3.5

m=m0*(1/sqrt(1-((vy0*Re)^2/(3*10^5)^2)));%kg relativistic mass 
A=q*M/m;
Energy=(m-m0)*c*c/(1.602E-19)/1.0E6; % Kinetic energy in MeV

%Prep fpr Phase shift using Canaocial angular momentum
r0_mag=sqrt(x0^2+y0^2); %Re
v0_mag=sqrt(vx0^2+vy0^2); %Re/s
r_st=r0_mag/(1-r0_mag^2*vy0/A);  %r star
v_st=A/r_st^2; %v star
r_max=2*r_st/(1+sqrt(1-4*v0_mag/v_st));%maximum bounds of particle 
r_min=2*r_st/(1+sqrt(1+4*v0_mag/v_st));%minimum bounds of particle 
p=m*r0_mag*vy0-q*M/r0_mag; %kgRe/s canonical angular momentum  

%Phase Distibution drift period stopping condition
iterations=1000; %number of phases 
colour=jet(iterations); %different color for each phase
for i=1:iterations 
angles=linspace(0,pi/8,iterations);%range of phases 
%inital conditions for each phase 
vx1=v0_mag*sin(angles(i)); 
vy1=v0_mag*cos(angles(i));
x1=(p+sqrt(p^2+4*m*vy1*q*M))/(2*m*vy1); 
y1=0.0;

%CALL ODE 45 to intgrate equation of motion 
r0=[x1, y1, 0];%Re
v0=[vx1,vy1, 0];%Re/s 
Tmin=0;Tmax=10000; %sufficiently high to reach stopping conditions first
%Determine integration time for velocity
optionsY1=odeset('Events',@eventsY1,'RelTol',1.0e-10,'AbsTol',1.0e-14); %integrate 1 quadrant 
optionsX=odeset('Events',@eventsX,'RelTol',1.0e-10,'AbsTol',1.0e-14); %integrate 1 quadrant
optionsY2=odeset('Events',@eventsY2,'RelTol',1.0e-10,'AbsTol',1.0e-14); %integrate 1 quadrant
optionsDP=odeset('Events',@eventsDP,'RelTol',1.0e-10,'AbsTol',1.0e-14); %integrate 1 quadrant
[t,y]=ode45(@Lorentz,[Tmin,Tmax],[r0, v0],optionsY1);
new_r0=[y(end,1), y(end,2), y(end,3)];
new_v0=[y(end,4), y(end,5), y(end,6)];
[t1,y1]=ode45(@Lorentz,[Tmin,Tmax],[new_r0,new_v0],optionsX);
new_r0=[y1(end,1), y1(end,2), y1(end,3)];
new_v0=[y1(end,4), y1(end,5), y1(end,6)];
[t2,y2]=ode45(@Lorentz,[Tmin,Tmax],[new_r0,new_v0],optionsY2);
new_r0=[y2(end,1), y2(end,2), y2(end,3)];
new_v0=[y2(end,4), y2(end,5), y2(end,6)];
[t3,y3]=ode45(@Lorentz,[Tmin,Tmax],[new_r0,new_v0],optionsDP);
%Combine data
t1=t1+t(end);
t2=t2+t1(end);
t3=t3+t2(end);
t_total=[t;t1;t2;t3]; 
y_total=[y;y1;y2;y3]; 

if i ==1
t_integrate=t_total(end);
end

%FIGURE 46-Zoom In
plot(y_total(:,1),y_total(:,2),'color',colour(i,:)) 
axis('equal')
grid on
hold on
end

%Phase Distibution intrgation time = Drift period for 1st phase 
iterations=20; %number of phases 
for i=1:iterations %first point below must be 0 in order to skip it w vx1 call
angles=linspace(0,2*pi,iterations);%range of phases 
%inital conditions for each phase 
%non-axissymmetric
vx=v0_mag*sin(angles(i)); 
vy=v0_mag*cos(angles(i));
x=(p+sqrt(p^2+4*m*vy*q*M))/(2*m*vy); 
y=0.0;

if v0_mag==6
t_integrate =  33.3993 ; %s %non-axissymmetric
t_integrate_d = 42.8352; %s %dipole
end
if v0_mag==3.5
t_integrate = 102.3964; %s %non-axissymmetric 
t_integrate_d = 134.6146; %s %dipole slower
end

%CALL ODE 45 to intgrate equation of motion 
r0=[x, y, 0];%Re
v0=[vx,vy, 0];%Re/s 
 
[t,y]=ode45(@Lorentz,[0,t_integrate],[r0,v0],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14)); %small phi gradient
[td,yd]=ode45(@LorentzDipole,[0,t_integrate_d],[r0,v0],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14)); %dipole field

%Phase Endpoints 
endX(i)= y(end,1);
endY(i)= y(end,2);
endXd(i)= yd(end,1);
endYd(i)= yd(end,2);

%Endpoint phases 
anglesF(i)=atan2(y(end,5),y(end,4))+pi; 
anglesF_d(i)=atan2(yd(end,5),yd(end,4))+pi; 
end

figure; %FIGURE 47 LEFT PANEL 
hold on
grid on
xlabel('X (Re)');
ylabel('Y (Re)');
plot(endXd(1),endYd(1),'*','Color','k') 
plot(endXd(iterations*0.25),endYd(iterations*0.25),'*','Color','r') 
plot(endXd(iterations*0.5),endYd(iterations*0.5),'*','Color','g') 
plot(endXd(iterations*0.75),endYd(iterations*0.75),'*','Color','m') 
legend('\delta = 0', '\delta = \pi/2', '\delta = \pi', '\delta = 3\pi/2','AutoUpdate','off')
plot(endXd,endYd,'b') 
axis equal

%Endpoints x vs y positions    
figure; %FIGURE 48 %CHANGE LINE 21 FROM vy0=6;TO vy0=3.5 FOR RIGHT PANEL OF FIGURE 47
hold on
grid on
xlabel('X (Re)');
ylabel('Y (Re)');
plot(endX(1),endY(1),'*','Color','k') 
plot(endX(iterations*0.25),endY(iterations*0.25),'*','Color','r') 
plot(endX(iterations*0.5),endY(iterations*0.5),'*','Color','g') 
plot(endX(iterations*0.75),endY(iterations*0.75),'*','Color','m') 
legend('\delta = 0', '\delta = \pi/2', '\delta = \pi', '\delta = 3\pi/2','AutoUpdate','off')
plot(endX,endY,'b') 
axis equal

%Endpoints intital phase vs final phase
figure; %FIGURE 49 LEFT PANEL %CHANGE LINE 21 FROM vy0=6;TO vy0=3.5 FOR RIGHT PANEL OF FIGURE 47
hold on
grid on
xlabel('\delta_{i} (Rad)');
xlim([0 2*pi])
xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('\delta_{f} (Rad)');
ylim([0 2*pi])
yticks([0 pi/2 pi 3/2*pi 2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(angles,anglesF,'b') 
figure; %FIGURE 49 RIGHT PANEL %CHANGE LINE 21 FROM vy0=6;TO vy0=3.5 FOR RIGHT PANEL OF FIGURE 47
hold on
grid on
xlabel('\delta (Rad)');
xlim([0 2*pi])
xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('\delta (Rad)');
ylim([0 2*pi])
yticks([0 pi/2 pi 3/2*pi 2*pi])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(angles,anglesF_d,'b') 
end

function [Bx, By, Bz] = smallphigrad(r)
global M;

time=0; %min %No time variation in magentic field
r_mag=sqrt(sum(r.^2,2)); %Re 
phi_d=acosd(r(:,1)./r_mag); 

%f(phi)
c0=(1+tanh(11/5))^2;
c1=(1+tanh(110/5))^2;
f_phi=(1/c0).*(1+tanh((55-(phi_d-180))/20)).*(1+tanh((55+(phi_d-180))/20));
f_phi=(f_phi/c1).*(1+tanh((110-(phi_d-180))/5)).*(1+tanh((110+(phi_d-180))/5));

%from reference paper
b1=648.9661-183.7327*tanh(1.6101*time-3.8109);
b2=0.2019-0.0859*tanh(1.6258*time-4.1687);
b3=3.3338+0.2906*tanh(1.6013*time-3.9748);
b4=1.3820-0.2966*tanh(1.8191*time-5.1973);
b5=93.1713+50.4476*tanh(1.6306*time-4.3549);
a1=(r_mag+b2)/b3;

%magnetic field strength
Btheta=(M./r_mag.^3)-f_phi.*10^-9.*(b1./a1.^2.884).*((exp(a1.^2.46)-b4)./(exp(a1.^2.46)+b5)).^2.065; %T 

%Bx=Btheta*cos(theta)*cos(phi); %T
%By=Btheta*cos(theta)*sin(phi); %T
%theta is pi/2 or 90 cos90=0  
Bx=zeros(size(Btheta)) ; %T
By=zeros(size(Btheta)) ; %T

if r_mag >= 3.0
Bz=Btheta;%T 
else
Bz=800.0.*10^-9;%T %pure dipole at r<3
end
end

function [Bx, By, Bz] = dipole(r)
global M;
%radial position 
r_mag=sqrt(sum(r.^2,2)); %Re 

%magnetic field strength
Bz=(M/r_mag.^3); %T 
Bx=zeros(size(Bz)) ; %T
By=zeros(size(Bz)) ; %T
end

function dxdt = Lorentz(t,x)
global m;
global ms;
global q;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %Re
v=[x(4), x(5), x(6)]; %Re

%Call Magnetic field Function
[Bx,By,Bz]=smallphigrad(r); %T
B=[Bx,By,Bz]; %T

%DIFFERENTIAL EQUATION 
%Lorentz Equation
a_L=(q/m)*(E + cross(v,B)); %Re/s^2
a_Lx=a_L(1); %Re/s^2
a_Ly=a_L(2); %Re/s^2
a_Lz=a_L(3); %Re/s^2

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; %Output new positon and velcoity
end

function dxdt = LorentzDipole(t,x)
global m;
global q;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %Re
v=[x(4), x(5), x(6)]; %Re

%Call Magnetic field Function
[Bx,By,Bz]=dipole(r); %T
B=[Bx,By,Bz]; %T

%DIFFERENTIAL EQUATION 
%Lorentz Equation
a_L=(q/m)*(E + cross(v,B)); %Re/s^2
a_Lx=a_L(1); %Re/s^2
a_Ly=a_L(2); %Re/s^2
a_Lz=a_L(3); %Re/s^2

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; %Output new positon and velcoity
end


function [value,isterminal,direction] = eventsY1(t,y)
%Stop particle when hits y axis, in our case negative y axis 
value=y(1);
isterminal = 1;% 1 to stop the integration
direction = 0;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end

function [value,isterminal,direction] = eventsX(t1,y1)
%Stop particle when hits x axis, in our case negative x axis 
value=y1(2); 
isterminal = 1;% 1 to stop the integration
direction = 0;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end

function [value,isterminal,direction] = eventsY2(t2,y2)
%Stop particle when hits y axis, in our case positive y axis 
value=y2(1); 
isterminal = 1;% 1 to stop the integration
direction = 0;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end

function [value,isterminal,direction] = eventsDP(t3,y3)
%Stop particle when hits x axis, in our case positve x axis, completing drift period  
value=y3(2); %crossing x axis
isterminal = 1;% 1 to stop the integration
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end
