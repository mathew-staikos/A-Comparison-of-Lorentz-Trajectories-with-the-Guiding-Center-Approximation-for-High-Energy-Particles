function Section_5dot2
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

%INTIAL CONSITIONS 
x0=7.0; %Re
y0=0.0; %Re
vx0=0.0; %Re
vy0=6; %Re/s

m=m0*(1/sqrt(1-((vy0*Re)^2/(3*10^5)^2)));%kg relativistic mass
A=q*M/m;
Energy=(m-m0)*c*c/(1.602E-19)/1.0E6; % Kinetic energy in MeV

%CALL ODE 45 to intgrate equation of motion 
r0=[x0, y0, 0];%Re
v0=[vx0,vy0, 0];%Re/s 
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
%Combine 
t1=t1+t(end);
t2=t2+t1(end);
t3=t3+t2(end);
t_total=[t;t1;t2;t3]; 
y_total=[y;y1;y2;y3]; 

%RESULTS
%Particle Motion Lorentz 
r=[y_total(:,1), y_total(:,2), y_total(:,3)];%Re
v=[y_total(:,4), y_total(:,5), y_total(:,6)];%Re/s
r_mag=sqrt(sum(r.^2,2));%Re 
v_mag=sqrt(sum(v.^2,2));%Re/s

%Call Magntic Field Function 
[Bx, By, Bz] = smallphigrad(r);
B=[Bx, By, Bz];
B_mag=sqrt(sum(B.^2,2));%Bz is B_mag 

%Radius of gyration/Larmor
r_gr=r+(m./(q.*B_mag.^2)).*cross(v,B);

%Plot Trajectory
figure; %FIGURE 39
xlabel('X (Re)');
ylabel('Y (Re)');
grid on;
axis equal;
hold on;
plot(y_total(:,1),y_total(:,2),'b') 
plot(r_gr(:,1),r_gr(:,2),'r') 

%CALL ODE 45 to intgrate equation of motion for unbound trajectory
r0=[x0, y0, 0];%Re
v0=[vx0,vy0, 0];%Re/s 
[tUB,yUB]=ode45(@Lorentz,[0,20],[[7,0,0],[0,7.98,0]],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14));

figure; %FIGURE 40
xlabel('X (Re)');
ylabel('Y (Re)');
grid on;
axis equal;
hold on;
plot(yUB(:,1),yUB(:,2),'b') 

%CALCULATE CANAOCNICAL ANGULAR MOMENTUM 
v_phi=(cross(r,v))./(r_mag); %Re/s
%Calculate Canonical Angular Momentum P
P=(m.*r_mag.*v_phi(:,3)-(q.*M./r_mag))*Re^2*1.0E6;%kgm2/s

%Plot Canonical Momentum
figure; %FIGURE 41
plot(t_total/t_total(end),P/P(1),'b');
xlabel('Time/Drift Period');
ylabel('P_{\phi}/P_{0}');
grid on;

%Calculate magnetic moment mu
mu=m.*(v_mag.^2./(Bz))*Re^2*1.0E6; %Am^2
%Magnetic moment at 1st gyrocenter
[B_grx, B_gry, B_grz] = smallphigrad(r_gr);%T
B_gr=[B_grx, B_gry, B_grz]; %T
mu_gr=m.*(v_mag.^2./(B_grz))*Re^2*1.0E6; %Am^2

%Exact analytical value of mu
k=(4*q*M*Re^3*1.0E9/m).*(v_mag*Re*1.0E3)./((P/m).^2); %unitless
[E_type1, E_type2] = ellipke(2.*k./(1+k)); %elliptic integral 
mu_exact=abs((q/pi)*P/m).*(pi-2.*sqrt(1+k).*E_type2);  %Am^2

%Plot Magnetic Moment
figure; %FIGURE 42 
xlabel('Time/Drift Period');
ylabel('\mu/\mu_{Exact}');
grid on;
hold on;
plot(t_total/t_total(end),mu/mu_exact(1),'b');
plot(t_total/t_total(end),mu_gr/mu_exact(1),'r');
plot(t_total/t_total(end),mu_exact/mu_exact(1), 'c');
%legend('\mu_{0}','Exact')
legend('Trajectory','1st Gyrocenter','Exact')

%PLOT TREJECTORIES IN LARGE PHI GRADIENT AND PIECEWISE LINEAR MAGENTIC FIELDS
%CALL ODE 45 FOR LARGE PHI GRADIENT 
[t,y]=ode45(@Lorentz2,[Tmin,Tmax],[r0, [0,5,0]],optionsY1);%INITIAL v 5Re/s not 6 Re/s
new_r0=[y(end,1), y(end,2), y(end,3)];
new_v0=[y(end,4), y(end,5), y(end,6)];
[t1,y1]=ode45(@Lorentz2,[Tmin,Tmax],[new_r0,new_v0],optionsX);
new_r0=[y1(end,1), y1(end,2), y1(end,3)];
new_v0=[y1(end,4), y1(end,5), y1(end,6)];
[t2,y2]=ode45(@Lorentz2,[Tmin,Tmax],[new_r0,new_v0],optionsY2);
new_r0=[y2(end,1), y2(end,2), y2(end,3)];
new_v0=[y2(end,4), y2(end,5), y2(end,6)];
[t3,y3]=ode45(@Lorentz2,[Tmin,Tmax],[new_r0,new_v0],optionsDP);
%Combine 
t1=t1+t(end);
t2=t2+t1(end);
t3=t3+t2(end);
t_total=[t;t1;t2;t3]; 
y_total=[y;y1;y2;y3]; 

%Plot Trajectory
figure; %FIGURE 43 LEFT PANEL 
xlabel('X (Re)');
ylabel('Y (Re)');
grid on;
axis equal;
hold on;
plot(y_total(:,1),y_total(:,2),'b') 

%CALL ODE 45 FOR SMALL PHI GRADIENT 
[t1,y1]=ode45(@Lorentz3,[Tmin,Tmax],[[0,-7,0],[6,0,0]],optionsX);%START ON NEGATIVE Y AXIS
new_r0=[y1(end,1), y1(end,2), y1(end,3)];
new_v0=[y1(end,4), y1(end,5), y1(end,6)];
[t2,y2]=ode45(@Lorentz3,[Tmin,Tmax],[new_r0,new_v0],optionsY2);
%Combine 
t2=t2+t1(end);
t_total=[t1;t2]; 
y_total=[y1;y2]; 

%Plot Trajectory
figure; %FIGURE 43 RIGHT PANEL 
xlabel('X (Re)');
ylabel('Y (Re)');
grid on;
axis equal;
hold on;
plot(y_total(:,1),y_total(:,2),'b') 

%COMPARE NUMERICAL INTGRATION THROUGH NON-AXISYMMETIC FIELD VS GC APPROXIAMTION THROUGH PURE DIPOPLE 
%Numerical Integration Time from Line 50 v0 set to 1 through 6 Re/s 0.5
%Re/s intevals
vrange2 = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5]; %Re/s
t_integration=[1.4017e+03,616.6957,342.8489,217.4550,148.8456,108.4525,81.9572,64.0200,51.6183,42.2683,34.6927,29.6307];%s
%GC approxiamtion
L=7.84; %rmax=rmin between rmin and rmax 
m=m0.*(1./sqrt(1-((vrange2.*6378).^2/(3*10^5)^2)));%relativistic mass at very velcoity
W=0.5*m.*vrange2.^2; %kgRe^2 Energy
alpha_eq=pi/2;%pitch angle
drift_period=(pi*q*M*1^2)./((3*L*W)*(0.35+0.15*sin(alpha_eq))); %s 

figure; %FIGURE 44
plot(vrange2, t_integration);
hold on
plot(vrange2, drift_period);
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Drift Period (s)');
legend('Numerical Integrator','Exact Drift Dipole')
grid on;
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

function dxdt = Lorentz2(t,x)
global m;
global q;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %Re
v=[x(4), x(5), x(6)]; %Re

%Call Magnetic field Function
[Bx,By,Bz]=largephigrad(r); %T
B=[Bx,By,Bz]; %T

%DIFFERENTIAL EQUATION 
%Lorentz Equation
a_L=(q/m)*(E + cross(v,B)); %Re/s^2
a_Lx=a_L(1); %Re/s^2
a_Ly=a_L(2); %Re/s^2
a_Lz=a_L(3); %Re/s^2

dxdt=[x(4);x(5);x(6);a_Lx;a_Ly;a_Lz]; %Output new positon and velcoity
end

function dxdt = Lorentz3(t,x)
global m;
global q;
global E;

%variables: x=x(1); y=x(2); z=x(3); vx=x(4); vy=x(5); vz=x(6)
%velocity equations: dx/dt=x(4); dy/dt=x(5); dz/dt=x(6);  
%acceleration equations: d(vx)/dt; d(vy)/dt=; d(vz)/dt=; see DIFFERENTIAL EQUATION
r=[x(1), x(2), x(3)]; %Re
v=[x(4), x(5), x(6)]; %Re

%Call Magnetic field Function
[Bx,By,Bz]=piecewiselinear(r); %T
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

function [Bx, By, Bz] = smallphigrad(r)
global M;

time=0; %min %No time variation in magentic field
r_mag=sqrt(sum(r.^2,2)); %Re 
phi_d=acosd(r(:,1)./r_mag); 

%f(phi)
c0=(1+tanh(11/4))^2;
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

function [Bx, By, Bz] = largephigrad(r)
global M;

time=0; %min %No time variation in magentic field
r_mag=sqrt(sum(r.^2,2)); %Re 
phi_d=acosd(r(:,1)./r_mag); 

%f(phi)
c0=(1+tanh(11))^2;
c1=(1+tanh(110/5))^2;
f_phi=(1/c0).*(1+tanh((55-(phi_d-180))/5)).*(1+tanh((55+(phi_d-180))/5));
f_phi=(f_phi/c1).*(1+tanh((110-(phi_d-180))/5)).*(1+tanh((110+(phi_d-180))/5));
%phi in terms of 0 to 360 not -180 to 180 to make f_phi work

%from reference paper
b1=648.9661-183.7327*tanh(1.6101*time-3.8109);
b2=0.2019-0.0859*tanh(1.6258*time-4.1687);
b3=3.3338+0.2906*tanh(1.6013*time-3.9748);
b4=1.3820-0.2966*tanh(1.8191*time-5.1973);
b5=93.1713+50.4476*tanh(1.6306*time-4.3549);
a1=(r_mag+b2)/b3;

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

function [Bx, By, Bz] = piecewiselinear(r)
global M;

time=0; %min %No time variation in magentic field
r_mag=sqrt(sum(r.^2,2)); %Re 
phi_d=acosd(r(:,1)./r_mag); 

%f(phi)
for i=1:361
if (-1<phi_d)&&(phi_d<=103) 
f_phi=0;
elseif (103<phi_d)&&(phi_d<146)  
f_phi=0.02369*phi_d-2.45216;
elseif (146<=phi_d)&&(phi_d<=214)
f_phi=1;
elseif (214<phi_d)&&(phi_d<256) 
f_phi=-0.02369*phi_d+6.07624;
elseif (256<=phi_d)&&(phi_d<361) 
f_phi=0;
end
end

%from reference paper
b1=648.9661-183.7327*tanh(1.6101*time-3.8109);
b2=0.2019-0.0859*tanh(1.6258*time-4.1687);
b3=3.3338+0.2906*tanh(1.6013*time-3.9748);
b4=1.3820-0.2966*tanh(1.8191*time-5.1973);
b5=93.1713+50.4476*tanh(1.6306*time-4.3549);
a1=(r_mag+b2)/b3;

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


