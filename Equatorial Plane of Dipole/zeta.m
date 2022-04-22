function zeta
clc

global m q;
global M n;
global r0 phi0 vr0 vphi0;
global A

m0=1;
q=1;
M=1000;
n=3;
c=3*10^8; %m/s

x0=7;
y0=0;
r0=sqrt(x0^2+y0^2);

vx0=0;
vy0=3;
v0=sqrt(vx0^2+vy0^2);

m=m0/sqrt(1-(v0^2/c^2));
A=q*M/m;

tmax=11;
optionsGP=odeset('Events',@eventsGP,'RelTol',1.0e-10,'AbsTol',1.0e-14);
[t1,y1]=ode45(@lorentz_force_Cartesian, [0 tmax], [x0 y0 vx0 vy0],optionsGP);
r=[y1(:,1), y1(:,2)]; 
v=[y1(:,3), y1(:,4)];
r_mag=sqrt(sum(r.^2,2)); 
v_mag=sqrt(sum(v.^2,2));

r_st=r0/(1-r0^2*vy0/A);
v_st=A/r_st^2;
r_max=2*r_st/(1+sqrt(1-4*v0/v_st));
r_min=2*r_st/(1+sqrt(1+4*v0/v_st));
p=(m*r0*vy0)-(q*M/r0);

ph_range=atan2(y1(:,2),y1(:,1));
max(ph_range);
min(ph_range);
ph=linspace(min(ph_range),max(ph_range));

start = linspace(r_min,r_max);

%figure 18
plot(y1(:,1),y1(:,2),'b',r_min*cos(ph),r_min*sin(ph),'k--',r_max*cos(ph),r_max*sin(ph),'k--')
hold on
plot(r_st*cos(ph),r_st*sin(ph),'r--')
grid on 
axis('equal')
xlabel('X')
ylabel('Y')
%axis([-0.3 0.9 0 0.9])
hold on

s=size(y1);
z=2100;
quiver(y1(z,1),y1(z,2), y1(z,3),y1(z,4), 0.75,'k');
quiver(0,0, y1(z,1),y1(z,2), 1,'k'); %r
quiver(0,0, r_min*cos(ph(100)),r_min*sin(ph(100)), 1,'k');% rmin
quiver(0,0, r_max*cos(ph(80)),r_max*sin(ph(80)), 1,'k');% rmax
quiver(0,0, r_st*cos(ph(90)),r_st*sin(ph(90)), 1,'k');% rstar

tet2=-0.1:0.1:2*pi;
plot(12+0.2*cos(tet2),-6+0.2*sin(tet2),'k' )
plot(12+0.00*cos(tet2),-6+0.00*sin(tet2),'k.' )
text(12.1, -6,'{\bf B}')
text(1.35, -6,'{\bf r}')
text(2.55, -5,'{\bf v}')
text(2, 0.75,'{r_{min}}')
text(6, 0.3,'{r_{*}}')
text(5.8, -0.75,'{r_{max}}')
text(1.6, -7,'{\zeta}')


%Calculate Canonical Angular Momentum P
v_phi=(y1(:,1).*y1(:,4)-y1(:,2).*y1(:,3))./r_mag;
P=m.*r_mag.*v_phi-(q.*M./r_mag);

%Plot Velocity 
figure;%figure 19
plot(t1/t1(end),(v_mag-v_mag(1))/v_mag(1));
xlabel('Time/Period');
ylabel('(v-v_{0})/v_{0}');
grid on;

%Plot Canonical Momentum
figure;%figure 19
plot(t1/t1(end),(P-P(1))/P(1));
xlabel('Time/Period');
ylabel('(P_{\phi}-P_{\phi0})/P_{\phi0}');
grid on;
%axis equal;


end

function dydt = lorentz_force_Cartesian( t,y )
global m q;
r=sqrt(y(1)^2+y(2)^2);
Bz=Bz_cyl(r);
dydt=[y(3);y(4); (q/m)*(y(4)*Bz); (q/m)*(-y(3)*Bz)];
end

function dydt = lorentz_force_Cyl( t,y )
global m q;
Bz=Bz_cyl(y(1));
dydt=[y(3);y(4)/y(1); y(4)^2/y(1)+(q/m)*(y(4)*Bz); -y(3)*y(4)/y(1)+(q/m)*(-y(3)*Bz)];
end

function Bz=Bz_cyl(r)
% returns the magnetic field (straight field lines, axisymmetric)
global QQ n;
Bz=QQ/r^n;
end

function [value,isterminal,direction] = eventsGP(t,y)
% Locate the time when height passes through zero in a decreasing direction
% and stop integration. 
tet=atan2(y(2),y(1));
vrad=y(3)*cos(tet)+y(4)*sin(tet);
value = vrad*sign(-10); %stop when radial velocity zero again 
isterminal = 0;% 1 to stop the integration
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both

end


