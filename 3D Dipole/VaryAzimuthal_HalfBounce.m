function VaryAzimuthal_HalfBounce
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

iterations=1000;
for i=1:iterations
%INITIAL CONDITIONS
%position
r0=[7, 0, 0];%Re
%Pitch Angle
Pangle=pi/4;  %Rad Pitch Angle 
%Azimuthal
Az=linspace(0,2*pi,iterations); %Rad
Aangle=Az(i); %Rad
%velocity
v0_mag = 5; %Re/s
v0=[v0_mag*sin(Pangle)*cos(Aangle), v0_mag*sin(Pangle)*sin(Aangle), v0_mag*cos(Pangle)];%Re/s
%Other
m=m0*(1/sqrt(1-((v0_mag*Re)^2/(3*10^5)^2)));%relativistic mass
Energy=(m-m0)*c*c/(1.602E-19)/1.0E6; % Kinetic energy in MeV

%Call ODE45 Functions
Tmax=1000; %Sufficiently high to reach stopping condition
options=odeset('Events',@events,'RelTol',1.0e-10,'AbsTol',1.0e-14); %Stop after half bounce
[t1,y1]=ode45(@Lorentz,[0,Tmax],[r0,v0],options); 

%RESULTS
r=[y1(:,1), y1(:,2), y1(:,3)];%Re
v=[y1(:,4), y1(:,5), y1(:,6)];%Re/s
r_mag=sqrt(sum(r.^2,2));%Re
v_mag=sqrt(sum(v.^2,2));%Re/s
size=[(length(r)),1]; %size of matrix
rho=[y1(:,1), y1(:,2), zeros(size)];%Re
rho_mag = sqrt(sum(rho.^2,2));%Re

%Call Magntic Field Function 
[Bx, By, Bz] = mag_vec(r); %T
B=[Bx, By, Bz]; %T
B_mag=sqrt(sum(B.^2,2)); %T 

%Pitch Angle Change 
alpha = acos((dot(v,B,2))./(v_mag.*B_mag));%Rad
AlphaI=alpha(1); %Rad
AlphaF=pi-alpha(end); %Rad Correcting for B direction
AlphaDiff(i)=(AlphaF-AlphaI); %Rad
AlplaDiffFIT(i)=0.19*sin(Aangle);%LEAST SQUARES FIT

%Magnetic Moment
v_perp = v_mag.*sin(alpha); %Re/s Perpendicular velocity at every position
mu=m.*((v_perp*Re*1.0E3).^2)./B_mag;  % Am^2
muindex(i)=(mu(end)-mu(1))/mu(1); 
muindexFit(i)=0.19+0.61*sin(Aangle);%least squares fit

%GC approximation of Magnetic Mirror Latitude
f= @(lam) ((cos(lam)).^6./(sqrt(1+3.*(sin(lam)).^2)))-(B_mag(1)/max(B_mag));
GC_lam_Num(i) = abs(fzero(f,0)); %Rad Numerical 
f2= @(lam) ((cos(lam)).^6./(sqrt(1+3.*(sin(lam)).^2)))-(sin(alpha(1))).^2;
GC_lam_An(i) = abs(fzero(f2,0)); %Rad Analytical 
f3= @(lam) ((cos(lam)).^6./(sqrt(1+3.*(sin(lam)).^2)))-(sin(alpha(1)+0.19*sin(Aangle))).^2;
GC_lam_AnMod(i) = abs(fzero(f3,0)); %Rad Modified Analytical 

%GC approximation of Bounce Period 
%Bounce Period 
Bounce_Num(i)=2*t1(end); %s Numerical 
L=r_mag(1)/1; %unitless shell pararamter 
W=Energy*1.602E-19*1.0E6; %Energy in J
Bounce_An(i)=L*Re*1.0E3*(3.6-1.6*sin(alpha(1)))/(sqrt(W/m)); %s Analytical
Bounce_AnMod(i)=L*Re*1.0E3*(3.6-1.6*sin(alpha(1)+0.19*sin(Aangle)))/(sqrt(W/m)); %s Modified Analytical 

%GC approximation of Drift Period 
%MEASURE THE CHANGE IN PHI 
x=[ones(size), zeros(size), zeros(size)];%Re %x-axis unit vector 
x_mag = sqrt(sum(x.^2,2)); %Re
phi = acos((dot(rho,x,2))./(rho_mag.*x_mag));%Rad
Drift_Num(i)=t1(end)*2*pi/phi(end); %s Numerical 
Drift_An(i)=(pi*q*M*(Re*1.0E3)^2)./(3*L*W*(0.35+0.15*sin(alpha(1)))); %s Analytical
Drift_AnMod(i)=(pi*q*M*(Re*1.0E3)^2)./(3*L*W*(0.35+0.15*sin(alpha(1)+0.19*sin(Aangle)))); %s Modified Analytical 
end

%Accuracy Comparison of GC approximation vs our modification
BounceACC_An = (sum(Bounce_Num)-sum(Bounce_An))/iterations;
BounceACC_AnMod  = (sum(Bounce_Num)-sum(Bounce_AnMod))/iterations;
DriftACC_An  = (sum(Drift_Num)-sum(Drift_An))/iterations;
DriftACC_AnMod  = (sum(Drift_Num)-sum(Drift_AnMod))/iterations;


%PLOTS
%Change in Pitch Angle as function of Azimuthal Angle 
figure; %FIGURE 81
hold on
xlim([0 2*pi])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(Az,AlphaDiff,'b');
xlabel('\xi (Rad)');
ylabel('\Delta\alpha_{eq} (Rad)');
grid on;
plot(Az,AlplaDiffFIT,'k')

%Change in Magnetic Moment as function of Azimuthal Angle 
figure; %FIGURE 83,FIGURE 82 LEFT+RIGHT PANEL CHANGE LINE 31 TO 0.1 AND 1 RESPECTIVELY 
hold on
xlim([0 2*pi])
xticks([0 pi/2 pi 3*pi/2 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(Az,muindex,'b'); 
plot(Az,muindexFit,'k');
xlabel('\xi (Rad)');
ylabel('\Delta\mu/\mu_{0}');
grid on;

%GC approximation of Magnetic Mirror Latitude as function of Azimuthal Angle 
figure;  %FIGURE 84
hold on
xlim([0 2*pi])
xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(Az,GC_lam_Num,'b');
plot(Az,GC_lam_An,'r');
plot(Az,GC_lam_AnMod,'g');
xlabel('\xi (Rad)');
ylabel('\lambda_{m} (Rad)');
grid on;
legend('Numerical', 'GC Theory','Modified GC Theory')

%GC approximation of Bounce period as function of Azimuthal Angle 
figure; % FIGURE 85
hold on
xlim([0 2*pi])
xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(Az,Bounce_Num,'b');
plot(Az,Bounce_An,'r')
plot(Az,Bounce_AnMod,'g')
xlabel('\xi (Rad)');
ylabel('Bounce (s)');
grid on;
legend('Numerical', 'GC Theory','Modified GC Theory')

%GC approximation of Drift period as function of Azimuthal Angle 
figure; %FIGURE 86
hold on
xlim([0 2*pi])
xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
plot(Az,Drift_Num,'b');
plot(Az,Drift_An,'r');
plot(Az,Drift_AnMod,'g');
xlabel('\xi (Rad)');
ylabel('Drift Period (s)');
grid on;
legend('Numerical', 'GC Theory','Modified GC Theory')

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
global q;
% Locate the time when height passes through zero in a decreasing direction

value = y(3); %stop when hit equitorial plane  
isterminal = 1;% 1 to stop the integration
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
end




