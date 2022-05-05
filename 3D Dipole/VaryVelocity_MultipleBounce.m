function VaryVelocity_MultipleBounce
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
r0_mag=sqrt(sum(r0.^2,2));%Re
rho0_mag=sqrt(r0(1)^2 + r0(2)^2);%Re
%velocity
Pangle=pi/4; %Rad Pitch Angle 
vv=linspace(0.01,10,iterations); 
v0_mag=vv(i); %Re/s
v0=[0, v0_mag*sin(Pangle), v0_mag*cos(Pangle)];%Re/s
%Other
m=m0*(1/sqrt(1-((v0_mag*Re)^2/(3*10^5)^2)));%relativistic mass
Energy=(m-m0)*c*c/(1.602E-19)/1.0E6; % Kinetic energy in MeV

%Call ODE45 Functions
Tmax=1300/v0_mag; 
[t1,y1]=ode45(@Lorentz,[0,Tmax],[r0,v0],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14));

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

%Calculate Pitch Angle 
alpha = acos((dot(v,B,2))./(v_mag.*B_mag));%Rad

%Equitorial crossing when z=0, Below determines Pitch angle at each
%crossing adjusted for magnetic field direction 
indices = (find(diff(sign(y1(:,3))))+1);
crossing=1:1:length(indices);
alphacross=[];%Initilization to prevent problems in loops 
for k=1:length(indices)
    if alpha(indices(k))<pi/2
    alphacross(k)=alpha(indices(k)); %Rad
    else 
    alphacross(k)=pi-alpha(indices(k)); %Rad
    end
end

%Pitch Angle at every equitorial crossing 
    if v0_mag == 1
    figure;%FIGURES 67 
    plot(crossing, alphacross,'b')
    xlabel('Crossing #');
    ylabel('\alpha (Rad)');
    xlim([1 max(crossing)])
    end
    if v0_mag == 3
    figure;%FIGURE 68 TOP LEFT PANEL(CHANGE LINE 36 TO Tmax=2600/v0_mag; FOR TOP RIGHT PANEL)
    plot(crossing, alphacross,'b')
    xlabel('Crossing #');
    ylabel('\alpha (Rad)');
    grid on;
    xlim([1 max(crossing)])
    end
    if v0_mag == 5 
    figure;%FIGURE 68 MIDDLE LEFT PANEL(CHANGE LINE 36 TO Tmax=2600/v0_mag; FOR MIDDLE RIGHT PANEL)
    plot(crossing, alphacross,'b')
    xlabel('Crossing #');
    ylabel('\alpha (Rad)');
    grid on;
    xlim([1 max(crossing)])
    end
    if v0_mag == 7
    figure;%FIGURE 68 BOTTOM LEFT PANEL(CHANGE LINE 36 TO Tmax=2600/v0_mag; FOR BOTTOM RIGHT PANEL)
    plot(crossing, alphacross,'b')
    xlabel('Crossing #');
    ylabel('\alpha (Rad)');
    grid on;
    xlim([1 max(crossing)])
    end
%Maximum range of alpha
alpharange(i) = (max(alphacross) - min(alphacross));   
end

%Maximum Range of Alpha over multiple bounces  
figure;%FIGURE 69 CHANGE LINE 36 TO Tmax=332/v0_mag TO REDUCE COMPUTATIONAL TIME
plot(vv,alpharange,'b');
hold on
xlabel('Velocity (Re/s)');
ylabel('(\alpha_{eq max}-\alpha_{eq min})');
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




