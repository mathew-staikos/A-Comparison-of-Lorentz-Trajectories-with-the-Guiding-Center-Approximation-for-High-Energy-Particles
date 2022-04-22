function Plasma_Motion_Constantgrad_Varyvelocity
clc

%GLOBALS 
global Re;
global m;
global q;
global E;
global Beta;
global B0;

Re=6378;%km Equatorial Radius of Earth 
m0=1;%kg rest mass 
q=1; %C charge 
E=0; %V/m or N/C (mkg/s^3A) Electric Field, Assumed neglible
Beta=1; %T/m gradient magentic field factor 
B0=1; %T uniform magentic field factor 
c=3*10^8; %m/s  speed of light 

%INITIAL CONDITIONS 
r0=[10, 0, 0];%m 
for i=1:3000  
v0=[0, i*-0.01, 0];%m/s
v_abs(i,:)=sqrt(v0(1)^2+v0(2)^2+v0(3)^2); %used for ode45 guiding center 
Tmin=0;Tmax=5;

%Energy MeV
m=m0/sqrt(1-(v0(2)^2/c^2)); %kg realativistic mass  
energy=((m-m0)*c^2)/(1.602E-19*1.0E6); %MeV energy 

%Call ODE45 Functions
options=odeset('Events',@events,'RelTol',1.0e-10,'AbsTol',1.0e-14);
[t,y]=ode45(@Lorentz,[Tmin,Tmax],[r0,v0],options);

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
%Call del Magntic Field Function 
[delBx,delBy,delBz]=position_vec(r);
delB_prep=[delBx,delBy,delBz];
delB=repmat(delB_prep,size(Bz));

%Radius of gyration/Larmor/Cyclotron Equation
%r_g=m*v_perp/(|q|*B)
%use vectorized form to get correct direction
%r_gr=r+m*v_perp/(q*B)[B/B] 
%r_gr=r+m/(q*B^2)[v_perp x B] 
r_gr=r+(m./(q.*B_mag.^2)).*cross(v,B);
r_grx=r_gr(:,1);
r_gry=r_gr(:,2);
r_gr_mag=sqrt(sum(r_gr.^2,2));

%Calculate Canonical Angular Momentum P
Ay=B0.*y(:,1)+Beta.*(y(:,1).^2)./2; %dont need seperate function as dont need as vecotor
Py=m.*y(:,5)+q*Ay;
vstar=(B0^2+2*Beta*Py(1)/q)/(2*abs(Beta/q));

%Iterative gyrocenter approximation 
r_gc=r_gr;%feed intial condition 
for iterations=1:11
[B_gcx, B_gcy, B_gcz] = mag_vec(r_gc);
r_gc=r+(m./(q.*B_gcz)).*cross(v,B./B_mag);
end

%Calculate magnetic moment mu
%M=m*v_perpendicular^2/(B);% Am^2 magnetic moment
mu=m.*(v_mag.^2./(Bz));
%Magnetic moment of gyrocenter/radius
%create new Magntic Field Function using r_gr
[B_grx, B_gry, B_grz] = mag_vec(r_gr);
B_gr=[B_grx, B_gry, B_grz];
B_gr_mag=sqrt(sum(B_gr.^2,2));
mu_gr=m.*(v_mag.^2./(B_grz));
%Magnetic moment for guding center 
mu_gc=m.*(v_mag.^2./(B_gcz));

%Exact anayltical value of mu
v_star=q*B_mag.^2/(Beta*m);
u_x=y(:,4)./v_star;
u_y=y(:,5)./v_star;
u=sqrt(u_x.^2+u_y.^2);
k=(4*u)./(1+2*(u+u_y)); %unitless
[E_type1, E_type2] = ellipke(k);
mu_exact=(mu/m).*(2/(3*pi)).*(sqrt(1+2*(u+u_y))./u.^2).*((1+2*u_y).*E_type2-(1-2*(u-u_y)).*E_type1);
%Second Order Approximation for Adiabatic Invarient %eq A1 from 1.5054594
I=(mu/m).*(1-u_y+(3/8).*(u.^2+4*u_y.^2)); %using Lorentz trajecory position for B
v_star2=q*B_gcz.^2/(Beta*m);
u_x2=y(:,4)./v_star2;
u_y2=y(:,5)./v_star2;
u2=sqrt(u_x2.^2+u_y2.^2);
I2=(mu/m).*(1-u_y2+(3/8).*(u2.^2+4*u_y2.^2));%using gc position for B


%Magnetic Moment Ecaxt Value - Mean for each velocity 
MUexact_pos(i,:)=max(abs(mu_exact(1)-max(mu)),abs(mu_exact(1)-min(mu)));
MUexact_gr(i,:)=max(abs(mu_exact(1)-max(mu_gr)),abs(mu_exact(1)-min(mu_gr)));
MUexact_gc(i,:)=max(abs(mu_exact(1)-max(mu_gc)),abs(mu_exact(1)-min(mu_gc)));
MUexact_I(i,:)=max(abs(mu_exact(1)-max(I)),abs(mu_exact(1)-min(I)));
MUexact_I2(i,:)=max(abs(mu_exact(1)-max(I2)),abs(mu_exact(1)-min(I2)));


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
distance_pos(i,:)=sqrt((r(end,1)-r_trace(end,1))^2+(r(end,2)-r_trace(end,2))^2+(r(end,3)-r_trace(end,3))^2);
distance_gr(i,:)=sqrt((r_gr(end,1)-r_gr_trace(end,1))^2+(r_gr(end,2)-r_gr_trace(end,2))^2+(r_gr(end,3)-r_gr_trace(end,3))^2);
distance_gc(i,:)=sqrt((r_gc(end,1)-r_gc_trace(end,1))^2+(r_gc(end,2)-r_gc_trace(end,2))^2+(r_gc(end,3)-r_gc_trace(end,3))^2);


%Where to start gyrocenter approcimation? compare:
%@ what position is local Magnetic Moment same as exact Value
%anaylically use mu equation, isolate for x, feed mu exact 
%-min(x))/(max(x)-min(x)) is normalization for linear 
x_mu_prep=((m*v_mag(1)^2/(mu_exact(1)))-B0)/Beta;
x_mu(i,:)=(x_mu_prep-min(y(:,1)))/(max(y(:,1))-min(y(:,1))); %analyitally
%@ what position is local gyroperiod same as exact Value
%inverse of gyrofrequency which is freq=abs(q)*B/m, dont forget 2pi factor!
x_gyroperiod(i,:)=((((m*2*pi/(t(end)*abs(q)))-B0)/Beta)-min(y(:,1)))/(max(y(:,1))-min(y(:,1))); %analyitally
%@ what position is local drift velocity same as exact Value
%drift velocity  is (distance travelled - distance start)/t(end)
%Vdrift=(m/q)*(v_abs.^2/(2*B_mag.^2))*Beta, cant do thi way because Bmag is function of x
Vdrift=(abs(r_gc(end,2)-r_gc(1,2)))/(t(end));
x_vdrift_prep=(sqrt((m/q)*(v_mag(1)^2./(2*Vdrift))*Beta)-B0)/Beta;
x_vdrift(i,:)=((x_vdrift_prep)-min(y(:,1)))/(max(y(:,1))-min(y(:,1)));%analyitally
%}
%Tracing error xmin and to xmax
xmax=max(y(:,1));
xmin=min(y(:,1));
xnum=250;
xstart=linspace(xmin,xmax,xnum);
if v_abs(i)==0.1 
for k=1:xnum
options2=odeset('RelTol',1.0e-10,'AbsTol',1.0e-14);
[t4,z4]=ode45(@GuidingCenter,[Tmin,t(end)],[xstart(k),0,0],options2,v0(2));
distance_x(k)=sqrt((r(end,2)-z4(end,2))^2);
end
end
if v_abs(i)==1 
for k=1:xnum
options2=odeset('RelTol',1.0e-10,'AbsTol',1.0e-14);
[t4,z4]=ode45(@GuidingCenter,[Tmin,t(end)],[xstart(k),0,0],options2,v0(2));
distance_x(k)=sqrt((r(end,2)-z4(end,2))^2);
end
end
%Trajecory Error Correlation GC approxiamtion 
%figure 11
if v_abs(i)==0.1
set(gca, 'YScale', 'log')
xlabel('(x_{0} - x_{min})/(x_{max} - x_{min})');
ylabel('Tracing Error');
grid on;
hold on
xlim([0,1])
semilogy((xstart-xmin)/(xmax-xmin),distance_x,'b');
semilogy((xstart-xmin)/(xmax-xmin),distance_x*10^3,'k');
hold on
end
if v_abs(i)==1
semilogy((xstart-xmin)/(xmax-xmin),distance_x,'r');
hold on
end
%end of function
end

%Plot Velocity vs Distance Error(GC Tracing)
figure; %figure 10
axes('XScale', 'log', 'YScale', 'log')
xlabel('Velocity');
ylabel('Tracing Error');
grid on;
hold on;
loglog(v_abs,distance_pos,'b');
loglog(v_abs,distance_gr,'r');
loglog(v_abs,distance_gc,'k');
xlim([0.01,30])
%loglog(v_abs,distance_pos(1)*(v_abs/v_abs(1)).^3,'m');
%loglog(v_abs,distance_gc(1)*(v_abs/v_abs(1)).^4,'m');


%Magnteic moment error
figure;%figure 14
axes('XScale', 'log', 'YScale', 'log')
xlabel('Velocity');
ylabel('Max(\mu-\mu_{Exact})');
grid on;
hold on;
loglog(v_abs,MUexact_pos,'b');
loglog(v_abs,MUexact_gr,'r');
loglog(v_abs,MUexact_gc,'k');
loglog(v_abs,MUexact_I, 'm');
loglog(v_abs,MUexact_I2, 'g');
legend('Trajectory','1st Gyrocenter','Converged Gyrocenter','2nd Order Taylor','2nd Order Taylor @ gc');
xlim([0.01,30])
%loglog(v_abs,MUexact_pos(1)*(v_abs/v_abs(1)).^3,'y');
%loglog(v_abs,MUexact_I(1)*(v_abs/v_abs(1)).^5,'y');
%loglog(v_abs,MUexact_gc(1)*(v_abs/v_abs(1)).^4,'y');

%Plot Velocity vs Position Exact Value Achieved
figure; %figure 15
xlabel('Velocity');
ylabel('(x_{\mu}x_{GP}x_{VD} - x_{min})/(x_{max} - x_{min})');
grid on;
hold on;
plot(v_abs,x_mu,'b');
plot(v_abs,x_gyroperiod,'r');
plot(v_abs,x_vdrift,'g');
legend('x_{\mu}','x_{GP}','x_{VD}');

end

function [Bx, By, Bz] = mag_vec(r)
global Beta;
global B0;
%From draft0_02
Bz=B0+Beta*r(:,1); %T
n=size(Bz);
Bx=zeros(n);
By=zeros(n);
end

function [delBx, delBy, delBz] = position_vec(r)
global Beta;
delBx=Beta;%T %delBx becasue partial der of B not Bz and B=,Bx,By,Bz
n=size(delBx);
delBy=zeros(n);%T
delBz=zeros(n);%T
end

function dxdt = Lorentz(t,x)
global m;
global q;
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
a_L=(q/m)*(E + cross(v,B));
a_Lx=a_L(1);
a_Ly=a_L(2);
a_Lz=a_L(3);

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
value = y(4); %stop when x velocity zero again 
isterminal = 1;% 1 to stop the integration
if q>0 %direction -1 for ion, 1 for electon 
direction = -1;%-1 to approach value from positve direction(decreasing),1 negative, 0 both
else
direction = 1;
end
end




