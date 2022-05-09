function Section_5dot4_start
clc
%GLOBALS 
global Re;
global E;
global m;
global q;
global M; 
global A;

Re=6378;%km Equatorial Radius of Earth
m0=1.6726231*10^-27; %kg rest mass of proton 
q=1.602*10^-19; %C q_p=1.602*10^-19;%C(As) Charge of proton
M=31000.*10^-9; %TRe^3 Magnetic Field Strength of Earth 
E=0; %V/m or N/C (mkg/s^3A) Electric Field
c=3.0E8; %speed of light m/s

v_iter=3; 
vrange=linspace(5.5,6.5,v_iter);%CHNAGE VALUES HERE SEE LINES 
for i_v=1:v_iter
%INITIAL CONDITIONS 
x0=7.0; %Re
y0=0.0; %Re
vx0=0.0; %Re/s
vy0=vrange(i_v);%Re/s

m=m0*(1/sqrt(1-((vy0*Re)^2/(3*10^5)^2)));%kg relativistic mass 
A=q*M/m;
Energy=(m-m0)*c*c/(1.602E-19)/1.0E6; % Kinetic energy in MeV

%CALL ODE 45 FOR FIRST PHASE TO GET INTEGRATION TIME OF 1 DP 
r0=[x0, y0, 0];%Re
v0=[vx0,vy0, 0];%Re/s 
Tmin=0;Tmax=10000; %sufficiently high to reach stopping conditions first
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

t_integrate=t_total(end);

%Prep for Phase shift using Canaocial angular momentum
r0_mag=sqrt(x0^2+y0^2); %Re
v0_mag=sqrt(vx0^2+vy0^2); %Re/s
r_st=r0_mag/(1-r0_mag^2*vy0/A);  %r star
v_st=A/r_st^2; %v star
r_max=2*r_st/(1+sqrt(1-4*v0_mag/v_st));%maximum bounds of particle 
r_min=2*r_st/(1+sqrt(1+4*v0_mag/v_st));%minimum bounds of particle 
p=m*r0_mag*vy0-q*M/r0_mag; %kgRe/s canonical angular momentum  

iterations=1000; %number of phases 
for i=1:iterations 
angles=linspace(0,2*pi,iterations);
%INITIAL CONDITIONS of each phase
vx1=v0_mag*sin(angles(i)); %Re/s
vy1=v0_mag*cos(angles(i)); %Re/s
x1=(p+sqrt(p^2+4*m*vy1*q*M))/(2*m*vy1);%Re
y1=0.0; %Re

%CALL ODE 45 
[t_phase,y_phase]=ode45(@Lorentz,[0,t_integrate],[[x1,y1,0], [vx1,vy1,0]],odeset('RelTol',1.0e-10,'AbsTol',1.0e-14));

%RESULTS
r=[y_phase(:,1), y_phase(:,2), y_phase(:,3)];%Re
v=[y_phase(:,4), y_phase(:,5), y_phase(:,6)];%Re/s
r_mag=sqrt(sum(r.^2,2));%Re 
v_mag=sqrt(sum(v.^2,2));%Re/s

%Calculate Canonical Angular Momentum P
v_phi=(cross(r,v))./(r_mag); 
P=(m.*r_mag.*v_phi(:,3)-(q.*M./r_mag))*Re^2*1.0E6; %kgm^2/s

%Exact anayltical value of mu
k=(4*q*M*Re^3*1.0E9/m).*(v_mag*Re*1.0E3)./((P/m).^2); %unitless
[E_type1, E_type2] = ellipke(2.*k./(1+k));
mu_exact=abs((q/pi)*P/m).*(pi-2.*sqrt(1+k).*E_type2); %Am^2

%Rmin and Rmax Difference, Last Value vs First, anayltically @ Different Phases 
r_st_a=r_mag./(1-r_mag.^2.*v_phi./A); %new r star for calculating rmin/rmax after DP
v_st_a=A./r_st_a.^2; %new v star for calculating rmin/rmax after DP
r_max_i=2*r_st_a(1)/(1+sqrt(1-4*v_mag(1)/v_st_a(1)));
r_min_i=2*r_st_a(1)/(1+sqrt(1+4*v_mag(1)/v_st_a(1)));
r_max_f=2*r_st_a(end)/(1+sqrt(1-4*v_mag(end)/v_st_a(end)));
r_min_f=2*r_st_a(end)/(1+sqrt(1+4*v_mag(end)/v_st_a(end)));
    
    %index normalized difference in P Mu and RMIN/RMAX as function of phase 
    if i_v==1
    P_diff(i)=(P(end)-P(1))/P(1); %index normalized difference in P 
    mu_exactDiff(i,:)=(mu_exact(end)-mu_exact(1))/mu_exact(1); %index normalized difference in Mu
    Rmin_diff(i)=(r_min_f - r_min_i)/r_min_i; %index normalized difference in RMIN
    Rmax_diff(i)=(r_max_f - r_max_i)/r_max_i; %index normalized difference in RMAX
    end
    if i_v==2
    P_diff2(i)=(P(end)-P(1))/P(1);
    mu_exactDiff2(i,:)=(mu_exact(end)-mu_exact(1))/mu_exact(1);
    Rmin_diff2(i)=(r_min_f - r_min_i)/r_min_i;
    Rmax_diff2(i)=(r_max_f - r_max_i)/r_max_i;
    end
    if i_v==3
    P_diff3(i)=(P(end)-P(1))/P(1);
    mu_exactDiff3(i,:)=(mu_exact(end)-mu_exact(1))/mu_exact(1);
    Rmin_diff3(i)=(r_min_f - r_min_i)/r_min_i;
    Rmax_diff3(i)=(r_max_f - r_max_i)/r_max_i;
    end    

end
    
    %index average and max amplitude across phases for P & Mu as function of velocity
    if i_v==1
    %Average P_diff and mu_exactDiff to see if net drift 
    AverageP(i_v)=sum(P_diff)/iterations;
    AverageMu(i_v)=sum(mu_exactDiff)/iterations;
    %Amplitude of P_diff and mu_exactDiff to see scaling with V 
    min_P_diff = min(P_diff, [], 'all');
    max_P_diff = max(P_diff, [], 'all');
    AmplitudeP(i_v)=max_P_diff-min_P_diff;
    min_mu_exactDiff = min(mu_exactDiff, [], 'all');
    max_mu_exactDiff = max(mu_exactDiff, [], 'all');
    AmplitudeMu(i_v)=max_mu_exactDiff-min_mu_exactDiff;
    end
    if i_v==2
    AverageP(i_v)=sum(P_diff2)/iterations;
    AverageMu(i_v)=sum(mu_exactDiff2)/iterations;
    min_P_diff = min(P_diff2, [], 'all');
    max_P_diff = max(P_diff2, [], 'all');
    AmplitudeP(i_v)=max_P_diff-min_P_diff;
    min_mu_exactDiff = min(mu_exactDiff2, [], 'all');
    max_mu_exactDiff = max(mu_exactDiff2, [], 'all');
    AmplitudeMu(i_v)=max_mu_exactDiff-min_mu_exactDiff;
    end
    if i_v==3
    AverageP(i_v)=sum(P_diff3)/iterations;
    AverageMu(i_v)=sum(mu_exactDiff3)/iterations;
    min_P_diff = min(P_diff3, [], 'all');
    max_P_diff = max(P_diff3, [], 'all');
    AmplitudeP(i_v)=max_P_diff-min_P_diff;
    min_mu_exactDiff = min(mu_exactDiff3, [], 'all');
    max_mu_exactDiff = max(mu_exactDiff3, [], 'all');
    AmplitudeMu(i_v)=max_mu_exactDiff-min_mu_exactDiff;
    end 
end

%FOR NO MAGNETOTAIL-COMMENT OUT LINE 229 TO REMOVE TAIL AND MAKE MAGENTIC
%FIELD PURE DIPOLE 
figure; %FIGURE 53
subplot(2,1,1);%rows columns, position
plot(angles,P_diff2,'k');
xticks([0 pi/2 pi 3/2*pi 2*pi]);
xlim([0 2*pi]);
set(gca,'xticklabel',[]);
ylabel('\DeltaP/(P_0)');
grid on;
subplot(2,1,2);
plot(angles,mu_exactDiff2,'k');
xlabel('\delta (Rad)');
xlim([0 2*pi]);
xticks([0 pi/2 pi 3/2*pi 2*pi]);
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
ylabel('\Delta\mu_{Exact}/(\mu_{0Exact})');
grid on;

%''''''''''''''''''SAVE DATA BY REMOVING % BEFORE matlab.io....., VELOCITES IN SETS OF 3''''''''''''''''''''''''''''''''
%SMALL PHI GRADIENT MAGNETIC FIELD 
%SET LINE 19 = vrange=linspace(5.5,6.5,v_iter);             
%matlab.io.saveVariablesToScript('velocities_velocities_SmallPhiGrad_h')
%SET LINE 19 = vrange=linspace(4,5,v_iter);
%matlab.io.saveVariablesToScript('velocities_SmallPhiGrad_m')
%SET LINE 19 = vrange=linspace(2.5,3.5,v_iter);
%matlab.io.saveVariablesToScript('velocities_SmallPhiGrad_l')
%SET LINE 19 = vrange=linspace(1,2,v_iter);
%matlab.io.saveVariablesToScript('velocities_SmallPhiGrad_l2')

%LARGE PHI GRADIENT MAGNETIC FIELD -CHANGE LINE 316 TO [Bx,By,Bz]=largephigrad(r);
%SET LINE 19 = vrange=linspace(4,5,v_iter);
%matlab.io.saveVariablesToScript('velocities_LargePhiGrad_h')
%SET LINE 19 = vrange=linspace(2.5,3.5,v_iter);
%matlab.io.saveVariablesToScript('velocities_LargePhiGrad_m')
%SET LINE 19 = vrange=linspace(1,2,v_iter);
%matlab.io.saveVariablesToScript('velocities_LargePhiGrad_l')

%PIECEWISE LINEAR MAGNETIC FIELD -CHANGE LINE 316 TO [Bx,By,Bz]=piecewiselinear(r);
%SET LINE 19 = vrange=linspace(5.5,6.5,v_iter);
%matlab.io.saveVariablesToScript('velocities_PiecewiseLin_h')
%SET LINE 19 = vrange=linspace(4,5,v_iter);
%matlab.io.saveVariablesToScript('velocities_PiecewiseLin_m')
%SET LINE 19 = vrange=linspace(2.5,3.5,v_iter);
%matlab.io.saveVariablesToScript('velocities_PiecewiseLin_l')
%SET LINE 19 = vrange=linspace(1,2,v_iter);
%matlab.io.saveVariablesToScript('velocities_PiecewiseLin_l2')
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
[Bx,By,Bz]=smallphigrad(r); %T %CHANGE FOR DIFFERENT MAGENTIC FIELD
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

