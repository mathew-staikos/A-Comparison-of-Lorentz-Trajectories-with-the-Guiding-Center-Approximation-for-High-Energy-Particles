function Magnetic_Fields
clc
%GLOBALS 
global Re;
global E;
global q;
global M; 

%INITIAL CONDITONS
Re=6378;%km Equatorial Radius of Earth
M=31000.*10^-9; %TRe^3 Magnetic Field Strength of Earth 
E=0; %V/m or N/C (mkg/s^3A) Electric Field

%PREP
phi=0:1:360;
%Plot f(phi)
%20 case 
%f(phi) small phi gradient
c0=(1+tanh(11/4))^2;
c1=(1+tanh(110/5))^2;
f_phi=(1/c0).*(1+tanh((55-(phi-180))/20)).*(1+tanh((55+(phi-180))/20));
f_phi=(f_phi/c1).*(1+tanh((110-(phi-180))/5)).*(1+tanh((110+(phi-180))/5));
%f(phi) large phi gradient
c0_L=(1+tanh(11))^2;
f_phi_L=(1/c0_L).*(1+tanh((55-(phi-180))/5)).*(1+tanh((55+(phi-180))/5));
f_phi_L=(f_phi_L/c1).*(1+tanh((110-(phi-180))/5)).*(1+tanh((110+(phi-180))/5));

%f(phi) functions
figure;%FIGURE 33
hold on
plot(phi,f_phi,'b');
plot(phi,f_phi_L,'r');
xlim([0 360])
xlabel('\phi (Deg)');
ylabel('f(\phi)');
grid on;
%piecewise linear f(phi) based off of figure 37
xx=0:1:360;
for i=1:361
if (-1<xx(i))&&(xx(i)<=103) 
yy(i)=0;
elseif (103<xx(i))&&(xx(i)<146) 
yy(i)=0.02369*xx(i)-2.45216;
elseif (146<=xx(i))&&(xx(i)<=214)
yy(i)=1;
elseif (214<xx(i))&&(xx(i)<256) 
yy(i)=-0.02369*xx(i)+6.07624;
elseif (256<=xx(i))&&(xx(i)<361) 
yy(i)=0;
end
end
plot(xx,yy,'g')
legend('small \nabla\phi','large \nabla\phi','piecewise linear')

%Contour Plot Magnetic Field Lines for small phi gradient
figure; % FIGURE 34 
x = linspace(-10,10,1000);
y = linspace(-10,10,1000);
[X,Y] = meshgrid(x,y); %X and Y are matrix's of rows=length(y) and columns=length(x)
for k1=1:length(x) %call magentic field for every position in meshgrid
for k2=1:length(y)
[Bxtmp, Bytmp, BZZ(k2,k1)] = smallphigrad([x(k1),y(k2),0]); 
end
end
[C,h]=contour(X,Y,BZZ,50,'LabelSpacing',300);
xlabel('X (Re)') 
ylabel('Y (Re)')
grid on;
axis equal;
hold on
%%plot earth with dayside and nightside radius 1 Re
t_circle=0:0.001:2*pi;
x_circle=sin(t_circle);
z_circle=cos(t_circle); 
%shad half
t_semi_circle=pi:0.001:2*pi;
x_semi_circle=sin(t_semi_circle); 
z_semi_circle=cos(t_semi_circle); 
plot(x_circle,z_circle,'k');
patch(x_semi_circle,z_semi_circle,'k')
hold on
%clabel(C,h,'manual')%Use to determine magentic field strength to be
%displayed 
clabel(C,h,[2.6e-08,4.8e-08,7.0e-08])%,4.8332e-08,7.0465e-08]) %display desired magentic field strength on contour lines 

%Plot MAGNTIC FIELD STRENGTH AS FUNCTION OF RADIUS AT VARIOUS 'SLICES'
%OF PHI
x = linspace(3,15); %RANGE OF INTEREST IN Re 
for i1=1:length(x)
[Bx, By, Bz1(i1)] = smallphigrad([x(i1),0,0]);
[Bx, By, Bz2(i1)] = smallphigrad([-x(i1),0,0]);
phi=2.5;
[Bx, By, Bz3(i1)] = smallphigrad([x(i1)*cos(phi),x(i1)*sin(phi),0]);
end
figure; %FIGURE 35
semilogy(x,Bz1,'k',x,Bz2,'r',x,Bz3,'g')
legend('\phi = 0','\phi = \pi','\phi = 2.2')
xlabel('r (Re)');
ylabel('B_{z} (T)');
grid on;

%Contour Plot Magnetic Field Lines for large phi gradient
figure; % FIGURE 36
x = linspace(-10,10,1000);
y = linspace(-10,10,1000);
[X,Y] = meshgrid(x,y); %X and Y are matrix's of rows=length(y) and columns=length(x)
for k1=1:length(x) 
for k2=1:length(y)
[Bxtmp, Bytmp, BZZ(k2,k1)] = largephigrad([x(k1),y(k2),0]); 
end
end
[C,h]=contour(X,Y,BZZ,50,'LabelSpacing',325);
xlabel('X (Re)') 
ylabel('Y (Re)')
grid on;
axis equal;
hold on
%%plot earth with dayside and nightside radius 1 Re
plot(x_circle,z_circle,'k');
patch(x_semi_circle,z_semi_circle,'k')
hold on
%clabel(C,h,'manual')%Use to determine magentic field strength to be
%displayed 
clabel(C,h,[2.6e-08,4.8e-08,7.0e-08])%,4.8332e-08,7.0465e-08]) %display desired magentic field strength on contour lines 

%DERIVATIVE OF f(phi) for SMALL AND LARGE PHI GRADIENTS
syms phii
c0=(1+tanh(11/4))^2;
c1=(1+tanh(110/5))^2;
f_phi=(1/c0).*(1+tanh((55-(phii-180))/20)).*(1+tanh((55+(phii-180))/20));
f=(f_phi/c1).*(1+tanh((110-(phii-180))/5)).*(1+tanh((110+(phii-180))/5));
ff = diff(f);
c02=(1+tanh(11))^2;
f_phi2=(1/c02).*(1+tanh((55-(phii-180))/5)).*(1+tanh((55+(phii-180))/5));
f2=(f_phi2/c1).*(1+tanh((110-(phii-180))/5)).*(1+tanh((110+(phii-180))/5));
ff2 = diff(f2);

figure; %FIGURE 37
hold on
grid on 
xlim([0 360])
xlabel('\phi (Deg)');
ylabel('df(\phi)/d\phi');
fplot(ff,[0 360],'b');
fplot(ff2,[0 360],'r');
legend('small \nabla\phi','large \nabla\phi')

%Contour Plot Magnetic Field Lines for piecewise linear
figure; % FIGURE 38
x = linspace(-10,10,1000);
y = linspace(-10,10,1000);
[X,Y] = meshgrid(x,y); %X and Y are matrix's of rows=length(y) and columns=length(x)
for k1=1:length(x) 
for k2=1:length(y)
[Bxtmp, Bytmp, BZZ(k2,k1)] = piecewiselinear([x(k1),y(k2),0]); 
end
end
[C,h]=contour(X,Y,BZZ,50,'LabelSpacing',300);
xlabel('X (Re)') 
ylabel('Y (Re)')
grid on;
axis equal;
hold on
%%plot earth with dayside and nightside radius 1 Re
plot(x_circle,z_circle,'k');
patch(x_semi_circle,z_semi_circle,'k')
hold on
%clabel(C,h,'manual')%Use to determine magentic field strength to be
%displayed 
clabel(C,h,[2.6e-08,4.8e-08,7.0e-08])%,4.8332e-08,7.0465e-08]) %display desired magentic field strength on contour lines

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


