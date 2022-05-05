%Magnetic Field Lines 3D Dipole

global Re;
global M;
Re= 6378; %km Radius of Earth @ Equator 
M= 31000.*10^-9; %TRe^3 Magnteic Field Strength of Earth 

%Inital Point of Integration for Field Lines 
%8 points at equal radial distances, r=10Re, evenly distributed around Earth 
Line1=[-10,0,0]; %Re
Line2=[10,0,0]; %Re 
Line3=[0,-10,0]; %Re 
Line4=[0,10,0]; %Re 
Line5=[-7.071,-7.071,0]; %Re 
Line6=[-7.071,7.071,0]; %Re 
Line7=[7.071,-7.071,0]; %Re  
Line8=[7.071,7.071,0]; %Re 

%Call ODE45 for 16 field lines
options3=odeset('Events',@events_FieldLines,'RelTol',1.0e-10,'AbsTol',1.0e-14); %Stopping Condition
arc_length=[0, 100]; % Sufficiently High to Reach Stopping Condition of r = 1 Re
[s1,w1]=ode45(@FieldLines,arc_length,Line1,options3);
[s2,w2]=ode45(@FieldLines,arc_length,Line2,options3);
[s3,w3]=ode45(@FieldLines,arc_length,Line3,options3);
[s4,w4]=ode45(@FieldLines,arc_length,Line4,options3);
[s5,w5]=ode45(@FieldLines,arc_length,Line5,options3);
[s6,w6]=ode45(@FieldLines,arc_length,Line6,options3);
[s7,w7]=ode45(@FieldLines,arc_length,Line7,options3);
[s8,w8]=ode45(@FieldLines,arc_length,Line8,options3);
%feed arc length end of each and add - to get opposite direction
[s9,w9]=ode45(@FieldLines,[0,-s1(end)],Line1,options3);
[s10,w10]=ode45(@FieldLines,[0,-s2(end)],Line2,options3);
[s11,w11]=ode45(@FieldLines,[0,-s3(end)],Line3,options3);
[s12,w12]=ode45(@FieldLines,[0,-s4(end)],Line4,options3);
[s13,w13]=ode45(@FieldLines,[0,-s5(end)],Line5,options3);
[s14,w14]=ode45(@FieldLines,[0,-s6(end)],Line6,options3);
[s15,w15]=ode45(@FieldLines,[0,-s7(end)],Line7,options3);
[s16,w16]=ode45(@FieldLines,[0,-s8(end)],Line8,options3);
%Generate Sphere
[x,y,z]=sphere;

%Plot Magnetic Field Lines
figure; %Figure 59 
surf(x,y,z) % Plot sphere
xlabel('X (Re)');
ylabel('Y (Re)');
zlabel('Z (Re)');
hold on;
grid on;
plot3(w1(:,1),w1(:,2),w1(:,3),'r'); 
plot3(w2(:,1),w2(:,2),w2(:,3),'r');
plot3(w3(:,1),w3(:,2),w3(:,3),'r');
plot3(w4(:,1),w4(:,2),w4(:,3),'r');
plot3(w5(:,1),w5(:,2),w5(:,3),'r');
plot3(w6(:,1),w6(:,2),w6(:,3),'r');
plot3(w7(:,1),w7(:,2),w7(:,3),'r');
plot3(w8(:,1),w8(:,2),w8(:,3),'r');
plot3(w9(:,1),w9(:,2),w9(:,3),'b');
plot3(w10(:,1),w10(:,2),w10(:,3),'b');
plot3(w11(:,1),w11(:,2),w11(:,3),'b');
plot3(w12(:,1),w12(:,2),w12(:,3),'b');
plot3(w13(:,1),w13(:,2),w13(:,3),'b');
plot3(w14(:,1),w14(:,2),w14(:,3),'b');
plot3(w15(:,1),w15(:,2),w15(:,3),'b');
plot3(w16(:,1),w16(:,2),w16(:,3),'b');
axis equal; 

function dxds = FieldLines(s,x)
%s represents arc length when ode45 outputs are Bx/B_mag;By/B_mag;Bz/B_mag
%variables: x=x(1); y=x(2); z=x(3);
r=[x(1), x(2), x(3)]; %Re  

%Call Magnetic field Function
[Bx,By,Bz]=mag_vec(r);%T
B=[Bx,By,Bz];%T
B_mag=sqrt(sum(B.^2,2));%T

dxds=[Bx/B_mag;By/B_mag;Bz/B_mag];  %Re 
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

function [value,isterminal,direction] = events_FieldLines(t,y)
value = 1-sqrt(y(1).^2+y(2).^2+y(3).^2); %stop when radius 1 Re around orgin
isterminal = 1;% 1 to stop the integration
direction = 0; %from either direction 
end

