%MagFieldLines Script 
%Call ode45 for magentic field lines 
%options3=odeset('RelTol',1.0e-10,'AbsTol',1.0e-14);
options3=odeset('Events',@events_FieldLines,'RelTol',1.0e-10,'AbsTol',1.0e-14);
arc_length=[0, 200];%just pick high number so it stops at events not this 
Line1=[-8,0,0];
Line2=[-6,0,0];
Line3=[-4,0,0];
Line4=[-2,0,0];
Line5=[2,0,0];
Line6=[4,0,0];
Line7=[6,0,0];
Line8=[8,0,0];
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
%circle 
t_circle=0:0.001:2*pi;
x_circle=1*sin(t_circle); %1 here is radius
z_circle=1*cos(t_circle); %1 here is radius

%Plot Magnetic Field Lines
figure;
plot(w1(:,1),w1(:,3),'b'); hold on; grid on;
plot(w2(:,1),w2(:,3),'b');
plot(w3(:,1),w3(:,3),'b');
plot(w4(:,1),w4(:,3),'b');
plot(w5(:,1),w5(:,3),'b');
plot(w6(:,1),w6(:,3),'b');
plot(w7(:,1),w7(:,3),'b');
plot(w8(:,1),w8(:,3),'b');
plot(w9(:,1),w9(:,3),'b');
plot(w10(:,1),w10(:,3),'b');
plot(w11(:,1),w11(:,3),'b');
plot(w12(:,1),w12(:,3),'b');
plot(w13(:,1),w13(:,3),'b');
plot(w14(:,1),w14(:,3),'b');
plot(w15(:,1),w15(:,3),'b');
plot(w16(:,1),w16(:,3),'b');
plot(x_circle,z_circle,'k');
xlabel('X');
ylabel('Z');
axis equal; 

function dxds = FieldLines(s,x)
%s represents arc length when ode45 outputs is Bx/B_mag;By/B_mag;Bz/B_mag
%variables: x=x(1); y=x(2); z=x(3);
r=[x(1), x(2), x(3)]; %m  

%Call Magnetic field Function
[Bx,By,Bz]=mag_vec(r);
B=[Bx,By,Bz];
B_mag=sqrt(sum(B.^2,2));

dxds=[Bx/B_mag;By/B_mag;Bz/B_mag];  
end

function [Bx, By, Bz] = mag_vec(r)
global M;
%Tarsagov?(fix spelling)model 
SPS=0;
CPS=1;
P=r(:,1).*r(:,1); %m^2
T=r(:,2).*r(:,2); %m^2
U=r(:,3).*r(:,3); %m^2
V=3*r(:,3).*r(:,1); %m^2
Q=M./sqrt(P+T+U).^5; %Tm^2

Bx=Q.*((T+U-2.*P).*SPS-V.*CPS); %T
By=-3.*r(:,2).*Q.*(r(:,1).*SPS+r(:,3).*CPS); %T
Bz=Q.*((P+T-2.*U).*CPS-V.*SPS); %T
end

function [value,isterminal,direction] = events_FieldLines(t,y)
value = 1-sqrt(y(1).^2+y(3).^2); %stop when radius 1 around orgin, graph in x and z 
isterminal = 1;% 1 to stop the integration
direction = 0; %from either direction 
end