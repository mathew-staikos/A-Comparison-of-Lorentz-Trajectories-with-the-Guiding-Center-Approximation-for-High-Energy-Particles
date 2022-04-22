%Markers 
%Markers w different colours at rmin, rmax and intervals between
%only works if events on 
%get first rmin and rmax occurance in terms of row indexes 
[rmin, index_min]=mink(r_mag,1);
[rmax, index_max]=maxk(r_mag,1);
%get 2nd rmin or rmax occurance in terms of row indexes 
indexEnd=size(r_mag,1);
%get evenly distributed points between rmin and rmax in terms of row indexes
%total points want is 11, 3 for rmin/rmax and 8 for inbetween
%distrbute by r_mag values not row as ode45 rows not even distribution 

%Mu exact,Gyroperiod exact, and driftperiod exact 
if index_max>index_min
index_r_mag_mu1=find(r_mag>(r_mag_mu(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1);%unormalize with *(max(r_mag)-min(r_mag))+min(r_mag))
index_r_mag_gp1=find(r_mag>(r_mag_gyroperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1);
index_r_mag_dp1=find(r_mag>(r_mag_driftperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1);
index_r_mag_mu2=find(r_mag>(r_mag_mu(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1,'last');
index_r_mag_gp2=find(r_mag>(r_mag_gyroperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1,'last');
index_r_mag_dp2=find(r_mag>(r_mag_driftperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1,'last');
else 
index_r_mag_mu1=find(r_mag<(r_mag_mu(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1);
index_r_mag_gp1=find(r_mag<(r_mag_gyroperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1);
index_r_mag_dp1=find(r_mag<(r_mag_driftperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1);
index_r_mag_mu2=find(r_mag<(r_mag_mu(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1,'last');
index_r_mag_gp2=find(r_mag<(r_mag_gyroperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1,'last');
index_r_mag_dp2=find(r_mag<(r_mag_driftperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag)),1,'last');
end

%Make circles 
start_circle_mu=abs(atan2(y(index_r_mag_mu1,2),y(index_r_mag_mu1,1)));
end_circle_mu=abs(atan2(y(index_r_mag_mu2,2),y(index_r_mag_mu2,1)));
t_circle_mu=0.5*pi-start_circle_mu:0.001:0.5*pi+end_circle_mu;%0 is starting at y axis positive, x axis 0 
x_circle_mu=(r_mag_mu(end)*(max(r_mag)-min(r_mag))+min(r_mag))*sin(t_circle_mu); %unormalize with *(max(r_mag)-min(r_mag))+min(r_mag))
y_circle_mu=(r_mag_mu(end)*(max(r_mag)-min(r_mag))+min(r_mag))*cos(t_circle_mu); 
start_circle_gp=abs(atan2(y(index_r_mag_gp1,2),y(index_r_mag_gp1,1)));
end_circle_gp=abs(atan2(y(index_r_mag_gp2,2),y(index_r_mag_gp2,1)));
t_circle_gp=0.5*pi-start_circle_gp:0.001:0.5*pi+end_circle_gp;
x_circle_gp=(r_mag_gyroperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag))*sin(t_circle_gp); 
y_circle_gp=(r_mag_gyroperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag))*cos(t_circle_gp); 
start_circle_dp=abs(atan2(y(index_r_mag_dp1,2),y(index_r_mag_dp1,1)));
end_circle_dp=abs(atan2(y(index_r_mag_dp2,2),y(index_r_mag_dp2,1)));
t_circle_dp=0.5*pi-start_circle_dp:0.001:0.5*pi+end_circle_dp;
x_circle_dp=(r_mag_driftperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag))*sin(t_circle_dp); 
y_circle_dp=(r_mag_driftperiod(end)*(max(r_mag)-min(r_mag))+min(r_mag))*cos(t_circle_dp); 

figure;
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
colour=jet(11);
hold on;
plot(y(:,1),y(:,2),'b')
plot(r_grx,r_gry,'r')
plot(r_gc(:,1),r_gc(:,2),'k')
plot(y(index_r_mag_mu1,1),y(index_r_mag_mu1,2),'*','MarkerSize',5,'MarkerEdgeColor','b')
plot(y(index_r_mag_mu2,1),y(index_r_mag_mu2,2),'*','MarkerSize',5,'MarkerEdgeColor','b')
plot(y(index_r_mag_gp1,1),y(index_r_mag_gp1,2),'*','MarkerSize',5,'MarkerEdgeColor','r')
plot(y(index_r_mag_gp2,1),y(index_r_mag_gp2,2),'*','MarkerSize',5,'MarkerEdgeColor','r')
plot(y(index_r_mag_dp1,1),y(index_r_mag_dp1,2),'*','MarkerSize',5,'MarkerEdgeColor','g')
plot(y(index_r_mag_dp2,1),y(index_r_mag_dp2,2),'*','MarkerSize',5,'MarkerEdgeColor','g')
plot(x_circle_mu,y_circle_mu,'b')
plot(x_circle_gp,y_circle_gp,'r')
plot(x_circle_dp,y_circle_dp,'g')


