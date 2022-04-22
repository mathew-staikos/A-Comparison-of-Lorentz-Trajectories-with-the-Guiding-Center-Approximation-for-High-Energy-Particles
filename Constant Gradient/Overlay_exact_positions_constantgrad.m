%OVERLAY EXACT POSITIONS constantgrad 
%Markers 
%Markers w different colours at rmin, rmax and intervals between
%only works if events on 
%get first rmin and rmax occurance in terms of row indexes
[xmin, index_min]=mink(y(:,1),1)
[xmax, index_max]=maxk(y(:,1),1)
%get 2nd rmin or rmax occurance in terms of row indexes 
indexEnd=size(y(:,1),1)
%get evenly distributed points between rmin and rmax in terms of row indexes
%total points want is 11, 3 for rmin/rmax and 8 for inbetween
%distrbute by r_mag values not row as ode45 rows not even distribution 

%Mu exact,Gyroperiod exact, and driftperiod exact 
if index_max>index_min
index_x_mu1=find(y(:,1)>(x_mu(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1);%unormalize with *(max(x)-min(x))+min(x))
index_x_gp1=find(y(:,1)>(x_gyroperiod(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1);
index_x_Vd1=find(y(:,1)>(x_vdrift(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1);
index_x_mu2=find(y(:,1)>(x_mu(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1,'last');
index_x_gp2=find(y(:,1)>(x_gyroperiod(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1,'last');
index_x_Vd2=find(y(:,1)>(x_vdrift(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1,'last');
else 
index_x_mu1=find(y(:,1)<(x_mu(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1);
index_x_gp1=find(y(:,1)<(x_gyroperiod(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1);
index_x_Vd1=find(y(:,1)<(x_vdrift(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1);
index_x_mu2=find(y(:,1)<(x_mu(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1,'last');
index_x_gp2=find(y(:,1)<(x_gyroperiod(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1,'last');
index_x_Vd2=find(y(:,1)<(x_vdrift(end)*(max(y(:,1))-min(y(:,1)))+min(y(:,1))),1,'last');
end

figure;
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
colour=jet(11);
hold on
plot(y(:,1),y(:,2),'-','Color','b')
plot(r_grx,r_gry,'r')
plot(r_gc(:,1),r_gc(:,2),'k')
plot(y(index_x_mu1,1),y(index_x_mu1,2),'*','MarkerSize',5,'MarkerEdgeColor','b')
plot(y(index_x_mu2,1),y(index_x_mu2,2),'*','MarkerSize',5,'MarkerEdgeColor','b')
plot(y(index_x_gp1,1),y(index_x_gp1,2),'*','MarkerSize',5,'MarkerEdgeColor','r')
plot(y(index_x_gp2,1),y(index_x_gp2,2),'*','MarkerSize',5,'MarkerEdgeColor','r')
plot(y(index_x_Vd1,1),y(index_x_Vd1,2),'*','MarkerSize',5,'MarkerEdgeColor','g')
plot(y(index_x_Vd2,1),y(index_x_Vd2,2),'*','MarkerSize',5,'MarkerEdgeColor','g')
%Make Lines
plot([y(index_x_mu1,1), y(index_x_mu2,1)],[y(index_x_mu1,2), y(index_x_mu2,2)],'b')
plot([y(index_x_gp1,1), y(index_x_gp2,1)],[y(index_x_gp1,2), y(index_x_gp2,2)],'r')
plot([y(index_x_Vd1,1), y(index_x_Vd2,1)],[y(index_x_Vd1,2), y(index_x_Vd2,2)],'g')
%}