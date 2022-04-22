%Markers 
%Markers w different colours at rmin, rmax and intervals between
%only works if events on 
%get first rmin and rmax occurance in terms of row indexes
[xmin, index_min]=mink(y(:,1),1);
[xmax, index_max]=maxk(y(:,1),1);
%get 2nd rmin or rmax occurance in terms of row indexes 
indexEnd=size(y(:,1),1);
%get evenly distributed points between rmin and rmax in terms of row indexes
%total points want is 11, 3 for rmin/rmax and 8 for inbetween
%distrbute by r_mag values not row as ode45 rows not even distribution 

x_integer=(xmax-xmin)/5;
if index_max>index_min
index1=find(y(:,1)>xmin+x_integer*1,1);
index2=find(y(:,1)>xmin+x_integer*2,1);
index3=find(y(:,1)>xmin+x_integer*3,1);
index4=find(y(:,1)>xmin+x_integer*4,1);
index11=find(y(:,1)>xmin+x_integer*1,1,'last');
index22=find(y(:,1)>xmin+x_integer*2,1,'last');
index33=find(y(:,1)>xmin+x_integer*3,1,'last');
index44=find(y(:,1)>xmin+x_integer*4,1,'last');
else 
index1=find(y(:,1)<xmax-x_integer*1,1);
index2=find(y(:,1)<xmax-x_integer*2,1);
index3=find(y(:,1)<xmax-x_integer*3,1);
index4=find(y(:,1)<xmax-x_integer*4,1);
index11=find(y(:,1)<xmin+x_integer*1,1,'last');
index22=find(y(:,1)<xmin+x_integer*2,1,'last');
index33=find(y(:,1)<xmin+x_integer*3,1,'last');
index44=find(y(:,1)<xmin+x_integer*4,1,'last');
end

markers=[index_min;index_max;indexEnd;index1;index2;index3;index4;index11;index22;index33;index44];
Markers=sort(markers);
figure;
xlabel('X');
ylabel('Y');
grid on;
axis equal;
hold on;
colour=jet(11);
plot(y(:,1),y(:,2),'-','Color','b')
hold on
for i=1:11
plot(y(Markers(i),1),y(Markers(i),2),'*','MarkerSize',5,'MarkerEdgeColor',colour(i,:))
hold on 
end
plot(r_grx,r_gry,'-','Color','r')
hold on
for i=1:11
plot(r_grx(Markers(i)),r_gry(Markers(i)),'*','MarkerSize',5,'MarkerEdgeColor',colour(i,:))
hold on 
end
plot(r_gc(:,1),r_gc(:,2),'k')
hold on
for i=1:11
plot(r_gc(Markers(i),1),r_gc(Markers(i),2),'*','MarkerSize',5,'MarkerEdgeColor',colour(i,:))
hold on 
end
