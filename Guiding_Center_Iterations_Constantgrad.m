%Guiding_Center_Iterations
%SECOND ITERATION
B2z=B0+Beta*r_gr(:,1);
B2x=zeros(size(Bz));
B2y=zeros(size(Bz));
B2=[B2x,B2y,B2z];
B2_mag=sqrt(sum(B2.^2,2));

r_gc2=r+(m./(q.*B2_mag)).*cross(v,B./B_mag);

%THIRD ITERATION
B3z=B0+Beta*r_gc2(:,1);
B3x=zeros(size(Bz));
B3y=zeros(size(Bz));
B3=[B3x,B3y,B3z];
B3_mag=sqrt(sum(B3.^2,2));

r_gc3=r+(m./(q.*B3_mag)).*cross(v,B./B_mag);

plot(r_gc2(:,1),r_gc2(:,2),'color',[1,0.4,0]) %2nd ITERATION 
plot(r_gc3(:,1),r_gc3(:,2),'color',[1,0.8,0]) %3rd ITERATION 
%[1,0,0] is red [0,1,0] is green and [0,0,1] is blue
