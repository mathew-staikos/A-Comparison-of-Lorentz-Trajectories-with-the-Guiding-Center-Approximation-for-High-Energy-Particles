%Guiding_Center_Iterations
B2z=M./(r_gr_mag.^3);%dipole
B2x=zeros(size(Bz));
B2y=zeros(size(Bz));
B2=[B2x,B2y,B2z];
B2_mag=sqrt(sum(B2.^2,2));

r_gc2=r+(m./(q.*B2_mag)).*cross(v,B./B_mag);
r_gc_mag=sqrt(sum(r_gc2.^2,2));

B3z=M./(r_gc_mag.^3);%dipole
B3x=zeros(size(Bz));
B3y=zeros(size(Bz));
B3=[B3x,B3y,B3z];
B3_mag=sqrt(sum(B3.^2,2));

r_gc3=r+(m./(q.*B3_mag)).*cross(v,B./B_mag);

plot(r_gc2(:,1),r_gc2(:,2),'color',[1,0.4,0])
plot(r_gc3(:,1),r_gc3(:,2),'color',[1,0.8,0])
%[1,0,0] is red [0,1,0] is green and [0,0,1] is blue. Mix is get inbetween
