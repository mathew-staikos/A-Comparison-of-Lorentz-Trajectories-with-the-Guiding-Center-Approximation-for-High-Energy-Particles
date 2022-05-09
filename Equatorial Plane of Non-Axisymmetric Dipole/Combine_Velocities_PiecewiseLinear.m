 


a0max;
a1max; 
a2max; 
a3max;
%}

%P and Mu
vrange = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5];
%vrange2 = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5];

AverageP = [5.9849670313911767E-5 5.27059625411964E-6 -7.6308813755468258E-5 0.00010383402185377382 0.0001330283769752262 0.0025602632252642045 0.0019099069634313972 -0.0029846118848782218 0.0070376226274608527 -0.016439188372654061 0.013713634773913085 -0.040452676058561593];
AverageMu = [-0.00017966108937730345 -1.58307394119046E-5 0.00023093434429490898 -0.00031284837431795268 -0.0003391698007212133 -0.0076228067599687054 -0.0053127327919030012 0.0097960712763967232 -0.021245535713385343 0.058301152645092169 -0.038461230699576535 0.17836627664619803];
AmplitudeP = [0.0004896368926141189 0.00037383296299185872 0.001172988140457778 0.0014503192940338832 0.010520481294563178 0.020659173682934424 0.03140945894353251 0.031658537079007204 0.028254489798567844 0.083473408287095208 0.10117340945954328 0.19006415078823216];
AmplitudeMu = [0.0014707308443465824 0.0011252393329076304 0.0035424427808471871 0.004396019535828256 0.03197122955185977 0.062735072246795073 0.095396097029223484 0.10150647442925584 0.087186438398155069 0.27348597444163003 0.35440200183139431 0.84340646299238053];


%PLOTS VELOCITY
%Plot Canonical Momentum Average as function of V 
figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange,abs(AverageP),'b','Marker','*');
xlabel('Velocity (Re/s)');
xlim([min(vrange) max(vrange)])
ylabel('<P> & \DeltaP');
grid on;
%Plot Canonical Momentum Amplitude as function of V 
loglog(vrange,abs(AmplitudeP),'r','Marker','*');
%Scaling
loglog(vrange,abs(AverageP(1))*(vrange/vrange(1)).^30,'k');
loglog(vrange,abs(AmplitudeP(1))*(vrange/vrange(1)).^30,'k');
legend('<P>','\DeltaP','','')

%Plot Magnetic Moment Average as function of V 
figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange,abs(AverageMu),'b','Marker','*');
xlabel('Velocity (Re/s)');
xlim([min(vrange) max(vrange)])
ylabel('<\mu> & \Delta\mu');
grid on;
%Plot Magnetic Moment Amplitude as function of V 
loglog(vrange,abs(AmplitudeMu),'r','Marker','*');
%Scaling
loglog(vrange,abs(AverageMu(1))*(vrange/vrange(1)).^30,'k');
loglog(vrange,abs(AmplitudeMu(1))*(vrange/vrange(1)).^30,'k');
legend('<\mu>','\Delta\mu','','')

%p4=polyfit(log(vrange),log(abs(AmplitudeMu)),1);
%m4=p4(1);
%b4=exp(p4(2));
%shift=abs(AmplitudeMu(1))*20.^(vrange(1)/vrange(1))-abs(AmplitudeMu(1));
%loglog(vrange,(abs(AmplitudeMu(1))*100.^(vrange/vrange(1))),'k');
%loglog(vrange,(abs(AmplitudeMu(1))*20.^(vrange/vrange(1))),'k');
%loglog(vrange,abs(AmplitudeMu(1))*exp((vrange/vrange(1))),'k');
%loglog(vrange,abs(AmplitudeMu(1))*(vrange/vrange(1)).^1,'k');
%loglog(vrange2,abs(AmplitudeMu(7))*(vrange2/vrange2(1)).^20,'k');

%{
%FOURIER
%Plot Fourier Series Expansion Coefficients RMIN 
figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a0min));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a0');
grid on;

figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a1min));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a1');
grid on;
loglog(vrange2,abs(a1min(1))*(vrange2/vrange2(1)).^1,'k');

figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a2min));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a2');
grid on;

figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a3min));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a3');
grid on;

%Plot Fourier Series Expansion Coefficients RMAX
figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
loglog(vrange2, abs(a0max));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a0');
grid on;

figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a1max));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a1');
grid on;
loglog(vrange2,abs(a1max(1))*(vrange2/vrange2(1)).^1,'k');

figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a2max));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a2');
grid on;

figure;
axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
hold on
loglog(vrange2, abs(a3max));
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Fourier Expansion a3');
grid on;

%Integration Time
t_integration=[1.4017e+03,616.6957,342.8489,217.4550,148.8456,108.4525,81.9572,64.0200,51.6183,42.2683,34.6927,29.6307];
alpha_eq= pi/2; %90degrees
q=1.602*10^-19;
L=7; %L=r/cos^2(lambda) where r in Earth Radii and lamba geoomahnetic lattitude, but unitless
Re=1; %6378km, to get units to cancel this is just 1Re
Beq= 3.11.*10^-5; %T from Baumjohann
m0=1.6726231*10^-27; %kg 
m=m0.*(1./sqrt(1-((vrange2.*6378).^2/(3*10^5)^2)));%6378 to convert to km as speed of light in km/s
W=0.5*m.*vrange2.^2; %kgRe^2
drift_period=(pi*q*Beq*Re^2)./((3*L*W)*(0.35+0.15*sin(alpha_eq))); %units s 
figure;
%axes('XScale', 'log', 'YScale', 'log')
%axes('YScale', 'log') %just y scale
%hold on
plot(vrange2, t_integration);
hold on
plot(vrange2, drift_period);
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('Integration Time');
legend('Experimental','Analytical')
grid on;


%Net Drift 
iterations=length(angles);
iterations2=length(angles2);
iterations3=length(angles3);
Rmaxdrift=[sum(Rmax_diff_low3)/iterations3,sum(Rmax_diff2_low3)/iterations3,sum(Rmax_diff3_low3)/iterations3,sum(Rmax_diff_low2)/iterations2,sum(Rmax_diff2_low2)/iterations2,sum(Rmax_diff3_low2)/iterations2,sum(Rmax_diff_low)/iterations,sum(Rmax_diff2_low)/iterations,sum(Rmax_diff3_low)/iterations,sum(Rmax_diff_high)/iterations,sum(Rmax_diff2_high)/iterations,sum(Rmax_diff3_high)/iterations];
Rmindrift=[sum(Rmin_diff_low3)/iterations3,sum(Rmin_diff2_low3)/iterations3,sum(Rmin_diff3_low3)/iterations3,sum(Rmin_diff_low2)/iterations2,sum(Rmin_diff2_low2)/iterations2,sum(Rmin_diff3_low2)/iterations2,sum(Rmin_diff_low)/iterations,sum(Rmin_diff2_low)/iterations,sum(Rmin_diff3_low)/iterations,sum(Rmin_diff_high)/iterations,sum(Rmin_diff2_high)/iterations,sum(Rmin_diff3_high)/iterations];
figure;
plot(vrange2, Rmaxdrift)
hold on 
grid on 
plot(vrange2, Rmindrift)
xlabel('Velocity (Re/s)');
xlim([min(vrange2) max(vrange2)])
ylabel('<(\DeltaR_{min}/R_{0min})> & <(\DeltaR_{max}/R_{0max})>');
legend('RmaxDrift','RminDrift')
%}


