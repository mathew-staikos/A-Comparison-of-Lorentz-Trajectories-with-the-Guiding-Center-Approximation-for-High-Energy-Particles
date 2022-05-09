%P and Mu all Cases 
%EXTRACTED FROM ALL VELOCITY FILES:
%velocities_SmallPhiGrad...
%velocities_LargePhiGrad...
%velocities_PiecewiseLin...
vrange_PiecewiseLinear = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5];
vrange_LargePhiGrad = [2 2.5 3 3.5 4 4.5 5];
vrange_SmallPhiGrad = [3.5 4 4.5 5 5.5 6 6.5];

AverageP_SmallPhiGrad = [-7.8714722388501237E-11 7.1564416368888425E-8 8.1852529405396747E-7 -3.541698078800592E-5 -0.00066718286013773729 -0.0020570071025261039 -0.0094056606905195076];
AverageMu_SmallPhiGrad = [-8.0599286371245833E-9 -2.2756534213009334E-7 -2.5604871911223993E-6 0.00011230193925542998 0.0021549281181754568 0.00691092333569264 0.035257363632096229];
AmplitudeP_SmallPhiGrad = [7.848950104437384E-10 8.3206569095415625E-7 3.9250364975385314E-5 0.0005986058431781241 0.0035771964735961869 0.012712934665934306 0.057128504642958204];
AmplitudeMu_SmallPhiGrad = [2.50006553983223E-9 2.57290537532003E-6 0.00012260263984453269 0.0018931686704975561 0.011490891021342714 0.041580082064949592 0.20715653156531033];
AverageP_LargePhiGrad = [-5.4915051849600712E-11 5.5724125088409071E-7 -2.5250784507797213E-5 0.00010692123870993447 -0.00044810813160166 -0.011495219191541896 -0.071763224900170189];
AverageMu_LargePhiGrad = [-2.738287379543717E-8 -1.7066969359147578E-6 7.6960007289829136E-5 -0.00032434016410978674 0.0014864817886361091 0.038843989404186909 0.3736864918661954];
AmplitudeP_LargePhiGrad = [2.9057342093124297E-9 5.7329482666601378E-6 0.00022029608226091426 0.002164469404420305 0.011426858930728541 0.051270186942394824 0.24484489068058904];
AmplitudeMu_LargePhiGrad = [9.1233221183606072E-9 1.7376574027472736E-5 0.00067118356699281262 0.0066360764279056244 0.035478424159733804 0.17032445100920615 1.2465067286174258];
AverageP_PiecewiseLinear = [5.9849670313911767E-5 5.27059625411964E-6 -7.6308813755468258E-5 0.00010383402185377382 0.0001330283769752262 0.0025602632252642045 0.0019099069634313972 -0.0029846118848782218 0.0070376226274608527 -0.016439188372654061 0.013713634773913085 -0.040452676058561593];
AverageMu_PiecewiseLinear = [-0.00017966108937730345 -1.58307394119046E-5 0.00023093434429490898 -0.00031284837431795268 -0.0003391698007212133 -0.0076228067599687054 -0.0053127327919030012 0.0097960712763967232 -0.021245535713385343 0.058301152645092169 -0.038461230699576535 0.17836627664619803];
AmplitudeP_PiecewiseLinear = [0.0004896368926141189 0.00037383296299185872 0.001172988140457778 0.0014503192940338832 0.010520481294563178 0.020659173682934424 0.03140945894353251 0.031658537079007204 0.028254489798567844 0.083473408287095208 0.10117340945954328 0.19006415078823216];
AmplitudeMu_PiecewiseLinear = [0.0014707308443465824 0.0011252393329076304 0.0035424427808471871 0.004396019535828256 0.03197122955185977 0.062735072246795073 0.095396097029223484 0.10150647442925584 0.087186438398155069 0.27348597444163003 0.35440200183139431 0.84340646299238053];

%Plot Canonical Momentum Average as function of V 
figure; %FIGURE 54
axes('XScale', 'log', 'YScale', 'log')
hold on
loglog(vrange_SmallPhiGrad,abs(AverageP_SmallPhiGrad),'b','Marker','*');
loglog(vrange_LargePhiGrad,abs(AverageP_LargePhiGrad),'Color','#A2142F','Marker','*');
loglog(vrange_PiecewiseLinear,abs(AverageP_PiecewiseLinear),'Color','#77AC30','Marker','*');
xlabel('Velocity (Re/s)');
xlim([min(vrange_PiecewiseLinear) max(vrange_PiecewiseLinear)])
ylabel('<P> & \DeltaP');
grid on;
%Plot Canonical Momentum Amplitude as function of V 
loglog(vrange_SmallPhiGrad,abs(AmplitudeP_SmallPhiGrad),'c','Marker','*');
loglog(vrange_LargePhiGrad,abs(AmplitudeP_LargePhiGrad),'r','Marker','*');
loglog(vrange_PiecewiseLinear,abs(AmplitudeP_PiecewiseLinear),'g','Marker','*');
%Scaling
loglog(vrange_PiecewiseLinear,abs(AmplitudeP_PiecewiseLinear(end))*(vrange_PiecewiseLinear/vrange_PiecewiseLinear(end)).^3,'k');
loglog(vrange_SmallPhiGrad,abs(AmplitudeP_SmallPhiGrad(end)*(exp(-960./vrange_SmallPhiGrad.^3))/(exp(-960./vrange_SmallPhiGrad(end).^3))),'k');
loglog(vrange_LargePhiGrad,abs(AmplitudeP_LargePhiGrad(end)*(exp(-80./vrange_LargePhiGrad.^2))/(exp(-80./vrange_LargePhiGrad(end).^2))),'k');
legend('<P_{small \nabla\phi}>','<P_{large \nabla\phi}>','<P_{piecewise}>','\DeltaP_{small \nabla\phi}','\DeltaP_{large \nabla\phi}','\DeltaP_{piecewise}')

%Plot Magnetic Moment Average as function of V 
figure; %FIGURE 55 
axes('XScale', 'log', 'YScale', 'log')
hold on
loglog(vrange_SmallPhiGrad,abs(AverageMu_SmallPhiGrad),'b','Marker','*');
loglog(vrange_LargePhiGrad,abs(AverageMu_LargePhiGrad),'Color','#A2142F','Marker','*');
loglog(vrange_PiecewiseLinear,abs(AverageMu_PiecewiseLinear),'Color','#77AC30','Marker','*');
xlabel('Velocity (Re/s)');
xlim([min(vrange_PiecewiseLinear) max(vrange_PiecewiseLinear)])
ylabel('<\mu> & \Delta\mu');
grid on;
%Plot Magnetic Moment Amplitude as function of V  
loglog(vrange_SmallPhiGrad,abs(AmplitudeMu_SmallPhiGrad),'c','Marker','*');
loglog(vrange_LargePhiGrad,abs(AmplitudeMu_LargePhiGrad),'r','Marker','*');
loglog(vrange_PiecewiseLinear,abs(AmplitudeMu_PiecewiseLinear),'g','Marker','*');
%Scaling
loglog(vrange_PiecewiseLinear,abs(AmplitudeMu_PiecewiseLinear(end))*(vrange_PiecewiseLinear/vrange_PiecewiseLinear(end)).^3,'k');
loglog(vrange_SmallPhiGrad,abs(AmplitudeMu_SmallPhiGrad(end)*(exp(-960./vrange_SmallPhiGrad.^3))/(exp(-960./vrange_SmallPhiGrad(end).^3))),'k');
loglog(vrange_LargePhiGrad,abs(AmplitudeMu_LargePhiGrad(end)*(exp(-80./vrange_LargePhiGrad.^2))/(exp(-80./vrange_LargePhiGrad(end).^2))),'k');
legend('<\mu_{small \nabla\phi}>','<\mu_{large \nabla\phi}>','<\mu_{piecewise}>','\Delta\mu_{small \nabla\phi}','\Delta\mu_{large \nabla\phi}','\Delta\mu_{piecewise}')

