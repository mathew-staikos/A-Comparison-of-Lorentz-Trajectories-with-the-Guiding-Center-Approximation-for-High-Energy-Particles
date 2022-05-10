% -------------------------------------------------------------------
%  Generated by MATLAB on 9-Dec-2021 23:25:45
%  MATLAB version: 9.11.0.1769968 (R2021b)
% -------------------------------------------------------------------
%DATA USED IN FILE Section_5dot5, NO PLOTS HERE 
saveVarsMat = load('velocities_SmallPhiGrad_h_l2.mat');

Re = 6378;

E = 0;

m = 1.6741371631242128E-27;

q = 1.602E-19;

A = 2966.4236057768776;

M = 3.1E-5;

m0 = 1.6726231E-27;

c = 3.0E+8;

v_iter = 3;

vrange = [1 1.5 2];

i_v = 3;

x0 = 7;

y0 = 0;

vx0 = 0;

vy0 = 2;

Energy = 0.85059726079368581;

r0_mag = 7;

v0_mag = 2;

r_st = 7.2391557504332438;

v_st = 56.60533008987192;

r_max = 7.514779907636072;

r_min = 7;

p = -6.86019222573404E-25;

iterations = 50;

i = 50;

angles = [0 0.12820786341182827 0.25641572682365654 0.38462359023548487 ...
          0.51283145364731308 0.64103931705914141 0.76924718047096974 0.897455043882798 ...
          1.0256629072946262 1.1538707707064546 1.2820786341182828 1.410286497530111 ...
          1.5384943609419395 1.6667022243537679 1.7949100877655959 1.9231179511774241 ...
          2.0513258145892523 2.1795336780010808 2.3077415414129092 2.4359494048247372 ...
          2.5641572682365656 2.6923651316483941 2.8205729950602221 2.9487808584720505 ...
          3.076988721883879 3.2051965852957074 3.3334044487075358 3.4616123121193634 ...
          3.5898201755311918 3.7180280389430203 3.8462359023548482 3.9744437657666767 ...
          4.1026516291785047 4.2308594925903336 4.3590673560021616 4.48727521941399 ...
          4.6154830828258184 4.7436909462376464 4.8718988096494744 5.0001066730613033 ...
          5.1283145364731313 5.2565223998849593 5.3847302632967882 5.5129381267086162 ...
          5.6411459901204442 5.769353853532273 5.897561716944101 6.02576958035593 ...
          6.1539774437677579 6.2821853071795859];

vx1 = -0.0019999996666678413;

vy1 = 1.9999990000000833;

x1 = 7.0000001119296877;

y1 = saveVarsMat.y1; % <67781x6 double> too many elements

r1 = [7.0000001119296877 0 0];

v1 = [-0.0019999996666678413 1.9999990000000833 0];

Tmin = 0;

Tmax = 10000;

optionsY1 = saveVarsMat.optionsY1; % <1x1 struct> unsupported class at optionsY1.Events

optionsX = saveVarsMat.optionsX; % <1x1 struct> unsupported class at optionsX.Events

optionsY2 = saveVarsMat.optionsY2; % <1x1 struct> unsupported class at optionsY2.Events

optionsDP = saveVarsMat.optionsDP; % <1x1 struct> unsupported class at optionsDP.Events

t = saveVarsMat.t; % <124313x1 double> too many elements

y = saveVarsMat.y; % <124313x6 double> too many elements

new_r1 = [3.1061308039537217E-15 7.3073484231912458 0];

new_v1 = [0.99779643169416365 -1.7333211497161383 0];

t1 = saveVarsMat.t1; % <67781x1 double> too many elements

t2 = saveVarsMat.t2; % <64737x1 double> too many elements

y2 = saveVarsMat.y2; % <64737x6 double> too many elements

t3 = saveVarsMat.t3; % <126725x1 double> too many elements

y3 = saveVarsMat.y3; % <126725x6 double> too many elements

t_combine = saveVarsMat.t_combine; % <383556x1 double> too many elements

y_combine = saveVarsMat.y_combine; % <383556x6 double> too many elements

r = saveVarsMat.r; % <383556x3 double> too many elements

v = saveVarsMat.v; % <383556x3 double> too many elements

r_mag = saveVarsMat.r_mag; % <383556x1 double> too many elements

v_mag = saveVarsMat.v_mag; % <383556x1 double> too many elements

v_phi = saveVarsMat.v_phi; % <383556x3 double> too many elements

P = saveVarsMat.P; % <383556x1 double> too many elements

k = saveVarsMat.k; % <383556x1 double> too many elements

E_type1 = saveVarsMat.E_type1; % <383556x1 double> too many elements

E_type2 = saveVarsMat.E_type2; % <383556x1 double> too many elements

mu_exact = saveVarsMat.mu_exact; % <383556x1 double> too many elements

r_st_a = saveVarsMat.r_st_a; % <383556x3 double> too many elements

v_st_a = saveVarsMat.v_st_a; % <383556x3 double> too many elements

r_max_i = 7.247926305095409;

r_min_i = 6.7828691519103721;

r_max_f = 7.5147799040839534;

r_min_f = 7.0000000033605483;

P_diff = [-4.3457753301922959E-11 -4.3410851327148189E-11 -4.3419275468248558E-11 ...
          -4.3464469846854288E-11 -4.3425422814456842E-11 -4.3401744147580271E-11 ...
          -4.3395027602648922E-11 -4.3346531871449721E-11 -4.3439652782531814E-11 ...
          -4.3335147896989957E-11 -4.3346873390683751E-11 -4.337624404479015E-11 ...
          -4.3350857781744752E-11 -4.3349947063787774E-11 -4.3384668185890376E-11 ...
          -4.3331277345673273E-11 -4.3358029685654074E-11 -4.3399922711666656E-11 ...
          -4.3366453826754507E-11 -4.33937753654582E-11 -4.340948525021278E-11 ...
          -4.3410282128425054E-11 -4.3356891288208131E-11 -4.3432480878621911E-11 ...
          -4.34767645392709E-11 -4.3515470052434408E-11 -4.3472324789231489E-11 ...
          -4.3476195340547792E-11 -4.3518657565283329E-11 -4.3563168905421497E-11 ...
          -4.3533342892336525E-11 -4.3577171194007033E-11 -4.3511713340862513E-11 ...
          -4.3516722289625242E-11 -4.3557818437425318E-11 -4.3547572860410805E-11 ...
          -4.3555769322022876E-11 -4.35879859697468E-11 -4.3519112924262025E-11 ...
          -4.3554744764321464E-11 -4.3551898770706231E-11 -4.35719345657557E-11 ...
          -4.3570568488820615E-11 -4.3544043828328978E-11 -4.349771105227717E-11 ...
          -4.349998784716921E-11 -4.3524691071747096E-11 -4.3493271302237841E-11 ...
          -4.3468112718681434E-11 -4.3465266725066491E-11];

mu_exactDiff = [-1.1268376663076024E-7; -1.1262941813864298E-7; -1.1257512496259472E-7; ...
                -1.1251214167710913E-7; -1.1251258741282838E-7; -1.1246601965456633E-7; ...
                -1.124582554888254E-7; -1.1241748095877559E-7; -1.1240884696703614E-7; ...
                -1.1236753130088733E-7; -1.1236511890350011E-7; -1.1236029998713382E-7; ...
                -1.1235542628937257E-7; -1.1235494046945089E-7; -1.1236176292500943E-7; ...
                -1.1236316399946604E-7; -1.1235737277546112E-7; -1.1240638299449407E-7; ...
                -1.1240489493732796E-7; -1.1245486036159333E-7; -1.124980265789963E-7; ...
                -1.1249657272707878E-7; -1.1254985123968759E-7; -1.1260180416948699E-7; ...
                -1.1265518048741553E-7; -1.1270952069542065E-7; -1.1272111329819093E-7; ...
                -1.1277832794046218E-7; -1.1282976298307696E-7; -1.1287441231404357E-7; ...
                -1.1287292773089169E-7; -1.1292484912777309E-7; -1.1291702630567229E-7; ...
                -1.1297036320741401E-7; -1.1297137386249341E-7; -1.1296215183630935E-7; ...
                -1.1301355254028277E-7; -1.1301552401216329E-7; -1.1301206140999888E-7; ...
                -1.1296458307341167E-7; -1.1296700455675646E-7; -1.1296169140301643E-7; ...
                -1.1291369130460834E-7; -1.1291172531068165E-7; -1.1286513497173784E-7; ...
                -1.1285883427816457E-7; -1.1281037548349568E-7; -1.1275846677981804E-7; ...
                -1.126944431770386E-7; -1.1268377411313442E-7];

Rmin_diff = [0.016243160129227219 0.016111960825937726 0.015720314998477394 ...
             0.015074058493872826 0.014182849504936593 0.013060065454286323 ...
             0.011722657132284044 0.0101909587285734 0.00848845226086677 ...
             0.0066414849063163115 0.0046789380191530677 0.0026318471277125876 ...
             0.00053297308767918591 -0.0015836742332494646 -0.0036833567630761826 ...
             -0.0057311686016671111 -0.0076926287809105411 -0.009534296759794678 ...
             -0.011224398476558603 -0.012733448383305727 -0.014034850890029045 ...
             -0.015105463228849121 -0.015926101236927868 -0.016481970073659188 ...
             -0.016763003558319424 -0.016764098652987298 -0.016485235455013167 ...
             -0.015931477664788446 -0.015112853515502107 -0.01404412214995988 ...
             -0.01274443505248371 -0.011236905984763899 -0.0095481057227678953 ...
             -0.007707499560170231 -0.0057468461019596993 -0.0036995753175351124 ...
             -0.0016001624629295582 0.000516487544840014 0.0026156331093110329 ...
             0.00466325661239012 0.0066265855536529941 0.0084745692241974657 ...
             0.010178307933902805 0.011711433395117652 0.013050440086441222 ...
             0.014174968293994721 0.01506804004409896 0.015716249421071176 ...
             0.016109908769377443 0.016243152137501903];

Rmax_diff = [0.017382678061540143 0.017242352260888397 0.01682345533490974 ...
             0.016132207417417145 0.015178904139379147 0.013977808821965171 ...
             0.012546999764361789 0.010908171032936754 0.0090863849214770521 ...
             0.00710977426320992 0.0050091930008319225 0.0028178139849163427 ...
             0.000570673872738173 -0.0016958336782770629 -0.0039445136025464556 ...
             -0.0061379687386619438 -0.0082392353061568942 -0.010212444252893658 ...
             -0.012023495355515074 -0.013640728289108325 -0.015035572584429129 ...
             -0.016183156818957798 -0.01706285673162037 -0.017658762477154435 ...
             -0.017960047061995017 -0.017961221081860969 -0.017662263125550765 ...
             -0.017068620286472855 -0.016191078755393881 -0.015045510018028316 ...
             -0.013652503255018069 -0.012036898790461489 -0.010227240531452372 ...
             -0.0082551671792299528 -0.0061547625060136765 -0.0039618843837962827 ...
             -0.001713490639650684 0.000553022446323274 0.0028004558600212837 ...
             0.0049924074832800461 0.0070938280463569946 0.00907152833283381 ...
             0.010894634687541153 0.012534991632260193 0.013967511705696956 ...
             0.015170473547862461 0.016125769824672672 0.016819106809203926 ...
             0.017240157447932811 0.017382669513926467];

P_diff2 = [-4.7264577216474316E-11 -4.7218311204450394E-11 -4.7204649528443246E-11 ...
           -4.7204993940443432E-11 -4.7162746068421566E-11 -4.714609948841299E-11 ...
           -4.7105573676392013E-11 -4.7121761040400454E-11 -4.7140703700410237E-11 ...
           -4.709478210038652E-11 -4.70334767643547E-11 -4.7032558332353772E-11 ...
           -4.7034395196356362E-11 -4.7028884604352476E-11 -4.7053797072365233E-11 ...
           -4.7095700532386655E-11 -4.7136455952407861E-11 -4.7148395568414119E-11 ...
           -4.7055519132366012E-11 -4.713932605240945E-11 -4.7125549572402204E-11 ...
           -4.7158727928419465E-11 -4.716136842042085E-11 -4.7229561996456165E-11 ...
           -4.72187704204505E-11 -4.7317961076501867E-11 -4.7261247900472475E-11 ...
           -4.7313139308499451E-11 -4.736319385252533E-11 -4.7377659156532804E-11 ...
           -4.7360438556523865E-11 -4.740371966454619E-11 -4.737260778053032E-11 ...
           -4.7424384384556952E-11 -4.7423695560556846E-11 -4.7457218328574079E-11 ...
           -4.7500155024595605E-11 -4.7484541680589065E-11 -4.7457218328574305E-11 ...
           -4.7452396560571559E-11 -4.7441834592566E-11 -4.7369393268528578E-11 ...
           -4.7406245352547548E-11 -4.7367900816527822E-11 -4.7348843352517928E-11 ...
           -4.7363308656525422E-11 -4.7355846396521516E-11 -4.7325997356506019E-11 ...
           -4.7241042396462094E-11 -4.7280305364482389E-11];

mu_exactDiff2 = [-4.8702921206793533E-8; -4.8696221675829208E-8; -4.8641942701437792E-8; ...
                 -4.8637784842247707E-8; -4.858576455364238E-8; -4.8584292517963765E-8; ...
                 -4.8533313543737467E-8; -4.853291403415581E-8; -4.8531269293825674E-8; ...
                 -4.8530391898586106E-8; -4.8479807925626887E-8; -4.8479599038486619E-8; ...
                 -4.8480432506147357E-8; -4.8477516409749379E-8; -4.8478372650390462E-8; ...
                 -4.852852162490258E-8; -4.8527730582019078E-8; -4.8528158413337613E-8; ...
                 -4.8526194504591479E-8; -4.8573682278156144E-8; -4.8572213016859706E-8; ...
                 -4.8619858120449927E-8; -4.8620484550580928E-8; -4.8670867497172017E-8; ...
                 -4.8675846801731272E-8; -4.8728339658407356E-8; -4.873306499232855E-8; ...
                 -4.8784263257560408E-8; -4.8833796668243059E-8; -4.8834018849219554E-8; ...
                 -4.8832754197875564E-8; -4.8881241323428884E-8; -4.888058610922687E-8; ...
                 -4.8878974661322809E-8; -4.8928873017563695E-8; -4.8929530312595833E-8; ...
                 -4.892978116209489E-8; -4.8930389443198916E-8; -4.8930569777383428E-8; ...
                 -4.8928693839366848E-8; -4.8926604043470278E-8; -4.8876008973012293E-8; ...
                 -4.88747984219268E-8; -4.8874344349671071E-8; -4.8823178683313056E-8; ...
                 -4.8820906010056909E-8; -4.88186117198069E-8; -4.87649402173261E-8; ...
                 -4.8759449734682485E-8; -4.8702936928226E-8];

Rmin_diff2 = [0.024182697658463204 0.023988873395122991 0.023410145566283635 ...
              0.022454727659403272 0.02113623688380303 0.019473584166920323 ...
              0.017490815969592791 0.015216903901858867 0.012685477141052252 ...
              0.0099344919638246452 0.0070058322620101636 0.0039448349581463161 ...
              0.00079973483743461565 -0.0023789752948966971 -0.0055392709638479769 ...
              -0.0086284073932209986 -0.011593810018388837 -0.014384032849699678 ...
              -0.016949767209047527 -0.019244876750381293 -0.021227428320436118 ...
              -0.0228606827961768 -0.024114006482697695 -0.024963662630260761 ...
              -0.02539344486276994 -0.025395119868289229 -0.024968655513159185 ...
              -0.024122221731197935 -0.022871964158710355 -0.021241563031481576 ...
              -0.0192616011310302 -0.016968773353949869 -0.01440497610057611 ...
              -0.011616316892511809 -0.00865208368467172 -0.00556370975135339 ...
              -0.0024037646420121809 0.00077500465236653583 0.0039205647911583911 ...
              0.0069824080992239206 0.009912279590668864 0.012664817484645203 ...
              0.015198109019201317 0.017474165642565275 0.019459323061067034 ...
              0.021124572225666775 0.022445827373266211 0.023404136843689374 ...
              0.023985841651365483 0.024182685852785369];

Rmax_diff2 = [0.026812649464246472 0.02659801799594759 0.025957137754548387 ...
              0.024899029176143023 0.023438652092775507 0.02159679184422264 ...
              0.019399894827096344 0.0168798485700084 0.014073700293151492 ...
              0.011023306929401448 0.0077749090213599971 0.0043786207943626813 ...
              0.00088782931311352938 -0.0026415029303533447 -0.0061516350897417409 ...
              -0.0095839440377070775 -0.012879913699939568 -0.01598220561095615 ...
              -0.018835791320193227 -0.021389119474328529 -0.023595282748977073 ...
              -0.025413143283315922 -0.026808370819111846 -0.027754346407336063 ...
              -0.028232886906337389 -0.028234751989190954 -0.027759905585750894 ...
              -0.026817516912417138 -0.025425701014370517 -0.02361101358667746 ...
              -0.021407728031828596 -0.018856932987035613 -0.016005495045382883 ...
              -0.012904933828179409 -0.0096102552064129924 -0.0061787840736578607 ...
              -0.0026690315558511274 0.00086037611054177543 0.0043516876037946193 ...
              0.0077489233465787982 0.010998673359874395 0.014050795442684781 ...
              0.0168590167508721 0.019381444407699098 0.021580992231378392 ...
              0.023425731272393772 0.024889171769484204 0.025950483522012797 ...
              0.026594660756062923 0.026812636391344725];

P_diff3 = [-5.4294595448767064E-11 -5.4274332944681008E-11 -5.420671407390249E-11 ...
           -5.4206945645377823E-11 -5.41571577781949E-11 -5.4116401198547658E-11 ...
           -5.4187493641455206E-11 -5.4104243696096105E-11 -5.4106096267898077E-11 ...
           -5.4065455473988521E-11 -5.4037203754005846E-11 -5.4034888039252742E-11 ...
           -5.4056308400715725E-11 -5.40307197526986E-11 -5.4049245470719639E-11 ...
           -5.4074718332999421E-11 -5.4072402618246575E-11 -5.4101001695442389E-11 ...
           -5.4117327484448773E-11 -5.4116169627072467E-11 -5.4176030853429473E-11 ...
           -5.4185756855390627E-11 -5.4182051711786431E-11 -5.4276532873696171E-11 ...
           -5.4248744296663841E-11 -5.4329562741532685E-11 -5.4385255681334919E-11 ...
           -5.4340446600870371E-11 -5.4404707685257523E-11 -5.442786483278442E-11 ...
           -5.4463411054238249E-11 -5.44848314157006E-11 -5.4515514636173589E-11 ...
           -5.4478347414392986E-11 -5.45119252783071E-11 -5.4582786149739407E-11 ...
           -5.4578617863184421E-11 -5.4566923503681629E-11 -5.4556965930246525E-11 ...
           -5.4479736843244523E-11 -5.4477305342754143E-11 -5.4469200341119912E-11 ...
           -5.4522924923382258E-11 -5.4487031344715471E-11 -5.4398918398375679E-11 ...
           -5.4412117972466063E-11 -5.4428791118685477E-11 -5.4348551602504835E-11 ...
           -5.430860552302077E-11 -5.4298184806633757E-11];

mu_exactDiff3 = [-2.680844626621074E-8; -2.6809891280577273E-8; -2.6759446016282221E-8; ...
                 -2.6746147281536921E-8; -2.6749140363003366E-8; -2.6703681207517053E-8; ...
                 -2.6694961681309964E-8; -2.6690595659172045E-8; -2.6695105449361804E-8; ...
                 -2.6648519288885829E-8; -2.6647927652357954E-8; -2.664510780064657E-8; ...
                 -2.6643100358577681E-8; -2.6643187605679276E-8; -2.6645572865647477E-8; ...
                 -2.6647627092405029E-8; -2.6642778429403215E-8; -2.6686648043806737E-8; ...
                 -2.6685311898656247E-8; -2.6692748892423407E-8; -2.6733832036616565E-8; ...
                 -2.6730009854766588E-8; -2.6739022227779844E-8; -2.6777886134390386E-8; ...
                 -2.6790255497804584E-8; -2.6837332644519583E-8; -2.6839980530926311E-8; ...
                 -2.6846359558716911E-8; -2.68869964782291E-8; -2.68961483201615E-8; ...
                 -2.694092618956187E-8; -2.6933734752135246E-8; -2.693827349824813E-8; ...
                 -2.6942969794932661E-8; -2.6941876297888988E-8; -2.69890563709044E-8; ...
                 -2.69865727366553E-8; -2.6987688108696863E-8; -2.698677631323551E-8; ...
                 -2.6939364592835356E-8; -2.6940263870240654E-8; -2.6938790658624583E-8; ...
                 -2.693061737515903E-8; -2.692596040289407E-8; -2.6881919076073223E-8; ...
                 -2.688339734549257E-8; -2.6871693083200303E-8; -2.687240168142036E-8; ...
                 -2.6823900130880563E-8; -2.6808449806675988E-8];

Rmin_diff3 = [0.032011668186722651 0.031757030747942783 0.030996550062461697 ...
              0.02974050263788754 0.028005953743231952 0.025816660738457085 ...
              0.023202931082515244 0.0202014283013274 0.016854917279927138 ...
              0.013211938346652557 0.0093263979822417256 0.0052570626995430015 ...
              0.0010669420959515481 -0.0031774524855216718 -0.007406985302074948 ...
              -0.011550917340275907 -0.015538072836559303 -0.019298144067967218 ...
              -0.022763116053676733 -0.025868779518429783 -0.028556286241158021 ...
              -0.030773687260928939 -0.03247738349061105 -0.033633412309629379 ...
              -0.034218494560159234 -0.0342207752540061 -0.033640208084134031 ...
              -0.032488556931862866 -0.030789014278951931 -0.028575463110537106 ...
              -0.025891432063595066 -0.022788810424318098 -0.019326398485380893 ...
              -0.01556836974410279 -0.011582715441855598 -0.0074397310685918481 ...
              -0.0032105908389797318 0.0010339581270010819 0.0052247634401980537 ...
              0.0092952897087984036 0.013182496985514096 0.016827583119992849 ...
              0.020176601632262495 0.023180968637171461 0.025797872840601178 ...
              0.027990602056913411 0.029728798427748598 0.030988652890554757 ...
              0.031753047548391915 0.032011652677837968];

Rmax_diff3 = [0.036817941676116551 0.03652574477465699 0.035653029963924353 ...
              0.034211407465214796 0.0322201683026969 0.029706191181349761 ...
              0.026703803850830418 0.023254590582031286 0.019407134849226456 ...
              0.015216683871084545 0.010744719368415468 0.0060584170760614757 ...
              0.0012299764519368291 -0.0036641977548253742 -0.0085444784510724653 ...
              -0.013329175643949125 -0.017935877022661875 -0.022282963657243915 ...
              -0.026291280802958752 -0.029885927085485112 -0.032998107109511635 ...
              -0.035566974948687685 -0.037541381603994957 -0.03888143129569635 ...
              -0.039559751986078959 -0.039562396267614208 -0.038889309646387024 ...
              -0.037554332388620167 -0.035584734779406282 -0.03302031940041221 ...
              -0.029912153190146853 -0.026321013013497641 -0.022315639330738886 ...
              -0.017970893099583458 -0.013365902933208272 -0.0085822751909215081 ...
              -0.0037024220769858261 0.0011919553586822228 0.0060212092115920473 ...
              0.010708905572958534 0.01518280881147495 0.019375701235383829 ...
              0.023226054461810752 0.026678570850835262 0.029684613543871013 ...
              0.032202542588232824 0.03419797281257745 0.035643966810943388 ...
              0.036521173963189726 0.036817923879957491];

AverageP = [-4.3460701751308074E-11 -4.7255739604549735E-11 -5.42947413387965E-11 ...
            ];

AverageMu = [-1.1267159047507115E-7 -4.870845517079459E-8 -2.6808568010062559E-8 ...
             ];

min_P_diff = -5.4582786149739407E-11;

max_P_diff = -5.40307197526986E-11;

AmplitudeP = [2.5670862407352867E-13 4.7127042024312837E-13 5.5206639704080461E-13 ...
              ];

min_mu_exactDiff = -2.69890563709044E-8;

max_mu_exactDiff = -2.6642778429403215E-8;

AmplitudeMu = [6.60583542712397E-10 4.5305336763404851E-10 3.4627794150118608E-10 ...
               ];
%Plot Rmin and Rmax Difference w First Value vs Phase 
figure;
hold on
plot(angles,Rmin_diff,'b');
plot(angles,Rmin_diff2,'r');
plot(angles,Rmin_diff3,'g');
xlabel('\delta (Rad)');
xlim([0 2*pi])
xticks([0 pi/2 pi 3/2*pi 2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ylabel('\Deltar_{min}/r_{0min} & \Deltar_{max}/r_{0max}');
grid on;
plot(angles,Rmax_diff,'b','LineWidth',1.3)
plot(angles,Rmax_diff2,'r','LineWidth',1.3)
plot(angles,Rmax_diff3,'g','LineWidth',1.3)
legend('r_{min} v=1','r_{min} v=1.5','r_{min} v=2','r_{max} v=1','r_{max} v=1.5','r_{max} v=2')

clear saveVarsMat;
