% -------------------------------------------------------------------
%  Generated by MATLAB on 1-Mar-2021 20:56:18
%  MATLAB version: 9.8.0.1451342 (R2020a) Update 5
% -------------------------------------------------------------------
                                 

tmp = [-8.8004110099528364E-5; -0.00021158681564218752];

t = 1;

X = [0 0.09801714032956077 0.19509032201612833 0.29028467725446233 0.38268343236508984 ...
     0.47139673682599781 0.55557023301960229 0.63439328416364549 0.70710678118654757 ...
     0.773010453362737 0.83146961230254524 0.88192126434835494 0.92387953251128674 ...
     0.95694033573220894 0.98078528040323043 0.99518472667219693 1 1 0.99518472667219682 ...
     0.98078528040323043 0.95694033573220894 0.92387953251128674 0.88192126434835494 ...
     0.83146961230254524 0.77301045336273688 0.70710678118654746 0.63439328416364549 ...
     0.55557023301960218 0.47139673682599775 0.38268343236508978 0.29028467725446228 ...
     0.19509032201612825 0.098017140329560715 0];

LUT = ...
  [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 0;
   0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 -1];

L0 = [0.0981751911023479; 0.09817519110234757; 0.0981751911023476; 0.09817519110234782; ...
      0.09817519110234782; 0.098175191102347542; 0.098175191102347736; ...
      0.09817519110234782; 0.098175191102347667; 0.098175191102347709; ...
      0.098175191102347639; 0.09817519110234782; 0.09817519110234775; 0.098175191102347653; ...
      0.09817519110234757; 0.098174618883555623];

R0 = [0.99880200643017625; 0.98918299710753721; 0.970037614780322; 0.94155023994627152; ...
      0.90399522159782142; 0.8577342350913284; 0.80321279901567755; 0.74095598460472611; ...
      0.67156335901428832; 0.59570321116266778; 0.51410611574305076; 0.42755789738983779; ...
      0.33689206275783895; 0.24298177339754734; 0.14673143673208888; 0.049008570164780357; ...
      ];

K = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10];

MU = [10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10; 10];

ext_verts = [1 17];

ext_force_status = 0;

P = 1;

N = 17;

rv = ...
  [0 0.09801714032956077 0.19509032201612833 0.29028467725446233 0.38268343236508984 ...
   0.47139673682599781 0.55557023301960229 0.63439328416364549 0.70710678118654757 ...
   0.773010453362737 0.83146961230254524 0.88192126434835494 0.92387953251128674 ...
   0.95694033573220894 0.98078528040323043 0.99518472667219693 1;
   1 0.99518472667219682 0.98078528040323043 0.95694033573220894 0.92387953251128674 ...
   0.88192126434835494 0.83146961230254524 0.77301045336273688 0.70710678118654746 ...
   0.63439328416364549 0.55557023301960218 0.47139673682599775 0.38268343236508978 ...
   0.29028467725446228 0.19509032201612825 0.098017140329560715 0];

cut = 6;

slope = -0.53451113595079158;

horizA = [-0.024553270825918756 -0.076982730606830613 -0.14067083423431992 ...
          -0.22556930854325458 -0.35297921715735775];

horizB = [-0.49880008229204753 -0.48338312514180143 -0.44610832720595112 ...
          -0.37217414856345504 -0.22590116512615277];

horizC = [0 -0.0015111260524505038 -0.0087830783848416583 -0.0302450375701524 ...
          -0.0862212849342211];

horizD = [1 1.0000493720847796 1.0005222679255168 1.0025989605573025 1.0097393547147022 ...
          ];

vertA = [11.773768309524469 5.1498680454640464 2.6718854357819986 1.4773688651722798 ...
         0.88090996383498876 0.54804994409593943 0.35123522927847045 0.22600956282620668 ...
         0.14055772729944982 0.07701273001089444 0.024543399325741413];

vertB = [26.543470512444607 10.020755158977467 4.2742357773715742 1.7402834754180803 ...
         0.6051149115540686 0.050333555366069538 -0.22799988760681311 -0.37176525118128767 ...
         -0.44618132667134569 -0.48337236862196048 -0.4988010498682548];

vertC = [21.216985291697565 7.4788495625649318 3.0367300101320049 1.2449551542174275 ...
         0.52481184090840838 0.21659183357611089 0.085386356809149022 0.030369734021254091 ...
         0.0087678875650788468 0.0015122752148181151 -1.3877787807814457E-17 ...
         ];

vertD = [6.6141511762550769 2.8065370467678581 1.6619354303953071 1.2396100467363911 ...
         1.0873253528701918 1.0302460657718893 1.0096294545713409 1.0026114712228051 ...
         1.0005212428805947 1.0000494096306491 1];

arclength = [0.098174652356210243; 0.098174646402923962; 0.098174636672935175; ...
             0.098174612264088038; 0.09817460546904079; 0.0981740754944422; ...
             0.0981740280081984; 0.098174337650390936; 0.098174406160202926; ...
             0.098174500719211355; 0.098174550044057965; 0.098174590251013466; ...
             0.098174616804665146; 0.098174635356304982; 0.098174646768666954; ...
             0.098174652280987054];

bisector = [0.049008570164780385; 0.14655373117284454; 0.24268749963529535; ...
            0.33648405480977606; 0.42704008459554382; -0.85774199031676657; ...
            -0.80321261558808887; -0.74095585227632277; -0.6715616838314068; ...
            -0.59498175859162383; -0.51348348492279994; -0.42704008459554377; ...
            -0.33648405480977606; -0.24268749963529526; -0.14655373117284448; ...
            -0.049008570164780357];

rmid = [0.99879907185176919; 0.98920349299929378; 0.97010548779425976; ...
        0.94169029804374915; 0.904234588376165; 0.85774199031676657; 0.80321261558808887; ...
        0.74095585227632277; 0.6715616838314068; 0.59498175859162383; 0.51348348492279994; ...
        0.42704008459554377; 0.33648405480977606; 0.24268749963529526; ...
        0.14655373117284448; 0.049008570164780357];

i = 16;

Ls = [0.99999451240041792; 0.9999944517610051; 0.99999435265257752; 0.99999410402716526; ...
      0.999994034813679; 0.99998863655988013; 0.99998815287104326; 0.99999130684700177; ...
      0.99999200467922778; 0.99999296784525094; 0.9999934702618607; 0.99999387980478971; ...
      0.99999415027690641; 0.99999433924154946; 0.99999445548641663; 1.0000003401839683; ...
      ];

Ltheta = [0.99999706190177018; 1.0000207200202758; 1.0000699694660327; ...
          1.0001487526544368; 1.0002647876588557; 1.000009041524893; 0.9999997716326372; ...
          0.99999982140855048; 0.9999975055475272; 0.99878890602312542; ...
          0.99878890602312531; 0.99878890602312531; 0.99878890602312531; ...
          0.99878890602312542; 0.99878890602312531; 1];

sigmaS = [-0.0001097521526338685; -0.00011095995163079575; -0.00011287794921710503; ...
          -0.00011759690408785772; -0.0001182687388656678; -0.00022727054054871854; ...
          -0.00023694465664747533; -0.00017386417753728622; -0.00015990708155833033; ...
          -0.00011852186016958477; -0.0001084795105055214; -0.0001002935341631872; ...
          -9.4887318937697529E-5; -9.1110281990047781E-5; -8.8786772941551462E-5; ...
          6.803677630728977E-6];

sigmaTheta = [-5.8761481142033389E-5; 0.00041439327805592452; 0.0013993124185845973; ...
              0.0029747129952961693; 0.0052946865994052406; 0.00018083018115511162; ...
              -4.5652156333009231E-6; -3.5706803669377507E-6; -4.9887984463481416E-5; ...
              -0.024243830440296987; -0.024243836627244608; -0.024243841664904897; ...
              -0.024243844989140229; -0.02424384731030349; -0.024243848737678375; ...
              1.7358336990014323E-12];

r_dot = ...
  [0 0.00959656644408233 0.0191007121153352 0.02841807832751396 0.037465782996258586 ...
   0 0 0 0 0 0 0 0 0 0 0 0;
   0.098033793149011964 0.0975390410025647 0.096055281570608261 0.093586315756413629 ...
   0.090130687943782539 0 0 0 0 0 0 0 0 0 0 0 0];

rA = 11.773768309524469;

rB = 26.543470512444607;

rC = 21.216985291697565;

rD = 6.6141511762550769;

rSigmaS = -0.00022727054054871854;

rSigmaTheta = 0.00018083018115511162;

rZmid = 0.42704008459554382;

index = 5;

lA = -0.35297921715735775;

lB = -0.22590116512615277;

lC = -0.0862212849342211;

lD = 1.0097393547147022;

lSigmaS = -0.0001182687388656678;

lSigmaTheta = 0.0052946865994052406;

lZmid = 0.42704008459554382;

rRmid = -0.85774199031676657;

ans = [0.47139673682599781; 0.88192126434835494];
