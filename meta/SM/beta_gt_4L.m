(* 4-loop beta function of Chetyrkin, Zoller [JHEP 1606 (2016) 175]
   [arXiv:1604.00853] Chetyrkin:2016ruf.
   Copied from https://www.ttp.kit.edu/Progdata/ttp16/ttp16-008/

   Note: Chetyrkin:2016ruf defines the beta function to be

       beta = d yt / d log(mu^2) = 1/2 d yt / d log(mu)

   In order to convert the beta function of Chetyrkin:2016ruf to the
   SARAH convention

       beta_SARAH = d yt / d log(mu)

   the expression has been multiplied by an overall factor 2.

   Note: In the 4-loop expression "betaytl4" in ttp16_008.m a global
   factor yt is missing, which is present in Eq.(3.5).  This factor
   has been added below.

   Note: In the 4-loop expression "betaytl4" g1 = gY denotes the
   non-GUT normalized gauge coupling of U(1)_Y.

   h = 1/(4Pi)^2;

   cf = 4/3; ca = 3; tr = 1/2; dR = 3; ng = 8; nf = 6;
   dFA4 = 15/2; dFF4 = 5/12;
   z3 = Zeta[3]; z4 = Zeta[4]; z5 = Zeta[5];

   toSARAHConvention = { lambda -> \[Lambda]/2, g1 -> Sqrt[3/5] g1 }

 *)

betayt = 2 (

+ h*yt/2*((-17*g1^2)/12 - (9*g2^2)/4 - 8*gs^2 + (3*yb^2)/2 + (9*yt^2)/2 + ytau^2)

+ h^2*yt/2*(g1^4/8 - (3*g1^2*g2^2)/4 - (35*g2^4)/4 + (19*g1^2*gs^2)/9 + 9*g2^2*gs^2 - (404*gs^4)/3 + 6*lambda^2 
+ (145*g1^4*nf)/162 + (g2^4*nf)/2 + (40*gs^4*nf)/9 + (7*g1^2*yb^2)/48 + (99*g2^2*yb^2)/16 + 4*gs^2*yb^2 
-  yb^4/4 + (131*g1^2*yt^2)/16 + (225*g2^2*yt^2)/16 + 36*gs^2*yt^2 - 12*lambda*yt^2 - (11*yb^2*yt^2)/4 - 12*yt^4 
+ (25*g1^2*ytau^2)/8 + (15*g2^2*ytau^2)/8 + (5*yb^2*ytau^2)/4 - (9*yt^2*ytau^2)/4 - (9*ytau^4)/4) 

+ h^3*yt/2*((18103*g1^6)/5184 + (1081*g1^4*g2^2)/192 + (309*g1^2*g2^4)/64 - (14677*g2^6)/576 - (1187*g1^4*gs^2)/108 
- (107*g1^2*g2^2*gs^2)/4 + 66*g2^4*gs^2 - (127*g1^2*gs^4)/36 + (531*g2^2*gs^4)/4 - 2498*gs^6 - (121*g1^4*lambda)/16 
+ (39*g1^2*g2^2*lambda)/8 - (171*g2^4*lambda)/16 + 15*g1^2*lambda^2 + 45*g2^2*lambda^2 - 36*lambda^3 + (267065*g1^6*nf)/23328 
+ (241*g1^4*g2^2*nf)/288 - (3*g1^2*g2^4*nf)/32 - (1139*g2^6*nf)/288 + (5281*g1^4*gs^2*nf)/648 + (57*g2^4*gs^2*nf)/8 
+ (220*g1^2*gs^4*nf)/27 + 19*g2^2*gs^4*nf + (4432*gs^6*nf)/27 + (9125*g1^6*nf^2)/4374 + (25*g2^6*nf^2)/18 
+ (280*gs^6*nf^2)/81 - (35153*g1^4*yb^2)/6912 + (1245*g1^2*g2^2*yb^2)/128 + (13653*g2^4*yb^2)/256 - (457*g1^2*gs^2*yb^2)/18 
- (27*g2^2*gs^2*yb^2)/2 - (277*gs^4*yb^2)/2 - (291*lambda^2*yb^2)/4 - (115*g1^4*nf*yb^2)/864 - (69*g2^4*nf*yb^2)/32 
- (7*gs^4*nf*yb^2)/3 - (959*g1^2*yb^4)/96 - (2283*g2^2*yb^4)/32 + 82*gs^2*yb^4 + 15*lambda*yb^4 + (477*yb^6)/16 
- (44179*g1^4*yt^2)/6912 + (2699*g1^2*g2^2*yt^2)/128 + (49239*g2^4*yt^2)/256 - 42*g1^2*gs^2*yt^2 - 168*g2^2*gs^2*yt^2 
+ (4799*gs^4*yt^2)/6 - (127*g1^2*lambda*yt^2)/6 - (135*g2^2*lambda*yt^2)/2 + 16*gs^2*lambda*yt^2 + (15*lambda^2*yt^2)/4 
- (2875*g1^4*nf*yt^2)/288 - (351*g2^4*nf*yt^2)/32 - 27*gs^4*nf*yt^2 - (461*g1^2*yb^2*yt^2)/32 - (2307*g2^2*yb^2*yt^2)/32 
+ 27*gs^2*yb^2*yt^2 + 93*lambda*yb^2*yt^2 + (825*yb^4*yt^2)/8 - (2437*g1^2*yt^4)/48 - (1593*g2^2*yt^4)/16 - 157*gs^2*yt^4 
+ 198*lambda*yt^4 + (739*yb^2*yt^4)/16 + (339*yt^6)/8 - (20215*g1^4*ytau^2)/1152 - (347*g1^2*g2^2*ytau^2)/64 
+ (2121*g2^4*ytau^2)/128 - (45*lambda^2*ytau^2)/2 - (65*g1^4*nf*ytau^2)/16 - (21*g2^4*nf*ytau^2)/16 + (491*g1^2*yb^2*ytau^2)/72 
- (153*g2^2*yb^2*ytau^2)/8 - (43*gs^2*yb^2*ytau^2)/6 + 22*yb^4*ytau^2 - 21*g1^2*yt^2*ytau^2 - (81*g2^2*yt^2*ytau^2)/4 
+ (5*gs^2*yt^2*ytau^2)/2 + 30*lambda*yt^2*ytau^2 + (7*yb^2*yt^2*ytau^2)/2 + (21*yt^4*ytau^2)/2 - (45*g1^2*ytau^4)/16 
- (315*g2^2*ytau^4)/16 + 15*lambda*ytau^4 + (53*yb^2*ytau^4)/4 + (207*yt^2*ytau^4)/8 + (71*ytau^6)/16 - (17*g1^6*z3)/24 
- (17*g1^4*g2^2*z3)/8 - (9*g1^2*g2^4*z3)/8 + (45*g2^6*z3)/8 - (1615*g1^6*nf*z3)/162 - (17*g1^4*g2^2*nf*z3)/6
- (3*g1^2*g2^4*nf*z3)/2 + (45*g2^6*nf*z3)/2 - (374*g1^4*gs^2*nf*z3)/27 - 18*g2^4*gs^2*nf*z3 - (88*g1^2*gs^4*nf*z3)/9
- 24*g2^2*gs^4*nf*z3 + (320*gs^6*nf*z3)/3 - (199*g1^4*yb^2*z3)/72 + (9*g1^2*g2^2*yb^2*z3)/2 - (225*g2^4*yb^2*z3)/8 
- (28*g1^2*gs^2*yb^2*z3)/3 - 108*g2^2*gs^2*yb^2*z3 - 44*gs^4*yb^2*z3 + (19*g1^2*yb^4*z3)/6 + (63*g2^2*yb^4*z3)/2 
- 64*gs^2*yb^4*z3 + (9*yb^6*z3)/2 - (31*g1^4*yt^2*z3)/24 + (123*g1^2*g2^2*yt^2*z3)/4 - (729*g2^4*yt^2*z3)/8 + 60*g1^2*gs^2*yt^2*z3 
+ 180*g2^2*gs^2*yt^2*z3 - 228*gs^4*yt^2*z3 + (5*g1^2*yb^2*yt^2*z3)/6 - (9*g2^2*yb^2*yt^2*z3)/2 - 32*gs^2*yb^2*yt^2*z3 
- 48*yb^4*yt^2*z3 + (27*yt^6*z3)/2 - (269*g1^4*ytau^2*z3)/12 - 3*g1^2*g2^2*ytau^2*z3 - (81*g2^4*ytau^2*z3)/4 
- 9*g1^2*yb^2*ytau^2*z3 + 9*g2^2*yb^2*ytau^2*z3 + 8*g1^2*yt^2*ytau^2*z3 - 9*g2^2*yt^2*ytau^2*z3 - 9*g1^2*ytau^4*z3
+ 9*g2^2*ytau^4*z3 + 3*ytau^6*z3)

       + h^4*gs^8 * yt (
          + 32*dFA4*dR^-1
          + 1261/8*cf^4
          - 15349/12*ca*cf^3
          + 34045/36*ca^2*cf^2
          - 70055/72*ca^3*cf
          - 64*nf*dFF4*dR^-1
          + 280/3*nf*tr*cf^3
          + 8819/27*nf*ca*tr*cf^2
          + 65459/162*nf*ca^2*tr*cf
          - 304/27*nf^2*tr^2*cf^2
          - 1342/81*nf^2*ca*tr^2*cf
          + 664/81*nf^3*tr^3*cf
          - 440*z5*ca^2*cf^2
          + 440*z5*ca^3*cf
          + 480*z5*nf*tr*cf^3
          - 80*z5*nf*ca*tr*cf^2
          - 400*z5*nf*ca^2*tr*cf
          + 264*z4*nf*ca*tr*cf^2
          - 264*z4*nf*ca^2*tr*cf
          - 96*z4*nf^2*tr^2*cf^2
          + 96*z4*nf^2*ca*tr^2*cf
          - 240*z3*dFA4*dR^-1
          + 336*z3*cf^4
          - 316*z3*ca*cf^3
          + 152*z3*ca^2*cf^2
          - 1418/9*z3*ca^3*cf
          + 480*z3*nf*dFF4*dR^-1
          - 552*z3*nf*tr*cf^3
          - 368*z3*nf*ca*tr*cf^2
          + 2684/3*z3*nf*ca^2*tr*cf
          + 160*z3*nf^2*tr^2*cf^2
          - 160*z3*nf^2*ca*tr^2*cf
          - 128/9*z3*nf^3*tr^3*cf
          )
);
