(* 4-loop beta function of Chetyrkin, Zoller [JHEP 1606 (2016) 175]
   [arXiv:1604.00853] Chetyrkin:2016ruf.
   Copied from https://www.ttp.kit.edu/Progdata/ttp16/ttp16-008/

   Note: Chetyrkin:2016ruf defines the beta function to be

       beta = d yt / d log(mu^2) = 1/2 d yt / d log(mu)

   In order to convert the beta function of Chetyrkin:2016ruf to the
   SARAH convention

       beta_SARAH = d yt / d log(mu)

   the expression has been multiplied by an overall factor 2.

   Note: In Chetyrkin:2016ruf the quartic Higgs coupling lambda is
   defined by the Lagrangian

       L = - lambda |\Phi|^4

   In order to convert this definition of lambda to the one used in
   SARAH:

       L = - lambda_SARAH/2 |\Phi|^4

   one has to multiply the beta function by an additional factor 2 and
   replace

       lambda -> lambdaSARAH/2

   This has not yet been done!

   Definitions:

   h = 1/(4Pi)^2;

   cf = 4/3; ca = 3; tr = 1/2; dR = 3; ng = 8; nf = 6;
   dFA4 = 15/2; dFF4 = 5/12;
   z3 = Zeta[3]; z4 = Zeta[4]; z5 = Zeta[5];

 *)

betalambda = 2 (

+ h*((3*g1^4)/16 + (3*g1^2*g2^2)/8 + (9*g2^4)/16 - (3*g1^2*lambda)/2 - (9*g2^2*lambda)/2 
+ 12*lambda^2 + 6*lambda*yb^2 - 3*yb^4 + 6*lambda*yt^2 - 3*yt^4 + 2*lambda*ytau^2 - ytau^4)

+  h^2*((39*g1^2*g2^2*lambda)/8 + 18*g1^2*lambda^2 + 54*g2^2*lambda^2 - 156*lambda^3 + g2^6*(497/32 - nf) 
+ g1^4*g2^2*(-239/96 - (5*nf)/9) + g1^6*(-59/96 - (5*nf)/9) + g1^2*g2^4*(-97/96 - nf/3)
+ g1^4*lambda*(229/48 + (25*nf)/18) + g2^4*lambda*(-313/16 + (5*nf)/2) + (5*g1^4*yb^2)/8 
+ (9*g1^2*g2^2*yb^2)/4 - (9*g2^4*yb^2)/8 + (25*g1^2*lambda*yb^2)/12 + (45*g2^2*lambda*yb^2)/4 
+ 40*gs^2*lambda*yb^2 -    72*lambda^2*yb^2 + (2*g1^2*yb^4)/3 - 16*gs^2*yb^4 - (3*lambda*yb^4)/2 
+ 15*yb^6 - (19*g1^4*yt^2)/8 + (21*g1^2*g2^2*yt^2)/4 - (9*g2^4*yt^2)/8 + (85*g1^2*lambda*yt^2)/12 
+ (45*g2^2*lambda*yt^2)/4 + 40*gs^2*lambda*yt^2 - 72*lambda^2*yt^2 - 21*lambda*yb^2*yt^2 - 3*yb^4*yt^2 
- (4*g1^2*yt^4)/3 - 16*gs^2*yt^4 - (3*lambda*yt^4)/2 - 3*yb^2*yt^4 + 15*yt^6 - (25*g1^4*ytau^2)/8 
+ (11*g1^2*g2^2*ytau^2)/4 - (3*g2^4*ytau^2)/8 + (25*g1^2*lambda*ytau^2)/4 + (15*g2^2*lambda*ytau^2)/4 
- 24*lambda^2*ytau^2 - 2*g1^2*ytau^4 - (lambda*ytau^4)/2 + 5*ytau^6) 

+  h^3*(873*lambda^3*yb^2 + 873*lambda^3*yt^2 + 192*gs^4*yb^2*yt^2 + (477*g2^2*yb^4*yt^2)/32 
+ (477*g2^2*yb^2*yt^4)/32 + 291*lambda^3*ytau^2 + (41*g1^4*yb^2*ytau^2)/24 - (5*g1^2*g2^2*yb^2*ytau^2)/4 
+ (9*g2^4*yb^2*ytau^2)/8 -  9*g1^2*lambda*yb^2*ytau^2 - 27*g2^2*lambda*yb^2*ytau^2 - 216*lambda^2*yb^2*ytau^2
+ 240*lambda*yb^4*ytau^2 - (297*yb^6*ytau^2)/8 + (701*g1^4*yt^2*ytau^2)/24 + (29*g1^2*g2^2*yt^2*ytau^2)/4 
+ (9*g2^4*yt^2*ytau^2)/8 -  9*g1^2*lambda*yt^2*ytau^2 - 27*g2^2*lambda*yt^2*ytau^2 - 216*lambda^2*yt^2*ytau^2 
+ 21*lambda*yb^2*yt^2*ytau^2 + (45*yb^4*yt^2*ytau^2)/8 + 240*lambda*yt^4*ytau^2 + (45*yb^2*yt^4*ytau^2)/8 
- (297*yt^6*ytau^2)/8 +  240*lambda*yb^2*ytau^4 - 72*yb^4*ytau^4 + 240*lambda*yt^2*ytau^4 + 12*yb^2*yt^2*ytau^4
- 72*yt^4*ytau^4 - (297*yb^2*ytau^6)/8 - (297*yt^2*ytau^6)/8 + gs^2*lambda*yb^4*(895 - 1296*z3) 
+ gs^2*lambda*yt^4*(895 - 1296*z3) + lambda^2*yb^2*yt^2*(117 - 864*z3) + g2^4*lambda^2*(1995/8 - (141*nf)/2 - 513*z3) 
+ g2^2*lambda^2*yb^2*(639/4 - 432*z3) + g2^2*lambda^2*yt^2*(639/4 - 432*z3) + lambda*yb^6*(117/8 - 198*z3) 
+ lambda*yt^6*(117/8 - 198*z3) + g1^2*lambda^2*yb^2*(417/4 - 192*z3) + g2^4*lambda*yb^2*(-3933/64 - (63*nf)/8 - (351*z3)/2)
+ g2^4*lambda*yt^2*(-3933/64 - (63*nf)/8 - (351*z3)/2) + g1^2*g2^2*lambda^2*(-333 - 162*z3) 
+ g2^2*lambda^2*ytau^2*(213/4 - 144*z3) + g1^2*lambda*ytau^4*(507/8 - 117*z3) + gs^2*lambda*yb^2*yt^2*(82 - 96*z3) 
+ g1^2*g2^2*yt^4*(-1079/192 - (743*z3)/8) + g1^4*lambda^2*(-183 - (235*nf)/6 - 81*z3) + g1^4*lambda*yt^2*(-112447/1728 
- (635*nf)/72 - (449*z3)/6) + lambda*ytau^6*(-1241/8 - 66*z3) + g1^4*lambda*ytau^2*(-1783/64 - (65*nf)/8 - (123*z3)/2) 
+ g2^4*lambda*ytau^2*(-1311/64 - (21*nf)/8 - (117*z3)/2) + g2^4*gs^2*yb^2*(651/8 - 54*z3) + g2^4*gs^2*yt^2*(651/8 - 54*z3) 
+ g2^4*yb^4*(13653/128 - (39*nf)/8 - (819*z3)/16) + g2^4*yt^4*(13653/128 - (39*nf)/8 - (819*z3)/16) 
+ g1^2*lambda^2*yt^2*(-195/4 - 48*z3) + gs^2*yb^4*yt^2*(-2 - 48*z3) + gs^2*yb^2*yt^4*(-2 - 48*z3) 
+ gs^4*lambda*yb^2*(1820/3 - 32*nf - 48*z3) + gs^4*lambda*yt^2*(1820/3 - 32*nf - 48*z3) + g1^2*g2^2*ytau^4*(-15/64 
- (381*z3)/8) + g1^2*g2^2*yb^4*(-3239/192 - (311*z3)/8) + yb^8*(-1599/8 - 36*z3) + yt^8*(-1599/8 - 36*z3) 
+ yb^6*yt^2*(-717/8 - 36*z3) + yb^2*yt^6*(-717/8 - 36*z3) + g1^2*g2^2*gs^2*yb^2*(233/4 - 36*z3)
+ g1^2*g2^2*gs^2*yt^2*(249/4 - 36*z3) + g1^2*yb^2*yt^4*(1337/96 - 28*z3) + g2^2*yb^6*(3411/32 - 27*z3) 
+ g2^2*yt^6*(3411/32 - 27*z3) + g1^2*yb^6*(5111/96 - 25*z3) + g1^2*gs^2*yt^4*(931/18 - (56*z3)/3) 
+ g1^4*gs^2*yt^2*(587/24 - 18*z3) + g1^4*gs^2*yb^2*(683/24 - 18*z3) + g2^4*ytau^4*(4503/128 - (13*nf)/8 - (273*z3)/16) 
+ g1^4*yb^4*(15137/3456 - (415*nf)/72 - (2035*z3)/144) + g1^2*g2^6*(-54053/3456 - (8341*nf)/1728 - (5*nf^2)/54 
- (405*z3)/32) + ytau^8*(-143/8 - 12*z3) + g2^2*ytau^6*(1137/32 - 9*z3) + g1^4*lambda*yb^2*(-127303/1728
- (155*nf)/72 - (47*z3)/6) + g1^4*g2^2*ytau^2*(6657/256 - (5*nf)/12 - (15*z3)/2) + g1^6*ytau^2*(3929/256 
+ (55*nf)/12 - (15*z3)/4) + g1^4*g2^2*yt^2*(23521/768 + (5*nf)/12 - 3*z3) + g1^6*yt^2*(42943/2304 
+ (215*nf)/36 - (5*z3)/2) + g1^2*lambda*yb^2*yt^2*(-929/12 - 2*z3) + g1^2*g2^4*ytau^2*(1833/256 - nf/4 
- (3*z3)/2) + g1^4*yb^2*yt^2*(-709/64 - z3) + 72*yb^4*yt^4*z3 + g1^6*yb^2*(12043/2304 + (95*nf)/36 + (5*z3)/4) 
+ g1^4*g2^2*yb^2*(4403/256 + (25*nf)/12 + (9*z3)/2) + g1^2*g2^4*yt^2*(3103/256 + nf/4 + (27*z3)/4) 
+ g1^2*g2^4*yb^2*(4179/256 + (5*nf)/4 + 9*z3) + g1^2*g2^2*lambda*yb^2*(-3009/32 + 12*z3) + g1^2*g2^2*yb^2*yt^2*(1001/96 
+ (31*z3)/2) + g1^2*yt^6*(3467/96 + 17*z3) + g1^4*yt^4*(100913/3456 - (115*nf)/72 + (2957*z3)/144) 
+ g1^4*ytau^4*(5697/128 + (65*nf)/24 + (375*z3)/16) + g1^2*lambda^3*(-158 + 24*z3) 
+ g2^2*gs^2*yb^4*(-31/2 + 24*z3) + g2^2*gs^2*yt^4*(-31/2 + 24*z3) + g2^6*ytau^2*(-5739/256 + (9*nf)/4 + (99*z3)/4) 
+ g1^2*yb^4*yt^2*(-2299/96 + 26*z3) + gs^4*yb^4*(-626/3 + 20*nf + 32*z3) + gs^4*yt^4*(-626/3 + 20*nf + 32*z3) 
+ g1^2*ytau^6*(135/32 + 33*z3) + g1^2*gs^2*lambda*yb^2*(-991/18 + 40*z3) + g1^2*gs^2*yb^4*(-641/18 + (136*z3)/3) 
+ g2^2*lambda*yb^2*yt^2*(-531/4 + 54*z3) + g1^2*lambda*yt^4*(-2485/24 + 57*z3) + g2^4*yb^2*yt^2*(-351/64 - 6*nf + (117*z3)/2) 
+ g2^2*lambda^3*(-474 + 72*z3) + g2^6*yb^2*(-17217/256 + (27*nf)/4 + (297*z3)/4) + g2^6*yt^2*(-17217/256 + (27*nf)/4 + (297*z3)/4) 
+ g1^2*lambda^2*ytau^2*(-541/4 + 96*z3) + g2^2*gs^2*yb^2*yt^2*(-8 + 96*z3) + g1^2*g2^2*lambda*ytau^2*(-3771/32 + 126*z3) 
+ g1^2*gs^2*lambda*yt^2*(-2419/18 + 136*z3) + lambda*yb^4*yt^2*(6399/8 + 144*z3) + lambda*yb^2*yt^4*(6399/8 + 144*z3) 
+ g2^2*lambda*ytau^4*(-1587/8 + 171*z3) + g1^2*g2^2*lambda*yt^2*(-6509/32 + 177*z3) + g2^2*gs^2*lambda*yb^2*(-489/2 + 216*z3) 
+ g2^2*gs^2*lambda*yt^2*(-489/2 + 216*z3) + gs^2*yb^6*(-38 + 240*z3) + gs^2*yt^6*(-38 + 240*z3) + g1^2*lambda*yb^4*(-5737/24 + 249*z3) 
+ lambda^2*ytau^4*(717/2 + 252*z3) + g2^2*lambda*yb^4*(-4977/8 + 513*z3) + g2^2*lambda*yt^4*(-4977/8 + 513*z3) 
+ lambda^2*yb^4*(1719/2 + 756*z3) + lambda^2*yt^4*(1719/2 + 756*z3) + gs^2*lambda^2*yb^2*(-1224 + 1152*z3) 
+ gs^2*lambda^2*yt^2*(-1224 + 1152*z3) + lambda^4*(3588 + 2016*z3) + g2^4*gs^2*lambda*((135*nf)/4 - 36*nf*z3) 
+ g2^8*(982291/3072 - (14749*nf)/384 - (5*nf^2)/12 - (2781*z3)/128 - (45*nf*z3)/2) + g1^4*gs^2*lambda*((55*nf)/4 - (44*nf*z3)/3) 
+ g1^6*lambda*(12679/432 + (5995*nf)/324 + (875*nf^2)/486 + (9*z3)/8 - (95*nf*z3)/9) 
+ g1^2*g2^4*lambda*(4553/32 + (33*nf)/4 - (249*z3)/8 - 3*nf*z3) + g1^4*g2^2*lambda*(979/8 + (95*nf)/8 - (3*z3)/8 - 3*nf*z3) 
+ g1^2*g2^4*gs^2*((-51*nf)/16 + 3*nf*z3) + g1^6*g2^2*(-29779/6912 - (18001*nf)/5184 - (125*nf^2)/486 + (75*z3)/32 + (61*nf*z3)/18) 
+ g1^4*g2^4*(-64693/3456 + (149*nf)/1296 - (25*nf^2)/162 + (873*z3)/64 + (7*nf*z3)/2) + g1^6*gs^2*((-187*nf)/48 + (11*nf*z3)/3) 
+ g1^4*g2^2*gs^2*((-187*nf)/48 + (11*nf*z3)/3) + g1^8*(-6845/9216 - (20735*nf)/3456 - (125*nf^2)/324 + (99*z3)/128 + (95*nf*z3)/18) 
+ g2^6*gs^2*((-153*nf)/16 + 9*nf*z3) + g2^6*lambda*(-46489/288 + (3515*nf)/72 + (35*nf^2)/18 + (2259*z3)/8 + 45*nf*z3))

     + h^4*gs^6*yt^4 * (
          - 2942/3*cf^3*dR
          - 64*tr*cf^2*dR
          + 3584/3*ca*cf^2*dR
          + 5888/9*ca*tr*cf*dR
          - 121547/243*ca^2*cf*dR
          + 562/3*nf*tr*cf^2*dR
          - 256/9*nf*tr^2*cf*dR
          - 2644/243*nf*ca*tr*cf*dR
          + 10912/243*nf^2*tr^2*cf*dR
          + 160*z5*cf^3*dR
          + 720*z5*ca*cf^2*dR
          - 160*z5*ca*tr*cf*dR
          - 520*z5*ca^2*cf*dR
          + 288*z4*cf^3*dR
          + 32*z4*ca*cf^2*dR
          - 88*z4*ca^2*cf*dR
          - 160*z4*nf*tr*cf^2*dR
          + 128*z4*nf*ca*tr*cf*dR
          + 48*z3*cf^3*dR
          - 3304/3*z3*ca*cf^2*dR
          + 352*z3*ca*tr*cf*dR
          + 1880/3*z3*ca^2*cf*dR
          + 32/3*z3*nf*tr*cf^2*dR
          + 16*z3*nf*ca*tr*cf*dR
          - 128/3*z3*nf^2*tr^2*cf*dR
          )
);
