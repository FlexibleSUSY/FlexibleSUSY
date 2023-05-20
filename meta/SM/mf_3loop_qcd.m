(*

  Relation between the quark pole mass and the running quark mass in
  QCD at 3-loop level [arxiv:hep-ph/9912391, arxiv:hep-ph/9911434].

  Extracted from meta/ThreeLoopQCD.m via:

Needs["SARAH`"];
Get["meta/ThreeLoopQCD.m"];

MfOvermf = Collect[
   GetMTopPoleOverMTopMSbar[{1,h^1,h^2,h^3}, TopQuark, Q, NH, NL] //. {
      FlexibleSUSY`M[Fu] -> mf,
      h -> (4Pi)^2 k,
      Log[mf^2/Q^2] -> L,
      Log[Q^2/mf^2] -> -L
   }
   , {k, g3, L}, FullSimplify
];

mfOverMf = Collect[
   GetMTopMSbarOverMTopPole[{1,h^1,h^2,h^3}, TopQuark, Q, NH, NL] //. {
      FlexibleSUSY`M[Fu] -> mf,
      h -> (4Pi)^2 k,
      Log[mf^2/Q^2] -> L,
      Log[Q^2/mf^2] -> -L
   }
   , {k, g3, L}, FullSimplify
];

 *)

(*
   k  = 1/(4Pi)^2;
   L  = Log[mf^2/Q^2];
   mf = quark MS-bar mass
   Mf = quark pole mass
   g3 = strong gauge coupling
   Q  = renormalization scale
   NL = number of light flavours (usually NL = 5)
   NH = number of heavy flavours (usually NH = 1)
 *)

(* ratio Mf / mf *)
MfOvermf = 1 + g3^2*k*(16/3 - 4*L) + g3^4*k^2*((-242*L)/3 + 22*L^2 + 
   (3049 - 2*NL*(71 + 8*Pi^2) + NH*(-286 + 32*Pi^2) + 32*Pi^2*(2 + Log[2]) - 
     48*Zeta[3])/18) + g3^6*k^3*((-4*L^3*(1231 + 4*(-22 + NL)*NL))/27 + 
   (2*L^2*(17395 + 2*NL*(-911 + 26*NL)))/27 + 
   (2*L*(-194453 + 22658*NL + NH*(5148 - 576*Pi^2) - 4*NL^2*(89 + 12*Pi^2) - 
      48*Pi^2*(99 + 37*Log[2]) + 24*NL*Pi^2*(49 + Log[16]) + 
      72*(67 + 28*NL)*Zeta[3]))/81 + 
   (-4*NH^2*(-47405 + 1152*Pi^2 + 23760*Zeta[3]) + 
     5*(12*NL*(-88715 + 244*Pi^4 + 96*Log[2]^4 + 
         12*Pi^2*(-965 + 8*Log[2]*(-11 + Log[4])) + 2304*PolyLog[4, 1/2] - 
         25452*Zeta[3]) + 4*NL^2*(2353 + 936*Pi^2 + 3024*Zeta[3]) - 
       9*(-1205357 + 1364*Pi^4 + 3648*Log[2]^4 + 87552*PolyLog[4, 1/2] + 
         11664*Zeta[3] + 4*Pi^2*(-25379 + 48*Log[2]*(-235 + 14*Log[2]) + 
           7986*Zeta[3]) - 81840*Zeta[5])) + 
     20*NH*(-2*NL*(-5917 + 468*Pi^2 + 864*Zeta[3]) + 
       3*(-162911 + 328*Pi^4 + 96*Log[2]^4 + 2304*PolyLog[4, 1/2] + 
         27036*Zeta[3] - 4*Pi^2*(-13627 + 24*Log[2]*(640 + Log[2]) + 
           486*Zeta[3]) + 9720*Zeta[5])))/7290);

(* ratio mt / Mt *)
(* Eq. (10) of arxiv:hep-ph/9912391 *)
mfOverMf = 1 + g3^2*k*(-16/3 + 4*L) + g3^4*k^2*(38*L - 6*L^2 + 
   (-2537 + NH*(286 - 32*Pi^2) + 2*NL*(71 + 8*Pi^2) - 32*Pi^2*(2 + Log[2]) + 
     48*Zeta[3])/18) + g3^6*k^3*((4*L^3*(475 + 4*(-22 + NL)*NL))/27 - 
   (2*L^2*(8971 + 2*NL*(-911 + 26*NL)))/27 + 
   (2*L*(118547 - 20102*NL + 4*NL^2*(89 + 12*Pi^2) + 1200*Pi^2*(3 + Log[2]) - 
      24*NL*Pi^2*(37 + Log[16]) - 72*(55 + 28*NL)*Zeta[3]))/81 + 
   (4*NH^2*(-47405 + 1152*Pi^2 + 23760*Zeta[3]) - 
     5*(12*NL*(-78491 + 96*Log[2]^4 + 4*Pi^2*(-2607 + 61*Pi^2 + 
           8*Log[2]*(-33 + Log[64])) + 2304*PolyLog[4, 1/2] - 25452*Zeta[3]) + 
       4*NL^2*(2353 + 936*Pi^2 + 3024*Zeta[3]) - 
       9*(-937229 + 1364*Pi^4 + 3648*Log[2]^4 + 87552*PolyLog[4, 1/2] + 
         7056*Zeta[3] + 4*Pi^2*(-23843 + 48*Log[2]*(-219 + 14*Log[2]) + 
           7986*Zeta[3]) - 81840*Zeta[5])) - 
     20*NH*(-2*NL*(-5917 + 468*Pi^2 + 864*Zeta[3]) + 
       3*(-142319 + 328*Pi^4 + 96*Log[2]^4 + 2304*PolyLog[4, 1/2] + 
         27036*Zeta[3] - 4*Pi^2*(-13051 + 24*Log[2]*(640 + Log[2]) + 
           486*Zeta[3]) + 9720*Zeta[5])))/7290);
