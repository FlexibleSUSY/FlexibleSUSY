(* 1-, 2-, 3- and 4-loop beta functions for the strong coupling in the
   Standard Model O(as^n at^n alambda^n).
   Provided by A. Bednyakov and A.F. Pikelner [1508.02680]
*)

(*
  To convert to beta(g3):

  kappa = 1/(4 Pi)^2;

  pref = (g3 kappa)^(-1) (-as);

  repl = {
     cA -> 3, Tf -> 1/2, NG -> nf/2, nf -> 6, h -> 1, cR -> 4/3,
     z3 -> Zeta[3], NA -> 8, R -> 3,
     dAA -> NA 3^2 (3^2 + 36)/24,
     dFA -> NA 3 (3^2 + 6)/48,
     dFF -> NA (3^4 - 6 3^2 + 18)/(96 3^2),
     as -> g3^2 kappa, ay -> gt^2 kappa, alam -> \[Lambda] kappa
  };

  toSUSYHDConvention = \[Lambda] -> \[Lambda]/2;

  pref {
     betaas1/kappa,
     betaas2/kappa^2,
     betaas3/kappa^3,
     betaas4/kappa^4
  } //. repl /. toSUSYHDConvention

*)

betaas1 = +as*h^2*(+11/3*cA - 8/3*Tf*NG);

betaas2 = +as^2*h^4*(+34/3*cA^2 - 40/3*Tf*cA*NG - 8*Tf*cR*NG) + 
   ay*as*h^4*(+4*Tf);

betaas3 = +as^3*
    h^6*(+2857/54*cA^3 - 2830/27*Tf*cA^2*NG - 410/9*Tf*cR*cA*NG + 
      4*Tf*cR^2*NG + 632/27*Tf^2*cA*NG^2 + 176/9*Tf^2*cR*NG^2) + 
   ay*as^2*h^6*(+24*Tf*cA + 6*Tf*cR) + ay^2*as*h^6*(-30*Tf);

betaas4 = +as^4*
    h^8*(-80/9*NA^-1*dAA + 1024/9*NA^-1*NG*dFA - 
      2816/9*NA^-1*NG^2*dFF + 704/3*z3*NA^-1*dAA - 
      3328/3*z3*NA^-1*NG*dFA + 2048/3*z3*NA^-1*NG^2*dFF + 
      150653/486*cA^4 - 44/9*cA^4*z3 - 78286/81*Tf*cA^3*NG + 
      272/3*Tf*cA^3*z3*NG + 14146/243*Tf*cR*cA^2*NG - 
      1312/9*Tf*cR*cA^2*z3*NG - 8408/27*Tf*cR^2*cA*NG + 
      704/9*Tf*cR^2*cA*z3*NG + 92*Tf*cR^3*NG + 
      31720/81*Tf^2*cA^2*NG^2 + 896/9*Tf^2*cA^2*z3*NG^2 + 
      68608/243*Tf^2*cR*cA*NG^2 + 1792/9*Tf^2*cR*cA*z3*NG^2 + 
      5408/27*Tf^2*cR^2*NG^2 - 2816/9*Tf^2*cR^2*z3*NG^2 + 
      3392/243*Tf^3*cA*NG^3 + 9856/243*Tf^3*cR*NG^3) + 
   ay*as^3*h^8*(+1970/9*Tf*cA^2 + 523/9*Tf*cR*cA - 72*Tf*cR*cA*z3 + 
      6*Tf*cR^2 - 144*Tf*cR^2*z3 - 872/9*Tf^2*cA*NG - 
      1288/9*Tf^2*cR*NG) + 
   ay^2*as^2*
    h^8*(-222*Tf*cA - 117*Tf*cR + 144*Tf*cR*z3 - 48*Tf^2 + 
      96*Tf^2*z3) + ay^2*as^2*R*h^8*(-16/3*Tf^2 - 32*Tf^2*z3) + 
   ay^3*as*h^8*(+423/2*Tf + 12*Tf*z3) + alam*ay^2*as*h^8*(+60*Tf) + 
   alam^2*ay*as*h^8*(-72*Tf);
