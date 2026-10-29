(*

CorrectionListVectorVector[EWSB][[1, 2]] gives data like
{bar[Fu], Fu, C[Fu, VP, bar[Fu]], C[Fu, VZ, bar[Fu]], FFV, 3, 1/2}

whereas we want data like ThreeParticleVertex2[VP] gives, along the lines of {{bar[Fd], Fd, Cp[VP, bar[Fd[{gI1}]], Fd[{gI2}]], FFV, 3, 1/2}

These normally feed into AddLoopScalar etc
*)

AddLoopVectorDiff[particles_, corrections_] :=
  Block[{particle1, particle2, temp, m1, m2, m0, mS, mV, v, fac, m11,
    m12, m21, m22, ind1, ind2, coupL1, coupR1, coupL2, coupR2, coup1,
    coup2, coup3, pi1Loop},

   {particle1, particle2} = getBlank /@ particles; (* remove indices*)

   If[getType[corrections[[1]]] === F,
    temp = OrderMasses[corrections[[1]], corrections[[2]], F];
    If[FreeQ[m11, corrections[[1]]] && FreeQ[m21, corrections[[2]]],
     m11 = m11 /. gI1 -> gI2;
     m12 = m12 /. gI1 -> gI2;
     m22 = m22 /. gI2 -> gI1;
     m21 = m21 /. gI2 -> gI1;];
    m11 = temp[[1]]; m12 = temp[[1]] /. Mass -> Mass2; m21 = temp[[2]];
     m22 = temp[[2]] /. Mass -> Mass2;
    coupL1 = corrections[[3]][PL];
    coupR1 = corrections[[3]][PR];
    coupL2 = corrections[[4]][PL];
    coupR2 = corrections[[4]][PR];

    (* original pi1Loop=(coupL1 conj[coupL1]+coupR1 conj[coupR1])*H0[p^2,m12,m22]+4 Re[conj[coupL1] coupR1]*B0[p^2,m12,m22] m11 m21;
    NB conj[coup[PL]] = -coup[PR] = coup[PL] for identical fields *)

    pi1Loop = (coupL1 coupL2 + coupR1 coupR2)*H0[p^2, m12, m22] + 2 (coupL1 coupR2 + coupR1 coupL2)*B0[p^2, m12, m22] m11 m21;
    If[AntiField[particle1] === particle1 && AntiField[particle2] === particle2,
     pi1Loop = 2 pi1Loop*corrections[[7]]*corrections[[6]];,
     pi1Loop = pi1Loop*corrections[[7]]*corrections[[6]];];
    Return[sum[gI1, 1, getGen[corrections[[1]]],sum[gI2, 1, getGen[corrections[[2]]], pi1Loop]]];
    , (* Not fermionic amplitude *)
    temp = OrderMasses[corrections[[1]], corrections[[2]], V];
    m11 = temp[[1]];
    m12 = temp[[1]] /. Mass -> Mass2; m21 = temp[[2]];
    m22 = temp[[2]] /. Mass -> Mass2;
    If[FreeQ[m11, corrections[[1]]] && FreeQ[m21, corrections[[2]]],
     m11 = m11 /. gI1 -> gI2;
     m12 = m12 /. gI1 -> gI2;
     m22 = m22 /. gI2 -> gI1;
     m21 = m21 /. gI2 -> gI1;];
    coup1 = corrections[[3]];
    coup2 = corrections[[4]];
    Switch[corrections[[5]],
     SSV,
     (*pi1Loop=-4 B22[p^2,m12,m22]-(A0[m12]+A0[m22]);*)
     pi1Loop = -4 B00[p^2, m12, m22];
     pi1Loop = pi1Loop*corrections[[7]]*corrections[[6]] coup1 coup2;
     If[conj[particle1] === particle1 && conj[particle2] === particle2,
      Return[sum[gI1, 1, getGen[corrections[[1]]],sum[gI2, 1, getGen[corrections[[2]]], 2 pi1Loop]]];,
      Return[sum[gI1, 1, getGen[corrections[[1]]],sum[gI2, 1, getGen[corrections[[2]]], pi1Loop]]];];,
     SVV,
     pi1Loop = coup1 coup2 B0[p^2, m12, m22];
     pi1Loop = pi1Loop*corrections[[7]]*corrections[[6]];
     If[conj[particle1] === particle1 && conj[particle2] === particle2, pi1Loop = 2 pi1Loop];
     Return[sum[gI1, 1, getGen[corrections[[1]]],sum[gI2, 1, getGen[corrections[[2]]], pi1Loop]]];,
     VVV,
     pi1Loop = 10 B00[p^2, m12, m22] + (m12 + m22 + 4*p^2)*B0[p^2, m12, m22] + A0[m12] + A0[m22] - 2 (m12 + m22 - 1/3 p^2)*rMS;
     pi1Loop = -pi1Loop*corrections[[7]]*corrections[[6]] coup1 coup2;
     If[conj[particle1] === particle1 && conj[particle2] === particle2, pi1Loop = 2 pi1Loop];
     Return[sum[gI1, 1, getGen[corrections[[1]]],
       sum[gI2, 1, getGen[corrections[[2]]], pi1Loop]]];, GGV,
     pi1Loop = coup1 coup2 B00[p^2, m12, m22];
     pi1Loop = pi1Loop*corrections[[7]]*corrections[[6]];
     If[conj[particle1] === particle1 && conj[particle2] === particle2, pi1Loop = 2 pi1Loop];
     Return[sum[gI1, 1, getGen[corrections[[1]]], sum[gI2, 1, getGen[corrections[[2]]], pi1Loop]]];,
     VVVV,
     coup1 = corrections[[3]][1];
     coup2 = corrections[[3]][2];
     coup3 = corrections[[3]][3];
     pi1Loop = 2*coup1*m12*rMS - (4*coup1 + coup2 + coup3)*A0[m12];
     pi1Loop = pi1Loop*corrections[[7]]*corrections[[6]];
     Return[sum[gI1, 1, getGen[corrections[[1]]], pi1Loop]];,
     SSVV,
     pi1Loop = coup1*A0[m12];
     pi1Loop = pi1Loop*corrections[[7]]*corrections[[6]];
     Return[sum[gI1, 1, getGen[corrections[[1]]], pi1Loop]];];
     ]; (* end switch *)

     ]; (* End block *)

GetVPVZSelfEnergy:=Module[{temp, temp2, particles, i, j, k, part, temploop,
  loopCorrection, correction},
 particles = {VectorP, VectorZ};
  (* First get the data for the 3-vertices *)
 temp = InsFields[{{C[particles[[1]], AntiField[FieldToInsert[1]],
      FieldToInsert[2]],
     C[AntiField[particles[[2]]], AntiField[FieldToInsert[2]],
      FieldToInsert[1]]}, {Internal[1] -> FieldToInsert[1],
     Internal[2] -> FieldToInsert[2], External[1] -> particles[[1]],
     External[2] -> particles[[2]], Index[1] -> gO1, Index[2] -> gO2,
     Index[3] -> gI1, Index[4] -> gI2}}];

 If[temp =!= {},
  temp2 = {};
  For[k = 1, k <= Length[temp], k++,
   temp2 =
     Join[temp2, {ReleaseHold[
        Hold[{AntiField[Internal[1]], Internal[2],
           C[External[1][{Index[1]}], AntiField[Internal[1][{Index[3]}]],
             Internal[2][{Index[4]}]],
           C[AntiField[External[2][{Index[2]}]],
            AntiField[Internal[2][{Index[4]}]], Internal[1][{Index[3]}]],
            VType[getType[Internal[1]], getType[Internal[2]],
            getType[External[1]]],
           CalculateColorFactor[getBlank[External[1]], Internal[1],
            Internal[2]],
           CalculateSymmetryFactor[Internal[1], Internal[2]]}] /.
         temp[[k, 2]]]}];
   ];

   ];(* temp =!= {}*)

 (* Then get the data for the 4-vertices*)
 temp = InsFields[{{C[particles[[1]], AntiField[particles[[2]]],
      AntiField[FieldToInsert[1]], FieldToInsert[1]]}, {Internal[1] -> FieldToInsert[1], External[1] -> particles[[1]], External[2] -> particles[[2]], Index[1] -> gO1, Index[2] -> gO2, Index[3] -> gI1}}];
 If[temp =!= {},
  For[k = 1, k <= Length[temp], k++,
   temp2 = Join[
      temp2, {ReleaseHold[
        Hold[{AntiField[Internal[1]], Internal[1],
           Cp[External[1][{Index[1]}], AntiField[External[2][{Index[2]}]],
             AntiField[Internal[1][{Index[3]}]], Internal[1][{Index[3]}]],
            Cp[External[1][{Index[1]}],
            AntiField[External[2][{Index[2]}]],
            AntiField[Internal[1][{Index[3]}]], Internal[1][{Index[3]}]],
            VType[getType[Internal[1]], getType[Internal[1]],
            getType[External[1]], getType[External[2]]],
           CalculateColorFactor[getBlank[External[1]], Internal[1],
            AntiField[Internal[1]]],
           2 CalculateSymmetryFactor[Internal[1],
            AntiField[Internal[1]]]}] /. temp[[k, 2]]]}];
   ];
  ];

 loopCorrection = 0;
 For[k = 1, k <= Length[temp2], k++,
  part = temp2[[k, 1]];
  correction = temp2[[k]];

  temploop = AddLoopVectorDiff[particles, correction];
  temploop = temploop /. {Mass[conj[x_]] -> Mass[x], Mass[bar[x_]] -> Mass[x], Mass2[conj[x_]] -> Mass2[x], Mass2[bar[x_]] -> Mass2[x]};
  loopCorrection += temploop;
  ];


 Return[loopCorrection];
 ];
