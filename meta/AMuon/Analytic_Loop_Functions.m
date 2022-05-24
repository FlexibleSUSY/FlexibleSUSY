(* ::Package:: *)

(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

(*One Loop Functions*)
OneLoopFunctionA[0]  := 2
OneLoopFunctionA[1]  := 3/2.
OneLoopFunctionA[Infinity] := 0
OneLoopFunctionA[x_] := (2 - 9*x + 18*x^2 - 11*x^3 + 6*x^3*Log[x])/(1 - x)^4

OneLoopFunctionB[0]  := 2
OneLoopFunctionB[1]  := 1
OneLoopFunctionB[Infinity] := 0
OneLoopFunctionB[x_] := (2 (1 - 6 x + 3 x^2 + 2 x^3 - 6 x^2 Log[x]))/(1 - x)^4

OneLoopFunctionC[0]  := 3
OneLoopFunctionC[1]  := 1
OneLoopFunctionC[Infinity] := 0
OneLoopFunctionC[x_] := (3 (1 - x^2 + 2 x Log[x]))/(1 - x)^3

OneLoopFunctionD[0]  := -Infinity
OneLoopFunctionD[1]  := -(9/2.)
OneLoopFunctionD[Infinity] := 0
OneLoopFunctionD[x_] := (16 - 45 x + 36 x^2 - 7 x^3 + 6 (2 - 3 x)*Log[x])/(1 - x)^4

OneLoopFunctionE[0]  := 4
OneLoopFunctionE[1]  := 1
OneLoopFunctionE[Infinity] := 0
OneLoopFunctionE[x_] := (2 (2 + 3 x - 6 x^2 + x^3 + 6 x Log[x]))/(1 - x)^4

OneLoopFunctionF[0]  := Infinity
OneLoopFunctionF[1]  := 1
OneLoopFunctionF[Infinity] := 0
OneLoopFunctionF[x_] := (3 (-3 + 4 x - x^2 - 2 Log[x]))/(2 (1 - x)^3)

OneLoopFunctionH[0]  := 5
OneLoopFunctionH[1]  := 3
OneLoopFunctionH[Infinity] := 0
OneLoopFunctionH[x_] := (5 - 27*x + 27*x^2 - 5*x^3 + 6*(x - 3)*x^2*Log[x])/(1 - x)^4

OneLoopFunctionI[0]  := 1
OneLoopFunctionI[1]  := 1/2.
OneLoopFunctionI[Infinity] := 0
OneLoopFunctionI[x_] := (1 - x + x*Log[x])/(1 - x)^2

OneLoopFunctionJ[0]  := 7
OneLoopFunctionJ[1]  := 5
OneLoopFunctionJ[Infinity] := 0
OneLoopFunctionJ[x_] := (7 - 33 x + 57 x^2 - 31 x^3 + 6 x^2 (3 x - 1) Log[x])/(1 - x)^4

OneLoopFunctionK[0]  := 1
OneLoopFunctionK[1]  := 2/3.
OneLoopFunctionK[Infinity] := 0
OneLoopFunctionK[x_] := (1 - 4 x + 3 x^2 - 2 x^2 Log[x])/(1 - x)^3

OneLoopFunctionL[0]  := Infinity
OneLoopFunctionL[1]  := 21/2.
OneLoopFunctionL[Infinity] := 0
OneLoopFunctionL[x_] := (2 + 27*x - 54*x^2 + 25*x^3 - 6*(2 - 9*x + 6*x^2)*Log[x])/(1 - x)^4

OneLoopFunctionD[0]  := 4
OneLoopFunctionD[1]  := 3/2.
OneLoopFunctionD[Infinity] := 0
OneLoopFunctionM[x_] := (4 - 9 x + 5 x^3 + 6 (1 - 2 x) x Log[x])/(1 - x)^4

OneLoopFunctionN[0,0] := 1/2.
OneLoopFunctionN[1,1] := 1/3.
OneLoopFunctionN[x_, y_] :=
 If[x == y,
   (1 - 4 x + 3 x^2 - 2 x^2 Log[x])/ ((1 - x)^3)]],
   y/(y - 1) (x (x - y)^2 Log[x] + (x - 1) ((x - y) (y - 1) - x (x - 1) Log[x/y]))/((1 - x)^2 (x - y)^2)]

(*SSF Diagram
  - mi mass of the outgoing fermion
  - mj mass of the incoming fermion
  - mF mass of the fermion in the loop
  - mS mass of the scalar in the loop
  - SFinleft left-chirality of the incoming fermion-fermion-scalar vertex
  - SFinright right-chirality of the incoming fermion-fermion-scalar vertex
  - SFoutleft left-chirality of the outcoming fermion-fermion-scalar vertex
  - SFoutright right-chirality of the outcoming fermion-fermion-scalar vertex
  - r ratio of loop masses which we expand over in the simplified case
  *)

SSFA2L[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r},
  r=Power[mF/mS,2];
  (*eq.16 of hep-ph/9510309 (possibly with different sign) *)
  Return[-SFoutleft[]*SFinright[]*OneLoopFunctionB[r]/12.-SFoutleft[]*SFinleft[]*mF/mj*OneLoopFunctionC[r]/3.-mi/mj*SFoutright[]*SFinleft[]*OneLoopFunctionB[r]/12.]
]
SSFA2LPV[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r,mS2},
  (*Massless emitter (in this case mS==0,mi==mj==mF) will fail *)
  (*Massless mediator (in this case mF==0) case can be calculated*)
  mS2=Power[mS,2];
  If[mi==mj&&mi==mS&&mF==0.,
    Return[-0.5*mS2*(SFinright[]*SFoutleft[]+SFinleft[]*SFoutright[])/Power[mi,2]],
    Return[2 mS^2/mj (mF (PVC[0,0,0,mj^2,0,mi^2,mF,mS,mS]+PVC[0,0,1,mj^2,0,mi^2,mF,mS,mS]+PVC[0,1,0,mj^2,0,mi^2,mF,mS,mS]) SFinleft[] SFoutleft[]-mj (PVC[0,1,0,mj^2,0,mi^2,mF,mS,mS]+PVC[0,1,1,mj^2,0,mi^2,mF,mS,mS]+PVC[0,2,0,mj^2,0,mi^2,mF,mS,mS]) SFinright[] SFoutleft[]-mi (PVC[0,0,1,mj^2,0,mi^2,mF,mS,mS]+PVC[0,0,2,mj^2,0,mi^2,mF,mS,mS]+PVC[0,1,1,mj^2,0,mi^2,mF,mS,mS]) SFinleft[] SFoutright[])]
  ]
]

SSFA2R[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r},
  r=Power[mF/mS,2];
  (*eq.15 of hep-ph/9510309 after replacement L<->R (possibly with different sign) *)
  Return[-SFoutright[]*SFinleft[]*OneLoopFunctionB[r]/12.-SFoutright[]*SFinright[]*mF/mj*OneLoopFunctionC[r]/3.-mi/mj*SFoutleft[]*SFinright[]*OneLoopFunctionB[r]/12.]
]
SSFA2RPV[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r,mS2},
  (*Massless emitter (in this case mS==0,mi==mj==mF) will fail *)
  (*Massless mediator (in this case mF==0) case can be calculated*)
  mS2=Power[mS,2];
  If[mi==mj&&mi==mS&&mF==0.,
    Return[-0.5*mS2*(SFinright[]*SFoutleft[]+SFinleft[]*SFoutright[])/Power[mi,2]],
    Return[mS^2/mj (-mi (PVC[0,0,1,mj^2,0,mi^2,mF,mS,mS]+PVC[0,0,2,mj^2,0,mi^2,mF,mS,mS]+PVC[0,1,1,mj^2,0,mi^2,mF,mS,mS]) SFinright[] SFoutleft[]-mj (PVC[0,1,0,mj^2,0,mi^2,mF,mS,mS]+PVC[0,1,1,mj^2,0,mi^2,mF,mS,mS]+PVC[0,2,0,mj^2,0,mi^2,mF,mS,mS]) SFinleft[] SFoutright[]+mF (PVC[0,0,0,mj^2,0,mi^2,mF,mS,mS]+PVC[0,0,1,mj^2,0,mi^2,mF,mS,mS]+PVC[0,1,0,mj^2,0,mi^2,mF,mS,mS]) SFinright[] SFoutright[])]
  ]
]

(*FFS Diagram
  - mi mass of the outgoing fermion
  - mj mass of the incoming fermion
  - mF mass of the fermion in the loop
  - mS mass of the scalar in the loop
  - SFinleft left-chirality of the incoming fermion-fermion-scalar vertex
  - SFinright right-chirality of the incoming fermion-fermion-scalar vertex
  - SFoutleft left-chirality of the outcoming fermion-fermion-scalar vertex
  - SFoutright right-chirality of the outcoming fermion-fermion-scalar vertex
  - r ratio of loop masses which we expand over in the simplified case
  *)

FSA2L[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r},
  r=Power[mF/mS,2];
  (*eq.19 of hep-ph/9510309 (possibly with different sign)*)
  Return[-SFoutleft[]*SFinright[]*OneLoopFunctionE[r]/12.-SFoutleft[]*SFinleft[]*mF/mj*OneLoopFunctionF[r]*2./3.-mi/mj*SFoutright[]*SFinleft[]*OneLoopFunctionE[r]/12.]
]
FFSA2LPV[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r,mS2},
  (*Massless emitter (in this case mF==0) will fail*)
  (*Massless mediator (in this case mS==0) case can be calculated*)*
  mS2=Power[mS,2];
  If[mi==mj&&mi==mF&&mS==0,
    Return[-0.5*mS2*(SFinright[]*SFoutleft[]+SFinleft[]*(4.*SFoutleft[]+SFoutright[]))/Power[mi,2]],
    Return[2 mS2/mj (mF (-PVC[0,0,1,mj^2,0,mi^2,mS,mF,mF]-PVC[0,1,0,mj^2,0,mi^2,mS,mF,mF]) SFinleft[] SFoutleft[]-mj (PVC[0,1,0,mj^2,0,mi^2,mS,mF,mF]+PVC[0,1,1,mj^2,0,mi^2,mS,mF,mF]+PVC[0,2,0,mj^2,0,mi^2,mS,mF,mF]) SFinright[] SFoutleft[]-mi (PVC[0,0,1,mj^2,0,mi^2,mS,mF,mF]+PVC[0,0,2,mj^2,0,mi^2,mS,mF,mF]+PVC[0,1,1,mj^2,0,mi^2,mS,mF,mF]) SFinleft[] SFoutright[])];
  ]
]

FFSA2R[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r},
  r=Power[mF/mS,2];
  (*eq.19 of hep-ph/9510309after replacement L<->R (possibly with different sign) *)
  Return[-SFoutright[]*SFinleft[]*OneLoopFunctionE[r]/12.-SFoutright[]*SFinright[]*mF/mj*OneLoopFunctionF[r]*2./3.-mi/mj*SFoutleft[]*SFinright[]*OneLoopFunctionE[r]/12.]
]
FFSA2RPV[mi_,mj_,mF_,mS_,SFinleft_,SFoutleft_,SFinright_,SFoutright_]:=Module[{r,mS2},
  (*Massless emitter (in this case mF==0) will fail*)
  (*Massless mediator (in this case mS==0) case can be calculated*)
  mS2=Power[mS,2];
  If[mi==mj&&mi==mF&&mS==0,
    Return[-0.5*mS2*(SFinleft[]*SFoutright[]+SFinright[]*(SFoutleft[]+4.*SFoutright[]))/Power[mi,2]],
    Return[2 mS2/mj(-mi (PVC[0,0,1,mj^2,0,mi^2,mS,mF,mF]+PVC[0,0,2,mj^2,0,mi^2,mS,mF,mF]+PVC[0,1,1,mj^2,0,mi^2,mS,mF,mF]) SFinright[] SFoutleft[]-mj (PVC[0,1,0,mj^2,0,mi^2,mS,mF,mF]+PVC[0,1,1,mj^2,0,mi^2,mS,mF,mF]+PVC[0,2,0,mj^2,0,mi^2,mS,mF,mF]) SFinleft[] SFoutright[]+mF (-PVC[0,0,1,mj^2,0,mi^2,mS,mF,mF]-PVC[0,1,0,mj^2,0,mi^2,mS,mF,mF]) SFinright[] SFoutright[])]
  ]
]

(*VVF Diagram
  - mi mass of the outgoing fermion
  - mj mass of the incoming fermion
  - mF mass of the fermion in the loop
  - mV mass of the vector in the loop
  - VFinleft left-chirality of the incoming fermion-fermion-vector vertex
  - VFinright right-chirality of the incoming fermion-fermion-vector vertex
  - VFoutleft left-chirality of the outcoming fermion-fermion-vector vertex
  - VFoutright right-chirality of the outcoming fermion-fermion-vector vertex
  - r ratio of loop masses which we expand over in the simplified case
  - dim the number of spacetime dimensions, used for dimensional regularization, equivalent to dim = 4 - 2 Eps
  *)

VVFA2L[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r},
  r=Power[mF/mV,2];
  Return/3.*mF/mj*VFinleft[]*VFoutright[]*OneLoopFunctionK[r]+mi/mj*VFinleft[]*VFoutleft[]*OneLoopFunctionJ[r]/6.+VFinright[]*VFoutright[]*OneLoopFunctionJ[r]/6.]
]
VVFA2LPV[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r,mV2,dim},
  (*Massless emitter (in this case mV==0,mi==mj==mF) will fail*)
  (*Massless mediator (in this case mF==0) case can be calculated*)
  mV2=Power[mV,2];
  If[mi==mj&&mi==mV&&mF==0,
    2.*mV2*(VFinleft[]*VFoutleft[]+VFinright[]*VFoutright[])/Power[mi,2],
  Return[2 mV2/mj (mi (-(-4+dim) PVC[0,0,1,mj^2,0,mi^2,mF,mV,mV]-(-2+dim) PVC[0,0,2,mj^2,0,mi^2,mF,mV,mV]+PVC[0,1,0,mj^2,0,mi^2,mF,mV,mV]+2 PVC[0,1,1,mj^2,0,mi^2,mF,mV,mV]-dim PVC[0,1,1,mj^2,0,mi^2,mF,mV,mV]) VFinleft[] VFoutleft[]+mF ((-4+dim) PVC[0,0,0,mj^2,0,mi^2,mF,mV,mV]+(-1+dim) (PVC[0,0,1,mj^2,0,mi^2,mF,mV,mV]+PVC[0,1,0,mj^2,0,mi^2,mF,mV,mV])) VFinleft[] VFoutright[]+mj (PVC[0,0,1,mj^2,0,mi^2,mF,mV,mV]-(-4+dim) PVC[0,1,0,mj^2,0,mi^2,mF,mV,mV]-(-2+dim) (PVC[0,1,1,mj^2,0,mi^2,mF,mV,mV]+PVC[0,2,0,mj^2,0,mi^2,mF,mV,mV])) VFinright[] VFoutright[])]
  ]
]

VVFA2R[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r},
  r=Power[mF/mV,2];
  Return[3.*mF/mj*VFinright[]*VFoutleft[]*OneLoopFunctionK[r]+mi/mj*VFinright[]*VFoutright[]*OneLoopFunctionJ[r]/6.+VFinleft[]*VFoutleft[]*OneLoopFunctionJ[r]/6.]
]
VVFA2RPV[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r,mV2,dim},
  (*Massless emitter (in this case mV==0,mi==mj==mF) will fail*)
  (*Massless mediator (in this case mF==0) case can be calculated*)
  mV2=Power[mV,2];
  If[mi==mj&&mi==mV&&mF==0,
    2.*mV2*(VFinleft[]*VFoutleft[]+VFinright[]*VFoutright[])/Power[mi,2],
    Return[2 mV2/mj (mj (PVC[0,0,1,mj^2,0,mi^2,mF,mV,mV]-(-4+dim) PVC[0,1,0,mj^2,0,mi^2,mF,mV,mV]-(-2+dim) (PVC[0,1,1,mj^2,0,mi^2,mF,mV,mV]+PVC[0,2,0,mj^2,0,mi^2,mF,mV,mV])) VFinleft[] VFoutleft[]+mF ((-4+dim) PVC[0,0,0,mj^2,0,mi^2,mF,mV,mV]+(-1+dim) (PVC[0,0,1,mj^2,0,mi^2,mF,mV,mV]+PVC[0,1,0,mj^2,0,mi^2,mF,mV,mV])) VFinright[] VFoutleft[]+mi (-(-4+dim) PVC[0,0,1,mj^2,0,mi^2,mF,mV,mV]-(-2+dim) PVC[0,0,2,mj^2,0,mi^2,mF,mV,mV]+PVC[0,1,0,mj^2,0,mi^2,mF,mV,mV]+2 PVC[0,1,1,mj^2,0,mi^2,mF,mV,mV]- dim PVC[0,1,1,mj^2,0,mi^2,mF,mV,mV]) VFinright[] VFoutright[])]
  ]
]

(*FFV Diagram
  - mi mass of the outgoing fermion
  - mj mass of the incoming fermion
  - mF mass of the fermion in the loop
  - mV mass of the vector in the loop
  - VFinleft left-chirality of the incoming fermion-fermion-vector vertex
  - VFinright right-chirality of the incoming fermion-fermion-vector vertex
  - VFoutleft left-chirality of the outcoming fermion-fermion-vector vertex
  - VFoutright right-chirality of the outcoming fermion-fermion-vector vertex
  - r ratio of loop masses which we expand over in the simplified case
  - dim the number of spacetime dimensions, used for dimensional regularization, equivalent to dim = 4 - 2 Eps
*)

FFVA2L[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r},
  r=Power[mF/mV,2];
  Return[-4.*mF/mj*VFinleft[]*VFoutright[]*OneLoopFunctionC[r]/3.+mi/mj*VFinleft[]*VFoutleft[]*OneLoopFunctionM[r]/3.+VFinright[]*VFoutright[]*OneLoopFunctionM[r]/3.]
]
FFVA2LPV[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r,mV2,dim},
  mV2=Power[mV,2];
  Return[2 mV2/mj (-mi (2 PVC[0,0,0,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,0,1,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,0,2,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,0,2,mj^2,0,mi^2,mV,mF,mF]+2 PVC[0,1,0,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]) VFinleft[] VFoutleft[]+mF (4 PVC[0,0,0,mj^2,0,mi^2,mV,mF,mF]+dim (PVC[0,0,1,mj^2,0,mi^2,mV,mF,mF]+PVC[0,1,0,mj^2,0,mi^2,mV,mF,mF])) VFinleft[] VFoutright[]-mj (2 PVC[0,0,0,mj^2,0,mi^2,mV,mF,mF]+2 PVC[0,0,1,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,1,0,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,2,0,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,2,0,mj^2,0,mi^2,mV,mF,mF]) VFinright[] VFoutright[])]
]

FFVA2R[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r},
  r=Power[mF/mV,2];
  Return[-4.*mF/mj*VFinright[]*VFoutleft[]*OneLoopFunctionC[r]/3.+mi/mj*VFinright[]*VFoutright[]*OneLoopFunctionM[r]/3.+VFinleft[]*VFoutleft[]*OneLoopFunctionM[r]/3.]
]
FFVA2RPV[mi_,mj_,mF_,mV_,VFinleft_,VFoutleft_,VFinright_,VFoutright_]:=Module[{r,mV2,dim},
  mV2=Power[mV,2];
  Return[2 mV2/mj (-mj (2 PVC[0,0,0,mj^2,0,mi^2,mV,mF,mF]+2 PVC[0,0,1,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,1,0,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,2,0,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,2,0,mj^2,0,mi^2,mV,mF,mF]) VFinleft[] VFoutleft[]+mF (4 PVC[0,0,0,mj^2,0,mi^2,mV,mF,mF]+dim (PVC[0,0,1,mj^2,0,mi^2,mV,mF,mF]+PVC[0,1,0,mj^2,0,mi^2,mV,mF,mF])) VFinright[] VFoutleft[]-mi (2 PVC[0,0,0,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,0,1,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,0,2,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,0,2,mj^2,0,mi^2,mV,mF,mF]+2 PVC[0,1,0,mj^2,0,mi^2,mV,mF,mF]-2 PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]+dim PVC[0,1,1,mj^2,0,mi^2,mV,mF,mF]) VFinright[] VFoutright[])]
]

(*VSF Diagram
  - mi mass of the outgoing fermion
  - mj mass of the incoming fermion
  - mF mass of the fermion in the loop
  - mV mass of the vector in the loop
  - mG mass of the Goldstone boson (corresponding to the vector) in the loop
  - VFinleft left-chirality of the incoming fermion-fermion-vector vertex
  - VFinright right-chirality of the incoming fermion-fermion-vector vertex
  - SFoutleft left-chirality of the outcoming fermion-fermion-scalar vertex
  - SFoutright right-chirality of the outcoming fermion-fermion-scalar vertex
  - r ratio of loop masses mF and mV which we expand over in the simplified case
  - s ratio of loop masses mF and mG which we expand over in the simplified case
  *)

VSFA2L[mi_,mj_,mF_,mV_,mG_,VFinleft_,SFoutleft_,VFinright_,SFoutright_]:=Module[{r,s},
  r=Power[mF/mV,2];
  s=Power[mF/mG,2];
  Return[1./mj*VFinleft[]*SFoutleft[]*OneLoopFunctionN[r,s]]
]
VSFA2LPV[mi_,mj_,mF_,mV_,mG_,VFinleft_,SFoutleft_,VFinright_,SFoutright_]:=Module[{r,mV2},
  mV2=Power[mV,2];
  If[mi==mj&&mi==mF&&mV==0&&mG==0,
    mV2,
    Return[(2 mV2 PVC[0,1,0,mj^2,0,mi^2,mF,mV,mG] SFoutleft[] VFinleft[])/mj]
  ]
]

VSFA2R[mi_,mj_,mF_,mV_,mG_,VFinleft_,SFoutleft_,VFinright_,SFoutright_]:=Module[{r,s},
  r=Power[mF/mV,2];
  s=Power[mF/mG,2];
  Return[1./mj*VFinright[]*SFoutright[]*OneLoopFunctionN[r,s]]
]
VSFA2RPV[mi_,mj_,mF_,mV_,mG_,VFinleft_,SFoutleft_,VFinright_,SFoutright_]:=Module[{r,mV2},
  mV2=Power[mV,2];
  If[mi==mj&&mi==mF&&mV==0&&mG==0,
    mV2,
    Return[(2 mV2 PVC[0,1,0,mj^2,0,mi^2,mF,mV,mG] SFoutright[] VFinright[])/mj]
  ]
]

(*SVF Diagram
  - mi mass of the outgoing fermion
  - mj mass of the incoming fermion
  - mF mass of the fermion in the loop
  - mV mass of the vector in the loop
  - mG mass of the Goldstone boson (corresponding to the vector) in the loop
  - SFinleft left-chirality of the incoming fermion-fermion-scalar vertex
  - SFinright right-chirality of the incoming fermion-fermion-scalar vertex
  - VFoutleft left-chirality of the outcoming fermion-fermion-vector vertex
  - VFoutright right-chirality of the outcoming fermion-fermion-vector vertex
  - r ratio of loop masses mF and mV which we expand over in the simplified case
  - s ratio of loop masses mF and mG which we expand over in the simplified case
  *)

SVFA2L[mi_,mj_,mF_,mV_,mG_,SFinleft_,VFoutleft_,SFinright_,VFoutright_]:=Module[{r,s},
  r=Power[mF/mV,2];
  s=Power[mF/mG,2];
  Return[1./mj*SFinleft[]*VFoutright[]*OneLoopFunctionN[r,s]]
]
SVFA2LPV[mi_,mj_,mF_,mV_,mG_,SFinleft_,VFoutleft_,SFinright_,VFoutright_]:=Module[{r,mV2},
  (*Massless emitter (in this case mV==0,mi==mj==mF) will fail*)
  (*Massless mediator (in this case mF==0) case can be calculated*)
  mV2=Power[mV,2];
  If[mi==mj&&mi==mV&&mF==0,
    mV2,
    Return[(2 mV2 PVC[0,0,1,mj^2,0,mi^2,mF,mG,mV] SFinleft[] VFoutright[])/mj]
  ]
]

SVFA2R[mi_,mj_,mF_,mV_,mG_,SFinleft_,VFoutleft_,SFinright_,VFoutright_]:=Module[{r,s},
  r=Power[mF/mV,2];
  s=Power[mF/mG,2];
  Return[1./mj*SFinright[]*VFoutleft[]*OneLoopFunctionN[r,s]]
]
SVFA2RPV[mi_,mj_,mF_,mV_,mG_,SFinleft_,VFoutleft_,SFinright_,VFoutright_]:=Module[{r,mV2},
  (*Massless emitter (in this case mV==0,mi==mj==mF) will fail*)
  (*Massless mediator (in this case mF==0) case can be calculated*)
  mV2=Power[mV,2];
  If[mi==mj&&mi==mV&&mF==0,
    mV2,
    Return[(2 mV^2 PVC[0,0,1,mj^2,0,mi^2,mF,mG,mV] SFinright[] VFoutleft[])/mj]
  ]
]
