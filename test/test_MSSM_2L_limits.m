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

(* Script for converting MSSM 2L Higgs corrections from Fortran to
   Mathematica form and derivation of the limits s2t -> 0 and mst1 ->
   mst2.
 *)

repl = {t -> mt^2, g -> mg^2, q -> Q^2, T1 -> mst1^2, T2 -> mst2^2, T -> T2, phi -> phiExpr};
A = Xt - mu/tanb; (* notice sign convention for mu! *)
Xt = (T1 - T2) s2t/(2 mt);
del = g^2 + t^2 + T1^2 - 2 (g t + g T1 + t T1);

(* phi[x,y,z] function from Pietro *)
phiExpr[x_, y_, z_] :=
    Module[{u = x/z, v = y/z, lambda, xp, xm},
           lambda = Sqrt[(1 - u - v)^2 - 4 u v];
           xp = 1/2 (1 + (u - v) - lambda);
           xm = 1/2 (1 - (u - v) - lambda);
           1/lambda (2 Log[xp] Log[xm] - Log[u] Log[v] - 
                     2 (PolyLog[2, xp] + PolyLog[2, xm]) + Pi^2/3)
          ];

(*****************************)
(* CP-even tadpoles O(at*as) *)
(*****************************)

(* Fortran expressions from Pietro Slavich *)

F2lStr = "4*mg*mt*(1+4*c2t**2)/s2t
       $     -(2*(T1-T2)+4*mg*mt/s2t)*Log(g/q)*Log(t/q)
       $     -2*(4-s2t**2)*(T1-T2)
       $     +(4*T1*T2-s2t**2*(T1+T2)**2)/(T1-T2)*Log(T1/q)*Log(T2/q)
       $     + strF2lc(t,mg,T1,T2,s2t,c2t,q)
       $     - strF2lc(t,mg,T2,T1,-s2t,c2t,q)";

G2lStr = "5*mg/mt*s2t*(T1-T2)-10*(T1+T2-2*t)-4*g
       $     + 12*t*(Log(t/q)**2-2*Log(t/q))
       $     +(4*g-s2t*mg/mt*(T1-T2))*Log(g/q)*Log(t/q)
       $     +s2t**2*(T1+T2)*Log(T1/q)*Log(T2/q)
       $     + strG2lc(t,mg,T1,T2,s2t,q)
       $     + strG2lc(t,mg,T2,T1,-s2t,q)";

strF2lcStr = "(4*(g+t+2*T1)-s2t**2*(3*T1+T2)-4*s2t*mg*mt
       $     -16*mg*mt*T1*c2t**2/s2t/(T1-T2))*Log(T1/q)       
       $     +T1/(T1-T2)*(s2t**2*(T1+T2)-2*(2*T1-T2))*Log(T1/q)**2
       $     +2*(T1-g-t+mg*mt*s2t
       $     +2*c2t**2*mg*mt*T1/s2t/(T1-T2))*Log(g*t/q**2)*Log(T1/q)
       $     +4*mg*mt*c2t**2*(t-g)/s2t/(T1-T2)*Log(t/g)*Log(T1/q)
       $     +((2*del+4*g*t)/T1-2*mg*mt*s2t/T1*(g+t-T1)
       $     +4*c2t**2*mg*mt/T1/(T1-T2)/s2t*del)*phi(g,t,T1)";

strG2lcStr = "(4*(g+t+2*T1)+s2t**2*(T1-T2)
       $     -4*mg/mt*s2t*(t+T1))*Log(T1/q) 
       $     +(mg/mt*s2t*(5*t-g+T1)-2*(g+2*t))*Log(t/q)*Log(T1/q)
       $     +(mg/mt*s2t*(g-t+T1)-2*g)*Log(g/q)*Log(T1/q)
       $     -(2+s2t**2)*T1*Log(T1/q)**2
       $     +(2*g/T1*(g+t-T1-2*mg*mt*s2t)
       $     +mg/mt*s2t*del/T1)*phi(g,t,T1)";

(**********************************************)
(* CP-even dMh O(at*as) in the limit s2t -> 0 *)
(**********************************************)

(* Fortran expressions from Pietro Slavich *)

S12Str = "ht**2*mt*mu*F2s";

S22Str = "2 * ht**2 * mt**2 * F1 +
       $        2 * ht**2 * mt * A * F2s";

F1Str = "strF1ab(t,T1,T2,s2t,c2t,q) 
       $     + strF1c(t,mg,T1,s2t,q)
       $     + strF1c(t,mg,T2,-s2t,q)";

strF1abStr = "-6*(1-Log(t/q))+5*Log(T1*T2/t**2)+Log(T1*T2/t**2)**2
       $     +8*Log(t/q)**2-4*Log(T1/q)**2-4*Log(T2/q)**2
       $     -c2t**2*(2-Log(T1/q)-Log(T2/q)-Log(T1/T2)**2)
       $     -s2t**2*(T1/T2*(1-Log(T1/q))+T2/T1*(1-Log(T2/q)))";

strF1cStr = "+4*(t+g-mg*mt*s2t)/T1*(1-Log(g/q))
       $     +4*Log(t/g) - 2*Log(T1/g)
       $     +2/del*(4*g**2*Log(T1/g)
       $     +(g**2-T1**2+t*(10*g+3*t+2*t*g/T1-2*t**2/T1))*Log(t/g))
       $     +2*mg/mt*s2t*(Log(T1/q)**2+2*Log(t/q)*Log(T1/q))
       $     +4*mg/mt*s2t/del*(g*(T1-t-g)
       $     *Log(T1/g)+t*(T1-3*g-2*t-(t*g-t**2)/T1)*Log(t/g))
       $     +(4*g*(t+g-T1-2*mg*mt*s2t)/del
       $     -4*mg/mt*s2t)*phi(t,T1,g)";

F2sStr = "-8*mg*mt/(T1-T2)*(
       $     (Log(T1/q)-Log(t/q)*Log(T1/q)+phi(t,T1,g))-
       $     (Log(T2/q)-Log(t/q)*Log(T2/q)+phi(t,T2,g)))";

exprStr = ("(" <> # <> ")") &/@ {F2lStr, G2lStr, strF2lcStr, strG2lcStr,
                                 S12Str, S22Str, F1Str, strF1abStr, strF1cStr, F2sStr};

NoBracket[c_String] := !MemberQ[{"(", ")"}, c];

MakeFunc[head_String] :=
    head <> "(" ~~ q__?NoBracket ~~ ")" :> head <> "[" <> q <> "]";

exprInputForm = \
    StringReplace[
        StringReplace[exprStr, {"$" -> "", "**" -> "^"}],
        {MakeFunc["Log"],
         MakeFunc["strF2lc"],
         MakeFunc["strG2lc"],
         MakeFunc["phi"],
         MakeFunc["strF1ab"],
         MakeFunc["strF1c"]}
    ];

expr = ToExpression /@ exprInputForm;

F2l = expr[[1]];
G2l = expr[[2]];
strF2lc[t_, mg_, T1_, T2_, s2t_, c2t_, q_] := Evaluate[expr[[3]]];
strG2lc[t_, mg_, T1_, T2_, s2t_, q_] := Evaluate[expr[[4]]];

S12 = expr[[5]];
S22 = expr[[6]];
F1 = expr[[7]];
strF1ab[t_, T1_, T2_, s2t_, c2t_, q_] := Evaluate[expr[[8]]];
strF1c[t_, mg_, T1_, s2t_, q_] := Evaluate[expr[[9]]];
F2s = expr[[10]];

k = 4 gs^2/(16 Pi^2)^2;

S1 = mt*mu/tanb*s2t*F2l;
S1 = S1/v1^2;

S2 = mt*A*s2t*F2l + 2*mt^2*G2l;
S2 = S2/v2^2;

S1 = k*S1;
S2 = k*S2;

S11 = 0;
S12 = k*S12;
S22 = k*S22;

(* CP-even tadpoles O(at*as) *)
tadpoles = {S1, S2};

(* CP-even dMh O(at*as) in the limit s2t -> 0 *)
selfEnergy = {S11, S12, S22} /. s2t -> 0 /. c2t -> 1;

(* Limits *)

(* s2t -> 0 *)
tadpolesS2t0 =
    Simplify[Limit[tadpoles, s2t -> 0] /. c2t -> 1 /. Xt -> 0];

(* s2t -> 0 and mst1 -> mst2 *)
tadpolesS2t0T1eqT2 =
    Simplify[Normal[
             Series[tadpolesS2t0 /. T1 -> T2 + eps,
                    {eps, 0, 0}]] /. T2 -> T];

(* s2t -> 0 and mst1 -> mst2 *)
selfEnergyS2t0T1eqT2 = 
    Simplify[Normal[
             Series[selfEnergy /. T1 -> T2 + eps,
                    {eps, 0, 0}]] /. T2 -> T];

gs /: gs^2 = gs2;
ht /: ht^2 = ht2;
mt /: mt^2 = mt2;
g  /: g^2  = g2;
t  /: t^2  = t2;
t  /: t^3  = t3;
T  /: T^2  = Tsqr;
T  /: T^3  = Tcub;

simp1 = {
    Log[t/g] -> ltg,
    Log[T/g] -> lTg,
    Log[T/q] -> lTq,
    Log[T1/q]-> lT1q,
    Log[T2/q]-> lT2q,
    Log[t/q] -> ltq,
    Log[g/q] -> lgq,
    Log[Tsqr/t2] -> lT2t2,
    Log[g t/q^2] -> lgtq2
};

simp2 = {
    a_^2 :> sqr[a],
    a_^-2 :> 1/sqr[a]
};

(* simplify CP-even tadpole s2t = 0 *)
pref = mt^2 gs^2 / (4 Pi)^4;
Print[CForm /@ FullSimplify[1/pref tadpolesS2t0 //. simp1] //. simp1 //. simp2];

(* simplify CP-even tadpole s2t = 0 and T1 = T2 *)
pref = mt^2 gs^2 / (4 Pi)^4;
Print[CForm /@ Simplify[1/pref tadpolesS2t0T1eqT2 //. simp1] //. simp1 //. simp2];

(* simplify CP-even self-energy *)
pref = ht^2 mt^2 gs^2 / (4 Pi)^4;
Print[CForm /@ Simplify[1/pref selfEnergyS2t0T1eqT2 //. simp1] //. simp1 //. simp2];
