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

BeginPackage["RGIntegrator`"];
EndPackage[];

(* options *)
{ loopOrder };

RGIntegrate::usage = "Integrates a given list of renormalization group
 equations from a starting scale to a destination scale.

 The first parameter of RGIntegrate[] is the list of beta functions.
 Each list element is again a list, where the first entry is the name
 of the parameter and the rest are the 1-loop, 2-loop, ... beta
 functions.

 The loopOrder of the integration can be given as an optional
 argument, see Example 3, (default: Automatic).

 Example 1
 =========

 betas = { {g, h a g^2} }

 RGIntegrate[betas, Q1, Q2]

 Result:

 {g[Q1] -> g[Q2] + a h g[Q2]^2 Log[Q1/Q2]}

 Example 2
 =========

 betas = { {g, h a g^2, h^2 b g^2 l^2},
           {l, h (c l^2 + d g^2), h^2 (e l^4 + f g^2 l^2)} }

 RGIntegrate[betas, Q1, Q2]

 Result:

 {g[Q1] -> g[Q2] + a*h*g[Q2]^2*Log[Q1/Q2] +
    h^2*(b*g[Q2]^2*l[Q2]^2*Log[Q1/Q2] + a^2*g[Q2]^3*Log[Q1/Q2]^2),
  l[Q1] -> l[Q2] + h*(d*g[Q2]^2*Log[Q1/Q2] + c*l[Q2]^2*Log[Q1/Q2]) +
    h^2*(f*g[Q2]^2*l[Q2]^2*Log[Q1/Q2] + e*l[Q2]^4*Log[Q1/Q2] +
      a*d*g[Q2]^3*Log[Q1/Q2]^2 + c*d*g[Q2]^2*l[Q2]*Log[Q1/Q2]^2 +
      c^2*l[Q2]^3*Log[Q1/Q2]^2)}

 Example 3
 =========

 Integrate the beta functions up to the 2nd order:

 betas = { {g, h a g^2} }

 RGIntegrate[betas, Q1, Q2, loopOrder -> 2]

 Result:

 {g[Q1] -> g[Q2] + a h g[Q2]^2 Log[Q1/Q2] + a^2 h^2 g[Q2]^3 Log[Q1/Q2]^2}
";

Begin["RGIntegrator`Private`"];

(* loop factor for expansion *)
{h};

AddScale[par_, Q_] := Rule[par, par[Q]];

MultiplyLoopFactor[{par_, betas___}, h_] :=
    Join[{par}, MapIndexed[#1 h^(First[#2])&, {betas}]];

PerformIntegrals[expr_] :=
    expr /.
    Integral[ex_, rest__] :> Distribute[Integral[Expand[ex], rest]] /.
    {
        (* Log[]^0 *)
        Integral[ex_ /; FreeQ[ex,Log], {Log[Qp_], Log[Q1_], Log[Q2_]}, ___] :>
        (ex Log[Q2/Q1]),
        (* Log[]^n *)
        Integral[ex_, {Log[Qp_], Log[Q1_], Log[Q2_]}, ___] :>
        (ex /. {Log[Qp/Q1]^n_ :> Log[Q2/Q1]^(n+1)/(n+1),
                Log[Qp/Q1] -> Log[Q2/Q1]^2/2} )
    };

IntegrateSingleRHS[{par_, betas___}, Q1_, Q2_, Qp_, addScales_, sol_] :=
    Module[{loopOrder = Length[{betas}], integrand},
           integrand = Total[{betas} /. addScales /. (sol /. Q1 -> Qp)];
           integrand = Normal @ Series[integrand, {RGIntegrator`Private`h,0,loopOrder}];
           par[Q1] -> par[Q2] + Integral[integrand, {Log[Qp], Log[Q2], Log[Q1]}]
          ];

IntegrateRHS[{}, Q1_, Q2_, sol_] := {};

IntegrateRHS[betas_List, Q1_, Q2_, sol_] :=
    Module[{Qp, addScales, ints},
           addScales = AddScale[First[#], Qp]& /@ betas;
           ints = IntegrateSingleRHS[#, Q1, Q2, Qp, addScales, sol]& /@ betas;
           PerformIntegrals /@ ints
          ];

(* integrate beta functions of order h^N,
   given a solution of the (N-1)th order beta functions *)
RGIntegrateRecursively[betas_List, Q1_, Q2_, sol_] :=
    IntegrateRHS[betas, Q1, Q2, sol];

(* solution of empty beta functions *)
RGIntegrateRecursively[{}, Q1_, Q2_] := {};

(* solution of trivial beta functions *)
RGIntegrateRecursively[betas : {{par_}, ___}, Q1_, Q2_] :=
    Rule[First[#][Q1], First[#][Q2]]& /@ betas;

RGIntegrateRecursively[betas_List, Q1_, Q2_] :=
    Module[{sol, betasLower},
           betasLower = Drop[#, -1]& /@ betas;
           sol = RGIntegrateRecursively[betasLower, Q1, Q2];
           RGIntegrateRecursively[betas, Q1, Q2, sol]
      ];

Options[RGIntegrate] = {
    loopOrder -> Automatic
};

RGIntegrate[beta_List, Q1_, Q2_, OptionsPattern[]] :=
    Module[{lbeta, lo = OptionValue[loopOrder]},
           lbeta = MultiplyLoopFactor[#, RGIntegrator`Private`h]& /@ beta;
           If[IntegerQ[lo],
              (* chop or fill beta functions with zeros *)
              lbeta = PadRight[#, lo + 1]& /@ lbeta;
             ];
           RGIntegrateRecursively[lbeta, Q1, Q2] /. RGIntegrator`Private`h -> 1
          ];

End[];
