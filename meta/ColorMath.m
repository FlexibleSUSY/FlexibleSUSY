(* ::Package:: *)

BeginPackage["ColorMath`", {"Utils`"}];

CM$VersionNumber = "1.0";
CM$VersionDate = "February 12 2013";
Utils`PrintHeadline["Loading ColorMath package"];

Print["Version: ", CM$VersionNumber, " (" CM$VersionDate, "), for Mathematica 7, 8 and 9."];
Print["Author: Malin Sjodahl"];
Print["For suggestions and bug reports contact malin.sjodahl@thep.lu.se."];
Print["If you use this package for research, please cite the ColorMath paper, arXiv:1211.2099."];

CSimplify::usage = "";
Nc::usage = "";
TR::usage = "";
CMf::usage = "";
CMt::usage = "";
CM\[Delta]::usage = "";
CM\[CapitalDelta]::usage = "";
CMdelta::usage = "";
CMDelta::usage = "";
CMo::usage = "";
CMd::usage = "";

Begin["Private`"];
Unprotect[Conjugate];
Conjugate[Nc]=Nc;
Conjugate[TR]=TR;
Conjugate[Nc^Expr_ ]=Nc^Expr;
Conjugate[TR^Expr_ ]=TR^Expr;
Conjugate[Sqrt[Nc] ]=Sqrt[Nc];
Conjugate[Sqrt[TR] ]=Sqrt[TR];
Conjugate[Superscript[CMf,{G1$1_,G2$1_,G3$1_}]]:=Superscript[CMf,{G1$1,G2$1,G3$1}]
Conjugate[Superscript[CMd,{G1$1_,G2$1_,G3$1_}]]:=Superscript[CMd,{G1$1,G2$1,G3$1}]
Conjugate[Subscript[(Superscript[(Superscript[CMt,G1$1s_]),Q2$1_]), Q1$1_]]:=Subscript[(Superscript[(Superscript[CMt,Reverse[G1$1s]]),Q1$1]), Q2$1];
Conjugate[Superscript[CMo,G1$1s_]]:=Superscript[CMo,Reverse[G1$1s]];
Conjugate[Superscript[CM\[CapitalDelta],{G1$1_,G2$1_}]]:=Superscript[CM\[CapitalDelta],{G1$1,G2$1}]
Conjugate[Subscript[(Superscript[CM\[Delta],Q1$1_]), Q2$1_]]:=Subscript[(Superscript[CM\[Delta],Q2$1]), Q1$1]
Protect[Conjugate ];

AllIndices[CMExpr$_] := Module[{AllInd$, LowerQuarks$, UpperQuarks$, GluInd$, 
      inters$}, AllInd$ = {}; LowerQuarks$ = LowerQuarkIndices[CMExpr$]; 
      UpperQuarks$ = UpperQuarkIndices[CMExpr$]; GluInd$ = GluonIndices[CMExpr$]; 
      If[Intersection[LowerQuarks$, GluInd$] != {}, 
       {inters$ = Intersection[LowerQuarks$, GluInd$]; AllIndices::indices = 
          "The set of indices `1` appears to be both quark and Gluon type \
indices in `2`."; Message[AllIndices::indices, inters$, CMExpr$]; }]; 
      If[Intersection[GluInd$, UpperQuarks$] != {}, 
       {inters$ = Intersection[GluInd$, UpperQuarks$]; AllIndices::indices = 
          "The set of indices `1` appears to be both quark and Gluon type \
indices in `2`."; Message[AllIndices::indices, inters$, CMExpr$]; }]; 
      AllInd$ = Union[LowerQuarks$, UpperQuarks$, GluInd$]]
 
AllIndices /: AllIndices::usage = "AllIndices[Expr] returns a list of all \
(external and dummy) indices in Expr."
 
LowerQuarkIndices[CMExpr$_] := Module[{AllInd$, NewInd$}, 
     AllInd$ = {}; NewInd$ = Reap[CMExpr$ /. 
         Subscript[Superscript[Superscript[CMt, Gs$_], Q1$1_], Q2$1_] :>
          Module[{}, Sow[Q2$1]; Subscript[Superscript[Superscript[CMt, {Gs$}],
              Q1$1], Q2$1]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      NewInd$ = Reap[CMExpr$ /. Subscript[Superscript[CM\[Delta], Q1$1_], 
           Q2$1_] :> Module[{}, Sow[Q2$1]; Subscript[Superscript[CM\[Delta], 
              Q1$1], Q2$1]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      AllInd$]
 
LowerQuarkIndices /: LowerQuarkIndices::usage = "LowerQuarkIndices[Expr] \
returns a list of all (external and dummy) quark-type indices placed \
downstairs, i.e. Q2 in \
\!\(\*SubscriptBox[TemplateBox[{\"\[Delta]\",\"Q1\"},\n\"Superscript\"], \
\"Q2\"]\), and \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", RowBox[{\"G1\", \",\", RowBox[{\"...\", \"GNg\"}]}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\)."
 
Attributes[Subscript] = {NHoldRest}
 
Attributes[Superscript] = {NHoldRest, ReadProtected}
 
CMt[Gs$_, Q1$_, Q2$_] := Subscript[Superscript[Superscript[CMt, Gs$], Q1$], Q2$]
 
CM\[Delta][Q1$_, Q2$_] := Subscript[Superscript[CM\[Delta], Q1$], Q2$]
 
UpperQuarkIndices[CMExpr$_] := Module[{AllInd$, NewInd$}, 
     AllInd$ = {}; NewInd$ = Reap[CMExpr$ /. 
         Subscript[Superscript[Superscript[CMt, Gs$_], Q1$1_], Q2$1_] :>
          Module[{}, Sow[Q1$1]; Subscript[Superscript[Superscript[CMt, {Gs$}],
              Q1$1], Q2$1]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      NewInd$ = Reap[CMExpr$ /. Subscript[Superscript[CM\[Delta], Q1$1_], 
           Q2$1_] :> Module[{}, Sow[Q1$1]; Subscript[Superscript[CM\[Delta], 
              Q1$1], Q2$1]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      AllInd$]
 
UpperQuarkIndices /: UpperQuarkIndices::usage = "UpperQuarkIndices[Expr] \
returns a list of all (external and dummy) quark-type indices placed \
upstairs, i.e. Q1 in \
\!\(\*SubscriptBox[TemplateBox[{\"\[Delta]\",\"Q1\"},\n\"Superscript\"], \
\"Q2\"]\), and \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", RowBox[{\"G1\", \",\", RowBox[{\"...\", \"GNg\"}]}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\)."
 
GluonIndices[CMExpr$_] := Module[{AllInd$, NewInd$}, 
     AllInd$ = {}; NewInd$ = Reap[CMExpr$ /. Superscript[CM\[CapitalDelta], 
           {G1$1_, G2$1_}] :> Module[{}, Sow[G1$1]; Sow[G2$1]; 
            Superscript[CM\[CapitalDelta], {G1$1, G2$1}]; ]]; 
      If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      NewInd$ = Reap[CMExpr$ /. Superscript[CMf, {G1$1_, G2$1_, G3$1_}] :>
          Module[{}, Sow[G1$1]; Sow[G2$1]; Sow[G3$1]; Superscript[CMf,
             {G1$1, G2$1, G3$1}]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      NewInd$ = Reap[CMExpr$ /. Superscript[CMd, {G1$1_, G2$1_, G3$1_}] :> 
          Module[{}, Sow[G1$1]; Sow[G2$1]; Sow[G3$1]; Superscript[CMd, 
             {G1$1, G2$1, G3$1}]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      NewInd$ = Reap[CMExpr$ /. Superscript[CMo, Gs$_] :> 
          Module[{}, Table[Sow[Gs$[[GsI$]]], {GsI$, 1, Length[Gs$]}]; 
            Superscript[CMo, Gs$]; ]]; If[Length[NewInd$[[2]]] != 0, 
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      NewInd$ = Reap[CMExpr$ /. Subscript[Superscript[Superscript[CMt, Gs$_],
            Q1$1_], Q2$1_] :> Module[{}, Table[Sow[Gs$[[gInd$]]], 
             {gInd$, 1, Length[Gs$]}]; Subscript[Superscript[Superscript[CMt, {
                Gs$}], Q1$1], Q2$1]; ]]; If[Length[NewInd$[[2]]] != 0,
       AllInd$ = DeleteDuplicates[Join[AllInd$, NewInd$[[2]][[1]]]]; ]; 
      AllInd$]
 
GluonIndices /: GluonIndices::usage = "GluonIndices[Expr] returns a list of \
all gluon indices (external and dummy) in the expression Expr."
 
CM\[CapitalDelta][G1$_, G2$_] := Superscript[CM\[CapitalDelta], {G1$, G2$}]
 
CMf[G1$_, G2$_, G3$_] := Superscript[CMf, {G1$, G2$, G3$}]
 
d[G1$_, G2$_, G3$_] := Superscript[CMd, {G1$, G2$, G3$}]
 
o[Gs$_] := Superscript[CMo, Gs$]
 
AllPairs[CMExpr$_] := Select[FindSymbols[CMExpr$], 
     Count[FindSymbols[CMExpr$], #1] == 2 & ]
 
AllPairs /: AllPairs::usage = "AllPairs[Expr] takes a list as argument and \
finds all the symbols that appear twice (including, for example, Times)."
 
FindSymbols[CMExpr$_] := Module[{symbs$ = {}}, 
     CMExpr$ /. cfac$_Symbol :> (symbs$ = Prepend[symbs$, cfac$]; cfac$); 
      symbs$]
 
FindSymbols /: FindSymbols::usage = "FindSymbols[Expr] returns a list of all \
symbols used in an expression. If they appear more than once, they will be \
listed more than once."
 
AllPermutations[Superscript[TCM\[CapitalDelta]$_, {G1$1_, G2$1_}]] := 
    Superscript[TCM\[CapitalDelta]$, {G1$1, G2$1}] | 
     Superscript[TCM\[CapitalDelta]$, {G2$1, G1$1}]
 
AllPermutations[Superscript[Tfd$_, {G1$1_, G2$1_, G3$1_}]] := 
    Superscript[Tfd$, {G1$1, G2$1, G3$1}] | Superscript[Tfd$, 
      {G2$1, G3$1, G1$1}] | Superscript[Tfd$, {G3$1, G1$1, G2$1}] | 
     Superscript[Tfd$, {G1$1, G3$1, G2$1}] | Superscript[Tfd$, 
      {G2$1, G1$1, G3$1}] | Superscript[Tfd$, {G3$1, G2$1, G1$1}]
 
AllSimpleRules := Union[SimpleRules, OTSimpleRules, Remove0To2ORules]
 
AllSimpleRules /: AllSimpleRules::usage = "AllSimpleRules is the set of all \
rules involving \
(\!\(\*TemplateBox[{\"\[Delta]\",\"Q2\"},\n\"Superscript\"]\)\!\(\*SubscriptB\
ox[\()\), \(Q1\)]\), \!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\"}], \"}\"}]},\n\"Superscript\"]\), \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \
\",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\) which do not increase \
the number of terms, i.e. the union of SimpleRules, Remove0To2ORules and \
OTSimpleRules."
 
SimpleRules = {Subscript[Superscript[CM\[Delta], Q1$1_], q2$1_]*
       Subscript[Superscript[Superscript[CMt, G1$1s_], q2$1_], Q3$1_] ->
      Subscript[Superscript[Superscript[CMt, G1$1s], Q1$1], Q3$1],
     Subscript[Superscript[CM\[Delta], Q1$1_], Q2$1_]*
       Subscript[Superscript[Superscript[CMt, G1$1s_], Q3$1_], Q1$1_] ->
      Subscript[Superscript[Superscript[CMt, G1$1s], Q3$1], Q2$1],
     Subscript[Superscript[CM\[Delta], Q1$1_], q2$1_]*
       Subscript[Superscript[CM\[Delta], q2$1_], Q3$1_] -> 
      Subscript[Superscript[CM\[Delta], Q1$1], Q3$1], 
     Subscript[Superscript[CM\[Delta], q1$1_], q1$1_] -> Nc, 
     Subscript[Superscript[Superscript[CMt, {G1$1_}], q1$1_], q2$1_]*
       Subscript[Superscript[Superscript[CMt, {G2$1_}], q2$1_], q1$1_] ->
      TR*Superscript[CM\[CapitalDelta], {G1$1, G2$1}], 
     Subscript[Superscript[Superscript[CMt, {G1$1_}], q1$1_], q1$1_] -> 0,
     Subscript[Superscript[Superscript[CMt, G1$1s_], q1$1_], q2$1_]*
       Subscript[Superscript[Superscript[CMt, G2$1s_], q2$1_], q1$1_] ->
      Superscript[CMo, Join[G1$1s, G2$1s]], 
     Subscript[Superscript[Superscript[CMt, G1$1s_], q1$1_], q1$1_] ->
      Superscript[CMo, G1$1s], 
     Subscript[Superscript[Superscript[CMt, G1$1s_], Q1$1_], q1$1_]*
       Subscript[Superscript[Superscript[CMt, G2$1s_], q1$1_], Q2$1_] ->
      Subscript[Superscript[Superscript[CMt, Join[G1$1s, G2$1s]], Q1$1], Q2$1],
     Subscript[Superscript[Superscript[CMt, {}], Q1$1_], Q2$1_] ->
      Subscript[Superscript[CM\[Delta], Q1$1], Q2$1], 
     Superscript[CMo, {G1$1_, G2$1_}] -> TR*Superscript[CM\[CapitalDelta], 
        {G1$1, G2$1}], (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
        Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
       (Superscript[CMd, {g2$1_, G3$1_, G4$1_}] | Superscript[CMd, 
         {G3$1_, G4$1_, g2$1_}] | Superscript[CMd, {G4$1_, g2$1_, G3$1_}] | 
        Superscript[CMd, {g2$1_, G4$1_, G3$1_}] | Superscript[CMd, 
         {G3$1_, g2$1_, G4$1_}] | Superscript[CMd, {G4$1_, G3$1_, g2$1_}]) -> 
      Superscript[CMd, {G1$1, G3$1, G4$1}], 
     (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
        Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
       (Superscript[CMf, {g2$1_, G3$1_, G4$1_}] | Superscript[CMf,
         {G3$1_, G4$1_, g2$1_}] | Superscript[CMf, {G4$1_, g2$1_, G3$1_}]) ->
      Superscript[CMf, {G1$1, G3$1, G4$1}],
     (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
        Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
       (Superscript[CM\[CapitalDelta], {g2$1_, G2$1_}] | 
        Superscript[CM\[CapitalDelta], {G2$1_, g2$1_}]) -> 
      Superscript[CM\[CapitalDelta], {G1$1, G2$1}], 
     Superscript[CM\[CapitalDelta], {G1$1_, G1$1_}] -> -1 + Nc^2, 
     (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
         Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
        Subscript[Superscript[Superscript[CMt, G3$1s_], Q1$1_], Q2$1_] /;
       Count[G3$1s, g2$1] == 1 :> Subscript[Superscript[
        Superscript[CMt, G3$1s /. g2$1 -> G1$1], Q1$1], Q2$1],
     (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
         Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
        Superscript[CMo, G3$1s_] /; Count[G3$1s, g2$1] == 1 :> 
      Superscript[CMo, G3$1s /. g2$1 -> G1$1], 
     Superscript[CM\[CapitalDelta], {G1$1_, G2$1_}]^2 -> -1 + Nc^2}
 
SimpleRules /: SimpleRules::usage = "Basic rules for quark and gluon \
contraction. These rules, which involve \
\!\(\*SubscriptBox[TemplateBox[{\"\[Delta]\",\"Q1\"},\n\"Superscript\"], \
\"Q2\"]\), \!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\"}], \"}\"}]},\n\"Superscript\"]\) or quark \
contraction, never increase the number of terms."
 
OTSimpleRules = {Superscript[CMo, Gs$_ /; NMaxGluonCheck[Gs$] == 2] :> 
      Module[{Place1$, Place2$, Newgs1$, CMres$, Multiplicities$, theg$, 
        muli$}, Multiplicities$ = Transpose[Tally[Gs$]][[2]]; 
        For[muli$ = 1, muli$ <= Length[Multiplicities$], muli$++, 
         If[Multiplicities$[[muli$]] == 2, {theg$ = Gs$[[muli$]]; 
            Break[]; }]]; Place1$ = Position[Gs$, theg$][[1]][[1]]; 
        Place2$ = Position[Gs$, theg$][[2]][[1]]; 
        If[Place2$ - Place1$ == 1 || (Place1$ == 1 && Place2$ == 
            Length[Gs$]), {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
           If[ !Length[Union[Gs$]] - Length[Union[Newgs1$]] == 1, 
            {Print["OTSimpleRules Warning "]; Print["OTSimpleRules Old ", 
               Superscript[CMo, Gs$]]; Print["OTSimpleRules New ", Superscript[
                CMo, Newgs1$]]; }]; CMres$ = (TR/Nc)*(Nc^2 - 1)*Superscript[CMo, 
               Newgs1$] /. Remove0To2ORules; }, 
         If[Place2$ - Place1$ == 2 || (Place1$ == 1 && Place2$ == 
              Length[Gs$] - 1) || (Place1$ == 2 && Place2$ == Length[Gs$]), 
           {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
             If[ !Length[Union[Gs$]] - Length[Union[Newgs1$]] == 1, 
              {Print["Warning "]; Print["Old ", Superscript[CMo, Gs$]]; 
                Print["New ", Superscript[CMo, Newgs1$]]; }]; 
             CMres$ = (-(TR/Nc))*Superscript[CMo, Newgs1$] /. 
               Remove0To2ORules; }, {CMres$ = Superscript[CMo, Gs$]; }]; ]; 
        CMres$], Subscript[Superscript[Superscript[CMt,
         Gs$_ /; NMaxGluonCheck[Gs$] == 2], Q1$1_], Q2$1_] :> 
      Module[{place$, Newgs1$, CMres$, Multiplicities$, NNei$, NNNei$}, 
       Multiplicities$ = Transpose[Tally[Gs$]][[2]]; 
        CMres$ = Subscript[Superscript[Superscript[CMt, Gs$], Q1$1], Q2$1];
        Newgs1$ = Gs$; NNei$ = 0; For[place$ = 1, place$ < Length[Newgs1$], 
         place$++, If[Newgs1$[[place$]] == Newgs1$[[place$ + 1]], 
           {Newgs1$ = Drop[Newgs1$, {place$, place$ + 1}]; NNei$++; 
             CMres$ = Subscript[Superscript[Superscript[CMt, Newgs1$], Q1$1],
                Q2$1] /. Remove0ORules; place$ = place$ - 1; }]; ]; 
        CMres$ = CMres$*((TR/Nc)*(Nc^2 - 1))^NNei$; NNNei$ = 0; 
        For[place$ = 1, place$ < Length[Newgs1$] - 1, place$++, 
         If[Newgs1$[[place$]] == Newgs1$[[place$ + 2]], 
           {NNNei$++; Newgs1$ = Newgs1$ = Drop[Drop[Newgs1$, {place$}], 
                {place$ + 1}]; CMres$ = Subscript[Superscript[Superscript[CMt,
                   Newgs1$], Q1$1], Q2$1]*((TR/Nc)*(Nc^2 - 1))^NNei$ /. 
               Remove0ORules; place$ = place$ - 1; }]; ]; 
        CMres$ = CMres$*(-(TR/Nc))^NNNei$; CMres$ /. Remove0ORules], 
     Superscript[CMo, Gs$_]^2 :> OSquare[Length[Gs$]]}
 
OTSimpleRules /: OTSimpleRules::usage = "Rules for contracting neighboring \
and next to neighboring gluons in closed and open quark-lines, \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \
\",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\)."
 
NMaxGluonCheck[Gs$_] := Module[{CMres$}, If[Length[Gs$] == 0, CMres$ = 0]; 
      If[Length[Gs$] > 0, CMres$ = Max[Transpose[Tally[Gs$]][[2]]]]; CMres$]
 
NMaxGluonCheck /: NMaxGluonCheck::usage = "Little helper function to see if a \
list contains two gluons of the same kind."
 
Remove0To2ORules = {Superscript[CMo, LIndV$_ /; Length[LIndV$] == 1] -> 0, 
     Superscript[CMo, LIndV$_ /; Length[LIndV$] == 0] -> Nc, 
     Superscript[CMo, {G1$1_, G2$1_}] -> TR*Superscript[CM\[CapitalDelta], 
        {G1$1, G2$1}]}
 
Remove0To2ORules /: Remove0To2ORules::usage = "Rules for simplifying closed \
quark-lines with 0 to 2 gluons, \!\(\*TemplateBox[{\"o\",RowBox[{\"{\", \
\"}\"}]},\n\"Superscript\"]\)=Nc, \!\(\*TemplateBox[{\"o\",RowBox[{\"{\", \
\"G1\", \"}\"}]},\n\"Superscript\"]\)=0, \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\"}], \
\"}\"}]},\n\"Superscript\"]\)=TR \
\!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \
\"G2\"}], \"}\"}]},\n\"Superscript\"]\)"
 
Remove0ORules = {Subscript[Superscript[Superscript[CMt,
         LIndV$_ /; Length[LIndV$] == 0], Q1$1_], Q2$1_] -> 
      Subscript[Superscript[CM\[Delta], Q1$1], Q2$1]}
 
Remove0ORules /: Remove0ORules::usage = "Rule for simplifying closed \
quark-lines with 0 gluons,\!\(\*TemplateBox[{RowBox[{\" \", \
\"o\"}],RowBox[{\"{\", \"}\"}]},\n\"Superscript\"]\)=Nc."
 
OSquare[Nglu$_] := OSquare[Nglu$] = Module[{Gs$, ii$}, 
      Gs$ = Table[Unique[CMd], {ii$, Nglu$}]; 
       Factor[Expand[(Superscript[CMo, Gs$] /. OTToTRules)*
            (Superscript[CMo, Gs$] /. OTToTRules) //. OTGluonRules] //. 
         SimpleRules]]
 
OSquare /: OSquare::usage = "OSquare[Ng] calculates \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"g1\", \",\", \
RowBox[{\"...\", \"gNg\"}]}], \"}\"}]},\n\"Superscript\"]\)^2."
 
OTToTRules = {Superscript[CMo, Gs_] :> Module[{qVecInd$, nInd$, qInd$}, 
       nInd$ = Length[Gs]; qVecInd$ = Table[Unique[qInd$], 
          {uind$, nInd$ + 1}]; Product[Subscript[Superscript[
            Superscript[CMt, {Gs[[gInd$]]}], qVecInd$[[gInd$]]],
           qVecInd$[[gInd$ + 1]]], {gInd$, Length[Gs]}] /. 
         qVecInd$[[nInd$ + 1]] -> qVecInd$[[1]]], 
     Subscript[Superscript[Superscript[CMt, Gs_], Q1$1_], Q2$1_] :>
      Module[{qVecInd$, nInd$, qInd$}, nInd$ = Length[Gs]; 
        qVecInd$ = Table[Unique[qInd$], {uind$, nInd$ + 1}]; 
        Product[Subscript[Superscript[Superscript[CMt, {Gs[[gInd$]]}],
            qVecInd$[[gInd$]]], qVecInd$[[gInd$ + 1]]], 
          {gInd$, Length[Gs]}] /. {qVecInd$[[1]] -> Q1$1, 
          qVecInd$[[nInd$ + 1]] -> Q2$1}]}
 
OTToTRules /: OTToTRules::usage = "Rules for replacing open and closed \
quark-lines (\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\) and \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \
\",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\)) with products of \
SU(\!\(\*SubscriptBox[\(N\), \(c\)]\)) generators \
(\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \"G1\", \
\"}\"}]}, \"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\))."
 
OTGluonRules = {Superscript[CMo, g1$1s_]*Superscript[CMo, g2$1s_] /; 
       Intersection[g1$1s, g2$1s] != {}*TextCell[""]*TextCell[""]*
         TextCell[""] :> Module[{theg$, Place1$, Place2$, g1$1s1, g1$1s2, 
        g2$1s1, g2$1s2, FirstL$, SecondL$, NewL$, NewL1$, NewL2$, Supterm$}, 
       If[Order[g1$1s[[1]], g2$1s[[1]]] == -1, 
         {FirstL$ = g2$1s; SecondL$ = g1$1s}, 
         {FirstL$ = g1$1s; SecondL$ = g2$1s}]; 
        theg$ = Intersection[g1$1s, g2$1s][[1]]; 
        Place1$ = Position[FirstL$, theg$][[1,1]]; 
        Place2$ = Position[SecondL$, theg$][[1,1]]; 
        g1$1s1 = Take[FirstL$, Place1$ - 1]; g2$1s1 = Take[SecondL$, 
          Place2$ - 1]; g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewL$ = {}; 
        NewL$ = Join[g1$1s1, g2$1s2, g2$1s1, g1$1s2]; 
        NewL1$ = Drop[g1$1s, {Position[g1$1s, theg$][[1,1]]}]; 
        NewL2$ = Drop[g2$1s, {Position[g2$1s, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Superscript[CMo, NewL1$]*Superscript[CMo, NewL2$]; 
        TR*(Superscript[CMo, NewL$] - Supterm$) /. Remove0To2ORules], 
     Subscript[Superscript[Superscript[CMt, g1$1s_], Q1$1_], Q2$1_]*
        Superscript[CMo, g2$1s_] /; Intersection[g1$1s, g2$1s] != 
        {}*TextCell[""]*TextCell[""]*TextCell[""] :> 
      Module[{theg$, Place1$, Place2$, g1$1s1, g1$1s2, g2$1s1, g2$1s2, NewL$, 
        NewL1$, NewL2$, SecondL$, FirstL$, Supterm$}, 
       FirstL$ = g1$1s; SecondL$ = g2$1s; theg$ = Intersection[g1$1s, g2$1s][[
          1]]; Place1$ = Position[FirstL$, theg$][[1]][[1]]; 
        Place2$ = Position[SecondL$, theg$][[1,1]]; 
        g1$1s1 = Take[FirstL$, Place1$ - 1]; g2$1s1 = Take[SecondL$, 
          Place2$ - 1]; g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewL$ = {}; 
        NewL$ = Join[g1$1s1, g2$1s2, g2$1s1, g1$1s2]; 
        NewL1$ = Drop[g1$1s, {Position[g1$1s, theg$][[1,1]]}]; 
        NewL2$ = Drop[g2$1s, {Position[g2$1s, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Subscript[Superscript[Superscript[CMt, NewL1$],
            Q1$1], Q2$1]*Superscript[CMo, NewL2$]; 
        TR*(Subscript[Superscript[Superscript[CMt, NewL$], Q1$1], Q2$1] -
           Supterm$) /. Union[Remove0To2ORules, Remove0ORules]], 
     Subscript[Superscript[Superscript[CMt, g1$1s_], Q1$1_], Q2$1_]*
        Subscript[Superscript[Superscript[CMt, g2$1s_], Q3$1_], Q4$1_] /;
       Intersection[g1$1s, g2$1s] != {}*TextCell[""]*TextCell[""]*
         TextCell[""] :> Module[{theg$, Place1$, Place2$, g1$1s1, g1$1s2, 
        g2$1s1, g2$1s2, NewLa$, NewLb$, NewL1$, NewL2$, SecondL$, FirstL$, 
        Supterm$}, FirstL$ = g1$1s; SecondL$ = g2$1s; 
        theg$ = Intersection[g1$1s, g2$1s][[1]]; 
        Place1$ = Position[FirstL$, theg$][[1,1]]; 
        Place2$ = Position[SecondL$, theg$][[1,1]]; 
        g1$1s1 = Take[FirstL$, Place1$ - 1]; g2$1s1 = Take[SecondL$, 
          Place2$ - 1]; g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewLa$ = {}; 
        NewLb$ = {}; NewLa$ = Join[g1$1s1, g2$1s2]; 
        NewLb$ = Join[g2$1s1, g1$1s2]; NewL1$ = Drop[g1$1s, 
          {Position[g1$1s, theg$][[1,1]]}]; NewL2$ = 
         Drop[g2$1s, {Position[g2$1s, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Subscript[Superscript[Superscript[CMt, NewL1$],
            Q1$1], Q2$1]*Subscript[Superscript[Superscript[CMt, NewL2$], Q3$1],
           Q4$1]; TR*(Subscript[Superscript[Superscript[CMt, NewLa$], Q1$1],
             Q4$1]*Subscript[Superscript[Superscript[CMt, NewLb$], Q3$1],
             Q2$1] - Supterm$) /. Union[Remove0To2ORules, Remove0ORules]], 
     Subscript[Superscript[Superscript[CMt, Gs$_ /;
          Max[Transpose[Tally[Gs$]][[2]]] == 2], Q1$1_], Q2$1_] :> 
      Module[{muli$, theg$, Place1$, Place2$, Newgs1$, Newgs2$, CMres$, rem$, 
        Multiplicities$}, Multiplicities$ = Transpose[Tally[Gs$]][[2]]; 
        For[muli$ = 1, muli$ <= Length[Multiplicities$], muli$++, 
         If[Multiplicities$[[muli$]] == 2, {theg$ = Gs$[[muli$]], Break[]}]]; 
        Place1$ = Position[Gs$, theg$][[1]][[1]]; 
        Place2$ = Position[Gs$, theg$][[2]][[1]]; If[Place2$ - Place1$ == 1, 
         {Newgs1$ = Drop[Gs$, {Place1$, Place2$}]; 
           CMres$ = (TR/Nc)*(Nc^2 - 1)*Subscript[Superscript[Superscript[CMt,
                 Newgs1$], Q1$1], Q2$1] /. Remove0ORules; }, 
         If[Place2$ - Place1$ == 2, {Newgs1$ = Drop[Drop[Gs$, {Place2$}], 
              {Place1$}]; If[ !Length[Union[Gs$]] - Length[Union[Newgs1$]] == 
               1, {OTGluonRules::check = 
                "The old quark-line was `1` and the new is `2`."; Message[
                OTGluonRules::check, Subscript[Superscript[Superscript[CMt,
                   Gs$], Q1$1], Q2$1], Subscript[Superscript[Superscript[CMt,
                   Newgs1$], Q1$1], Q2$1]]; }]; CMres$ = 
             (-(TR/Nc))*Subscript[Superscript[Superscript[CMt, Newgs1$], Q1$1],
                Q2$1] /. Remove0ORules; }, 
          {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
            Newgs2$ = Drop[Gs$, {Place1$, Place2$}]; 
            rem$ = Take[Gs$, {Place1$ + 1, Place2$ - 1}]; 
            CMres$ = TR*(Subscript[Superscript[Superscript[CMt, Newgs2$], Q1$1],
                  Q2$1]*Superscript[CMo, rem$] - (1/Nc)*Subscript[Superscript[
                   Superscript[CMt, Newgs1$], Q1$1], Q2$1]) /.
              Remove0ORules; }]]; CMres$], 
     Superscript[CMo, Gs$_ /; Max[Transpose[Tally[Gs$]][[2]]] == 2] :> 
      Module[{muli$, theg$, Place1$, Place2$, Newgs1$, Newgs2$, CMres$, rem$, 
        Multiplicities$}, Multiplicities$ = Transpose[Tally[Gs$]][[2]]; 
        For[muli$ = 1, muli$ <= Length[Multiplicities$], muli$++, 
         If[Multiplicities$[[muli$]] == 2, {theg$ = Gs$[[muli$]], Break[]}]]; 
        Place1$ = Position[Gs$, theg$][[1]][[1]]; 
        Place2$ = Position[Gs$, theg$][[2]][[1]]; 
        If[Place2$ - Place1$ == 1 || (Place1$ == 1 && Place2$ == 
            Length[Gs$]), {Newgs1$ = Delete[Delete[Gs$, Place2$], Place1$]; 
           CMres$ = (TR/Nc)*(Nc^2 - 1)*Superscript[CMo, Newgs1$] /. 
             {Superscript[CMo, Gs1$_ /; Length[Gs1$] == 1] -> 0, 
              Superscript[CMo, Gs1$_ /; Length[Gs1$] == 0] -> Nc}; }, 
         If[Place2$ - Place1$ == 2, {Newgs1$ = Drop[Drop[Gs$, {Place2$}], 
              {Place1$}]; If[ !Length[Union[Gs$]] - Length[Union[Newgs1$]] == 
               1, {Print["Warning "]; Print["Old ", Superscript[CMo, Gs$]]; 
               Print["New ", Superscript[CMo, Newgs1$]]; }]; 
            CMres$ = (-(TR/Nc))*Superscript[CMo, Newgs1$] /. 
              {Superscript[CMo, Gs1$_ /; Length[Gs1$] == 1] -> 0, Superscript[
                 o, Gs1$_ /; Length[Gs1$] == 0] -> Nc}; }, 
          {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
            Newgs2$ = Drop[Gs$, {Place1$, Place2$}]; 
            rem$ = Take[Gs$, {Place1$ + 1, Place2$ - 1}]; 
            CMres$ = TR*(Superscript[CMo, Newgs2$]*Superscript[CMo, rem$] - 
                (1/Nc)*Superscript[CMo, Newgs1$]) /. Remove0To2ORules; }]]; 
        CMres$], Subscript[Superscript[Superscript[CMt, {g1$1_}], Q1$1_], Q2$1_]*
       Subscript[Superscript[Superscript[CMt, {g1$1_}], Q3$1_], Q4$1_] ->
      TR*(Subscript[Superscript[CM\[Delta], Q1$1], Q4$1]*
         Subscript[Superscript[CM\[Delta], Q3$1], Q2$1] - 
        (Subscript[Superscript[CM\[Delta], Q1$1], Q2$1]*
          Subscript[Superscript[CM\[Delta], Q3$1], Q4$1])/Nc)}
 
OTGluonRules /: OTGluonRules::usage = "Rules for contracting repeated gluon \
indices in \!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \
\"...\", \",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\), including the Fierz \
identity."
 
BasisType /: BasisType::usage = "Option for CGamma with default \
BasisType\[Rule] GeneralBasis, should be BasisType\[Rule]OrthonormalBasis, \
BasisType\[Rule]OrthogonalBasis or BasisType\[Rule]TraceBasis."
 
BasisWarning = "You are using the general versions of CGamma. This involves \
calculating and inverting the scalar product matrix (unless supplied) which \
is tedious for large matrices. If your basis is either:\n1) Orthonormal\n2) \
Orthogonal\n3) A trace basis where each vector is one product of closed and \
open quark-lines\nthe computation may be speeded up significantly if CGamma \
is used with one of the options:\n1) BasisType-> OrthonormalBasis\n2) \
BasisType-> OrthogonalBasis\n3) BasisType-> TraceBasis "
 
CalledByOther /: CalledByOther::usage = 
     "Option for CGamma, not intended for user."
 
CDot[Cv1$_, OptionsPattern[]] := Module[{CVectors$, viDotvi$, CMres$, CBInd$}, 
     If[Head[Cv1$] == List, {CVectors$ = Cv1$; CMres$ = {}; 
         For[CBInd$ = 1, CBInd$ <= Length[Cv1$], CBInd$++, 
          {viDotvi$ = PowerExpand[CDot[CVectors$[[CBInd$]], CVectors$[[
                CBInd$]], NcMin -> OptionValue[NcMin]]]; 
            If[OptionValue[Verbose], Print[CBInd$, " ", CVectors$[[CBInd$]], 
              " ", viDotvi$]]; CMres$ = Append[CMres$, viDotvi$]; }]; CMres$}]; 
      If[Head[Cv1$] =!= List, {CMres$ = PowerExpand[CDot[Cv1$, Cv1$, 
           NcMin -> OptionValue[NcMin]]]}]; CMres$]
 
CDot[Cv1$_, Cv2$_, OptionsPattern[]] := 
    Module[{NInd$, CDotIndV$, CMExpr$, Bothfd$, C1fdC2o$, C1oC2fd$}, 
     NInd$ = NTensorIndices[Cv1$]; CDotIndV$ = Table[Unique[d], 
        {uInd$, NInd$}]; CMExpr$ = Conjugate[Subscript[Cv1$, CDotIndV$]]*
        Subscript[Cv2$, CDotIndV$]; Bothfd$ = 
       ContainsFD[Subscript[Cv1$, CDotIndV$]] && ContainsFD[
         Subscript[Cv2$, CDotIndV$]]; C1fdC2o$ = 
       ContainsFD[Subscript[Cv1$, CDotIndV$]] && 
         !ContainsFD[Subscript[Cv2$, CDotIndV$]]; 
      C1oC2fd$ =  !ContainsFD[Subscript[Cv1$, CDotIndV$]] && 
        ContainsFD[Subscript[Cv2$, CDotIndV$]]; 
      If[Bothfd$, {CMExpr$ = (Expand //. AllSimpleRules)[
           Expand[CMExpr$] //. FDRules //. FDRules //. FDToORules]; }, 
       {If[C1fdC2o$, {CMExpr$ = CSimplify[Expand[CSimplify[Conjugate[
                   Subscript[Cv1$, CDotIndV$]] //. FDToORules]*Subscript[
                 Cv2$, CDotIndV$]]]; }]; If[C1oC2fd$, 
          {CMExpr$ = CSimplify[Expand[CSimplify[Subscript[Cv2$, CDotIndV$] //. 
                  FDToORules]*Conjugate[Subscript[Cv1$, CDotIndV$]]]]; }]}]; 
      CMExpr$ = CSimplify[CSimplify[CMExpr$]]; 
      CMExpr$ = PowerExpand[CMExpr$, Assumptions -> Nc >= OptionValue[NcMin]]; 
      Simplify[Simplify[PowerExpand[CMExpr$, Assumptions -> 
          Nc >= OptionValue[NcMin]], Assumptions -> 
         {Nc >= OptionValue[NcMin]}], Assumptions -> {Element[Nc, Reals], 
         Nc >= OptionValue[NcMin]}]]
 
Options[CDot] = {Verbose -> False, NcMin -> 3}
 
CDot /: CDot::usage = "CDot[C1,C2] calculates the scalar product between two \
color tensors C1 and C2. Before using CDot the color tensors \
\!\(\*SubscriptBox[\(C1\), \({I1, \(\(...\) \(Ik\)\)}\)]\) and  \
\!\(\*SubscriptBox[\(C2\), \({I1, \(\(...\) \(Ik\)\)}\)]\) should have been \
defined using pattern matching underscores and SetDelayed,:=. For improved \
simplification of roots the option NcMin\[Rule]6 (or any other value) may be \
used to assume that Nc is at least 6 when using Simplify. If CDot is called \
with one vector as argument, it returns the square of the vector. There is \
also a version of CDot which takes a List of basis vectors as argument and \
returns a List containing the corresponding squared norms, CDot[{C1,...Ck}]. \
For this version, the option Verbose-> True may be used to turn on progress \
information."
 
NTensorIndices[Tens$_] := Module[{NIndCand$, NInd$, MaxInd$, IndTab$}, 
     MaxInd$ = 20; NInd$ = 0; For[NIndCand$ = 1, NIndCand$ <= MaxInd$, 
       NIndCand$++, {IndTab$ = Table[iii$, {iii$, 1, NIndCand$}]; 
         If[ContainsColor[Subscript[Tens$, IndTab$]], 
          {If[NInd$ != 0, NTensorIndices::check = "Warning: There seems to be \
more than one tensor named `1` defined, one with `2` indices and one with \
`3`, indices."; Message[NTensorIndices::check, Tens$, NInd$, NIndCand$]; ]; 
            If[NInd$ == 0, NInd$ = NIndCand$]; }]; }]; 
      If[NInd$ == 0, {NTensorIndices::check = "There seems to be no tensor \
\!\(\*SubscriptBox[\(Tens\), \({\(I1 ... \) IN}\)]\) named `1` defined, or, \
the tensor or has more than `2` indices. In the latter case, consider \
increasing MaxInd$ in the function NTensorIndices."; 
         Message[NTensorIndices::check, Tens$, MaxInd$]; }]; NInd$]
 
NTensorIndices /: NTensorIndices::usage = "NTensorIndices[Tens] finds the \
number of free color indices in the tensor \!\(\*SubscriptBox[\(Tens\), \
\({\(I1 ... \) Ik}\)]\)."
 
ContainsColor[CMExpr$_] := Module[{TrueOrFalse$}, TrueOrFalse$ = False; 
      If[ContainsFD[CMExpr$] || ContainsO[CMExpr$] || ContainsT[CMExpr$] || 
        ContainsQuarkDelta[CMExpr$] || ContainsGluonDelta[CMExpr$], 
       TrueOrFalse$ = True]; TrueOrFalse$]
 
ContainsColor /: ContainsColor::usage = "ContainsColor[Expr] returns True if \
the expression Expr contains any of the color structure objects used in the \
ColorMath package, \
\!\(\*SubscriptBox[TemplateBox[{\"\[Delta]\",\"Q1\"},\n\"Superscript\"], \
\"Q2\"]\), \!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\"}], \"}\"}]},\n\"Superscript\"]\), \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \
\",\", \"GNg\"}], \"}\"}]},\n\"Superscript\"]\), \
\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"...\", \",\", \"GNg\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\), \
\!\(\*TemplateBox[{\"f\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) or \
\!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\), and False otherwise."
 
ContainsFD[CMExpr$_] := Module[{TrueOrFalse$}, TrueOrFalse$ = False; 
      If[ !(FreeQ[CMExpr$, Superscript[CMf, {g1$1_, g2$1_, g3$1_}]] &&
         FreeQ[CMExpr$, Superscript[CMd, {g1$1_, g2$1_, g3$1_}]]), 
       TrueOrFalse$ = True]; TrueOrFalse$]
 
ContainsFD /: ContainsFD::usage = "ContainsFD[Expr] returns True if the \
expression Expr contains structure constants, \
\!\(\*TemplateBox[{\"f\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) or \
\!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\), and False otherwise."
 
ContainsO[CMExpr$_] := Module[{TrueOrFalse$}, TrueOrFalse$ = False; 
      If[ !FreeQ[CMExpr$, Superscript[CMo, Gs$_]], TrueOrFalse$ = True]; 
      TrueOrFalse$]
 
ContainsO /: ContainsO::usage = "ContainsO[Expr] returns True if the \
expression Expr contains closed quark-lines \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \
\",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\), and False otherwise."
 
ContainsT[CMExpr$_] := Module[{TrueOrFalse$}, TrueOrFalse$ = False; 
      If[ !FreeQ[CMExpr$, Subscript[Superscript[Superscript[CMt, Gs$_], Q1$1_],
          Q2$1_]], TrueOrFalse$ = True]; TrueOrFalse$]
 
ContainsT /: ContainsT::usage = "ContainsT[Expr] returns True if the \
expression Expr contains open quark-lines, \
\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1$1\"},\n\"Superscript\"], \"Q2$1\"]\), and False \
otherwise."
 
ContainsQuarkDelta[CMExpr$_] := Module[{TrueOrFalse$}, 
     TrueOrFalse$ = False; 
      If[ !FreeQ[CMExpr$, Subscript[Superscript[CM\[Delta], Q1$1_], Q2$1_]], 
       TrueOrFalse$ = True]; TrueOrFalse$]
 
ContainsQuarkDelta /: ContainsQuarkDelta::usage = "ContainsQuarkDelta[Expr] \
returns True if the expression Expr contains a quark delta function, \
\!\(\*SubscriptBox[TemplateBox[{\"\[Delta]\",\"Q1\"},\n\"Superscript\"], \
\"Q2\"]\), and False otherwise."
 
ContainsGluonDelta[CMExpr$_] := Module[{TrueOrFalse$}, 
     TrueOrFalse$ = False; If[ !FreeQ[CMExpr$, Superscript[CM\[CapitalDelta], 
          {G1$1_, G2$1_}]], TrueOrFalse$ = True]; TrueOrFalse$]
 
ContainsGluonDelta /: ContainsGluonDelta::usage = "ContainsGluonDelta[Expr] \
returns True if the expression Expr contains a gluon delta function, \
\!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \
\"G2\"}], \"}\"}]},\n\"Superscript\"]\), and False otherwise."
 
FDRules = {Superscript[CMd, {g1$1_, g2$1_, g3$1_}]^2 -> 
      (2*(4 - 5*Nc^2 + Nc^4)*TR)/Nc, Superscript[CMf, {g1$1_, g2$1_, g3$1_}]^
       2 -> 2*Nc*(-1 + Nc^2)*TR, (Superscript[CMd, {g1$1_, g2$1_, G1$1_}] | 
        Superscript[CMd, {g2$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g2$1_}] | Superscript[CMd, {g1$1_, G1$1_, g2$1_}] | 
        Superscript[CMd, {g2$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g2$1_, g1$1_}])*(Superscript[CMd, {g1$1_, g2$1_, G2$1_}] | 
        Superscript[CMd, {g2$1_, G2$1_, g1$1_}] | Superscript[CMd, 
         {G2$1_, g1$1_, g2$1_}] | Superscript[CMd, {g1$1_, G2$1_, g2$1_}] | 
        Superscript[CMd, {g2$1_, g1$1_, G2$1_}] | Superscript[CMd, 
         {G2$1_, g2$1_, g1$1_}]) -> 
      (2*(-4 + Nc^2)*TR*Superscript[CM\[CapitalDelta], {G1$1, G2$1}])/Nc, 
     (Superscript[CMf, {g1$1_, g2$1_, G1$1_}] | Superscript[CMf,
         {g2$1_, G1$1_, g1$1_}] | Superscript[CMf, {G1$1_, g1$1_, g2$1_}])*
       (Superscript[CMf, {g1$1_, g2$1_, G2$1_}] | Superscript[CMf,
         {g2$1_, G2$1_, g1$1_}] | Superscript[CMf, {G2$1_, g1$1_, g2$1_}]) ->
      2*Nc*TR*Superscript[CM\[CapitalDelta], {G1$1, G2$1}], 
     (Superscript[CMf, {g1$1_, g2$1_, G1$1_}] | Superscript[CMf,
         {g2$1_, G1$1_, g1$1_}] | Superscript[CMf, {G1$1_, g1$1_, g2$1_}])*
       (Superscript[CMf, {g2$1_, g1$1_, G2$1_}] | Superscript[CMf,
         {g1$1_, G2$1_, g2$1_}] | Superscript[CMf, {G2$1_, g2$1_, g1$1_}]) ->
      -2*Nc*TR*Superscript[CM\[CapitalDelta], {G1$1, G2$1}], 
     (Superscript[CMd, {g1$1_, g2$1_, G2$1_}] | Superscript[CMd, 
         {g2$1_, G2$1_, g1$1_}] | Superscript[CMd, {G2$1_, g1$1_, g2$1_}] | 
        Superscript[CMd, {g1$1_, G2$1_, g2$1_}] | Superscript[CMd, 
         {g2$1_, g1$1_, G2$1_}] | Superscript[CMd, {G2$1_, g2$1_, g1$1_}])*
       (Superscript[CMf, {g1$1_, g2$1_, G1$1_}] | Superscript[CMf,
         {g2$1_, G1$1_, g1$1_}] | Superscript[CMf, {G1$1_, g1$1_, g2$1_}] |
        Superscript[CMf, {g1$1_, G1$1_, g2$1_}] | Superscript[CMf,
         {g2$1_, g1$1_, G1$1_}] | Superscript[CMf, {G1$1_, g2$1_, g1$1_}]) ->
      0, (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
        Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
       (Superscript[CMd, {g2$1_, G3$1_, G4$1_}] | Superscript[CMd, 
         {G3$1_, G4$1_, g2$1_}] | Superscript[CMd, {G4$1_, g2$1_, G3$1_}] | 
        Superscript[CMd, {g2$1_, G4$1_, G3$1_}] | Superscript[CMd, 
         {G3$1_, g2$1_, G4$1_}] | Superscript[CMd, {G4$1_, G3$1_, g2$1_}]) -> 
      Superscript[CMd, {G1$1, G3$1, G4$1}], 
     (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
        Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
       (Superscript[CMf, {g2$1_, G3$1_, G4$1_}] | Superscript[CMf,
         {G3$1_, G4$1_, g2$1_}] | Superscript[CMf, {G4$1_, g2$1_, G3$1_}]) ->
      Superscript[CMf, {G1$1, G3$1, G4$1}],
     (Superscript[CM\[CapitalDelta], {G1$1_, g2$1_}] | 
        Superscript[CM\[CapitalDelta], {g2$1_, G1$1_}])*
       (Superscript[CM\[CapitalDelta], {g2$1_, G2$1_}] | 
        Superscript[CM\[CapitalDelta], {G2$1_, g2$1_}]) -> 
      Superscript[CM\[CapitalDelta], {G1$1, G2$1}], 
     (Superscript[CMf, {g1$1_, G2$1_, g2$1_}] | Superscript[CMf,
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMf, {g2$1_, g1$1_, G2$1_}])*
       (Superscript[CMf, {g2$1_, G3$1_, g3$1_}] | Superscript[CMf,
         {G3$1_, g3$1_, g2$1_}] | Superscript[CMf, {g3$1_, g2$1_, G3$1_}])*
       (Superscript[CMf, {g3$1_, G1$1_, g1$1_}] | Superscript[CMf,
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMf, {g1$1_, g3$1_, G1$1_}]) ->
      -(Nc*TR*Superscript[CMf, {G1$1, G2$1, G3$1}]),
     (Superscript[CMf, {g1$1_, G2$1_, g2$1_}] | Superscript[CMf,
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMf, {g2$1_, g1$1_, G2$1_}])*
       (Superscript[CMf, {g2$1_, G3$1_, g3$1_}] | Superscript[CMf,
         {G3$1_, g3$1_, g2$1_}] | Superscript[CMf, {g3$1_, g2$1_, G3$1_}])*
       (Superscript[CMf, {g3$1_, g1$1_, G1$1_}] | Superscript[CMf,
         {g1$1_, G1$1_, g3$1_}] | Superscript[CMf, {G1$1_, g3$1_, g1$1_}]) ->
      Nc*TR*Superscript[CMf, {G1$1, G2$1, G3$1}],
     (Superscript[CMf, {g1$1_, G2$1_, g2$1_}] | Superscript[CMf,
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMf, {g2$1_, g1$1_, G2$1_}])*
       (Superscript[CMf, {g2$1_, G3$1_, g3$1_}] | Superscript[CMf,
         {G3$1_, g3$1_, g2$1_}] | Superscript[CMf, {g3$1_, g2$1_, G3$1_}])*
       (Superscript[CMd, {g3$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMd, {g1$1_, g3$1_, G1$1_}] | 
        Superscript[CMd, {g3$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g3$1_, g1$1_}] | Superscript[CMd, {g1$1_, G1$1_, g3$1_}]) -> 
      -(Nc*TR*Superscript[CMd, {G1$1, G2$1, G3$1}]), 
     (Superscript[CMf, {g1$1_, g2$1_, G2$1_}] | Superscript[CMf,
         {g2$1_, G2$1_, g1$1_}] | Superscript[CMf, {G2$1_, g1$1_, g2$1_}])*
       (Superscript[CMf, {g2$1_, G3$1_, g3$1_}] | Superscript[CMf,
         {G3$1_, g3$1_, g2$1_}] | Superscript[CMf, {g3$1_, g2$1_, G3$1_}])*
       (Superscript[CMd, {g3$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMd, {g1$1_, g3$1_, G1$1_}] | 
        Superscript[CMd, {g3$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g3$1_, g1$1_}] | Superscript[CMd, {g1$1_, G1$1_, g3$1_}]) -> 
      Nc*TR*Superscript[CMd, {G1$1, G2$1, G3$1}], 
     (Superscript[CMf, {g1$1_, G2$1_, g2$1_}] | Superscript[CMf,
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMf, {g2$1_, g1$1_, G2$1_}])*
       (Superscript[CMf, {g2$1_, g3$1_, G3$1_}] | Superscript[CMf,
         {g3$1_, G3$1_, g2$1_}] | Superscript[CMf, {G3$1_, g2$1_, g3$1_}])*
       (Superscript[CMd, {g3$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMd, {g1$1_, g3$1_, G1$1_}] | 
        Superscript[CMd, {g3$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g3$1_, g1$1_}] | Superscript[CMd, {g1$1_, G1$1_, g3$1_}]) -> 
      Nc*TR*Superscript[CMd, {G1$1, G2$1, G3$1}], 
     (Superscript[CMf, {g2$1_, G3$1_, g3$1_}] | Superscript[CMf,
         {G3$1_, g3$1_, g2$1_}] | Superscript[CMf, {g3$1_, g2$1_, G3$1_}])*
       (Superscript[CMd, {g1$1_, G2$1_, g2$1_}] | Superscript[CMd, 
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMd, {g2$1_, g1$1_, G2$1_}] | 
        Superscript[CMd, {g1$1_, g2$1_, G2$1_}] | Superscript[CMd, 
         {G2$1_, g1$1_, g2$1_}] | Superscript[CMd, {g2$1_, G2$1_, g1$1_}])*
       (Superscript[CMd, {g3$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMd, {g1$1_, g3$1_, G1$1_}] | 
        Superscript[CMd, {g3$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g3$1_, g1$1_}] | Superscript[CMd, {g1$1_, G1$1_, g3$1_}]) -> 
      ((-4 + Nc^2)*TR*Superscript[CMf, {G1$1, G2$1, G3$1}])/Nc,
     (Superscript[CMf, {g2$1_, g3$1_, G3$1_}] | Superscript[CMf,
         {g3$1_, G3$1_, g2$1_}] | Superscript[CMf, {G3$1_, g2$1_, g3$1_}])*
       (Superscript[CMd, {g1$1_, G2$1_, g2$1_}] | Superscript[CMd, 
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMd, {g2$1_, g1$1_, G2$1_}] | 
        Superscript[CMd, {g1$1_, g2$1_, G2$1_}] | Superscript[CMd, 
         {G2$1_, g1$1_, g2$1_}] | Superscript[CMd, {g2$1_, G2$1_, g1$1_}])*
       (Superscript[CMd, {g3$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMd, {g1$1_, g3$1_, G1$1_}] | 
        Superscript[CMd, {g3$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g3$1_, g1$1_}] | Superscript[CMd, {g1$1_, G1$1_, g3$1_}]) -> 
      -(((-4 + Nc^2)*TR*Superscript[CMf, {G1$1, G2$1, G3$1}])/Nc),
     (Superscript[CMd, {g1$1_, G2$1_, g2$1_}] | Superscript[CMd, 
         {G2$1_, g2$1_, g1$1_}] | Superscript[CMd, {g2$1_, g1$1_, G2$1_}] | 
        Superscript[CMd, {g1$1_, g2$1_, G2$1_}] | Superscript[CMd, 
         {G2$1_, g1$1_, g2$1_}] | Superscript[CMd, {g2$1_, G2$1_, g1$1_}])*
       (Superscript[CMd, {g2$1_, G3$1_, g3$1_}] | Superscript[CMd, 
         {G3$1_, g3$1_, g2$1_}] | Superscript[CMd, {g3$1_, g2$1_, G3$1_}] | 
        Superscript[CMd, {g2$1_, g3$1_, G3$1_}] | Superscript[CMd, 
         {G3$1_, g2$1_, g3$1_}] | Superscript[CMd, {g3$1_, G3$1_, g2$1_}])*
       (Superscript[CMd, {g3$1_, G1$1_, g1$1_}] | Superscript[CMd, 
         {G1$1_, g1$1_, g3$1_}] | Superscript[CMd, {g1$1_, g3$1_, G1$1_}] | 
        Superscript[CMd, {g3$1_, g1$1_, G1$1_}] | Superscript[CMd, 
         {G1$1_, g3$1_, g1$1_}] | Superscript[CMd, {g1$1_, G1$1_, g3$1_}]) -> 
      ((-12 + Nc^2)*TR*Superscript[CMd, {G1$1, G2$1, G3$1}])/Nc}
 
FDRules /: FDRules::usage = "Rules for gluon contraction for terms involving \
up to three structure constants, \!\(\*TemplateBox[{\"f\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\", \",\", \"G3\"}], \"}\"}]},\n\"Superscript\"]\) \
or \!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \
\",\", \"G3\"}], \"}\"}]},\n\"Superscript\"]\). The result after contraction \
is expressed in terms of structure constants and gluon deltas, \
\!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \
\"G2\"}], \"}\"}]},\n\"Superscript\"]\)."
 
FDToORules = {Superscript[CMf, {G1$1_, G2$1_, G3$1_}] :>
      (1/(TR*I))*(Superscript[CMo, {G1$1, G2$1, G3$1}] - 
        Superscript[CMo, {G1$1, G3$1, G2$1}]), 
     Superscript[CMd, {G1$1_, G2$1_, G3$1_}] :> 
      (1/TR)*(Superscript[CMo, {G1$1, G2$1, G3$1}] + Superscript[CMo, 
         {G1$1, G3$1, G2$1}])}
 
FDToORules /: FDToORules::usage = "Rules for replacing structure constants, \
\!\(\*TemplateBox[{\"f\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\), with sums of closed quark-lines \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\)."
 
CSimplify[CMExpr$_, OptionsPattern[]] := Module[{CMres$ = CMExpr$}, 
     CMres$ = Expand[CMres$]; If[ContainsFD[CMres$], 
       {CMres$ = CMres$ //. FDRules; If[OptionValue[RemoveFD] && 
            ContainsFD[CMres$], {CMres$ = CExpand[Expand[SortIndices[CExpand[
                   CSimplify[Expand[CMres$ //. FDToORules]]]]]]; 
              CMres$ = CMres$ //. FDToORules; }; ]; }; ]; 
      CMres$ = CMres$ //. AllSimpleRules //. OTThenAllSimpleRules; 
      If[AllIndices[CMres$] != {}, CMres$ = CMres$ //. ExpandThenRules]; 
      If[AllIndices[CMres$] != {}, CMres$ = CMres$ //. ExpandThenRules]; 
      CMres$ = CExpand[SortIndices[CMres$]]]
 
Options[CSimplify] = {RemoveFD -> True}
 
CSimplify /: CSimplify::usage = "CSimplify[Expr] is the most general function \
for simplifying color structure. If Expr contains structure constants, \
\!\(\*TemplateBox[{\"f\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) or \
\!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\), FDRules are first applied. If, after \
this, the expression still contains structure constants, they are -- by \
default -- removed using FDToORules, and all repeated indices are \
subsequently contracted by repeatedly using AllSimpleRules, then \
OTThenAllSimpleRules, and finally ExpandThenRules. For color structure \
containing structure constants, it could happen that it is desired not to \
replace the structure constants with quark-lines, as doing so could expand \
the expression. This can be achieved by setting the option RemoveFD -> \
False."
 
CExpand[CMExpr$_] := Collect[Expand[CMExpr$], 
     {Superscript[CM\[CapitalDelta], G1$1s_]*Superscript[CM\[CapitalDelta], 
        G2$1s_]*Superscript[CM\[CapitalDelta], G3$1s_], 
      Superscript[CM\[CapitalDelta], G1$1s_]*Superscript[CM\[CapitalDelta], 
        G2$1s_], Superscript[CMo, G1$1s_]*Superscript[CM\[CapitalDelta], G2$1s_], 
      Superscript[CMo, G1$1s_]*Superscript[CMo, G2$1s_], Superscript[CMo, G1$1s_]}, 
     Simplify]
 
SortIndices[CMExpr$_] := Module[{}, CMExpr$ /. 
      {Superscript[CMd, {G1$1_, G2$1_, G3$1_}] :> Superscript[CMd, 
         Sort[{G1$1, G2$1, G3$1}]], Superscript[CMf, {G1$1_, G2$1_, G3$1_}] :>
        Signature[{G1$1, G2$1, G3$1}]*Superscript[CMf,
          Sort[{G1$1, G2$1, G3$1}]], Superscript[CM\[CapitalDelta], 
         {G1$1_, G2$1_}] :> Superscript[CM\[CapitalDelta], Sort[{G1$1, G2$1}]], 
       Superscript[CMo, Gs$_] :> Module[{FirstInd$, FirstPlace$, GsNew$}, 
         FirstInd$ = Sort[Gs$][[1]]; Table[If[Gs$[[GsI$]] == FirstInd$, 
            FirstPlace$ = GsI$], {GsI$, 1, Length[Gs$]}]; 
          GsNew$ = Join[Take[Gs$, -(Length[Gs$] - FirstPlace$ + 1)], 
            Take[Gs$, FirstPlace$ - 1]]; Superscript[CMo, GsNew$]]}]
 
SortIndices /: SortIndices::usage = "SortIndices[Expr] sorts the gluon \
indices appearing in \!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\"}], \"}\"}]},\n\"Superscript\"]\), \
\!\(\*TemplateBox[{\"f\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \
\",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\), such that they stand in \
Mathematica default order. This is needed to ensure uniqueness."
 
OTThenAllSimpleRules = {CMExpr$_ :> Module[{Expr2$}, 
       Expr2$ = CMExpr$ /. OTGluonRules //. AllSimpleRules; Expr2$]}
 
OTThenAllSimpleRules /: OTThenAllSimpleRules::usage = "The rules \
OTThenAllSimpleRules first applies OTGluonRules and then repeatedly \
AllSimpleRules."
 
ExpandThenRules = {CMExpr$_ :> (Expand[CMExpr$] //. OTThenAllSimpleRules)}
 
ExpandThenRules /: ExpandThenRules::usage = "The rules ExpandThenRules first \
uses Expand, then applies OTThenAllSimpleRules and finally restores default \
index order with SortIndices."
 
CDotMatrix[CB$_, OptionsPattern[]] := Module[{VecI$, VecJ$, VIDotVJ$, 
      SPMrow$, SPMatr$}, SPMatr$ = {}; For[VecI$ = 1, VecI$ <= Length[CB$], 
       VecI$++, {SPMrow$ = {}; For[VecJ$ = 1, VecJ$ <= Length[CB$], VecJ$++, 
          {VIDotVJ$ = CDot[CB$[[VecI$]], CB$[[VecJ$]], NcMin -> OptionValue[
                NcMin]]; SPMrow$ = Append[SPMrow$, VIDotVJ$]; 
            If[OptionValue[Verbose], Print["Vector ", VecI$, " ", 
              CB$[[VecI$]], " dot vector ", VecJ$, " ", CB$[[VecJ$]], " ", 
              VIDotVJ$]]; }]; SPMatr$ = Append[SPMatr$, SPMrow$]; }]; SPMatr$]
 
Options[CDotMatrix] = {Verbose -> False, NcMin -> 3}
 
CDotMatrix /: CDotMatrix::usage = "CDotMatrix[{C1,...,Cn}] calculates the \
scalar product matrix between a list of color vectors (such as a basis) using \
CDot. I.e., element ij is CDot[Ci,Cj]. Before using CDotMatrix the color \
tensors \!\(\*SubscriptBox[\(C1\), \({I1, \(\(...\) \(Ik\)\)}\)]\) etc. \
should have been defined. CDotMatrix may be used with the Option Verbose, to \
write out progress information, CDotMatrix[{C1,...,Ck},Verbose\[Rule] True], \
and with the Option NcMin\[Rule] 100 (or any other value) for improved \
simplification of square roots."
 
CGamma[CB$_, p1$_, p2$_, OptionsPattern[]] := Module[{SPMMatr$, SPMInv$}, 
     If[p1$ == p2$, {Message[CGamma::partoncheck, p1$, p2$]; Abort[]}; ]; 
      If[OptionValue[BasisType] == GeneralBasis && 
         !OptionValue[CalledByOther], 
       {If[OptionValue[Verbose], Print[BasisWarning]]; }]; 
      If[OptionValue[NcMin] < 2, {Message[CGamma::NcMincheck, 
          OptionValue[NcMin]]; }, Null, 
       {Message[CGamma::NcMincheck, OptionValue[NcMin]]; }]; 
      If[OptionValue[Verbose], Null, Null, 
       {Message[CGamma::Verbosecheck, OptionValue[Verbose]]; }]; 
      {SPMMatr$, SPMInv$} = SPMAndSPMInvOfBasisType[CB$, 
        BasisType -> OptionValue[BasisType], NcMin -> OptionValue[NcMin], 
        Verbose -> OptionValue[Verbose]]; CGamma[CB$, p1$, p2$, SPMMatr$, 
       SPMInv$, MakeChecks -> OptionValue[MakeChecks], 
       BasisType -> OptionValue[BasisType], NcMin -> OptionValue[NcMin], 
       Verbose -> OptionValue[Verbose], CalledByOther -> True]]
 
CGamma[CB$_, p1$_, p2$_, spm$_, SPMInv$_, OptionsPattern[]] := 
    Module[{parton1$, parton2$, NInd$, Gam$, AEdecomp$, StartVec$, FinalVec$, 
      BasisLength$, IndV$, IndV2$, AfterExchange$, AE$, FinalTensor2$, 
      SquareDiff$, Decompv22$, GamCol$, cs$, ratio$, AEDotCB$, CMExpr$}, 
     If[p1$ == p2$, {Message[CGamma::partoncheck, p1$, p2$]; Abort[]}; ]; 
      If[OptionValue[BasisType] == GeneralBasis && 
         !OptionValue[CalledByOther], 
       {If[OptionValue[Verbose], Print[BasisWarning]]; }]; 
      If[OptionValue[NcMin] < 2, {Message[CGamma::NcMincheck, 
          OptionValue[NcMin]]; }, Null, 
       {Message[CGamma::NcMincheck, OptionValue[NcMin]]; }]; 
      If[OptionValue[Verbose], Null, Null, 
       {Message[CGamma::Verbosecheck, OptionValue[Verbose]]; }]; 
      If[ !OptionValue[BasisType] == GeneralBasis && 
         !OptionValue[BasisType] == TraceBasis && 
         !OptionValue[BasisType] == OrthogonalBasis && 
         !OptionValue[BasisType] == OrthonormalBasis, 
       {Message[CGamma::basischeck, OptionValue[BasisType]]; Abort[]}; , 
       Null, {Message[CGamma::basischeck, OptionValue[BasisType]]; Abort[]}]; 
      BasisLength$ = Length[CB$]; Gam$ = IdentityMatrix[BasisLength$] - 
        IdentityMatrix[BasisLength$]; NInd$ = NTensorIndices[CB$[[1]]]; 
      IndV$ = Table[Unique[d], {ii$, NInd$}]; 
      parton1$ = IdentifyParton[Subscript[CB$[[1]], IndV$], IndV$[[p1$]]]; 
      parton2$ = IdentifyParton[Subscript[CB$[[1]], IndV$], IndV$[[p2$]]]; 
      For[StartVec$ = 1, StartVec$ <= BasisLength$, 
       AfterExchange$ := Module[{i1$1, i2$1, Vertices$, sign$}, 
          IndV2$ = IndV$ /. {IndV$[[p1$]] -> i1$1, IndV$[[p2$]] -> i2$1}; 
           Vertices$ := Module[{gexch$}, CVertex[parton1$, IndV$[[p1$]], 
               i1$1, gexch$]*CVertex[parton2$, IndV$[[p2$]], i2$1, gexch$]]; 
           sign$ = 1; CSimplify[sign$*(Vertices$*Subscript[CB$[[StartVec$]], 
                IndV2$] /. FDToORules)]]; AE$ = AfterExchange$; 
        If[OptionValue[BasisType] == TraceBasis, 
         AEdecomp$ = SplitConstAndColor[Expand[AE$]]; ]; 
        If[OptionValue[BasisType] == GeneralBasis, 
         {AEDotCB$ = Table[Simplify[CSimplify[AE$*Conjugate[Subscript[
                  CB$[[cb$]], IndV$]]]], {cb$, 1, BasisLength$}]; 
           Table[Gam$[[FinalVec$,StartVec$]] = Simplify[Sum[Simplify[
                SPMInv$[[FinalVec$,cb$]]*AEDotCB$[[cb$]]], {cb$, 1, 
                BasisLength$}]], {FinalVec$, 1, BasisLength$}]; }]; 
        If[OptionValue[BasisType] == OrthogonalBasis, 
         {Table[{Gam$[[FinalVec$,StartVec$]] = Simplify[CSimplify[
                  AE$*Conjugate[Subscript[CB$[[FinalVec$]], IndV$]]]/
                 spm$[[FinalVec$,FinalVec$]]]; }; , {FinalVec$, 1, 
            BasisLength$}]}]; If[OptionValue[BasisType] == OrthonormalBasis, 
         {Table[Gam$[[FinalVec$,StartVec$]] = Simplify[CSimplify[
                AE$*Conjugate[Subscript[CB$[[FinalVec$]], IndV$]]]]; , 
            {FinalVec$, 1, BasisLength$}]; }]; 
        If[OptionValue[BasisType] == TraceBasis, 
         For[FinalVec$ = 1, FinalVec$ < BasisLength$ + 1, FinalVec$++, 
           {For[cs$ = 1, cs$ < Length[AEdecomp$] + 1, cs$++, 
              If[ !ContainsColor[ratio$ = AEdecomp$[[cs$,2]]/CSimplify[
                     Subscript[CB$[[FinalVec$]], IndV$]]], 
                {Gam$[[FinalVec$,StartVec$]] = Gam$[[FinalVec$,StartVec$]] + 
                   ratio$*AEdecomp$[[cs$,1]], Null}]; ]; }]; ]; 
        Table[{Gam$[[fVec$,StartVec$]] = Simplify[Simplify[Gam$[[fVec$,
               StartVec$]], Assumptions -> Nc >= OptionValue[NcMin]]]; 
           Gam$[[fVec$,StartVec$]] = PowerExpand[Gam$[[fVec$,StartVec$]], 
             Assumptions -> Nc >= OptionValue[NcMin]]; }, 
         {fVec$, 1, BasisLength$}]; If[OptionValue[MakeChecks], 
         GamCol$ = Table[Gam$[[fVec$,StartVec$]], {fVec$, 1, BasisLength$}]; 
          If[OptionValue[BasisType] == GeneralBasis || 
            OptionValue[BasisType] == TraceBasis, 
           {Decompv22$ = PowerExpand[Sum[Sum[Conjugate[GamCol$[[iVec$]]]*
                  GamCol$[[jVec$]]*spm$[[iVec$,jVec$]], {iVec$, 1, 
                  Length[GamCol$]}], {jVec$, 1, Length[GamCol$]}], 
               Assumptions -> Nc >= 2]; }]; If[OptionValue[BasisType] == 
            OrthonormalBasis, {Decompv22$ = Simplify[Sum[GamCol$[[iVec$]]*
                 Conjugate[GamCol$[[iVec$]]], {iVec$, BasisLength$}]]; }]; 
          If[OptionValue[BasisType] == OrthogonalBasis, 
           {Decompv22$ = Simplify[Sum[GamCol$[[iVec$]]*Conjugate[GamCol$[[
                   iVec$]]]*spm$[[iVec$,iVec$]], {iVec$, BasisLength$}]]; }]; 
          FinalTensor2$ = CSimplify[AE$*ReplaceDummyIndices[Conjugate[
                AE$]] /. FDToORules]; SquareDiff$ = 
           Simplify[Simplify[FinalTensor2$ - Decompv22$, Assumptions -> 
              Nc >= OptionValue[NcMin]]]; SquareDiff$ = 
           Simplify[PowerExpand[SquareDiff$, Assumptions -> Nc >= OptionValue[
                NcMin]]]; SquareDiff$ = Simplify[FullSimplify[SquareDiff$, 
             Assumptions -> Nc >= OptionValue[NcMin]]]; 
          If[SquareDiff$ == 0 && OptionValue[Verbose], Print["Completeness of \
basis for basis decomposed result verified for initial vector ", StartVec$, 
            "."]]; If[ContainsColor[SquareDiff$], 
           {Message[CGamma::colorcheck, StartVec$, SquareDiff$, AE$, 
              GamCol$]; }]; If[SquareDiff$ =!= 0 && 
             !ContainsColor[SquareDiff$], 
           {If[ !FreeQ[SquareDiff$, Sqrt[CMExpr$_]], Message[
               CGamma::squarerootcheck, StartVec$, SquareDiff$]]; 
             If[FreeQ[SquareDiff$, Sqrt[CMExpr$_]], Message[CGamma::check, 
                StartVec$, SquareDiff$, SquareDiff$ /. Nc -> 3]; ]}]; ]; 
        StartVec$++]; If[OptionValue[MakeChecks] && OptionValue[BasisType] == 
         OrthonormalBasis, {If[OptionValue[Verbose], 
          Print["Checking symmetry..."]]; CGammaSymmetryTest[Gam$, p1$, p2$, 
          NcMin -> OptionValue[NcMin]]; If[OptionValue[Verbose], 
          Print["Symmetry check done."]]; }]; Simplify[Gam$, 
       Assumptions -> {TR > 0, Nc >= OptionValue[NcMin]}]]
 
Options[CGamma] = {BasisType -> GeneralBasis, NcMin -> 3, Verbose -> True, 
     MakeChecks -> True, CalledByOther -> False}
 
CGamma /: CGamma::basischeck = "Error: the BasisType option should be \
BasisType\[Rule] GeneralBasis or BasisType\[Rule]OrthonormalBasis,\n or \
BasisType\[Rule]OrthogonalBasis or BasisType\[Rule]TraceBasis, but it was \
`1`."
 
CGamma /: CGamma::check = "Cannot verify completeness of basis for after \
decomposing initial vector `1`.\nPlease consider/check: \n1) the completeness \
of the basis \n2) simplifying square roots using the Option NcMin\[Rule] 100 \
(or any other value)\n3) that OrthogonalBasis/OrthonormalBasis should only be \
used if applicable.\n4) if mixed conventions of TR or Nc are used (e.g. TR is \
explicit 1/2 in basis but not generally defined) the completeness check which \
use the general parameters will fail.\nTensor^2-(decomposed \
tensor)^2\n`2`\nFor Nc=3: `3`\n"
 
CGamma /: CGamma::colorcheck = "Cannot determine completeness for initial \
vector `1` as color structure remains in: `2`.\nThe result after gluon \
exchange was: `3`\nThe decomposition into columns was: `4`"
 
CGamma /: CGamma::NcMincheck = 
     "NcMin should be \[GreaterEqual]2, but it was `1`,"
 
CGamma /: CGamma::partoncheck = "CGamma expects the partons to be different \
but they had number `1` and `2`. The color structure is trivial for exchange \
between a parton and itself."
 
CGamma /: CGamma::squarerootcheck = "Cannot determine completeness of basis \
for initial vector `1` due to square roots, consider using the option NcMin \
to set a minimal value for Nc\nTensor^2-(decomposed tensor)^2:\n`2`"
 
CGamma /: CGamma::usage = "CGamma[CB, p1, p2] calculates the effect of \
exchanging a gluon between parton p1 and parton p2 in the color basis \
CB={C1,..Cn}, where each basis vector, \!\(\*SubscriptBox[\(C1\), \({I1, \
\(\(...\) \(Ik\)\)}\)]\) etc., has been defined. The result is contained in a \
matrix, where column j contains the basis decomposed version of the vector \
resulting after a gluon exchange between p1 and p2 in the jth vector. The \
sign convention of the triple gluon vertex is given by ordering the partons \
as external index, internal dummy index, and index of gluon to be exchaged. \
The quark-gluon vertex comes without additional signs. Two more arguments, \
containing the scalar product matrix and its inverse, can be supplied, \
CGamma[CB,p1,p2,spm,Inverse[spm]]. As a self-consistency check, it is checked \
that the vector after gluon exchange has the same norm before and after being \
basis decomposed, unless the option MakeChecks\[Rule]False is used. For \
orthonormal bases it is also by default checked that the resulting matrix is \
symmetric. For more efficient calculations the option BasisType (with default \
value GeneralBasis) may be set as: BasisType \[Rule] OrthonormalBasis (for an \
orthonormal basis), BasisType \[Rule] OrthogonalBasis (for an orthogonal \
basis) or BasisType \[Rule] TraceBasis if the basis a is trace basis, i.e., a \
basis where each basis vector is proportional to one product of open and \
closed quark-lines (as opposed to a sum). For simplifying roots the option \
NcMin can be used, and for turning off progress messages Verbose\[Rule] False \
may be used."
 
CGamma /: CGamma::Verbosecheck = 
     "CGamma expects Verbose to be True or False, but it was `1`."
 
SPMAndSPMInvOfBasisType[CB$_, OptionsPattern[]] := 
    Module[{SPMatr$, SPMInv$}, 
     If[ !(OptionValue[BasisType] == OrthonormalBasis || 
         OptionValue[BasisType] == OrthogonalBasis || 
         OptionValue[BasisType] == TraceBasis || OptionValue[BasisType] == 
          GeneralBasis), Message[SPMAndSPMInvOfBasisType::basischeck, 
        OptionValue[BasisType]], Null, 
       Message[SPMAndSPMInvOfBasisType::basischeck, OptionValue[BasisType]]]; 
      If[OptionValue[BasisType] == GeneralBasis || OptionValue[BasisType] == 
         TraceBasis, {If[OptionValue[Verbose], 
          Print["Calculating general scalar product matrix..."]; ]; 
         SPMatr$ = CDotMatrix[CB$]; SPMatr$ = Simplify[PowerExpand[SPMatr$, 
            Assumptions -> Nc >= OptionValue[NcMin]]]; 
         If[OptionValue[Verbose], Print["The scalar product matrix is:"]]; 
         If[OptionValue[Verbose], Print[MatrixForm[SPMatr$]]; ]; 
         If[OptionValue[Verbose], 
          Print["Inverting scalar product matrix..."]; ]; 
         SPMInv$ = Simplify[Inverse[SPMatr$, Method -> 
             "OneStepRowReduction"]]; SPMInv$ = Simplify[PowerExpand[SPMInv$, 
            Assumptions -> Nc >= OptionValue[NcMin]]]; 
         If[OptionValue[Verbose], Print[
           "The inverse scalar product matrix is:"]]; 
         If[OptionValue[Verbose], Print[MatrixForm[SPMInv$]]; ]; }]; 
      If[OptionValue[BasisType] == OrthogonalBasis, 
       If[OptionValue[Verbose], 
         Print["Calculating diagonal scalar products ..."]; ]; 
        SPMatr$ = DiagonalMatrix[CDot[CB$, NcMin -> OptionValue[NcMin]]]; 
        If[OptionValue[Verbose], Print["The scalar product matrix is:"]; ]; 
        If[OptionValue[Verbose], Print[MatrixForm[SPMatr$]]; ]; 
        If[OptionValue[Verbose], 
         Print["Inverting scalar product matrix..."]; ]; 
        SPMInv$ = Simplify[Inverse[SPMatr$, Method -> 
            "OneStepRowReduction"]]; SPMInv$ = PowerExpand[SPMInv$, 
          Assumptions -> Nc >= OptionValue[NcMin]]; If[OptionValue[Verbose], 
         Print["The inverse scalar product matrix is:"]; ]; 
        If[OptionValue[Verbose], Print[MatrixForm[SPMInv$]]; ]; ]; 
      If[OptionValue[BasisType] == OrthonormalBasis, 
       SPMatr$ = DiagonalMatrix[Table[1, {uind$, 1, Length[CB$]}]]; 
        SPMInv$ = SPMatr$; ]; {SPMatr$, SPMInv$}]
 
Options[SPMAndSPMInvOfBasisType] = {BasisType -> GeneralBasis, NcMin -> 3, 
     Verbose -> True}
 
SPMAndSPMInvOfBasisType /: SPMAndSPMInvOfBasisType::basischeck = "Error: the \
BasisType option should be BasisType\[Rule] GeneralBasis or \
BasisType\[Rule]OrthonormalBasis,\n or BasisType\[Rule]OrthogonalBasis or \
BasisType\[Rule]TraceBasis, but it was `1`."
 
SPMAndSPMInvOfBasisType /: SPMAndSPMInvOfBasisType::usage = "Calculates the \
matrix of scalar products and its inverse using information about the basis \
supplied via the option BasisType. Used by CGamma."
 
IdentifyParton[CMExpr$_, parton$_] := Module[{LowerQuarks$, UpperQuarks$, 
      AllGluons$, PKind$}, LowerQuarks$ = LowerQuarkIndices[CMExpr$]; 
      UpperQuarks$ = UpperQuarkIndices[CMExpr$]; AllGluons$ = 
       GluonIndices[CMExpr$]; If[ !Intersection[LowerQuarks$, AllGluons$] == 
          {} ||  !Intersection[UpperQuarks$, AllGluons$] == {}, 
       {IdentifyParton::check = "The indices `1` seem to be of more than one \
kind. The UpperQuark indices are `2`, the LowerQuark indices are `3` and the \
Gluon indices are `4`."; Message[IdentifyParton::check, 
          Union[Intersection[LowerQuarks$, AllGluons$], Intersection[
            UpperQuarks$, AllGluons$]], UpperQuarks$, LowerQuarks$, 
          AllGluons$]; }]; If[ !FreeQ[LowerQuarks$, parton$], 
       {PKind$ = LowerQuark}]; If[ !FreeQ[UpperQuarks$, parton$], 
       PKind$ = UpperQuark]; If[ !FreeQ[AllGluons$, parton$], 
       PKind$ = Gluon]; If[PKind$ =!= LowerQuark && PKind$ =!= UpperQuark && 
        PKind$ =!= Gluon, {IdentifyParton::check = 
          "Cannot identify index `1` in `2`."; Message[IdentifyParton::check, 
          parton$, CMExpr$]}]; PKind$]
 
IdentifyParton /: IdentifyParton::usage = "IdentifyParton[Tens, parton] finds \
out if the parton is a \"LowerQuark\", i.e. a quark-type index sitting \
downstairs, an \"UpperQuark\", i.e. a quark-type index sitting upstairs or a \
Gluon."
 
CVertex[type$_, I1$_, I2$_, G3$1_] := Module[{CMres$}, 
     If[type$ =!= UpperQuark && type$ =!= LowerQuark && type$ =!= Gluon, 
       {CVertex::check = "CVertex[type,I1,I2,G3] expects its first argument \
to be UpperQuark, LowerQuark or Gluon, got `1`."; Message[CVertex::check, 
          type$]; }]; CMres$ = Switch[type$, UpperQuark, 
        Subscript[Superscript[Superscript[CMt, {G3$1}], I1$], I2$], LowerQuark,
        Subscript[Superscript[Superscript[CMt, {G3$1}], I2$], I1$], Gluon,
        I*Superscript[CMf, {I1$, I2$, G3$1}]]; CMres$]
 
CVertex /: CVertex::usage = "CVertex[type, I1, I2, G3] returns a vertex with \
indices I1, I2, G3. Depending on the type of particle to which a gluon G3 is \
attached either \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", \"G3\", \"}\"}]}, \"Superscript\"],\"I1\"},\n\"Superscript\"], \
\"I2\"]\) for type=UpperQuark, \
\!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \"G3\", \
\"}\"}]}, \"Superscript\"],\"I2\"},\n\"Superscript\"], \"I1\"]\) for \
type=LowerQuark,  or I \!\(\*TemplateBox[{\"f\",RowBox[{\"{\", \
RowBox[{\"I1\", \",\", \"I2\", \",\", \"G3\"}], \"}\"}]},\n\"Superscript\"]\) \
for type=Gluon is returned."
 
SplitConstAndColor[CMExpr$_] := Module[{Expr2$, ListExpr$, ColorDecomposed$, 
      ColorDecomposedRes$, UniqueColStrs$, ListExprI$, UniColStrInd$}, 
     If[ !ContainsColor[CMExpr$], {SplitConstAndColor::check = 
          "Cannot find color in color in `1`"; 
         Message[SplitConstAndColor::check, CMExpr$]; }]; 
      If[ !(Head[CMExpr$] =!= Times || Head[CMExpr$] =!= Superscript || 
         Head[CMExpr$] =!= Subscript), {SplitConstAndColor::check = "Expected a \
term with Head=Plus, Times, Superscript or Subscript, got `1` with Head=`2`"; 
         Message[SplitConstAndColor::check, CMExpr$, Head[CMExpr$]]; }]; 
      Expr2$ = Expand[CSimplify[CMExpr$]]; If[Head[Expr2$] == Plus, 
       {ListExpr$ = List @@ Expr2$; ColorDecomposed$ = 
          Table[SplitConstAndColorTerm[ListExpr$[[ListExprI$]]], 
           {ListExprI$, 1, Length[ListExpr$]}]; UniqueColStrs$ = 
          Sort[Union[Table[ColorDecomposed$[[ListExprI$,2]], 
             {ListExprI$, 1, Length[ColorDecomposed$]}]]]; 
         ColorDecomposedRes$ = Table[{0, UniqueColStrs$[[ListExprI$]]}, 
           {ListExprI$, 1, Length[UniqueColStrs$]}]; For[ListExprI$ = 1, 
          ListExprI$ <= Length[ColorDecomposed$], ListExprI$++, 
          For[UniColStrInd$ = 1, UniColStrInd$ <= Length[UniqueColStrs$], 
            UniColStrInd$++, If[ColorDecomposed$[[ListExprI$,2]] == 
               UniqueColStrs$[[UniColStrInd$]], {ColorDecomposedRes$[[
                  UniColStrInd$,1]] = ColorDecomposedRes$[[UniColStrInd$,
                   1]] + ColorDecomposed$[[ListExprI$,1]]; }]; ]; ]; }]; 
      If[Head[Expr2$] == Times || Head[Expr2$] == Superscript || 
        Head[Expr2$] == Subscript, ColorDecomposedRes$ = 
        {SplitConstAndColorTerm[Expr2$]}]; Simplify[ColorDecomposedRes$]]
 
SplitConstAndColor /: SplitConstAndColor::usage = "SplitConstAndColor[Expr] \
splits an expression Expr into a list of lists of color structures and \
corresponding multiplicative factors {{constants 1, color structure \
1},{constants 2, color structure 2},...}. This is done by first expanding the \
expression and then splitting the terms using SplitConstAndColorTerm[Term]."
 
SplitConstAndColorTerm[Term$_] := Module[{Termv2$, TermList$, ColorPart$, 
      ConstantPart$, ListLength$}, Termv2$ = Expand[Term$]; 
      If[Head[Termv2$] =!= Times && Head[Termv2$] =!= Superscript && 
        Head[Termv2$] =!= Subscript, {SplitConstAndColorTerm::check = "Expect\
ed a term with Head =Times, Superscript or Subscript got `1`, with Head= `2`  \
"; Message[SplitConstAndColorTerm::check, Termv2$, Head[Termv2$]]; }]; 
      ColorPart$ = 1; ConstantPart$ = 1; If[Head[Termv2$] == Times, 
       {TermList$ = List @@ Termv2$; For[ListLength$ = 1, 
          ListLength$ <= Length[TermList$], ListLength$++, 
          If[ContainsColor[TermList$[[ListLength$]]], ColorPart$ *= 
             TermList$[[ListLength$]], ConstantPart$ *= TermList$[[
              ListLength$]], {SplitConstAndColorTerm::check = 
               "Cannot decide on color in `1`"; Message[
               SplitConstAndColorTerm::check, TermList$[[
                ListLength$]]]; }]; ]; }]; If[Head[Termv2$] == Superscript || 
        Head[Termv2$] == Subscript, {If[ContainsColor[Termv2$], 
          ColorPart$ = Termv2$]; If[ !ContainsColor[Termv2$], 
          {SplitConstAndColorTerm::check = 
             "Cannot find color in color in `1`"; Message[
             SplitConstAndColorTerm::check, Termv2$]; }]; }]; 
      {ConstantPart$, ColorPart$}]
 
SplitConstAndColorTerm /: SplitConstAndColorTerm::usage = "SplitConstAndColor\
Term[Term] splits a single term consisting of one color structure, as opposed \
to a sum of color structures, into a list of term containing {constants, \
color structure}."
 
ReplaceDummyIndices[CMExpr$_] := Module[{ExpandExpr$, NewDummyIndex$, 
      ReplaceTerm$, DummySymbols$, IndRules$, ListExpr$, Expr2$}, 
     ExpandExpr$ = Expand[CMExpr$]; ListExpr$ = List @@ ExpandExpr$; 
      NewDummyIndex$ = Table[Unique[d], {uInd$, 100}]; 
      ReplaceTerm$[Term$_] := (DummySymbols$ = DummyIndicesTerm[Term$]; 
        IndRules$ = Table[DummySymbols$[[dInd$]] -> NewDummyIndex$[[dInd$]], 
          {dInd$, Length[DummySymbols$]}]; Term$ /. IndRules$); 
      ListExpr$ = ReplaceTerm$ /@ ListExpr$; If[Head[ExpandExpr$] == Plus, 
       Expr2$ = Plus @@ ListExpr$, Expr2$ = CMExpr$, Expr2$ = CMExpr$]; 
      If[Head[ExpandExpr$] == Times, Expr2$ = ReplaceTerm$[CMExpr$]]; 
      If[Head[ExpandExpr$] == Superscript || Head[ExpandExpr$] == Subscript, 
       {Expr2$ = ReplaceTerm$[CMExpr$]}]; Expr2$]
 
ReplaceDummyIndices /: ReplaceDummyIndices::usage = "ReplaceDummyIndices[Expr\
] replaces the dummy indices in Expr with a new set of unique dummy indices."
 
DummyIndicesTerm[Term$_] := Module[{AllInd$, DummySymbols$, PairSymbols$}, 
     If[Head[Term$] == Plus, DummyIndicesTerm::argx = "DummyIndicesTerm \
should be used on a single term, not on a sum of terms, the agument was `1`."\
; Message[DummyIndicesTerm::argx, Term$]; ]; AllInd$ = AllIndices[Term$]; 
      PairSymbols$ = AllPairs[Term$]; DummySymbols$ = 
       Intersection[PairSymbols$, AllInd$]; Union[DummySymbols$, 
       SquareIndices[Term$]]]
 
DummyIndicesTerm /: DummyIndicesTerm::usage = "DummyIndicesTerm[Term] \
replaces the dummy indices in one term (as opposed to a sum of terms) with a \
new set of unique dummy indices."
 
SquareIndices[Term$_] := Module[{SqInd$}, 
     If[Head[Term$] == Plus, {SquareIndices::argx = "SquareIndices should be \
used on a single term, not on a sum of terms, the agument was `1`."; 
         Message[SquareIndices::argx, Term$]; }]; SqInd$ = {}; 
      If[ !FreeQ[Term$, Superscript[CM\[CapitalDelta], Gs$_]^2], 
       Term$ /. Superscript[CM\[CapitalDelta], Gs$_]^2 :> 
          Module[{}, {Gs$; SqInd$ = Append[SqInd$, Gs$]; }]; ]; 
      If[ !FreeQ[Term$, Superscript[CMo, Gs$_]^2], 
       Term$ /. Superscript[CMo, Gs$_]^2 :> Module[{}, 
           {Gs$; SqInd$ = Append[SqInd$, Gs$]; }]; ]; 
      If[ !FreeQ[Term$, Superscript[CMf, Gs$_]^2],
       Term$ /. Superscript[CMf, Gs$_]^2 :> Module[{},
           {Gs$; SqInd$ = Append[SqInd$, Gs$]; }]; ]; 
      If[ !FreeQ[Term$, Superscript[CMd, Gs$_]^2], 
       Term$ /. Superscript[CMd, Gs$_]^2 :> Module[{}, 
           {Gs$; SqInd$ = Append[SqInd$, Gs$]; }]; ]; Flatten[SqInd$]]
 
SquareIndices /: SquareIndices::usage = "SquareIndices[Term] identifies dummy \
indices in squares, such as \!\(\*TemplateBox[{\"f\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\", \",\", \"G3\"}], \
\"}\"}]},\n\"Superscript\"]\)^2, used by DummyIndicesTerm[Term]."
 
CGammaSymmetryTest[Matr$_, parton1$_, parton2$_, OptionsPattern[]] := 
    Module[{MatrAsym$, IsSym$}, IsSym$ = True; 
      MatrAsym$ = Simplify[PowerExpand[Matr$ - Transpose[Matr$]], 
        Assumptions -> Nc >= OptionValue[NcMin]]; 
      If[Union[Flatten[MatrAsym$]] =!= {0}, 
       {IsSym$ = False; CGamma::symcheck = "Can not verify symmetry for \
excahnge between `1` and  `2`.\nGamma-Transpose[Gamma]:\n"; 
         Message[CGamma::symcheck, parton1$, parton2$]; 
         Print[MatrixForm[MatrAsym$]]; }]; IsSym$]
 
Options[CGammaSymmetryTest] = {NcMin -> 3}
 
CNorm[CBorV$_, OptionsPattern[]] := Module[{nInd$, CMres$, NormI$, VecI$, 
      CVectors$}, If[Head[CBorV$] == List, {CVectors$ = CBorV$; CMres$ = {}; 
         For[VecI$ = 1, VecI$ <= Length[CVectors$], VecI$++, 
          {NormI$ = Simplify[Sqrt[Simplify[CDot[CBorV$[[VecI$]], 
                 CBorV$[[VecI$]], NcMin -> OptionValue[NcMin]], 
                Assumptions -> {Nc >= OptionValue[NcMin]}]], Assumptions -> {
                Element[Nc, Reals], Nc >= OptionValue[NcMin]}]; 
            If[OptionValue[Verbose], Print[VecI$, " ", CVectors$[[VecI$]], 
              " ", NormI$]]; CMres$ = Append[CMres$, NormI$]; }]; }]; 
      If[Head[CBorV$] =!= List, {nInd$ = NTensorIndices[CBorV$]; 
         CMres$ = Simplify[Sqrt[Simplify[CDot[CBorV$, CBorV$, 
              NcMin -> OptionValue[NcMin]], Assumptions -> 
              {Nc >= OptionValue[NcMin]}]], Assumptions -> 
            {Element[Nc, Reals], Nc >= OptionValue[NcMin]}]; }]; CMres$]
 
Options[CNorm] = {Verbose -> False, NcMin -> 3}
 
CNorm /: CNorm::usage = "CNorm[C] calculates the norm of a color tensors \
using CDot[C]. By default CNorm assumes Nc\[GreaterSlantEqual]3. For improved \
simplification of roots the Option NcMin\[Rule]100 (or any other value) may \
be used to assume that Nc is at least 100 (or any other value) when \
simplifying. Before using CNorms the color tensor \!\(\*SubscriptBox[\(C\), \
\({I1, \(\(...\) \(Ik\)\)}\)]\) should have been defined."
 
CyclicPermutations[Superscript[Tfd$_, {G1$1_, G2$1_, G3$1_}]] := 
    Superscript[Tfd$, {G1$1, G2$1, G3$1}] | Superscript[Tfd$, 
      {G2$1, G3$1, G1$1}] | Superscript[Tfd$, {G3$1, G1$1, G2$1}]
 
CMdelta[Q1$_, Q2$_] := Subscript[Superscript[CM\[Delta], Q1$], Q2$]
 
CMDelta[G1$_, G2$_] := Superscript[CM\[CapitalDelta], {G1$, G2$}]
 
DummyIndices[CMExpr$_] := Module[{Expr2$, DInd$, ListExpr$, IndexList$}, 
     Expr2$ = Expand[CMExpr$]; DInd$ = {}; If[Head[Expr2$] != Plus, 
       DInd$ = Flatten[Append[DInd$, DummyIndicesTerm[Expr2$]]], Null, 
       DInd$ = Flatten[Append[DInd$, DummyIndicesTerm[Expr2$]]]]; 
      If[Head[Expr2$] == Plus, {ListExpr$ = List @@ Expr2$; 
         IndexList$ = DummyIndicesTerm /@ ListExpr$; 
         DInd$ = Union[Flatten[IndexList$]]; }]; DInd$]
 
DummyIndices /: DummyIndices::usage = "DummyIndices[Expr] finds the dummy \
indices in an expression Expr, by first expanding it and then finding all \
dummy indices in all terms."
 
FDToTRules := {Superscript[CMf, {G1$1_, G2$1_, G3$1_}] :>
      (1/(I*TR))*Module[{q1$1, q2$1, q3$1}, 
        Subscript[Superscript[Superscript[CMt, {G1$1}], q1$1], q2$1]*
          Subscript[Superscript[Superscript[CMt, {G2$1}], q2$1], q3$1]*
          Subscript[Superscript[Superscript[CMt, {G3$1}], q3$1], q1$1] -
         Subscript[Superscript[Superscript[CMt, {G2$1}], q1$1], q2$1]*
          Subscript[Superscript[Superscript[CMt, {G1$1}], q2$1], q3$1]*
          Subscript[Superscript[Superscript[CMt, {G3$1}], q3$1], q1$1]],
     Superscript[CMd, {G1$1_, G2$1_, G3$1_}] :> 
      (1/TR)*Module[{q1$1, q2$1, q3$1}, 
        Subscript[Superscript[Superscript[CMt, {G1$1}], q1$1], q2$1]*
          Subscript[Superscript[Superscript[CMt, {G2$1}], q2$1], q3$1]*
          Subscript[Superscript[Superscript[CMt, {G3$1}], q3$1], q1$1] +
         Subscript[Superscript[Superscript[CMt, {G2$1}], q1$1], q2$1]*
          Subscript[Superscript[Superscript[CMt, {G1$1}], q2$1], q3$1]*
          Subscript[Superscript[Superscript[CMt, {G3$1}], q3$1], q1$1]]}
 
FDToTRules /: FDToTRules::usage = "Replaces structure constants, \
\!\(\*TemplateBox[{\"f\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\) and \
\!\(\*TemplateBox[{\"d\",RowBox[{\"{\", RowBox[{\"G1\", \",\", \"G2\", \",\", \
\"G3\"}], \"}\"}]},\n\"Superscript\"]\), with a sum of products of SU(Nc) \
generators, \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", \"G1\", \"}\"}]}, \"Superscript\"],\"Q1\"},\n\"Superscript\"], \
\"Q2\"]\)."
 
FreeIndices[CMExpr$_] := Complement[AllIndices[CMExpr$], DummyIndices[CMExpr$]]

GluonContract[CMExpr$_, Gs$_] := Module[{ContrInd$, CMres$, theg$, Glist$}, 
     If[ !Head[Gs$] == List && Intersection[AllIndices[CMExpr$], {Gs$}] == 
         {Gs$}, {Glist$ = {Gs$}}, Glist$ = Gs$; , {Glist$ = {Gs$}}]; 
      If[ContainsFD[CMExpr$], {GluonContract::fdcheck = "The first argument to \
GluonContract contains structure constants, f and/or d which will be replaced \
using FDToORules."; Message[GluonContract::fdcheck]; 
         CMres$ = CMExpr$ //. FDToORules; }, CMres$ = Expand[CMExpr$]; ]; 
      CMres$ = Expand[CMExpr$]; For[ContrInd$ = 1, ContrInd$ <= Length[Glist$], 
       ContrInd$++, {theg$ = Glist$[[ContrInd$]]; 
         CMres$ = OGluonContract[CMres$, theg$]; CMres$ = OOGluonContract[CMres$, 
           theg$]; CMres$ = TGluonContract[CMres$, theg$]; 
         CMres$ = TTGluonContract[CMres$, theg$]; CMres$ = OTGluonContract[CMres$, 
           theg$]; CMres$ = ODeltaGluonContract[CMres$, theg$]; 
         CMres$ = O2GluonContract[CMres$, theg$]; 
         CMres$ = CMres$ /. AllPermutations[Superscript[CMo, {theg$, 
                G1$1_}]]->TR Superscript[CM\[CapitalDelta], {theg$, G1$1}]; 
         CMres$ = CMres$ /. AllPermutations[Superscript[CM\[CapitalDelta], {theg$, 
                G1$1_}]]*AllPermutations[Superscript[CM\[CapitalDelta], {theg$, 
                G2$1_}]] -> Superscript[CM\[CapitalDelta], {G1$1, G2$1}]; 
         CMres$ = CMres$ /. Superscript[CM\[CapitalDelta], {theg$, theg$}] -> 
            -1 + Nc^2; CMres$ = CMres$ /. AllPermutations[Superscript[
                CM\[CapitalDelta], {G1$1, theg$}]]*Superscript[CMo, G3$1s_] /; 
             Count[G3$1s, theg$] == 1 :> Superscript[CMo, 
             G3$1s /. theg$ -> G1$1]; CMres$ = CMres$ //. 
           Superscript[CM\[CapitalDelta], {theg$, G1$1_}]^2 -> 
            Superscript[CM\[CapitalDelta], {G1$1, G1$1}]; 
         CMres$ = CMres$ //. Superscript[CM\[CapitalDelta], {G1$1_, theg$}]^2 -> 
            Superscript[CM\[CapitalDelta], {G1$1, G1$1}]; 
         CMres$ = Expand[CMres$]}]; CMres$]
 
GluonContract /: GluonContract::usage = "GluonContract[Expr, Gs] contracts a \
set (List) of gluons Gs={g1,...,gk} or a single gluon Gs=g1, in the \
expression Expr, while leaving other indices uncontracted. This function is \
intended for quark-lines (os and ts) and will replace structure constants (fs \
and ds) with quark-lines."
 
OGluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. Superscript[CMo, Gs$_ /; Count[Gs$, theg$] == 2] :> 
      Module[{Place1$, Place2$, Newgs1$, Newgs2$, CMres$, rem$, 
        Multiplicities$}, Multiplicities$ = Transpose[Tally[Gs$]][[2]]; 
        Place1$ = Position[Gs$, theg$][[1]][[1]]; 
        Place2$ = Position[Gs$, theg$][[2]][[1]]; 
        If[Place2$ - Place1$ == 1 || (Place1$ == 1 && Place2$ == 
            Length[Gs$]), {Newgs1$ = Drop[Gs$, {Place1$, Place2$}]; 
           Newgs1$ = Delete[Delete[Gs$, Place2$], Place1$]; 
           CMres$ = (TR/Nc)*(Nc^2 - 1)*Superscript[CMo, Newgs1$] /. 
             Remove0To1ORules; }, If[Place2$ - Place1$ == 2, 
          {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
            If[ !Length[Union[Gs$]] - Length[Union[Newgs1$]] == 1, 
             {Print["Warning "]; Print["Old ", Superscript[CMo, Gs$]]; Print[
                "New ", Superscript[CMo, Newgs1$]]; }]; 
            CMres$ = (-(TR/Nc))*Superscript[CMo, Newgs1$] /. Remove0To1ORules; }, 
          {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
            Newgs2$ = Drop[Gs$, {Place1$, Place2$}]; 
            rem$ = Take[Gs$, {Place1$ + 1, Place2$ - 1}]; 
            CMres$ = TR*(Superscript[CMo, Newgs2$]*Superscript[CMo, rem$] - 
                (1/Nc)*Superscript[CMo, Newgs1$]) /. Remove0To1ORules; }]]; 
        CMres$]
 
OGluonContract /: OGluonContract::usage = "OGluonContract[Expr,theg] \
contracts the gluon theg in the expression Expr if theg appears twice in one \
closed quark-line \!\(\*TemplateBox[{\"o\",RowBox[{\"{\", RowBox[{\"G1\", \
\",\", \"...\", \",\", \"theg\", \",\", \"...\", \",\", \"theg\", \",\", \
\"...\", \",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\)."
 
Remove0To1ORules = {Superscript[CMo, LIndV$_ /; Length[LIndV$] == 1] -> 0, 
     Superscript[CMo, LIndV$_ /; Length[LIndV$] == 0] -> Nc}
 
Remove0To1ORules /: Remove0To1ORules::usage = "Rules for simplifying closed \
quark-lines with 0 or 1 gluons,\!\(\*TemplateBox[{RowBox[{\" \", \
\"o\"}],RowBox[{\"{\", \"}\"}]},\n\"Superscript\"]\)=Nc, \
\!\(\*TemplateBox[{\"o\",RowBox[{\"{\", \"G1\", \
\"}\"}]},\n\"Superscript\"]\)=0."
 
OOGluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. Superscript[CMo, g1$1s_]*Superscript[CMo, g2$1s_] /; 
       Count[g1$1s, theg$] == 1 && Count[g2$1s, theg$] == 1 :> 
      Module[{Place1$, Place2$, g1$1s1, g1$1s2, g2$1s1, g2$1s2, NewL$, 
        NewL1$, NewL2$, SecondL$, FirstL$, Supterm$}, 
       If[Order[g1$1s[[1]], g2$1s[[1]]] == -1, 
         {FirstL$ = g2$1s; SecondL$ = g1$1s}, 
         {FirstL$ = g1$1s; SecondL$ = g2$1s}]; Place1$ = 
         Position[FirstL$, theg$][[1,1]]; Place2$ = 
         Position[SecondL$, theg$][[1,1]]; g1$1s1 = Take[FirstL$, 
          Place1$ - 1]; g2$1s1 = Take[SecondL$, Place2$ - 1]; 
        g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewL$ = {}; 
        NewL$ = Join[g1$1s1, g2$1s2, g2$1s1, g1$1s2]; 
        NewL1$ = Drop[g1$1s, {Position[g1$1s, theg$][[1,1]]}]; 
        NewL2$ = Drop[g2$1s, {Position[g2$1s, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Superscript[CMo, NewL1$]*Superscript[CMo, NewL2$]; 
        TR*(Superscript[CMo, NewL$] - Supterm$) /. Remove0To1ORules]
 
OOGluonContract /: OOGluonContract::usage = "OOGluonContract[Expr,theg] \
contracts the gluon theg in the expression Expr if theg appears once in a \
closed quark-line (\!\(\*TemplateBox[{RowBox[{\"(\", TemplateBox[{\"t\", \
RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \",\", \"theg\", \",\", \
\"...\", \",\", \"Gk\"}], \"}\"}]}, \"Superscript\"], \
\")\"}],\"Q1\"},\n\"Superscript\"]\)\!\(\*SubscriptBox[\()\), \(Q2\)]\) an \
once in another closed quark-line."
 
TGluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. Subscript[Superscript[Superscript[CMt,
         Gs$_ /; Count[Gs$, theg$] == 2], Q1$1_], Q2$1_] :> 
      Module[{Place1$, Place2$, Newgs1$, Newgs2$, CMres$, rem$, 
        Multiplicities$}, Multiplicities$ = Transpose[Tally[Gs$]][[2]]; 
        Place1$ = Position[Gs$, theg$][[1]][[1]]; 
        Place2$ = Position[Gs$, theg$][[2]][[1]]; If[Place2$ - Place1$ == 1, 
         {Newgs1$ = Drop[Gs$, {Place1$, Place2$}]; 
           CMres$ = (TR/Nc)*(Nc^2 - 1)*Subscript[Superscript[Superscript[CMt,
                 Newgs1$], Q1$1], Q2$1] /. Remove0ORules; }, 
         If[Place2$ - Place1$ == 2, {Newgs1$ = Drop[Drop[Gs$, {Place2$}], 
              {Place1$}]; If[ !Length[Union[Gs$]] - Length[Union[Newgs1$]] == 
               1, {Print["TGluonContract: Warning "]; Print[
                "TGluonContract: Old ", Subscript[Superscript[Superscript[CMt,
                   Gs$], Q1$1], Q2$1]]; Print["TGluonContract: New ", 
                Subscript[Superscript[Superscript[CMt, Newgs1$], Q1$1],
                 Q2$1]]; }]; CMres$ = (-(TR/Nc))*Subscript[Superscript[
                 Superscript[CMt, Newgs1$], Q1$1], Q2$1] /. Remove0ORules; },
          {Newgs1$ = Drop[Drop[Gs$, {Place2$}], {Place1$}]; 
            Newgs2$ = Drop[Gs$, {Place1$, Place2$}]; 
            rem$ = Take[Gs$, {Place1$ + 1, Place2$ - 1}]; 
            CMres$ = TR*(Subscript[Superscript[Superscript[CMt, Newgs2$], Q1$1],
                  Q2$1]*Superscript[CMo, rem$] - (1/Nc)*Subscript[Superscript[
                   Superscript[CMt, Newgs1$], Q1$1], Q2$1]) /.
              Remove0ORules; }]]; CMres$]
 
TGluonContract /: TGluonContract::usage = "TGluonContract[Expr,theg] \
contracts the gluon theg in the expression Expr if theg appears twice in one \
open quark-line \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \",\", \"theg\", \",\", \
\"...\", \",\", \"theg\", \",\", \"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\)."
 
TTGluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. Subscript[Superscript[Superscript[CMt, g1$1s_], Q1$1_], Q2$1_]*
        Subscript[Superscript[Superscript[CMt, g2$1s_], Q3$1_], Q4$1_] /;
       Count[g1$1s, theg$] == 1 && Count[g2$1s, theg$] == 1 :> 
      Module[{Place1$, Place2$, g1$1s1, g1$1s2, g2$1s1, g2$1s2, NewLa$, 
        NewLb$, NewL1$, NewL2$, SecondL$, FirstL$, Supterm$}, 
       FirstL$ = g1$1s; SecondL$ = g2$1s; Place1$ = Position[FirstL$, theg$][[
          1,1]]; Place2$ = Position[SecondL$, theg$][[1,1]]; 
        g1$1s1 = Take[FirstL$, Place1$ - 1]; g2$1s1 = Take[SecondL$, 
          Place2$ - 1]; g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewLa$ = {}; 
        NewLb$ = {}; NewLa$ = Join[g1$1s1, g2$1s2]; 
        NewLb$ = Join[g2$1s1, g1$1s2]; NewL1$ = Drop[g1$1s, 
          {Position[g1$1s, theg$][[1,1]]}]; NewL2$ = 
         Drop[g2$1s, {Position[g2$1s, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Subscript[Superscript[Superscript[CMt, NewL1$],
            Q1$1], Q2$1]*Subscript[Superscript[Superscript[CMt, NewL2$], Q3$1],
           Q4$1]; TR*(Subscript[Superscript[Superscript[CMt, NewLa$], Q1$1],
             Q4$1]*Subscript[Superscript[Superscript[CMt, NewLb$], Q3$1],
             Q2$1] - Supterm$) /. Union[Remove0To1ORules, Remove0ORules]]
 
TTGluonContract /: TTGluonContract::usage = "TTGluonContract[Expr,theg] \
contracts the gluon theg in the expression Expr if theg appears once one open \
quark-line, \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", \
RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \",\", \"theg\", \",\", \
\"...\", \",\", \"Gk\"}], \"}\"}]}, \
\"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\), and once in another \
open quark-line."
 
OTGluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. Subscript[Superscript[Superscript[CMt, g1$1s_], Q1$1_], Q2$1_]*
        Superscript[CMo, g2$1s_] /; Count[g1$1s, theg$] == 1 && 
        Count[g2$1s, theg$] == 1 :> Module[{Place1$, Place2$, g1$1s1, g1$1s2, 
        g2$1s1, g2$1s2, NewL$, NewL1$, NewL2$, SecondL$, FirstL$, Supterm$}, 
       FirstL$ = g1$1s; SecondL$ = g2$1s; Place1$ = 
         Position[FirstL$, theg$][[1]][[1]]; Place2$ = 
         Position[SecondL$, theg$][[1,1]]; g1$1s1 = Take[FirstL$, 
          Place1$ - 1]; g2$1s1 = Take[SecondL$, Place2$ - 1]; 
        g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewL$ = {}; 
        NewL$ = Join[g1$1s1, g2$1s2, g2$1s1, g1$1s2]; 
        NewL1$ = Drop[g1$1s, {Position[g1$1s, theg$][[1,1]]}]; 
        NewL2$ = Drop[g2$1s, {Position[g2$1s, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Subscript[Superscript[Superscript[CMt, NewL1$],
            Q1$1], Q2$1]*Superscript[CMo, NewL2$]; 
        TR*(Subscript[Superscript[Superscript[CMt, NewL$], Q1$1], Q2$1] -
           Supterm$) /. Union[Remove0To1ORules, Remove0ORules]]
 
OTGluonContract /: OTGluonContract::usage = "OTGluonContract[Expr,theg] \
contracts the gluon theg in the expression Expr if theg appears once one open \
quark-line \!\(\*SubscriptBox[TemplateBox[{TemplateBox[{\"t\", RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"...\", \",\", \"theg\", \",\", \"...\", \",\", \
\"Gk\"}], \"}\"}]}, \"Superscript\"],\"Q1\"},\n\"Superscript\"], \"Q2\"]\), \
and once a closed quark-line, \!\(\*TemplateBox[{\"o\",RowBox[{\"{\", \
RowBox[{\"Ga\", \",\", \"...\", \",\", \"theg\", \",\", \"...\", \",\", \
\"GNg\"}], \"}\"}]},\n\"Superscript\"]\)."
 
ODeltaGluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. AllPermutations[Superscript[CM\[CapitalDelta], {G1$1_, theg$}]]*
        Superscript[CMo, G2$1s_] /; Count[G2$1s, theg$] == 1 :> 
      Superscript[CMo, G2$1s /. theg$ -> G1$1]
 
ODeltaGluonContract /: ODeltaGluonContract::usage = "ODeltaGluonContract[Expr\
,theg] contracts the gluon theg in the expression Expr if theg appears once \
in a closed quark-line \!\(\*TemplateBox[{RowBox[{\" \", RowBox[{\"(\", \
\"o\"}]}],RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \",\", \"theg\", \
\",\", \"...\", \",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\) an once in a \
gluon delta function \!\(\*TemplateBox[{\"\[CapitalDelta]\",RowBox[{\"{\", \
RowBox[{\"G1\", \",\", \"G2\"}], \"}\"}]},\n\"Superscript\"]\)."
 
O2GluonContract[CMExpr$_, theg$_] := 
    CMExpr$ /. Superscript[CMo, Gs$_]^2 /; Count[Gs$, theg$] == 1 :> 
      Module[{Place1$, Place2$, g1$1s1, g1$1s2, g2$1s1, g2$1s2, NewL$, 
        NewL1$, NewL2$, SecondL$, FirstL$, Supterm$}, 
       FirstL$ = Gs$; SecondL$ = Gs$; Place1$ = Position[FirstL$, theg$][[1,
          1]]; Place2$ = Position[SecondL$, theg$][[1,1]]; 
        g1$1s1 = Take[FirstL$, Place1$ - 1]; g2$1s1 = Take[SecondL$, 
          Place2$ - 1]; g1$1s2 = Take[FirstL$, -(Length[FirstL$] - Place1$)]; 
        g2$1s2 = Take[SecondL$, -(Length[SecondL$] - Place2$)]; NewL$ = {}; 
        NewL$ = Join[g1$1s1, g2$1s2, g2$1s1, g1$1s2]; 
        NewL1$ = Drop[Gs$, {Position[Gs$, theg$][[1,1]]}]; 
        NewL2$ = Drop[Gs$, {Position[Gs$, theg$][[1,1]]}]; 
        Supterm$ = (1/Nc)*Superscript[CMo, NewL1$]*Superscript[CMo, NewL2$]; 
        TR*(Superscript[CMo, NewL$] - Supterm$) /. Remove0To1ORules]
 
O2GluonContract /: O2GluonContract::usage = "O2GluonContract[Expr,theg] \
contracts the gluon theg in the expression Expr if theg appears once a \
squared quark-line\!\(\*TemplateBox[{RowBox[{\" \", RowBox[{\"(\", \
\"o\"}]}],RowBox[{\"{\", RowBox[{\"G1\", \",\", \"...\", \",\", \"theg\", \
\",\", \"...\", \",\", \"Gk\"}], \"}\"}]},\n\"Superscript\"]\))^2."
 

WhatIsWrong[CMExpr$_] := Module[{ExpandExpr$, ListExpr$, FreeInds$}, 
     ExpandExpr$ = Expand[CMExpr$]; ListExpr$ = {}; 
      If[Head[ExpandExpr$] == Plus, {ListExpr$ = List @@ ExpandExpr$}, 
       ListExpr$ = Append[ListExpr$, CMExpr$]; ]; 
      If[Head[ExpandExpr$] == Times || Head[ExpandExpr$] == Subscript || 
        Head[ExpandExpr$] == Superscript || Head[ExpandExpr$] == Power || 
        Head[ExpandExpr$] == List, {ListExpr$ = Append[ListExpr$, CMExpr$]}]; 
      FreeInds$ = Table[WhatIsWrongTerm[ListExpr$[[fInd$]]], 
        {fInd$, 1, Length[ListExpr$]}]; If[ !Length[Union[FreeInds$]] == 1, 
       {WhatIsWrong::indices = "It seems that the number of free indices in \
the expanded expression differs from term to term."; 
         Message[WhatIsWrong::indices]; }]; ]
 
WhatIsWrong /: WhatIsWrong::usage = "WhatIsWrong[Expr] checks if anything is \
obviously wrong with an expression, for example if Power is used instead of \
Superscript. The check is done by first expanding the expression, and then \
checking each term using WhatIsWrongTerm."
 
WhatIsWrongTerm[Term$_] := Module[{Gs$, G1$1, G2$1, Q1$1, Q2$1, AllInd$, 
      AllSymbols$, symb$, nIndi$, ConsTerm$, fdComb$}, 
     If[Head[Term$] == Plus, WhatIsWrongTerm::indices = "The function \
WhatIsWrongTerm expects a single term and got a sum of terms."; 
        Message[WhatIsWrongTerm::indices]]; WhatIsWrongTerm::indices = 
       "The indices in `1` should sit upstairs in Superscript."; 
      If[ !FreeQ[Term$, Subscript[CM\[CapitalDelta], Gs$_]], 
       Term$ /. Subscript[CM\[CapitalDelta], Gs$_] :> Module[{}, 
           Message[WhatIsWrongTerm::indices, Subscript[CM\[CapitalDelta], 
              Gs$]]; Term$]; ]; If[ !FreeQ[Term$, Subscript[CMf, Gs$_]],
       Term$ /. Subscript[CMf, Gs$_] :> Module[{},
           Message[WhatIsWrongTerm::indices, Subscript[CMf, Gs$]]; Term$]; ];
      If[ !FreeQ[Term$, Subscript[CMd, Gs$_]], 
       Term$ /. Subscript[CMd, Gs$_] :> Module[{}, 
           Message[WhatIsWrongTerm::indices, Subscript[CMd, Gs$]]; Term$]; ]; 
      WhatIsWrongTerm::indices = "The form of `1` should be `2`."; 
      If[ !FreeQ[Term$, Superscript[Subscript[CM\[Delta], Q1$1_], Q2$1_]], 
       {Term$ /. Superscript[Subscript[CM\[Delta], Q1$1_], Q2$1_] :> 
           Module[{}, Message[WhatIsWrongTerm::indices, Superscript[Subscript[
                CM\[Delta], Q1$1], Q2$1], Subscript[Superscript[CM\[Delta], 
                Q1$1], Q2$1]]; Term$]; }]; 
      If[ !FreeQ[Term$, Subscript[CM\[Delta], Q1$1_, Q2$1_]], 
       {Term$ /. Subscript[CM\[Delta], Q1$1_, Q2$1_] :> Module[{}, 
            Message[WhatIsWrongTerm::indices, Subscript[CM\[Delta], Q1$1, 
               Q2$1], Subscript[Superscript[CM\[Delta], Q1$1], Q2$1]]; 
             Term$]; }]; If[ !FreeQ[Term$, Subscript[CM\[Delta], 
          {Q1$1_, Q2$1_}]], {Term$ /. Subscript[CM\[Delta], {Q1$1_, Q2$1_}] :> 
           Module[{}, Message[WhatIsWrongTerm::indices, Subscript[CM\[Delta], {
                Q1$1, Q2$1}], Subscript[Superscript[CM\[Delta], Q1$1], Q2$1]]; 
             Term$]; }]; If[ !FreeQ[Term$, Superscript[
          Subscript[Superscript[CMt, Gs$_], Q1$1_], Q2$1_]],
       {Term$ /. Superscript[Subscript[Superscript[CMt, Gs$_], Q1$1_],
            Q2$1_] :> Module[{}, Message[WhatIsWrongTerm::indices, 
              Superscript[Subscript[Superscript[CMt, Gs$], Q1$1], Q2$1],
              Subscript[Superscript[Superscript[CMt, Gs$], Q1$1], Q2$1]];
             Term$]; }]; If[ !FreeQ[Term$, Subscript[CM\[CapitalDelta], G1$1_, 
          G2$1_]], {Term$ /. Subscript[CM\[CapitalDelta], G1$1_, G2$1_] :> 
           Module[{}, Message[WhatIsWrongTerm::indices, Subscript[
               CM\[CapitalDelta], G1$1, G2$1], Superscript[CM\[CapitalDelta], {
                G1$1, G2$1}]]; Term$]; }]; 
      If[ !FreeQ[Term$, Subscript[CMo, Gs$_]], 
       {Term$ /. Subscript[CMo, Gs$_] :> Module[{}, 
            Message[WhatIsWrongTerm::indices, Subscript[CMo, Gs$], 
              Superscript[CMo, Gs$]]; Term$]; }]; WhatIsWrongTerm::indices2 = 
       "The form of `1` should be `2` or `3`."; 
      If[ !FreeQ[Term$, Superscript[CM\[Delta], {G1$1_, G2$1_}]], 
       {Term$ /. Superscript[CM\[Delta], {G1$1_, G2$1_}] :> 
           Module[{}, Message[WhatIsWrongTerm::indices2, Superscript[
               CM\[Delta], {G1$1, G2$1}], Subscript[Superscript[CM\[Delta], 
                G1$1], G2$1], Superscript[CM\[CapitalDelta], {G1$1, G2$1}]]; 
             Term$]; }]; WhatIsWrongTerm::indices = 
       "The index `1` in `2` should sit inside Superscript, not Power."; 
      If[ !FreeQ[Term$, Subscript[CM\[Delta]^(Q1$1_), Q2$1_]], 
       {Term$ /. Subscript[CM\[Delta]^(Q1$1_), Q2$1_] :> 
           Module[{}, Message[WhatIsWrongTerm::indices, Q1$1, 
              Subscript[Superscript[CM\[Delta], Q1$1], Q2$1]]; Term$]; }]; 
      WhatIsWrongTerm::indices = 
       "The number of indices in `1` should be `2`."; 
      If[ !FreeQ[Term$, Superscript[CM\[CapitalDelta], Gs$_]], 
       Term$ /. Superscript[CM\[CapitalDelta], Gs$_] :> 
          Module[{}, If[Length[Gs$] != 2, {Message[WhatIsWrongTerm::indices, 
                Superscript[CM\[CapitalDelta], Gs$], 2]; }]; Term$]; ]; 
      If[ !FreeQ[Term$, Superscript[CMd, Gs$_]], 
       Term$ /. Superscript[CMd, Gs$_] :> Module[{}, 
           If[Length[Gs$] != 3, {Message[WhatIsWrongTerm::indices, 
                Superscript[CMd, Gs$], 3]; }]; Term$]; ]; 
      If[ !FreeQ[Term$, Superscript[CMf, Gs$_]],
       Term$ /. Superscript[CMf, Gs$_] :> Module[{},
           If[Length[Gs$] != 3, {Message[WhatIsWrongTerm::indices, 
                Superscript[CMf, Gs$], 3]; }]; Term$]; ];
      If[ !FreeQ[Term$, Superscript[If, Gs$_]], 
       WhatIsWrongTerm::term = "Encountered `1`, did you mean `2`?"; 
        Term$ /. Superscript[If, Gs$_] :> Module[{}, 
           Message[WhatIsWrongTerm::term, Superscript["(If)", Gs$], 
             Superscript["I*"*CMf, Gs$]]; Term$]; ];
      If[ !FreeQ[Term$, Superscript[Id, Gs$_]], 
       WhatIsWrongTerm::term = "Encountered `1`, did you mean `2`?"; 
        Term$ /. Superscript[Id, Gs$_] :> Module[{}, 
           Message[WhatIsWrongTerm::term, Superscript["(Id)", Gs$], 
             Superscript["I*"*d, Gs$]]; Term$]; ]; 
      If[ !FreeQ[Term$, Superscript[(ConsTerm$_)*(fdComb$_), Gs$_]], 
       Term$ /. Superscript[fdComb$_, Gs$_] :> 
          If[ !StringFreeQ[ToString[fdComb$], "f"], 
           {WhatIsWrongTerm::term = "Encountered `1`, did you mean `2`?"; 
             Message[WhatIsWrongTerm::term, Superscript["(const*f)", Gs$], 
              Superscript["const*"*f, Gs$]]; }]; 
        Term$ /. Superscript[fdComb$_, Gs$_] :> 
          If[ !StringFreeQ[ToString[fdComb$], "d"], 
           {WhatIsWrongTerm::term = "Encountered `1`, did you mean `2`?"; 
             Message[WhatIsWrongTerm::term, Superscript["(const*d)", Gs$], 
              Superscript["const*"*d, Gs$]]; }]; 
        Term$ /. Superscript[fdComb$_, Gs$_] :> 
          If[ !StringFreeQ[ToString[fdComb$], "\[CapitalDelta]"], 
           {WhatIsWrongTerm::term = "Encountered `1`, did you mean `2`?"; 
             Message[WhatIsWrongTerm::term, Superscript[
               "(const*\[CapitalDelta])", Gs$], Superscript["const*"*
                CM\[CapitalDelta], Gs$]]; }]; ]; If[Head[Term$] == List, 
       If[ !FreeQ[Term$, CMf^(Gs$_)], WhatIsWrongTerm::term = "The structure \
constant f should be defined using Superscript, not Power."; 
          Message[WhatIsWrongTerm::term]; ]; If[ !FreeQ[Term$, d^(Gs$_)], 
         WhatIsWrongTerm::term = "The structure constant d should be defined \
using Superscript, not Power."; Message[WhatIsWrongTerm::term]; ]; 
        If[ !FreeQ[Term$, o^(Gs$_)], WhatIsWrongTerm::term = "The closed \
quark-line o should be defined using Superscript, not Power."; 
          Message[WhatIsWrongTerm::term]; ]; 
        If[ !FreeQ[Term$, CM\[CapitalDelta]^(Gs$_)], 
         WhatIsWrongTerm::term = "The gluon delta function \[CapitalDelta] \
should be defined using Superscript, not Power."; 
          Message[WhatIsWrongTerm::term]; ]; ]; AllSymbols$ = 
       FindSymbols[Term$]; AllInd$ = AllIndices[Term$]; 
      For[symb$ = 1, symb$ <= Length[AllInd$], symb$++, 
       {nIndi$ = Count[AllSymbols$, AllInd$[[symb$]]]; If[nIndi$ > 2, 
           {WhatIsWrongTerm::indices = "The index `1` appeares `2` times in \
term `3`, should at most appeare twice."; Message[WhatIsWrongTerm::indices, 
              AllInd$[[symb$]], nIndi$, Term$]; }]; }; ]; FreeIndices[Term$]]
 
WhatIsWrongTerm /: WhatIsWrongTerm::usage = "WhatIsWrongTerm[Term] checks if \
anything is obviously wrong with one Term, as opposed to a sum of terms, for \
example if Power is used instead of Superscript."

Protect[CMt,o,CMf,d,CM\[Delta],CM\[CapitalDelta],CMdelta,CMDelta]
Protect[AllPermutations, CyclicPermutations,SimpleRules,FDRules,Remove0To2ORules,Remove0To1ORules,Remove0ORules,OTSimpleRules,AllSimpleRules,OTGluonRules,TTGluonContract,OTGluonContract,TGluonContract,OGluonContract,OOGluonContract,O2GluonContract,ODeltaGluonContract,GluonContract,FDToTRules,OTToTRules,OTToTRules,FDToORules,GluonIndices,UpperQuarkIndices,LowerQuarkIndices,SortIndices,FindSymbols,AllPairs,AllIndices,DummyIndicesTerm,DummyIndices,SquareIndices,ReplaceDummyIndices,NTensorIndices,IdentifyParton,WhatIsWrong,WhatIsWrongTerm,SplitConstAndColor,SplitConstAndColorTerm,ContainsFD,ContainsT,ContainsO,ContainsQuarkDelta,ContainsGluonDelta,ContainsColor,OTThenAllSimpleRules,ExpandThenRules,CSimplify,CDot,CNorm,CDotMatrix,CVertex,CGamma,CGammaSymmetryTest,SPMAndSPMInvOfBasisType,GeneralBasis,OrthonormalBasis,OrthogonalBasis,TraceBasis,RemoveFD,NcMin,Verbose,MakeChecks,NMaxGluonCheck,CalledByOther];
End[]
EndPackage[]
(* Private note *)
(* This document was created Nov 7 2012 using NewColorMathRules90.nb. *)
