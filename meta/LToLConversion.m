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

BeginPackage@"LToLConversion`";

create::usage = "";

Begin["`internal`"];

`type`lepton = _Symbol?TreeMasses`IsLepton;
`type`lepton ~ SetAttributes ~ {Protected, Locked};

`type`contribution = Alternatives[
   All,NPointFunctions`FourFermionScalarPenguins,
   NPointFunctions`FourFermionMassiveVectorPenguins,
   NPointFunctions`FourFermionFlavourChangingBoxes
];
`type`contribution ~ SetAttributes ~ {Protected, Locked};

`type`observable = FlexibleSUSYObservable`LToLConversion[
   (lIn:`type`lepton)[_Integer]->(lOut:`type`lepton)[_Integer],
   _Symbol,
   con:`type`contribution,
   CoefficientList->True
];
`type`observable ~ SetAttributes ~ {Protected, Locked};

`cxx`in = "";
`cxx`in // Protect;

setIn[lIn:`type`lepton] := (
   Unprotect@`cxx`in;
   `cxx`in = CConversion`ToValidCSymbolString@lIn;
   Protect@`cxx`in;
);
setIn // Utils`MakeUnknownInputDefinition;
setIn ~ SetAttributes ~ {Protected, Locked};

`cxx`out = "";
`cxx`out // Protect;

setOut[lOut:`type`lepton] := (
   Unprotect@`cxx`out;
   `cxx`out = CConversion`ToValidCSymbolString@lOut;
   Protect@`cxx`out;
);
setOut // Utils`MakeUnknownInputDefinition;
setOut ~ SetAttributes ~ {Protected, Locked};

`cxx`con = "";
`cxx`con // Protect;

setCon[con:`type`contribution] := (
   Unprotect@`cxx`con;
   `cxx`con = CConversion`ToValidCSymbolString@con;
   Protect@`cxx`con;
);
setCon // Utils`MakeUnknownInputDefinition;
setCon ~ SetAttributes ~ {Protected, Locked};

{`cxx`up, `cxx`down} := {
   CConversion`ToValidCSymbolString@SARAH`UpQuark,
   CConversion`ToValidCSymbolString@SARAH`DownQuark
};
{`cxx`up, `cxx`down} ~ SetAttributes ~ {Protected, Locked};

`cxx`classU = "";
`cxx`classD = "";
{`cxx`classU, `cxx`classD} ~ SetAttributes ~ {Protected};

setClass[] := (
   Unprotect@`cxx`classU;
   `cxx`classU = "conversion_"<>`cxx`in<>`cxx`up<>"_to_"<>
         `cxx`out<>`cxx`up<>"_for_"<>`cxx`con;
   Protect@`cxx`classU;

   Unprotect@`cxx`classD;
   `cxx`classD = "conversion_"<>`cxx`in<>`cxx`up<>"_to_"<>
         `cxx`out<>`cxx`down<>"_for_"<>`cxx`con;
   Protect@`cxx`classD;
);
setClass // Utils`MakeUnknownInputDefinition;
setClass ~ SetAttributes ~ {Protected, Locked};

`cxx`prototype :=
   "Eigen::Array<std::complex<double>,10,1>"<>
   " calculate_"<>`cxx`in<>"_to_"<>`cxx`out<>"_for_"<>`cxx`con<>"(\n"<>
   "   int generationIndex1,\n"<>
   "   int generationIndex2,\n"<>
   "   const " <> FlexibleSUSY`FSModelName <>
      "_l_to_l_conversion::Nucleus nucleus,\n" <>
   "   const " <> FlexibleSUSY`FSModelName <>
      "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)";
`cxx`prototype ~ SetAttributes ~ {Protected, Locked};

create::usage =
"@brief Main entrance point for the canculation.";
create[obs:`type`observable] :=
Module[
   {
      npfVertices, npfHeader, npfDefinition,
      calculatePrototype = getPrototype[lIn->lOut,contribution],
      calculateDefinition
   },
   setIn@lIn;
   setOut@lOut;
   setCon@con;
   setClass[];

   {npfVertices, npfHeader, npfDefinition} = `npf`create[obs];

   calculateDefinition = `cxx`prototype <> StringReplace[" {

   context_base context {model};
   // get Fermi constant from Les Houches input file
   const auto GF = qedqcd.displayFermiConstant();
   constexpr bool discard_SM_contributions = false;
   const auto photon_penguin = @photon_penguin_name@ (generationIndex1,
      generationIndex2,
      model,
      discard_SM_contributions);

   // translate from the convention of Hisano, Moroi & Tobe to Kitano, Koike & Okada
   // Hisano defines form factors A2 through a matrix element in eq. 14
   // Kitano uses a lagrangian with F_munu. There is a factor of 2 from translation
   // because Fmunu = qeps - eps q

   const auto A2L = -0.5 * photon_penguin[2]/(4.*GF/sqrt(2.));
   const auto A2R = -0.5 * photon_penguin[3]/(4.*GF/sqrt(2.));

   // penguins
   // 2 up and 1 down quark in proton (gp couplings)
   // 1 up and 2 down in neutron (gn couplings)

   // mediator: massless vector
   // construct 4-fermion operators from A1 form factors
   // i q^2 A1 * (- i gmunu/q^2) * (-i Qq e) = GF/sqrt2 * gpV
   // VP

   const auto vcU = vectorCurrent<
      typename fields::Fu::lorentz_conjugate,
      fields::Fu,
      typename fields::VP
   >(model);
   const auto vcD = vectorCurrent<
      typename fields::Fd::lorentz_conjugate,
      fields::Fd,
      typename fields::VP
   >(model);

   // the A1 term if the factor in front of q^2, the photon propagator is -1/q^2, we need only factor -1
   auto gpLV = -sqrt(2.0)/GF * photon_penguin[0] * (2.*vcU + vcD);
   auto gpRV = -sqrt(2.0)/GF * photon_penguin[1] * (2.*vcU + vcD);
   auto gnLV = -sqrt(2.0)/GF * photon_penguin[0] * (vcU + 2.*vcD);
   auto gnRV = -sqrt(2.0)/GF * photon_penguin[1] * (vcU + 2.*vcD);

   // all contributions

   auto npfU = @namespace_npf@@class_U@(
      model,
      std::array<int,4>{generationIndex1, 0, generationIndex2, 0},
      std::array<Eigen::Vector4d, 0>{}
   );

   auto npfD = @namespace_npf@@class_D@(
      model,
      std::array<int,4>{generationIndex1, 0, generationIndex2, 0},
      std::array<Eigen::Vector4d, 0>{}
   );

   auto npfS = @namespace_npf@@class_D@(
      model,
      std::array<int,4>{generationIndex1, 1, generationIndex2, 1},
      std::array<Eigen::Vector4d, 0>{}
   );

   // PDG 2018 data
   double m_p = 0.938272081, m_n = 0.939565413;
   //data from my notes
   double m_init = context.mass<fields::"<>`cxx`in<>">({generationIndex1});
   double m_u = context.mass<fields::"<>`cxx`up<>">({0});
   double m_d = context.mass<fields::"<>`cxx`down<>">({0});
   double m_s = context.mass<fields::"<>`cxx`down<>">({1});

   double GSpu = 0.021*m_p/m_u, GSpd = 0.041*m_p/m_d, GSps = 0.043*m_p/m_s;
   double GSnu = 0.019*m_n/m_u, GSnd = 0.045*m_n/m_d, GSns = 0.043*m_n/m_s;

   double GVpu = 2.,            GVpd = 1.;
   double GVnu = 1.,            GVnd = 2.;

   double GTpu = 0.77,          GTpd = -0.23,         GTps = 0.008;
   double GTnu = 0.77,          GTnd = -0.23,         GTns = 0.008;

   //minus because of descending order in FormCalc spinor chains
   auto CSLu = -( npfU.at(0)+npfU.at(1) )/2.;
   auto CSRu = -( npfU.at(2)+npfU.at(3) )/2.;
   auto CSLd = -( npfD.at(0)+npfD.at(1) )/2.;
   auto CSRd = -( npfD.at(2)+npfD.at(3) )/2.;
   auto CSLs = -( npfS.at(0)+npfS.at(1) )/2.;
   auto CSRs = -( npfS.at(2)+npfS.at(3) )/2.;

   //minus because of descending order in FormCalc spinor chains
   auto CVLu = -( npfU.at(4)+npfU.at(5) )/2.;
   auto CVRu = -( npfU.at(6)+npfU.at(7) )/2.;
   auto CVLd = -( npfD.at(4)+npfD.at(5) )/2.;
   auto CVRd = -( npfD.at(6)+npfD.at(7) )/2.;

   //plus because of descending order in FormCalc spinor chains and definition of tensor operators
   auto CTLu = npfU.at(8);
   auto CTRu = npfU.at(9);
   auto CTLd = npfD.at(8);
   auto CTRd = npfD.at(8);
   auto CTLs = npfS.at(8);
   auto CTRs = npfS.at(8);

   gpLV += (-sqrt(2.0)/GF)*( GVpu*CVLu + GVpd*CVLd );
   gpRV += (-sqrt(2.0)/GF)*( GVpu*CVRu + GVpd*CVRd );
   gnLV += (-sqrt(2.0)/GF)*( GVnu*CVLu + GVnd*CVLd );
   gnRV += (-sqrt(2.0)/GF)*( GVnu*CVRu + GVnd*CVRd );

   //scalar contribution from scalar coefficients
   auto gpLS = (-sqrt(2.0)/GF)*( GSpu*CSLu + GSpd*CSLd + GSps*CSLs );
   auto gpRS = (-sqrt(2.0)/GF)*( GSpu*CSRu + GSpd*CSRd + GSps*CSRs );
   auto gnLS = (-sqrt(2.0)/GF)*( GSnu*CSLu + GSnd*CSLd + GSns*CSLs );
   auto gnRS = (-sqrt(2.0)/GF)*( GSnu*CSRu + GSnd*CSRd + GSns*CSRs );

   //scalar contribution from tensor coefficients
   gpLS += (-sqrt(2.0)/GF)*(2*m_init/m_p)*( GTpu*CTLu + GTpd*CTLd + GTps*CTLs );
   gpRS += (-sqrt(2.0)/GF)*(2*m_init/m_p)*( GTpu*CTRu + GTpd*CTRd + GTps*CTRs );
   gnLS += (-sqrt(2.0)/GF)*(2*m_init/m_n)*( GTnu*CTLu + GTnd*CTLd + GTns*CTLs );
   gnRS += (-sqrt(2.0)/GF)*(2*m_init/m_n)*( GTnu*CTRu + GTnd*CTRd + GTns*CTRs );

   const auto nuclear_form_factors = get_overlap_integrals(nucleus, qedqcd);

   const auto left = A2R*nuclear_form_factors.D
      + gpLV*nuclear_form_factors.Vp + gnLV*nuclear_form_factors.Vn
      + gpRS*nuclear_form_factors.Sp + gnRS*nuclear_form_factors.Sn;

   const auto right = A2L*nuclear_form_factors.D
      + gpRV*nuclear_form_factors.Vp + gnRV*nuclear_form_factors.Vn
      + gpLS*nuclear_form_factors.Sp + gnLS*nuclear_form_factors.Sn;

   // eq. 14 of Kitano, Koike and Okada
   const double conversion_rate = 2.*pow(GF,2)*(std::norm(left) + std::norm(right));

   // normalize to capture
   const double capture_rate = get_capture_rate(nucleus);

   Eigen::Array<std::complex<double>,10,1> res =
      Eigen::Array<std::complex<double>,10,1>::Zero();

   res(0) = conversion_rate/capture_rate;
   return res;
}
",
   {
      "@photon_penguin_name@"->"calculate_"<>`cxx`in<>"_"<>`cxx`out<>"_"<>
         CXXDiagrams`CXXNameOfField@SARAH`Photon<>"_form_factors",
      "@namespace_npf@"->FlexibleSUSY`FSModelName<>"_cxx_diagrams::npointfunctions::",
      "@class_U@"->`cxx`classU,
      "@class_D@"->`cxx`classD
   }];

   {
      npfVertices,
      {npfHeader,npfDefinition},
      {`cxx`prototype <> ";", calculateDefinition}
   }
];

create[list:{__}] :=
   {
      DeleteDuplicates[ Join@@#[[All,1]] ],
      {
         #[[1,2,1]],
         StringJoin@Riffle[#[[All,2,2]], "\n\n"]
      },
      {
         StringJoin@Riffle[#[[All,3,1]], "\n\n"],
         StringJoin@Riffle[#[[All,3,2]], "\n\n"]
      }
   } & [create /@ list];

create // Utils`MakeUnknownInputDefinition;

`npf`create[`type`observable] :=
Module[
   {
      parsedCon = Switch[con,
         All, {
               NPointFunctions`FourFermionMassiveVectorPenguins,
               NPointFunctions`FourFermionScalarPenguins,
               NPointFunctions`FourFermionFlavourChangingBoxes
            },
         _, con
      ],
      npfU, npfD, (*@note objects, generated by NPointFunction*)

      l=SARAH`Lorentz, p=SARAH`Mom, m=SARAH`Mass, (*@note synonyms*)
      dc = NPointFunctions`internal`dc, (*@note current name for dirac chain*)

      fiG, foG, uiG, uoG, diG, doG, (*@note particle | inc/out | generation*)
      regulator, (*@note arbitrary sqr(3-momenta) of quarks*)
      inner, sp, dim6,

      codeU,codeD
   },

   Print["<<npf<< calculation for ",`cxx`in," to ",`cxx`out," conversion started ..."];

   {npfU, npfD} = NPointFunctions`NPointFunction[
      {lIn,#},{lOut,#},
      NPointFunctions`OnShellFlag -> True,
      NPointFunctions`UseCache -> False,
      NPointFunctions`ZeroExternalMomenta -> NPointFunctions`OperatorsOnly,
      NPointFunctions`KeepProcesses -> parsedCon] &/@ {SARAH`UpQuark,SARAH`DownQuark};

   {fiG, uiG, foG, uoG} = Flatten@NPointFunctions`internal`getProcess@npfU;
   {fiG, diG, foG, doG} = Flatten@NPointFunctions`internal`getProcess@npfD;
   regulator = m@fiG^2;

   inner = SARAH`sum[i_,1,4,SARAH`g[i_,i_]*p[#1,i_]*p[#2,i_]]&;
   {npfU, npfD} = {npfU, npfD} //. {
         inner[fiG,foG] :> m@fiG^2,
         inner[uiG,fiG] :> m@fiG*Sqrt[m@uiG^2+regulator],
         inner[uoG,fiG] :> inner[uiG,fiG],
         inner[uoG,foG] :> m@fiG^2/2+inner[uiG,fiG],
         inner[diG,fiG] :> m@fiG*Sqrt[m@diG^2+regulator],
         inner[doG,fiG] :> inner[diG,fiG],
         inner[doG,foG] :> m@fiG^2/2+inner[diG,fiG]
      };

   sp[particle:_,num:_Integer] := SARAH`DiracSpinor[#,p@num,m@#] &@
      particle@{Symbol["SARAH`gt"<>ToString@num]};

   dim6[i_,o_,q_] := {
      (*@note 6 means PR, 7 means PL.*)
      "S_LL" -> dc[o~sp~3,7,i~sp~1] dc[q~sp~4,7,q~sp~2],
      "S_LR" -> dc[o~sp~3,7,i~sp~1] dc[q~sp~4,6,q~sp~2],
      "S_RL" -> dc[o~sp~3,6,i~sp~1] dc[q~sp~4,7,q~sp~2],
      "S_RR" -> dc[o~sp~3,6,i~sp~1] dc[q~sp~4,6,q~sp~2],
      (*@note names are correct, one just need to commute projectors with
       *Dirac matrices. It changes 6 to 7 or 7 to 6.*)
      "V_LL" -> dc[o~sp~3,6,l@1,i~sp~1] dc[q~sp~4,6,l@1,q~sp~2],
      "V_LR" -> dc[o~sp~3,6,l@1,i~sp~1] dc[q~sp~4,7,l@1,q~sp~2],
      "V_RL" -> dc[o~sp~3,7,l@1,i~sp~1] dc[q~sp~4,6,l@1,q~sp~2],
      "V_RR" -> dc[o~sp~3,7,l@1,i~sp~1] dc[q~sp~4,7,l@1,q~sp~2],
      (*@note Minus, because FormCalc`s -6,Lor[1],Lor[2] is ours
       *-I*sigma[1,2] (according to FC definition of antisymmetrization), when
       *taking this twice we get I*I=-1. FC cites [Ni05] for Fierz identities,
       *where our conventions are used, but in FC manual on the page 20
       *weird convention for sigma_munu is shown.*)
      "minus_T_LL" -> dc[o~sp~3,-7,l@1,l@2,i~sp~1] dc[q~sp~4,-7,l@1,l@2,q~sp~2],
      "minus_T_RR" -> dc[o~sp~3,-6,l@1,l@2,i~sp~1] dc[q~sp~4,-6,l@1,l@2,q~sp~2]
   };
   npfU = npfU~WilsonCoeffs`InterfaceToMatching~dim6[lIn,lOut,SARAH`UpQuark];
   npfD = npfD~WilsonCoeffs`InterfaceToMatching~dim6[lIn,lOut,SARAH`DownQuark];

   Print[">>npf>> calculation for ",`cxx`in," to ",`cxx`out," conversion done."];

   Print["<<npf<< c++ code calculation for ",`cxx`in," to ",`cxx`out," conversion started ..."];

   codeU = NPointFunctions`CreateCXXFunctions[
      npfU,
      `cxx`classU,
      SARAH`Delta,
      dim6[lIn,lOut,SARAH`UpQuark]
   ][[2]];
   codeD = NPointFunctions`CreateCXXFunctions[
      npfD,
      `cxx`classD,
      SARAH`Delta,
      dim6[lIn,lOut,SARAH`DownQuark]
   ][[2]];

   Print[">>npf>> c++ code calculation for ",`cxx`in," to ",`cxx`out," conversion done."];
   {
      DeleteDuplicates@Join[
         NPointFunctions`VerticesForNPointFunction@npfU,
         NPointFunctions`VerticesForNPointFunction@npfD
      ],
      NPointFunctions`CreateCXXHeaders[],
      codeU<>"\n\n"<>codeD
   }
];
`npf`create // Utils`MakeUnknownInputDefinition;
`npf`create ~ SetAttributes ~ {Locked,Protected};

End[];
EndPackage[];
