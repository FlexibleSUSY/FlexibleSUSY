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

BeginPackage["FToFConversionInNucleus`",
    {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "CXXDiagrams`", "Utils`"}
];

FToFConversionInNucleusCreateInterface::usage = "";

Begin["`Private`"];

FToFConversionInNucleusCreateInterface[inFermion_ -> outFermion_] :=
Module[
   {
      inName = CXXNameOfField@inFermion,
      outName = CXXNameOfField@outFermion,
      prototype, definition, n=2
   },

   Utils`AssertWithMessage[
      Count[GetVectorBosons[], el_ /; IsMassless[el]] === 2,
      "We assume that there are only 2 massless vector bosons, i.e. photon and gluon"
   ];

   Utils`AssertWithMessage[
      TreeMasses`GetDimension[SARAH`UpQuark] === 3,
      "We assume that up quarks have three generations"
   ];

   Utils`AssertWithMessage[
      TreeMasses`GetDimension[SARAH`DownQuark] === 3,
      "We assume that down quarks have three generations"
   ];

   If[TreeMasses`GetDimension[inFermion] =!= 1,n=n+1];
   If[TreeMasses`GetDimension[outFermion] =!= 1,n=n+1];

   prototype =
      "double calculate_"<>inName<>"_to_"<>outName<>"_in_nucleus (\n"<>
         If[TreeMasses`GetDimension[inFermion] =!= 1,"int generationIndex1, "," "]<>
         If[TreeMasses`GetDimension[outFermion] =!= 1," int generationIndex2, "," "]<>
         "const " <> FlexibleSUSY`FSModelName <> "_f_to_f_conversion::Nucleus nucleus,\n" <>
         "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)";

   definition =
      prototype<>
      " {\n"<>
      IndentText[
         "\n" <>
         "context_base context {model};\n" <>
         "// get Fermi constant from Les Houches input file\n" <>
         "const auto GF = qedqcd.displayFermiConstant();\n" <>
         "constexpr bool discard_SM_contributions = false;\n" <>

         "const auto photon_penguin = calculate_"<>inName<>"_"<>outName<>"_" <> CXXNameOfField[SARAH`Photon]<>"_form_factors (" <>
         If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1, ", " "]<>
         If[TreeMasses`GetDimension[outFermion] =!= 1, " generationIndex2, ", " "]<>"model, discard_SM_contributions);\n" <>

         "\n// translate from the convention of Hisano, Moroi & Tobe to Kitano, Koike & Okada\n" <>
         (* TODO: check the statement below *)
         "// Hisano defines form factors A2 through a matrix element in eq. 14\n" <>
         "// Kitano uses a lagrangian with F_munu. There is a factor of 2 from translation\n" <>
         "// because Fmunu = qeps - eps q\n" <>
         "const auto A2L = -0.5 * photon_penguin[2]/(4.*GF/sqrt(2.));\n" <>
         "const auto A2R = -0.5 * photon_penguin[3]/(4.*GF/sqrt(2.));\n" <>

         "\n// ------ penguins ------\n" <>
         "// 2 up and 1 down quark in proton (gp couplings)\n" <>
         "// 1 up and 2 down in neutron (gn couplings)\n" <>

         "\n// mediator: massless vector\n" <>
         "\n// construct 4-fermion operators from A1 form factors\n" <>
         "// i q^2 A1 * (- i gmunu/q^2) * (-i Qq e) = GF/sqrt2 * gpV\n" <>
         "// VP\n" <>

         "const auto uEMVectorCurrent =\n" <>
             IndentText["vectorCurrent<typename Fu::lorentz_conjugate, Fu, typename VP::lorentz_conjugate>(model);\n"] <>
         "const auto dEMVectorCurrent =\n" <>
             IndentText["vectorCurrent<typename Fd::lorentz_conjugate, Fd, typename VP::lorentz_conjugate>(model);\n"] <>
         "\n" <>

         "// the A1 term if the factor in front of q^2, the photon propagator is -1/q^2, we need only factor -1\n" <>
         "auto gpLV = -sqrt(2.0)/GF * photon_penguin[0] * (2.*uEMVectorCurrent + dEMVectorCurrent);\n" <>
         "auto gpRV = -sqrt(2.0)/GF * photon_penguin[1] * (2.*uEMVectorCurrent + dEMVectorCurrent);\n" <>
         "auto gnLV = -sqrt(2.0)/GF * photon_penguin[0] * (uEMVectorCurrent + 2.*dEMVectorCurrent);\n" <>
         "auto gnRV = -sqrt(2.0)/GF * photon_penguin[1] * (uEMVectorCurrent + 2.*dEMVectorCurrent);\n" <>

         "\n// all contributions\n" <>

         "auto "<>FlexibleSUSY`FSModelName<>"_npf_up = "<>FlexibleSUSY`FSModelName<>
         "_cxx_diagrams::npointfunctions::zpinguins_u"<>ToString@inFermion<>ToString@outFermion<>"_1loop("<>
         "model, std::array<int,"<>ToString@n<>">{"<>
            If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1,", ""]<>"0,"<>
            If[TreeMasses`GetDimension[outFermion] =!= 1, "generationIndex2,", ""]<>"0}, std::array<Eigen::Vector4d, 0>{});\n"<>

         "auto "<>FlexibleSUSY`FSModelName<>"_npf_down = "<>FlexibleSUSY`FSModelName<>
         "_cxx_diagrams::npointfunctions::zpinguins_d"<>inName<>outName<>"_1loop("<>
         "model, std::array<int,"<>ToString@n<>">{"<>
            If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1,", ""]<>"0,"<>
            If[TreeMasses`GetDimension[outFermion] =!= 1, "generationIndex2,", ""]<>"0}, std::array<Eigen::Vector4d, 0>{});\n"<>

         "auto "<>FlexibleSUSY`FSModelName<>"_npf_strange = "<>FlexibleSUSY`FSModelName<>
         "_cxx_diagrams::npointfunctions::zpinguins_d"<>inName<>outName<>"_1loop("<>
         "model, std::array<int,"<>ToString@n<>">{"<>
            If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1,", ""]<>"1,"<>
            If[TreeMasses`GetDimension[outFermion] =!= 1, "generationIndex2,", ""]<>"1}, std::array<Eigen::Vector4d, 0>{});\n"<>
         "\n"<>
         "// PDG 2018 data\n"<>
         "double m_p = 0.938272081, m_n = 0.939565413;\n"<>
         "//data from my notes\n"<>
         "double m_init = context.mass<"<>ToString@inFermion<>">({generationIndex1});\n"<>
         "double m_u = context.mass<"<>ToString@SARAH`UpQuark<>">({0});\n"<>
         "double m_d = context.mass<"<>ToString@SARAH`DownQuark<>">({0});\n"<>
         "double m_s = context.mass<"<>ToString@SARAH`DownQuark<>">({1});\n"<>
         "\n"<>
         "double GSpu = 0.021*m_p/m_u, GSpd = 0.041*m_p/m_d, GSps = 0.043*m_p/m_s;\n"<>
         "double GSnu = 0.019*m_n/m_u, GSnd = 0.045*m_n/m_d, GSns = 0.043*m_n/m_s;\n"<>
         "\n"<>
         "double GVpu = 2.,            GVpd = 1.;\n"<>
         "double GVnu = 1.,            GVnd = 2.;\n"<>
         "\n"<>
         "double GTpu = 0.77,          GTpd = -0.23,         GTps = 0.008;\n"<>
         "double GTnu = 0.77,          GTnd = -0.23,         GTns = 0.008;\n"<>
         "\n"<>
         "//minus because of descending order in FormCalc spinor chains\n"<>
         "std::complex<double> CSLu = -( "<>FlexibleSUSY`FSModelName<>"_npf_up.at(0)+"<>FlexibleSUSY`FSModelName<>"_npf_up.at(1) )/2.;\n"<>
         "std::complex<double> CSRu = -( "<>FlexibleSUSY`FSModelName<>"_npf_up.at(2)+"<>FlexibleSUSY`FSModelName<>"_npf_up.at(3) )/2.;\n"<>
         "std::complex<double> CSLd = -( "<>FlexibleSUSY`FSModelName<>"_npf_down.at(0)+"<>FlexibleSUSY`FSModelName<>"_npf_down.at(1) )/2.;\n"<>
         "std::complex<double> CSRd = -( "<>FlexibleSUSY`FSModelName<>"_npf_down.at(2)+"<>FlexibleSUSY`FSModelName<>"_npf_down.at(3) )/2.;\n"<>
         "std::complex<double> CSLs = -( "<>FlexibleSUSY`FSModelName<>"_npf_strange.at(0)+"<>FlexibleSUSY`FSModelName<>"_npf_strange.at(1) )/2.;\n"<>
         "std::complex<double> CSRs = -( "<>FlexibleSUSY`FSModelName<>"_npf_strange.at(2)+"<>FlexibleSUSY`FSModelName<>"_npf_strange.at(3) )/2.;\n"<>
         "\n"<>
         "//minus because of descending order in FormCalc spinor chains\n"<>
         "std::complex<double> CVLu = -( "<>FlexibleSUSY`FSModelName<>"_npf_up.at(4)+"<>FlexibleSUSY`FSModelName<>"_npf_up.at(5) )/2.;\n"<>
         "std::complex<double> CVRu = -( "<>FlexibleSUSY`FSModelName<>"_npf_up.at(6)+"<>FlexibleSUSY`FSModelName<>"_npf_up.at(7) )/2.;\n"<>
         "std::complex<double> CVLd = -( "<>FlexibleSUSY`FSModelName<>"_npf_down.at(4)+"<>FlexibleSUSY`FSModelName<>"_npf_down.at(5) )/2.;\n"<>
         "std::complex<double> CVRd = -( "<>FlexibleSUSY`FSModelName<>"_npf_down.at(6)+"<>FlexibleSUSY`FSModelName<>"_npf_down.at(7) )/2.;\n"<>
         "\n"<>
         "//plus because of descending order in FormCalc spinor chains and definition of tensor operators\n"<>
         "std::complex<double> CTLu = +"<>FlexibleSUSY`FSModelName<>"_npf_up.at(8);\n"<>
         "std::complex<double> CTRu = +"<>FlexibleSUSY`FSModelName<>"_npf_up.at(9);\n"<>
         "std::complex<double> CTLd = +"<>FlexibleSUSY`FSModelName<>"_npf_down.at(8);\n"<>
         "std::complex<double> CTRd = +"<>FlexibleSUSY`FSModelName<>"_npf_down.at(8);\n"<>
         "std::complex<double> CTLs = +"<>FlexibleSUSY`FSModelName<>"_npf_strange.at(8);\n"<>
         "std::complex<double> CTRs = +"<>FlexibleSUSY`FSModelName<>"_npf_strange.at(8);\n"<>
         "\n"<>
         "gpLV += (-sqrt(2.0)/GF)*( GVpu*CVLu + GVpd*CVLd );\n" <>
         "gpRV += (-sqrt(2.0)/GF)*( GVpu*CVRu + GVpd*CVRd );\n" <>
         "gnLV += (-sqrt(2.0)/GF)*( GVnu*CVLu + GVnd*CVLd );\n" <>
         "gnRV += (-sqrt(2.0)/GF)*( GVnu*CVRu + GVnd*CVRd );\n" <>
         "\n//scalar contribution from scalar coefficients\n"<>
         "std::complex<double> gpLS = (-sqrt(2.0)/GF)*( GSpu*CSLu + GSpd*CSLd + GSps*CSLs );\n" <>
         "std::complex<double> gpRS = (-sqrt(2.0)/GF)*( GSpu*CSRu + GSpd*CSRd + GSps*CSRs );\n" <>
         "std::complex<double> gnLS = (-sqrt(2.0)/GF)*( GSnu*CSLu + GSnd*CSLd + GSns*CSLs );\n" <>
         "std::complex<double> gnRS = (-sqrt(2.0)/GF)*( GSnu*CSRu + GSnd*CSRd + GSns*CSRs );\n" <>
         "\n//scalar contribution from tensor coefficients\n"<>
         "gpLS += (-sqrt(2.0)/GF)*(2*m_init/m_p)*( GTpu*CTLu + GTpd*CTLd + GTps*CTLs );\n" <>
         "gpRS += (-sqrt(2.0)/GF)*(2*m_init/m_p)*( GTpu*CTRu + GTpd*CTRd + GTps*CTRs );\n" <>
         "gnLS += (-sqrt(2.0)/GF)*(2*m_init/m_n)*( GTnu*CTLu + GTnd*CTLd + GTns*CTLs );\n" <>
         "gnRS += (-sqrt(2.0)/GF)*(2*m_init/m_n)*( GTnu*CTRu + GTnd*CTRd + GTns*CTRs );\n" <>
         "\nconst auto nuclear_form_factors =\n" <>
                   IndentText@"get_overlap_integrals(nucleus, qedqcd);\n"<>

                "\nconst auto left =\n" <> IndentText[
                   "A2R*nuclear_form_factors.D\n" <>
                   "+ gpLV*nuclear_form_factors.Vp\n" <>
                   "+ gnLV*nuclear_form_factors.Vn\n" <>
                   "+ gpRS*nuclear_form_factors.Sp\n" <>
                   "+ gnRS*nuclear_form_factors.Sn"
                ] <> ";\n" <>
                "\nconst auto right =\n" <> IndentText[
                   "A2L*nuclear_form_factors.D\n" <>
                   "+ gpRV*nuclear_form_factors.Vp\n" <>
                   "+ gnRV*nuclear_form_factors.Vn\n" <>
                   "+ gpLS*nuclear_form_factors.Sp\n" <>
                   "+ gnLS*nuclear_form_factors.Sn"
                ] <> ";\n" <>

                "\n// eq. 14 of Kitano, Koike and Okada\n" <>
                "const double conversion_rate = 2.*pow(GF,2)*(std::norm(left) + std::norm(right));\n" <>

                 "\n// normalize to capture\n" <>
                 "const double capture_rate = get_capture_rate(nucleus);\n\n" <>

                 "return conversion_rate/capture_rate;\n"
            ] <>
            "}\n";

        {prototype <> ";", definition}
];

End[];
EndPackage[];
