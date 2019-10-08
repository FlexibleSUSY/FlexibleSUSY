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

Begin["Private`"];

FToFConversionInNucleusCreateInterface[inFermion_ -> outFermion_] :=
    Module[{prototype, definition},

       Utils`AssertWithMessage[
          Count[GetVectorBosons[], el_ /; IsMassless[el]] === 2,
          "We assume that there are only 2 massless vector bosons, i.e. photon and gluon"
       ];

        prototype =
            "double calculate_" <> CXXNameOfField[inFermion] <> "_to_" <>
                CXXNameOfField[outFermion] <> "_in_nucleus (\n" <>
                If[TreeMasses`GetDimension[inFermion] =!= 1,
                    "int generationIndex1, ",
                    " "
                ] <>
                If[TreeMasses`GetDimension[outFermion] =!= 1,
                    " int generationIndex2, ",
                    " "
                ] <>
                "const " <> FlexibleSUSY`FSModelName <> "_f_to_f_conversion::Nucleus nucleus,\n" <>
                "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)";

        definition =
            prototype <> " {\n" <>
            IndentText[
                "\n" <>
                "context_base context {model};\n" <>

                "// get Fermi constant from Les Houches input file\n" <>
                "const auto GF = qedqcd.displayFermiConstant();\n" <>
                "constexpr bool discard_SM_contributions = false;\n" <>

                "const auto photon_penguin = calculate_" <> CXXNameOfField[inFermion] <> "_" <>
                    CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[SARAH`Photon] <> "_form_factors (" <>
                    If[TreeMasses`GetDimension[inFermion] =!= 1, "generationIndex1, ", " "] <>
                    If[TreeMasses`GetDimension[outFermion] =!= 1, " generationIndex2, ", " "] <>
                    "model, discard_SM_contributions);\n" <>

                "\n// translate from the convention of Hisano, Moroi & Tobe to Kitano, Koike & Okada\n" <>
                (* TODO: check the statement below *)
                "// Hisano defines form factors A2 through a matrix element in eq. 14\n" <>
                "// Kitano uses a lagrangian with F_munu. There is a factor of 2 from translation\n" <>
                "// because Fmunu = qeps - eps q\n" <>
                "const auto A2L = -0.5 * photon_penguin[2]/(4.*GF/sqrt(2.));\n" <>
                "const auto A2R = -0.5 * photon_penguin[3]/(4.*GF/sqrt(2.));\n" <>

                (* TODO: remove *)
                "std::cout << std::setprecision(15);\n" <>
                (*"std::cout << \"A1 \" << photon_penguin[0] << ' ' << photon_penguin[1] << '\\n';\n" <>
                "std::cout << \"A2 \" << photon_penguin[2] << ' ' << photon_penguin[3] << '\\n';\n" <>*)

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
                (* TODO: remove *)
                (*"std::cout << \"A 4-fermion \" << -sqrt(2.0)/GF * photon_penguin[0] * uEMVectorCurrent << ' ' << -sqrt(2.0)/GF * photon_penguin[1] * uEMVectorCurrent << '\\n';\n" <>*)
                "auto gnLV = -sqrt(2.0)/GF * photon_penguin[0] * (uEMVectorCurrent + 2.*dEMVectorCurrent);\n" <>
                "auto gnRV = -sqrt(2.0)/GF * photon_penguin[1] * (uEMVectorCurrent + 2.*dEMVectorCurrent);\n" <>

                "\n// mediator: massive vector\n" <>
                StringJoin @ Map[
                    ("\n// " <> CXXNameOfField[#] <> "\n" <>

                        "const auto " <> CXXNameOfField[#] <> "_FF = " <>
                           "calculate_" <> CXXNameOfField[inFermion] <> "_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[#] <>
                           "_form_factors (generationIndex1,  generationIndex2, model);\n" <>

                        "const auto " <> CXXNameOfField[#] <> "_penguin = " <>
                           "create_massive_penguin_amp<" <> CXXNameOfField[#] <> ">(" <>
                           CXXNameOfField[#] <> "_FF, " <>
                           "model, qedqcd);\n" <>

                        (* TODO" remove *)
                        (*"std::cout << \"Z 4-fermion \" << VZ_penguin[0] / (-sqrt(2.0)/GF) * 16*Pi*Pi << ' ' << VZ_penguin[1]/ (-sqrt(2.0)/GF) * 16*Pi*Pi << '\\n';\n" <>*)
                        "std::cout << \"wkotlarski: gLVu \" << "<>CXXNameOfField[#]<>"_penguin[0] / (-sqrt(2.0)/GF) << \"\\n\";\n" <>
                        "std::cout << \"wkotlarski: gRVu \" << "<>CXXNameOfField[#]<>"_penguin[1] / (-sqrt(2.0)/GF) << \"\\n\";\n" <>
                        "std::cout << \"wkotlarski: gLVd \" << "<>CXXNameOfField[#]<>"_penguin[2] / (-sqrt(2.0)/GF) << \"\\n\";\n" <>
                        "std::cout << \"wkotlarski: gRVd \" << "<>CXXNameOfField[#]<>"_penguin[3] / (-sqrt(2.0)/GF) << \"\\n\";\n" <>
                        (*initialize my stuff*)
                        "auto "<>FlexibleSUSY`FSModelName<>"_npf_up = "<>FlexibleSUSY`FSModelName<>
                        "_cxx_diagrams::npointfunctions::zpinguins_u"<>ToString@inFermion<>ToString@outFermion<>"_1loop("<>
                        "model, std::array<int, 4>{generationIndex1, 0, generationIndex2, 0}, std::array<Eigen::Vector4d, 0>{});\n"<>
                        "auto "<>FlexibleSUSY`FSModelName<>"_npf_down = "<>FlexibleSUSY`FSModelName<>
                        "_cxx_diagrams::npointfunctions::zpinguins_d"<>ToString@inFermion<>ToString@outFermion<>"_1loop("<>
                        "model, std::array<int, 4>{generationIndex1, 0, generationIndex2, 0}, std::array<Eigen::Vector4d, 0>{});\n"<>
                        (*print my stuff*)
                        "std::cout << \"uukhas:     gLVu \" << ( "<>
                        FlexibleSUSY`FSModelName<>"_npf_up.at(4)+"<>FlexibleSUSY`FSModelName<>"_npf_up.at(5) )/2. << \"\\n\";\n" <>
                        "std::cout << \"uukhas:     gRVu \" << ( "<>
                        FlexibleSUSY`FSModelName<>"_npf_up.at(6)+"<>FlexibleSUSY`FSModelName<>"_npf_up.at(7) )/2.<< \"\\n\";\n" <>
                        "std::cout << \"uukhas:     gLVd \" << ( "<>
                        FlexibleSUSY`FSModelName<>"_npf_down.at(4)+"<>FlexibleSUSY`FSModelName<>"_npf_down.at(5) )/2. << \"\\n\";\n" <>
                        "std::cout << \"uukhas:     gRVd \" << ( "<>
                        FlexibleSUSY`FSModelName<>"_npf_down.at(6)+"<>FlexibleSUSY`FSModelName<>"_npf_down.at(7) )/2. << \"\\n\";\n" <>
                        "gpLV += 2.*" <> CXXNameOfField[#] <> "_penguin[0] + "    <> CXXNameOfField[#] <> "_penguin[2];\n" <>
                        "gpRV += 2.*" <> CXXNameOfField[#] <> "_penguin[1] + "    <> CXXNameOfField[#] <> "_penguin[3];\n" <>
                        "gnLV += "    <> CXXNameOfField[#] <> "_penguin[0] + 2.*" <> CXXNameOfField[#] <> "_penguin[2];\n" <>
                        "gnRV += "    <> CXXNameOfField[#] <> "_penguin[1] + 2.*" <> CXXNameOfField[#] <> "_penguin[3];\n")&,

                    (* create a list of massive, electrically neutral gauge bosons *)
                    Select[GetVectorBosons[],
                       !(IsMassless[#] || IsElectricallyCharged[#] || ColorChargedQ[#])&
                    ]
                ] <>

                (* TODO: add contributions from scalar penguins *)
                "\n// mediator: massive scalar\n" <>

                "std::complex<double> gpLS {0.0};\n" <>
                "std::complex<double> gpRS {0.0};\n" <>
                "std::complex<double> gnLS {0.0};\n" <>
                "std::complex<double> gnRS {0.0};\n" <>

                    StringJoin @ Map[
                      ("\n// " <> CXXNameOfField[#] <> "\n")& (*<>

                          "const auto " <> CXXNameOfField[#] <> "_FF = " <>
                          "calculate_" <> CXXNameOfField[inFermion] <> "_" <> CXXNameOfField[outFermion] <> "_" <> CXXNameOfField[#] <>
                          "_form_factors (generationIndex1,  generationIndex2, model);\n" <>

                          "const auto " <> CXXNameOfField[#] <> "_penguin = " <>
                          "create_massive_penguin_amp<" <> CXXNameOfField[#] <> ">(" <>
                          CXXNameOfField[#] <> "_FF, " <>
                          "model, qedqcd);\n" <>

                          "gpLV += 2.*" <> CXXNameOfField[#] <> "_penguin[0] + "    <> CXXNameOfField[#] <> "_penguin[2];\n" <>
                          "gpRV += 2.*" <> CXXNameOfField[#] <> "_penguin[1] + "    <> CXXNameOfField[#] <> "_penguin[3];\n" <>
                          "gnLV += "    <> CXXNameOfField[#] <> "_penguin[0] + 2.*" <> CXXNameOfField[#] <> "_penguin[2];\n" <>
                          "gnRV += "    <> CXXNameOfField[#] <> "_penguin[1] + 2.*" <> CXXNameOfField[#] <> "_penguin[3];\n")&*),

                      (* create a list of massive, electrically neutral scalars *)
                      Select[GetParticles[],
                         (IsScalar[#] && !IsMassless[#] && !IsElectricallyCharged[#] && !ColorChargedQ[#])&
                      ]
                    ] <>

                "\n// ------ boxes ------\n\n" <>

                "gpLV += 0.;\n" <>
                "gpRV += 0.;\n" <>
                "gnLV += 0.;\n" <>
                "gnRV += 0.;\n" <>

                "\nconst auto nuclear_form_factors =\n" <>
                   IndentText[
                      "get_overlap_integrals(nucleus, qedqcd);\n"
                   ] <>

                "\nconst auto left =\n" <> IndentText[
                   "A2R*nuclear_form_factors.D\n" <>
                      "+ gpLV*nuclear_form_factors.Vp\n" <>
                      "+ gnLV*nuclear_form_factors.Vn\n" <>
                      "+ gpLS*nuclear_form_factors.Sp\n" <>
                      "+ gnLS*nuclear_form_factors.Sn"
                ] <> ";\n" <>
                "\nconst auto right =\n" <> IndentText[
                   "A2L*nuclear_form_factors.D\n" <>
                      "+ gpRV*nuclear_form_factors.Vp\n" <>
                      "+ gnRV*nuclear_form_factors.Vn\n" <>
                      "+ gpRS*nuclear_form_factors.Sp\n" <>
                      "+ gnRS*nuclear_form_factors.Sn"
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
