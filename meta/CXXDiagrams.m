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

(*** This module generates c++ code capable of representing fields, vertices and diagrams. ***)
BeginPackage["CXXDiagrams`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`", "CConversion`", "ColorMath`", "Utils`"}];

CreateFieldTraitsDefinitions::usage="";

(*** Public interfaces that model fields ***)
CreateFields::usage="Creates c++ code that makes all fields and their properties \
available as c++ types.";
CXXNameOfField::usage="Returns the name of the c++ type corresponding to a \
(possibly conjugate) field.";
CXXNameOfVertex::usage="Returns the name of the c++ type corresponding to a \
given vertex.";
LorentzConjugateOperation::usage="Returns the appropriate c++ typename \
to conjugate a given field as it would be used by SARAH`AntiField[].";
LorentzConjugate::usage="Conjugate a given field as it would be by SARAH`AntiField[].";
RemoveLorentzConjugation::usage="Unconjugate a field if it is conjugate. \
Leave it untouched otherwise.";
NumberOfFieldIndices::usage="Return the number of indices a field has as would be \
determined by inspecting the result of TreeMasses`FieldInfo[].";
CreateMassFunctions::usage="Creates c++ code that makes functions available that \
return tree-level masses of given fields.";
CreatePhysicalMassFunctions::usage="Creates c++ code that makes functions available that \
return physical masses of given fields.";
LorentzIndexOfField::usage="Returns the Lorentz index of a given indexed field.";
ColourIndexOfField::usage="Returns the colour index of a given indexed field.";

(*** Public interfaces that model vertices ***)
(** We need to encode the vertex structure for every unbroken gauge group:
- Identity component of Poincare group
- SU(3)) **)
(* A basis for the currently implemented Lorentz structures of vertices *)
ScalarVertex::usage="A Lorentz scalar (invariant) vertex";
ChiralVertex::usage="A vertex decomposing into left/right projectors";
MomentumVertex::usage="A vertex proportional to a single momentum term";
TripleVectorVertex::usage="A vertex with the Lorentz structure of a tree-level VVV vertex";
QuadrupleVectorVertex::usage="A vertex with the Lorentz structure of a tree-level VVVV vertex";
MomentumDifferenceVertex::usage="A vertex with a momentum difference (p^\\mu - q^\\nv)";
InverseMetricVertex::usage="A vertex proportional to g^{\\mu \\nu}";

(* A basis for the currently implemented colour structures of vertices *)
ZeroColouredVertex::usage="A vertex that is identically zero";
UncolouredVertex::usage="A SU(3) invariant vertex";
KroneckerDeltaColourVertex::usage="A vertex proportional to \\delta( ct1, ct2 )";
GellMannVertex::usage="A vertex proportional to \\Lambda^t_{a b}";
AdjointlyColouredVertex::usage="A vertex proportional to f^{a b c}";

CreateVertices::usage="Creates c++ code that makes functions available that \
numerically evaluates any of the given vertices.";
VertexTypes::usage="Returns a list of the different Lorentz structures that \
are currently implemented by CXXDiagrams.";
VertexTypeForFields::usage="Returns the Lorentz structure of a given vertex.";

(*** Public interfaces that model diagrams ***)
FeynmanDiagramsOfType::usage="Obtain all instantiations of Feynman \
diagrams of a given topology with given external fields.";
VerticesForDiagram::usage="Returns a list of all vertices present in a diagram.";
IndexDiagramFromGraph::usage="Fully index all fields in a given diagram such that \
contracted fields share the same indices.";
ColourFactorForIndexedDiagramFromGraph::usage="Calculate the colour factor of \
a given diagram with a given topology.";
ContractionsBetweenVerticesForDiagramFromGraph::usage="Returns a list of the positions \
of the contracted fields in a given diagram in the given vertices.";

(*** Public utility functions ***)
AtomHead::usage="Repeatedly evaluate Head[] on the argument until it satisfies \
AtomQ[] and return the result.";
CreateUnitCharge::usage="Creates the c++ code for a function that returns the \
numerical value of the electrical charge of the electron.";

ColorFactorForDiagram::usage = "Given topology and diagram returns the color factors for the diagram";
ExtractColourFactor::usage = "Drop colour generator from the colour factor";

Begin["`Private`"];

LeftChiralVertex::usage="A left projector part of a vertex";
RightChiralVertex::usage="A right projector part of a vertex";
TwoMetricVertex::usage="A g[l1,l2] * g[l3,l4] part of a vertex";

(** \brief Returns a list of the different Lorentz structures that
 * are currently implemented by CXXDiagrams.
 * \returns a list of the different Lorentz structures that
 * are currently implemented by CXXDiagrams.
 **)
VertexTypes[] := {
	ScalarVertex,
	ChiralVertex,
	MomentumVertex,
	TripleVectorVertex,
	QuadrupleVectorVertex,
	MomentumDifferenceVertex,
	InverseMetricVertex
};

CreateSelfConjugateFieldDefinition[field_, namespacePrefix_] := "";
IsLorentzSelfConjugate[field_] := field === LorentzConjugate[field];
CreateSelfConjugateFieldDefinition[field_?IsLorentzSelfConjugate, namespacePrefix_] :=
    Module[{fieldName},
           fieldName = TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix];
           "template<> struct " <> FlexibleSUSY`FSModelName <>  "_cxx_diagrams::fields::" <> LorentzConjugateOperation[field] <> "<" <> fieldName <> ">" <>
           " { using type = " <> fieldName <> "; };"
          ];

CreateFieldTypeTraitDefinition[field_?TreeMasses`IsScalar, namespacePrefix_] :=
    "template<>\nstruct is_scalar<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTypeTraitDefinition[field_?TreeMasses`IsFermion, namespacePrefix_] :=
    "template<>\nstruct is_fermion<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTypeTraitDefinition[field_?TreeMasses`IsVector, namespacePrefix_] :=
    "template<>\nstruct is_vector<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTypeTraitDefinition[field_?TreeMasses`IsGhost, namespacePrefix_] :=
    "template<>\nstruct is_ghost<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTraitsDefinition[field_, namespacePrefix_] :=
    Module[{fieldTypeTraitDefinition = ""},
           fieldTypeTraitDefinition = CreateFieldTypeTraitDefinition[field, namespacePrefix];
           fieldTypeTraitDefinition <> "\n"
          ];

CreateFieldTraitsDefinitions[fields_, namespacePrefix_:""] :=
    StringJoin[CreateFieldTraitsDefinition[#, namespacePrefix]& /@ fields];

(* Return a string corresponding to the c++ class name of the field.
 Note that "bar" and "conj" get turned into bar<...>::type and
 conj<...>::type respectively! *)
(** \brief Returns the name of the c++ type corresponding to a
 * (possibly conjugate) field.
 * Note that `bar` and `conj` get turned into `typename bar<...>::type`
 * and `typename conj<...>::type` respectively.
 * \param p the given (possibly conjugate) field
 * \param prefixNamespace Either `False` to apply no prefix namespace.
 * Otherwise the string `OptionValue[prefixNamespace]` is prepended to
 * the field name.
 * \returns the name of the c++ type corresponding to a
 * (possibly conjugate) field.
 **)
CXXNameOfField[p_, OptionsPattern[{prefixNamespace -> False}]] :=
  If[StringQ[OptionValue[prefixNamespace]],
     OptionValue[prefixNamespace] <> "::",
     ""] <> SymbolName[p];
CXXNameOfField[SARAH`bar[p_], OptionsPattern[{prefixNamespace -> False}]] :=
  "typename " <> If[OptionValue[prefixNamespace] === False, "", OptionValue[prefixNamespace]<>"::"]  <> "bar<" <> CXXNameOfField[p, prefixNamespace -> OptionValue[prefixNamespace]] <>
  ">::type";
CXXNameOfField[Susyno`LieGroups`conj[p_],
               OptionsPattern[{prefixNamespace -> False}]] :=
  "typename " <> If[OptionValue[prefixNamespace] === False, "", OptionValue[prefixNamespace]<>"::"] <> "conj<" <> CXXNameOfField[p, prefixNamespace -> OptionValue[prefixNamespace]] <>
  ">::type";

(** \brief Returns the name of the c++ type corresponding to a
 * given vertex. The c++ field names will automatically have `fields`
 * prepended to them as determined by `CXXNameOfField[]`.
 * \param fields the list of fields in the given vertex
 * \returns the name of the c++ type corresponding to a
 * given vertex.
 **)
CXXNameOfVertex[fields_List] := "Vertex<" <> StringJoin[Riffle[
		CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields,
	", "]] <> ">"

(** \brief Returns the appropriate c++ typename to conjugate a
 * given field as it would be used by ``SARAH`AntiField[]``.
 * \param field the given field
 * \returns the appropriate c++ typename to conjugate a
 * given field as it would be used by ``SARAH`AntiField[]``.
 **)
LorentzConjugateOperation[field_] :=
	If[TreeMasses`IsFermion[field] || TreeMasses`IsGhost[field],
		"bar", "conj"]

(** \brief Conjugate a given field as it would be by
 * ``SARAH`AntiField[]``.
 * \param field the given field
 * \returns the conjugate a given field as it would be by
 * ``SARAH`AntiField[]``.
 **)
LorentzConjugate[field_] :=
	If[TreeMasses`IsFermion[field] || TreeMasses`IsGhost[field],
		SARAH`bar[field], Susyno`LieGroups`conj[field]]

(** \brief Unconjugate a field if it is conjugate. Leave it untouched
 * otherwise.
 * \param p the given field
 * \returns the unconjugated field.
 **)
RemoveLorentzConjugation[p_] := p
RemoveLorentzConjugation[SARAH`bar[p_]] := p
RemoveLorentzConjugation[Susyno`LieGroups`conj[p_]] := p

(** \brief Repeatedly evaluate `Head[]` on the argument until it
 * satisfies `AtomQ[]` and return the result.
 * \param x the given argument
 * \returns the result of repeatedly evaluating `Head[]` on the
 * argument until it satisfies `AtomQ[]` and return the result.
 **)
AtomHead[x_] := If[AtomQ[x], x, AtomHead[Head[x]]]

ParticleTypeAsString::argx = "Unknown type of particle `1`. Supported types are scalar, fermion, vector and ghost.";
ParticleTypeAsString[part_] := Module[
   {},
   If[TreeMasses`IsScalar[part],    Return["scalar"]];
   If[TreeMasses`IsVector[part],    Return["vector"]];
   If[TreeMasses`IsFermion[part],   Return["fermion"]];
   If[TreeMasses`IsGhost[part],   Return["ghost"]];

   Message[ParticleTypeAsString::argx, part]; Abort[];
];

(** \brief Creates c++ code that makes all fields and their properties
 * available as c++ types. Also creates two using declarations for
 * the following fields:
 * - Photon
 * - Electron
 * Furthermore create the necessary boilerplate code to conjugate any
 * given c++ field as well as convenience `boost::mpl::vector<>`
 * containers for the following field types:
 * - scalars
 * - fermions
 * - vectors
 * - ghosts
 * \returns the corresponding c++ code.
 **)
ParticleColorRepAsString::argx = "Unknown color representation `1` for particle `2`. Supported representations are singlet, (anti-)triplet and octet.";
ParticleColorRepAsString[part_] :=
   Module[{rep = TreeMasses`GetColorRepresentation[part]},
      Switch[rep,
         S, "singlet",
         T, "triplet",
         -T, "anti_triplet",
         O, "octet",
         _, Message[ParticleColorRepAsString::argx, rep, part]; Abort[];
      ]
   ];

CreateFields[] :=
  Module[{fields, scalars, fermions, vectors, ghosts},
       fields = TreeMasses`GetParticles[];
       scalars = Select[fields, TreeMasses`IsScalar];
       fermions = Select[fields, TreeMasses`IsFermion];
       vectors = Select[fields, TreeMasses`IsVector];
       ghosts = Select[fields, TreeMasses`IsGhost];

       StringJoin @ Riffle[
         ("struct " <> CXXNameOfField[#] <> " {\n" <>
            TextFormatting`IndentText[
              "static constexpr auto particle_type = ParticleType::" <> ParticleTypeAsString[#] <> ";\n" <>
              "static constexpr auto color_rep = ParticleColorRep::" <> ParticleColorRepAsString[#] <> ";\n" <>
              "static constexpr auto massless = " <> CConversion`CreateCBoolValue @ TreeMasses`IsMassless[#] <> ";\n" <>
              "using index_bounds = boost::mpl::pair<\n" <>
              "  boost::mpl::vector_c<int" <>
                   StringJoin[", " <> ToString[#] & /@
                     (IndexBoundsForField[#][[1]] - 1)] <> ">,\n" <>
              "  boost::mpl::vector_c<int" <>
                   StringJoin[", " <> ToString[#] & /@
                     IndexBoundsForField[#][[2]]] <> ">\n" <>
              ">;\n" <>
              "static constexpr int numberOfGenerations = " <>
                   ToString @ TreeMasses`GetDimension[#] <> ";\n" <>
                     "using sm_flags = boost::mpl::vector_c<bool, " <>
                        If[TreeMasses`GetDimension[#] === 1,
                           CConversion`CreateCBoolValue @ TreeMasses`IsSMParticle[#],
                           StringJoin @ Riffle[CConversion`CreateCBoolValue /@
                             TreeMasses`IsSMParticleElementwise[#],
                                               ", "]] <>
                        ">;\n" <>
                "static constexpr int numberOfFieldIndices = " <>
                   ToString @ NumberOfFieldIndices[#] <> ";\n" <>
                "static constexpr double electric_charge = " <>
                   CConversion`RValueToCFormString[TreeMasses`GetElectricCharge[#]] <> ";\n" <>
                "using lorentz_conjugate = " <>
                   CXXNameOfField[LorentzConjugate[#]] <> ";\n"] <>
              "};" &) /@ fields, "\n\n"] <> "\n\n" <>

       "// Named fields\n" <>
       "using Photon = " <> CXXNameOfField[SARAH`Photon] <> ";\n" <>
       "using Electron = " <> CXXNameOfField[AtomHead @ TreeMasses`GetSMElectronLepton[]] <> ";\n\n" <>

       "// Fields that are their own Lorentz conjugates.\n" <>
       StringJoin @ Riffle[
         ("template<> struct " <> LorentzConjugateOperation[#] <> "<" <> CXXNameOfField[#] <> ">" <>
            " { using type = " <> CXXNameOfField[#] <> "; };"
            &) /@ Select[fields, (# == LorentzConjugate[#] &)],
          "\n"] <> "\n\n" <>

       "using scalars = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ scalars, ", "]] <> ">;\n" <>
       "using fermions = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ fermions, ", "]] <> ">;\n" <>
       "using vectors = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ vectors, ", "]] <> ">;\n" <>
       "using ghosts = boost::mpl::vector<" <>
         StringJoin[Riffle[CXXNameOfField /@ ghosts, ", "]] <> ">;"
  ]

(** \brief Get the lorentz index of a given indexed field
 * \param field The indexed field
 * \returns the Lorentz index of a given indexed field.
 **)
LorentzIndexOfField[field_]:=
	Module[{lorentzIndices = Select[Vertices`FieldIndexList[field], Vertices`SarahLorentzIndexQ]},
		Utils`AssertWithMessage[Length[lorentzIndices] === 1,
			"CXXDiagrams`LorentzIndexOfField[]: Argument " <>
			ToString[field] <> " does not have exactly one Lorentz index."];
    lorentzIndices[[1]]
  ]

(** \brief Get the lorentz index of a given indexed field
 * \param field The indexed field
 * \returns the Lorentz index of a given indexed field.
 **)
ColourIndexOfField[field_ /; TreeMasses`ColorChargedQ[field]]:=
	Module[{colorIndices = Select[Vertices`FieldIndexList[field], Vertices`SarahColorIndexQ]},
		Utils`AssertWithMessage[Length[colorIndices] === 1,
			"CXXDiagrams`ColourIndexOfField[]: Argument " <>
			ToString[field] <> " does not have exactly one color index."];
    colorIndices[[1]]
  ]

(** \brief Obtain all instantiations of Feynman diagrams of a given
 * topology with given external fields.
 * \param adjacencyMatrix the symmetric adjacency matrix determining the
 * topology of the diagrams
 * \param externalFields a list of rules of the form
 * `{index1 -> externalField1, ... }`
 * where outgoing fields should be specified as conjugate fields and
 * incoming fields should be specified as unconjugate fields.
 * \returns a list of Feynman diagrams of the given topology with
 * given external fields.
 * All returned diagrams are of the form:
 * `{v1, v2, ..., vn}`
 * where `n` is the dimension of the given adjacency matrix and
 * `vk` for \f$k \in \{ 1, ..., n \} \f$ corresponds to the vertex
 * with index `k` as specified by the adjacency matrix.
 * A vertex in the returned diagram is either external or internal,
 * where external vertices are by definition taken to be those
 * specified by `externalFields`. Thus an external vertex `vk` is simply
 * represented as the external field `externalFieldk` as specified by
 * `externalFields`.
 * The internal vertices are given as lists of fields corresponding to
 * terms in the interaction Lagrangian. They are obtained in
 * correspondence with the standard perturbative QFT approach when
 * calculating
 * \f[ \langle \left( \text{outgoing} \right)^* \vert
 * 	\text{incoming} \rangle \f]
 * where the asterisk denotes Lorentz conjugation. Of course only
 * diagrams of the given topology are taken into account.
 * The order of fields in the internal vertices is given by the
 * adjacency matrix as follows:
 * Let `k` be the index of an internal vertex and `l` be a
 * non-zero entry in `adjacencyMatrix[[k, l]]`. Then, if `k =!= l`
 * `vk` will contain exactly `adjacencyMatrix[[k, l]]` entries of fields
 * that get contracted with `adjacencyMatrix[[k, l]]` fields in the
 * vertex `vl`. The order of those `adjacencyMatrix[[k, l]]` fields is
 * unspecified. If `k === l` `vk` will contain exactly
 * `2 * adjacencyMatrix[[k, l]]` entries of fields where
 * the first one gets contracted with the second one, the third one with
 * the fourth one and so on. Otherwise the order of these
 * `2 * adjacencyMatrix[[k, l]]` fields is unspecified.
 * \note Here, 'contracted' means that contracted fields will be
 * Lorentz conjugates to each other - unless one of the fields is
 * external, then they will be the same.
 * Finally, the vertex `vk` is given by concatenating the thusly
 * constructed lists of fields for all `l` where
 * `adjacencyMatrix[[k, l]]` is non-zero in the order given by
 * the index `l`.
 * Furthermore, two diagrams are considered equal if they only differ by
 * a reordering of fields in internal vertices and thus any appearing
 * duplicate diagrams are removed.
 **)
FeynmanDiagramsOfType[adjacencyMatrix_List,externalFields_List] :=
	Module[{externalVertices = externalFields[[All,1]],
			internalVertices,externalRules, internalFieldCouplings,
			unspecifiedEdgesLess,unspecifiedEdgesEqual,
			insertFieldRulesLess,insertFieldRulesGreater,insertFieldRulesEqual,
			fieldsToInsert,
			unresolvedFieldCouplings,resolvedFields,resolvedFieldCouplings,
			diagrams},
	Utils`AssertWithMessage[adjacencyMatrix === Transpose[adjacencyMatrix],
		"CXXDiagrams`FeynmanDiagramsOfType[]: The given adjacency matrix is \
not symmetric"];
	LoadVerticesIfNecessary[];

	internalVertices = Complement[Table[k,{k,Length[adjacencyMatrix]}],externalVertices];
	externalRules = Flatten @ ({{_,#,_} :> (# /. externalFields)} & /@
		externalVertices);

	internalFieldCouplings = (Flatten[(Flatten @ Position[adjacencyMatrix[[#]],Except[0],{1},Heads -> False]
		/. {i_Integer :> Table[{#,i,k},{k,adjacencyMatrix[[#,i]]}]}),1] &
		/@ internalVertices) /. externalRules;

	unspecifiedEdgesLess = Cases[internalFieldCouplings,{i_,j_,_} /; i < j,{2}];
	unspecifiedEdgesEqual = Cases[internalFieldCouplings,{i_,i_,_},{2}];

	insertFieldRulesLess = MapIndexed[#1 -> SARAH`FieldToInsert[#2[[1]]] &,unspecifiedEdgesLess];
	insertFieldRulesGreater = (insertFieldRulesLess /. {Rule[{i_,j_,k_},field_] :> Rule[{j,i,k},SARAH`AntiField[field]]});
	insertFieldRulesEqual = MapIndexed[#1 -> Sequence @@ {SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]],
		SARAH`AntiField[SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]]]} &,
			unspecifiedEdgesEqual];
	fieldsToInsert = Table[SARAH`FieldToInsert[k],
		{k,Length[insertFieldRulesLess] + Length[insertFieldRulesEqual]}];

	unresolvedFieldCouplings = internalFieldCouplings
		/. insertFieldRulesLess /. insertFieldRulesGreater /. insertFieldRulesEqual;

	resolvedFields = If[fieldsToInsert === {}, {{}},
		SARAH`InsFields[{C @@@ unresolvedFieldCouplings,
			fieldsToInsert}][[All,2]]];

	If[resolvedFields === {}, Return[{}]];

	resolvedFieldCouplings = unresolvedFieldCouplings /.
		((Rule @@@ Transpose[{fieldsToInsert,#}]) & /@ resolvedFields);

	diagrams = Table[k,{k,Length[adjacencyMatrix]}] /. externalFields /.
		((Rule @@@ Transpose[{internalVertices,#}]) & /@ resolvedFieldCouplings);

	(* Prevent overcounting of diagrams by removing diagrams that only
	 * differ by permutations of internal fields within the internal
	 * vertices. This is automatically performed by SARAH if
	 * SA`CheckSameVertices === True, but it is better to not depend on
	 * some internal SARAH state. *)
	DeleteDuplicates[diagrams,
		(And @@ ((Sort[#[[1]]] === Sort[#[[2]]] &) /@
			Cases[Transpose[{#1, #2}],{{___},{___}}] (* Only check internal vertices *)
		) &)]
  ]

(** \brief Returns a list of all vertices present in a diagram.
 * \param diagram The given Feynman diagram
 * \returns a list of all vertices present in a diagram.
 **)
VerticesForDiagram[diagram_] := Select[diagram,Length[#] > 1 &]

(** \brief Returns a list of the positions of the contracted fields in
 * a given diagram in the given vertices.
 * \param v1 the index of the first vertex in question
 * \param v2 the index of the second vertex in question
 * \param diagram the given Feynman diagram
 * \param graph the topology given as an adjacency matrix
 * \returns a list of the positions of the contracted fields in
 * a given diagram in the given vertices.
 **)
ContractionsBetweenVerticesForDiagramFromGraph[v1_Integer, v2_Integer,
		diagram_List, graph_List] :=
	Module[{fields1 = diagram[[v1]], fields2 = diagram[[v2]],
			preceedingNumberOfFields1 = Total[graph[[v1, ;;v2]]] - graph[[v1,v2]],
			preceedingNumberOfFields2 = Total[graph[[v2, ;;v1]]] - graph[[v2,v1]],
			contractedFieldIndices1, contractedFieldIndices2,
			stepSize = If[v1 === v2, 2, 1]},
		If[v1 === v2,
			preceedingNumberOfFields2 = preceedingNumberOfFields2 + 1];
		If[v1 < v2,
			preceedingNumberOfFields1 =
				preceedingNumberOfFields1 + graph[[v1,v1]]];
		If[v2 < v1,
			preceedingNumberOfFields2 =
				preceedingNumberOfFields2 + graph[[v2,v2]]];

		contractedFieldIndices1 = Table[k, {k, preceedingNumberOfFields1 + 1,
			preceedingNumberOfFields1 + stepSize * graph[[v1,v2]], stepSize}];
		contractedFieldIndices2 = Table[k, {k, preceedingNumberOfFields2 + 1,
			preceedingNumberOfFields2 + stepSize * graph[[v2,v1]], stepSize}];

		Transpose[{contractedFieldIndices1, contractedFieldIndices2}]
	]

(** \brief Convert color structures to the ColorMath convention **)
ConvertColourStructureToColorMathConvention[fields_List,
	ZeroColouredVertex] := 0;

ConvertColourStructureToColorMathConvention[fields_List,
	UncolouredVertex] := 1;

ConvertColourStructureToColorMathConvention[fields_List,
	GellMannVertex[cIndex1_, cIndex2_, cIndex3_]] :=
	(* FIXME: Factor of two? *)
	2 ColorMath`CMt[{cIndex1}, cIndex2, cIndex3];

ConvertColourStructureToColorMathConvention[fields_List,
	GellMannVertex[cIndex1_, cIndex2_, cIndex3_] GellMannVertex[cIndex4_, cIndex5_, cIndex6_] +
			GellMannVertex[cIndex7_, cIndex8_, cIndex9_] GellMannVertex[cIndex10_, cIndex11_, cIndex12_]] :=
       (4 ColorMath`CMt[{cIndex1}, cIndex2, cIndex3] ColorMath`CMt[{cIndex4}, cIndex5, cIndex6] +
				4 ColorMath`CMt[{cIndex7}, cIndex8, cIndex9] ColorMath`CMt[{cIndex10}, cIndex11, cIndex12]);

ConvertColourStructureToColorMathConvention[fields_List,
	AdjointlyColouredVertex[cIndex1_, cIndex2_, cIndex3_]] :=
	ColorMath`CMf[cIndex1, cIndex2, cIndex3];

ConvertColourStructureToColorMathConvention[indexedFields_List,
	KroneckerDeltaColourVertex[cIndex1_, cIndex2_]] :=
	Module[{colouredField1, colouredField2, colourRep1, colourRep2},
		(* If the result has a color Delta, we need to find out if it's adj. or fundamental
			because they are represented by different symbols in ColorMath.
			Also, in ColorMath the order of indices matters. As stated in  ColorMath tutorial notebook:

			"Note that lower and upper indices are different. One position is taken to define incoming
			quarks/outgoing anti-quarks, and the other position is used for outgoing quarks and incoming
			anti-quarks. What is taken to be what is a question of definition."
		*)

    colouredField1 = Select[indexedFields, !FreeQ[#, cIndex1] &][[1]];
    colouredField2 = Select[indexedFields, !FreeQ[#, cIndex2] &][[1]];

    colourRep1 = SARAH`getColorRep[colouredField1];
    colourRep2 = SARAH`getColorRep[colouredField2];

    Utils`AssertWithMessage[colourRep1 == colourRep2,
			"CXXDiagrams`ConvertColourStructureToColorMathConvention[]: " <>
			"Two colour indices in Kronecker delta that come from fields " <>
			"of incompatible representations."];

		(* FIXME: Are these orderings correct? *)
		Switch[{colourRep1},
			{SARAH`T}, If[RemoveLorentzConjugation[colouredField1] === colouredField1,
				ColorMath`CMdelta[cIndex2, cIndex1],
				ColorMath`CMdelta[cIndex1, cIndex2]
			],
			{O}, ColorMath`CMDelta[cIndex2, cIndex1],
			_,
			Print["CXXDiagrams`ConvertColourStructureToColorMathConvention[]: " <>
				"Two colour indices in Kronecker delta that come from fields " <>
				"or representations we cannot handle: " <>
				ToString[{colourRep1, colourRep2}]];
			Quit[1];
		]
	];

ConvertColourStructureToColorMathConvention[fields_List,
   AdjointlyColouredVertex[cIndex1_, cIndex2_, cIndex3_]  GellMannVertex[cIndex3_, cIndex4_, cIndex5_]] :=
   ColorMath`CMf[cIndex1, cIndex2, cIndex3] * 2 ColorMath`CMt[{cIndex3}, cIndex4, cIndex5];

(* FIXME: Are these correct? *)
ColorMathToSARAHConvention[expr_] :=
	expr /. {
		Subscript[Superscript[CM\[Delta], cIndex1_], cIndex2_] :>
			SARAH`Delta[cIndex1, cIndex2],
		Superscript[CM\[Delta], {indices__}] :> SARAH`Delta[indices],
		Superscript[CM\[CapitalDelta], {indices__}] :> SARAH`Delta[indices],
		Subscript[Superscript[Superscript[ColorMath`CMt, {cIndex1_}], cIndex2_], cIndex3_] :>
			1 / 2 * SARAH`Lam[cIndex1, cIndex2, cIndex3],
		ColorMath`Nc -> 3,
		ColorMath`TR -> 1/2
	};

(** \brief Fully index all fields in a given diagram such that
 * contracted fields share the same indices.
 * \param diagram The diagram to index
 * \param graph The graph giving the topology of \a diagram
 * \returns An indexed version of \a diagram in which
 * contracted fields share the same indices.
 **)
IndexDiagramFromGraph[diagram_, graph_] :=
	Module[{diagramWithUIDs, fields, indexedFields, indexedDiagram,
		applyUserIndices, vIndex1, vIndex2, contractIndices},
		fields = Vertices`StripFieldIndices /@ Flatten[diagram];

		indexedFields = IndexField /@ fields;

		applyUserIndices = Join[
			Sequence @@ (Rule @@@ Transpose[{
				Take[#[[1]], Length[#[[2]]]],
				#[[2]]
			}] & /@ Transpose[{
				Vertices`FieldIndexList /@ indexedFields,
				Vertices`FieldIndexList /@ Flatten[diagram]
			}])
		];

		diagramWithUIDs = If[Head[#] === List,
			(#[Unique[]] &) /@ #, #[Unique[]]] & /@ diagram;

		indexedDiagram = diagramWithUIDs /. Rule @@@ Transpose[
			{Flatten[diagramWithUIDs], indexedFields}];

		contractIndices = Module[{
				vIndex1 = #[[1]], vIndex2 = #[[2]], connectedFields},
			connectedFields = ContractionsBetweenVerticesForDiagramFromGraph[
				vIndex1, vIndex2, diagram, graph];

			Module[{field1, field2},
				field1 = If[Head[diagram[[vIndex1]]] === List,
					indexedDiagram[[vIndex1, #[[1]]]],
					indexedDiagram[[vIndex1]]];
				field2 = If[Head[diagram[[vIndex2]]] === List,
					indexedDiagram[[vIndex2, #[[2]]]],
					indexedDiagram[[vIndex2]]];

				Thread[Rule[Vertices`FieldIndexList[field1],
					Vertices`FieldIndexList[field2]]]
			] & /@ connectedFields
		] & /@ Flatten[
			Table[{v1, v2}, {v1, Length[diagram]}, {v2, v1, Length[diagram]}],
		1];

		indexedDiagram /. Flatten[contractIndices] /. applyUserIndices /.
			(Reverse /@ Flatten[contractIndices]) /. applyUserIndices
	]

(** \brief Calculate the colour factor of a given diagram with a given
 * topology.
 * \param diagram The indexed diagram specifiying all occurring
 * fields. All indices need to be properly connected!
 * \param graph The topology of the diagram as adjacency matrix.
 * \returns The color factor for \a diagram.
 **)
ColourFactorForIndexedDiagramFromGraph[indexedDiagram_, graph_] :=
	Module[{indexedVertexFields = Cases[indexedDiagram, {__}],
			vertexFields, sortedVertices, sortedColourStructures,
			verticesOrderings, inverseVOrderings, sortedVerticesOrderings,
			fOrderingWRTSortedVs, defaultIndexedVertexFields,
			defaultIndexedColorMathExpressions, colorMathExpressions},
		vertexFields = Vertices`StripFieldIndices /@ indexedVertexFields;
		sortedVertices = SortedVertex /@ vertexFields;

		sortedColourStructures = Transpose[{
			sortedVertices[[All,1]],
			(GaugeStructureOfVertex /@ sortedVertices)[[All,2]]
		}];

    (* Mathematica 7 does not know about permutations... :'-( *)
    verticesOrderings = Ordering /@ vertexFields;
    sortedVerticesOrderings = Ordering /@
			Map[Vertices`StripFieldIndices, sortedVertices[[All,1]], {2}];

    inverseVOrderings = Ordering /@ verticesOrderings;

    fOrderingWRTSortedVs = Part @@@ Transpose[
			{sortedVerticesOrderings, inverseVOrderings}];

    defaultIndexedVertexFields = Part @@@ Transpose[
			{sortedVertices[[All,1]], fOrderingWRTSortedVs}];
		defaultIndexedColorMathExpressions =
			ConvertColourStructureToColorMathConvention @@@ sortedColourStructures;

		colorMathExpressions = (#[[1]] /. (Flatten[Thread /@ Rule @@@
			Map[Vertices`FieldIndexList,
				Transpose[{#[[2]], #[[3]]}], {2}]]) &) /@
			Transpose[{defaultIndexedColorMathExpressions,
				defaultIndexedVertexFields, indexedVertexFields}];

		If[(Times @@ colorMathExpressions) === 1, 1,
			ColorMathToSARAHConvention[
				ColorMath`CSimplify[Times @@ colorMathExpressions]]
		]
	]
(** \brief Returns the index prefix string for an index of a given type **)
IndexPrefixForType[SARAH`generation] = "gt";
IndexPrefixForType[SARAH`lorentz] = "lt";
IndexPrefixForType[SARAH`color] = "ct";
IndexPrefixForType[indexType_] :=
	(Print["Unknown index type: " <> ToString[indexType] <> " encountered."]; Quit[1])

(** \brief Fully index a given field
 * \param field the given unindexed field
 * \returns a fully indexed version of the field with indices determined
 * by `Unique[]`.
 **)
IndexField[field_] :=
	Module[{indexSpecification, indexTypes, indexNames,
		conjugation, uniqueNumberString},
		indexSpecification =
			TreeMasses`GetParticleIndices[SARAH`getParticleName[field]];
		indexSpecification =
			Cases[indexSpecification, Except[{SARAH`generation, 1}]];

		indexTypes = First /@ indexSpecification;

		uniqueNumberString = StringDrop[SymbolName[Unique[]], 1];
		indexNames = "SARAH`" <> ToString[#]<> uniqueNumberString & /@
			(IndexPrefixForType /@ indexTypes);
		field[Symbol /@ indexNames] /. {
			field[{}] -> field,
			conjugation_[field][indices_List] :> conjugation[field[indices]]
		}
	]

(** \brief Sort a set of fields into SARAH canonical order and return
 * the corresponding vertex as determined by ``SARAH`Vertex[]`` with
 * appropriate FlexibleSUSY conventions applied. These are:
 * - GUT normalization
 * - an extra factor of `-I`
 * \param fields The given set of fields
 * \param ApplyGUTNormalization True if GUT normalization should be
 * applied and False otherwise.
 * \returns the sorted vertex as determined by ``SARAH`Vertex[]`` with
 * appropriate FlexibleSUSY conventions applied.
 **)
SortedVertex[fields_List, OptionsPattern[{ApplyGUTNormalization -> False}]] :=
	Module[{sortedFields, indexedSortedFields, vertexList, vertex, types,
			similarVertexList, similarVertex, fieldReplacementRules,
			indexReplacementRules},
    LoadVerticesIfNecessary[];
		sortedFields = Vertices`SortFieldsInCp[fields];

    vertexList = Switch[Length[fields],
			3, SARAH`VertexList3,
			4, SARAH`VertexList4,
			_, Print["CXXDiagrams`SortedVertex[]: Only three- and four-point vertices are supported"];
				Quit[1]];

    vertex = Select[vertexList, Vertices`StripFieldIndices[#[[1]]] === sortedFields &, 1];
    vertex = If[vertex =!= {}, vertex[[1]],
			types = SARAH`getType[#, False, FlexibleSUSY`FSEigenstates] & /@ sortedFields;
			similarVertexList = SA`VertexList[Symbol[StringJoin[SymbolName /@ types]]];

			Utils`AssertWithMessage[similarVertexList =!= {},
				"CXXDiagrams`SortedVertex[]: Cannot determine Lorentz structure of " <> ToString[fields] <> " vertex."];
			similarVertex = similarVertexList[[1]];

			indexedSortedFields = IndexField /@ sortedFields;

			fieldReplacementRules = Rule @@@ Transpose[{similarVertex[[1]], indexedSortedFields}];
			indexReplacementRules = Rule @@@ Transpose[{
				LorentzIndexOfField /@ Select[similarVertex[[1]], TreeMasses`IsVector],
				LorentzIndexOfField /@ Select[indexedSortedFields, TreeMasses`IsVector]
			}];

			{indexedSortedFields, Sequence @@ Transpose[{
				Table[0, {Length[similarVertex] - 1}],
					Rest[similarVertex][[All,2]] /.fieldReplacementRules /.
						indexReplacementRules}]}
		];

		(* Be consistent with the conventions in Vertices.m *)
		vertex = {vertex[[1]], Sequence @@ Transpose[
			{-I * Rest[vertex][[All,1]], Rest[vertex][[All,2]]}]};

		If[OptionValue[ApplyGUTNormalization],
			vertex /. Parameters`ApplyGUTNormalization[],
			vertex]
	]

(** \brief Turn SARAH-style Lorentz structures into CXXDiagrams ones. **)
LabelLorentzPart[{scalar_, 1}] := {scalar, ScalarVertex}

LabelLorentzPart[{scalar_,
		SARAH`PL | SARAH`LorentzProduct[
			SARAH`gamma[lIndex_] /; Vertices`SarahLorentzIndexQ[lIndex], SARAH`PL]}] :=
	{scalar, LeftChiralVertex}

LabelLorentzPart[{scalar_,
		SARAH`PR | SARAH`LorentzProduct[
			SARAH`gamma[lIndex_] /; Vertices`SarahLorentzIndexQ[lIndex], SARAH`PR]}] :=
	{scalar, RightChiralVertex}

LabelLorentzPart[{scalar_, SARAH`Mom[in_, lIndex_]}] :=
	{scalar, MomentumVertex[in]} /; Vertices`SarahLorentzIndexQ[lIndex]

LabelLorentzPart[{scalar_,
	SARAH`g[lIndex1_, lIndex2_] *
		(SARAH`Mom[field1_, lIndex3_] - SARAH`Mom[field2_, lIndex3_])
	+ SARAH`g[lIndex2_, lIndex3_] *
		(SARAH`Mom[field2_, lIndex1_] - SARAH`Mom[field3_, lIndex1_])
	+ SARAH`g[lIndex1_, lIndex3_] *
		(SARAH`Mom[field3_, lIndex2_] - SARAH`Mom[field1_, lIndex2_])}] :=
	{scalar, TripleVectorVertex[field1, field2, field3]} /;
	SARAH`g[lIndex1, lIndex2] ===
		SARAH`g[LorentzIndexOfField[field1], LorentzIndexOfField[field2]] &&
	SARAH`g[lIndex2, lIndex3] ===
		SARAH`g[LorentzIndexOfField[field2], LorentzIndexOfField[field3]] &&
	SARAH`g[lIndex1, lIndex3] ===
		SARAH`g[LorentzIndexOfField[field1], LorentzIndexOfField[field3]]

LabelLorentzPart[{scalar_,
	- SARAH`g[lIndex1_, lIndex2_] *
		(SARAH`Mom[field1_, lIndex3_] - SARAH`Mom[field2_, lIndex3_])
	- SARAH`g[lIndex2_, lIndex3_] *
		(SARAH`Mom[field2_, lIndex1_] - SARAH`Mom[field3_, lIndex1_])
	- SARAH`g[lIndex1_, lIndex3_] *
		(SARAH`Mom[field3_, lIndex2_] - SARAH`Mom[field1_, lIndex2_])}] :=
	{scalar, TripleVectorVertex[field3, field2, field1]} /;
	SARAH`g[lIndex1, lIndex2] ===
		SARAH`g[LorentzIndexOfField[field1], LorentzIndexOfField[field2]] &&
	SARAH`g[lIndex2, lIndex3] ===
		SARAH`g[LorentzIndexOfField[field2], LorentzIndexOfField[field3]] &&
	SARAH`g[lIndex1, lIndex3] ===
		SARAH`g[LorentzIndexOfField[field1], LorentzIndexOfField[field3]]

LabelLorentzPart[{scalar_,
		SARAH`g[lIndex1_, lIndex2_] * SARAH`g[lIndex3_, lIndex4_]}] :=
	{scalar, TwoMetricVertex[lIndex1, lIndex2, lIndex3, lIndex4]} /;
	Vertices`SarahLorentzIndexQ[lIndex1] && Vertices`SarahLorentzIndexQ[lIndex2] &&
	Vertices`SarahLorentzIndexQ[lIndex3] && Vertices`SarahLorentzIndexQ[lIndex4]

LabelLorentzPart[{scalar_, SARAH`Mom[in_, lIndex1_] - SARAH`Mom[out_, lIndex2_]}] :=
	{scalar, MomentumDifferenceVertex[in, out]} /;
	Vertices`SarahLorentzIndexQ[lIndex1] && Vertices`SarahLorentzIndexQ[lIndex2]

LabelLorentzPart[{scalar_, SARAH`g[lIndex1_, lIndex2_]}] :=
	{scalar, InverseMetricVertex} /;
	Vertices`SarahLorentzIndexQ[lIndex1] && Vertices`SarahLorentzIndexQ[lIndex2]

LabelLorentzPart[vertexPart_] :=
	(Print["Unknown Lorentz structure in vertex ", vertexPart]; Quit[1])

(** \brief Turn SARAH-style colour structures into CXXDiagrams ones. **)
GaugeStructureOfVertexLorentzPart[{0, lorentzStructure_}] :=
	{0, ZeroColouredVertex, lorentzStructure};

GaugeStructureOfVertexLorentzPart[{scalar_, lorentzStructure_}] :=
	{scalar, UncolouredVertex, lorentzStructure} /;
	FreeQ[scalar, atom_ /; Vertices`SarahColorIndexQ[atom], -1];

GaugeStructureOfVertexLorentzPart[
	{scalar_ * SARAH`Delta[cIndex1_, cIndex2_], lorentzStructure_}] :=
	{scalar, KroneckerDeltaColourVertex[cIndex1, cIndex2], lorentzStructure} /;
	FreeQ[scalar, atom_ /; Vertices`SarahColorIndexQ[atom], -1];

GaugeStructureOfVertexLorentzPart[
	{scalar_ * SARAH`Lam[cIndex1_, cIndex2_, cIndex3_], lorentzStructure_}] :=
	{scalar, GellMannVertex[cIndex1, cIndex2, cIndex3], lorentzStructure} /;
	FreeQ[scalar, atom_ /; Vertices`SarahColorIndexQ[atom], -1];

GaugeStructureOfVertexLorentzPart[
	{scalar_ * SARAH`fSU3[cIndex1_, cIndex2_, cIndex3_], lorentzStructure_}] :=
	{scalar, AdjointlyColouredVertex[cIndex1, cIndex2, cIndex3], lorentzStructure} /;
	FreeQ[scalar, atom_ /; Vertices`SarahColorIndexQ[atom], -1];

GaugeStructureOfVertexLorentzPart[
	{scalar_ * sum[j_, 1, 8, SARAH`fSU3[cIndex1_, cIndex3_, j_] SARAH`Lam[j_, cIndex4_, cIndex2_]], lorentzStructure_}
   ] /; Vertices`SarahDummyIndexQ[j] :=
	{scalar, AdjointlyColouredVertex[cIndex1, cIndex3, cIndex5]  GellMannVertex[cIndex5, cIndex4, cIndex2], lorentzStructure} /;
	FreeQ[scalar, atom_ /; Vertices`SarahColorIndexQ[atom], -1];

GaugeStructureOfVertexLorentzPart[
	{fullExpr_, lorentzStructure_}] /; MatchQ[
	Expand @ fullExpr,
	scalar_ sum[j1_, 1, 3, SARAH`Lam[c1__] SARAH`Lam[c2__]] +
     scalar_ sum[j2_, 1, 3, SARAH`Lam[c3__] SARAH`Lam[c4__]] /;
         MemberQ[{c1}, j1] && MemberQ[{c2}, j2] && MemberQ[{c3}, j2] && MemberQ[{c4}, j2]
] := Expand @ fullExpr /. scalar_ sum[j1_, 1, 3, SARAH`Lam[c1__] SARAH`Lam[c2__]] +
		scalar_ sum[j2_, 1, 3, SARAH`Lam[c3__] SARAH`Lam[c4__]] :>
{scalar, GellMannVertex[c1] GellMannVertex[c2] + GellMannVertex[c3] GellMannVertex[c4], lorentzStructure} /;
		FreeQ[scalar, atom_ /; Vertices`SarahColorIndexQ[atom], -1];

(* @todo:
      This case catches vertices with sum of products of 2 Kronecker deltas.
      Currently, if scalar1 and scalar2 coefficients would be equal in principle
      we could handle this vertex. The case of scalar1 != scalar2 we cannot.
      For the moment therefore we only print a warning message, similar to
      the generic catch all case *)
GaugeStructureOfVertexLorentzPart[
	{fullExpr_, lorentzStructure_}] /; MatchQ[
      Collect[ExpandAll@fullExpr, {SARAH`Delta[ct1, ct4] SARAH`Delta[ct2, ct3],
  SARAH`Delta[ct1, ct3] SARAH`Delta[ct2, ct4]}],
	scalar1_ SARAH`Delta[c1_, c4_] SARAH`Delta[c2_, c3_] +
     scalar2_ SARAH`Delta[c1_, c3_] SARAH`Delta[c2_, c4_] ] :=
Collect[ExpandAll@fullExpr, {SARAH`Delta[ct1, ct4] SARAH`Delta[ct2, ct3],
  SARAH`Delta[ct1, ct3] SARAH`Delta[ct2, ct4]}] /.
	scalar1_ SARAH`Delta[c1_, c4_] SARAH`Delta[c2_, c3_] +
     scalar2_ SARAH`Delta[c1_, c3_] SARAH`Delta[c2_, c4_]  :>
Parameters`DebugPrint["Vertices with sum of 2 deltas are currently not supported"];

GaugeStructureOfVertexLorentzPart[vertexPart_] :=
   (Print["Unknown colour structure in vertex ", vertexPart]; Quit[1]);

(** \brief Given a list of gauge (colour and Lorentz) structures
 * combine all left and right projector parts to chiral parts.
 * Pad with zeros if necessary.
 * \param gaugeParts a list of gauge structures of the format
 * `{{scalar1, colourStructure1, LorentzStructure1}, ...}`
 * \returns the list of gauge structures where all left and right
 * projector parts are combined to chiral parts.
 * \note There may be no two gauge parts with the same gauge structures
 * in the given list.
 **)
CombineChiralParts[gaugeParts_List] :=
	Module[{leftChiralParts, rightChiralParts, chiralParts, allParts},
		leftChiralParts = Cases[gaugeParts,
			{_, Except[ZeroColouredVertex], LeftChiralVertex}];
		rightChiralParts = Cases[gaugeParts,
			{_, Except[ZeroColouredVertex], RightChiralVertex}];

		(* First we pair up the nonzero chiral parts with matching
		 * colour structures *)
		chiralParts = Flatten[Cases[rightChiralParts,
			{rightScalar_, #[[2]], RightChiralVertex} :>
			{{#[[1]], rightScalar}, #[[2]], ChiralVertex}] &
			/@ leftChiralParts, 1];

		leftChiralParts = Complement[leftChiralParts,
			chiralParts, SameTest -> (#1[[2]] == #2[[2]] &)];
		rightChiralParts = Complement[rightChiralParts,
			chiralParts, SameTest -> (#1[[2]] == #2[[2]] &)];

		(* The rest is paired up with zeros *)
		chiralParts = Join[chiralParts,
			{{#[[1]], 0}, #[[2]], ChiralVertex} & /@ leftChiralParts,
			{{0, #[[1]]}, #[[2]], ChiralVertex} & /@ rightChiralParts
		];

		allParts = Join[
			Cases[gaugeParts,
				Except[{_, _, LeftChiralVertex | RightChiralVertex}]],
			chiralParts
		];

		If[Length[allParts] === 0 && Length[gaugeParts] =!= 0,
			allParts = {{{0, 0}, ZeroColouredVertex, ChiralVertex}}];
		allParts
	]

(** \brief Given a list of gauge (colour and Lorentz) structures
 * with `TwoMetricVertex[]` parts with the same set of Lorentz indices,
 * combine them to `QuadrupleVectorVertex[]` parts.
 * Pad with zeros if necessary.
 * \param gaugeParts a list of gauge structures of the format
 * `{{scalar1, colourStructure1, LorentzStructure1}, ...}`
 * \returns the list of gauge structures where all `TwoMetricVertex[]`
 * parts are combined to `QuadrupleVectorVertex[]` parts.
 * \note There may be no two gauge parts with the same gauge structures
 * in the given list.
 **)
CombineQuadrupleVectorPartsInGroup[gaugeParts_List] :=
	Module[{g12g34Parts, g13g24Parts, g14g23Parts,
			parts1And2, parts1And3, parts2And3, vvvvParts},
		g12g34Parts = Cases[gaugeParts,
			{_, _, TwoMetricVertex[lIndex1_, lIndex2_, lIndex3_, lIndex4_]} /;
			OrderedQ[{lIndex1, lIndex2, lIndex3, lIndex4}]
		];
		g13g24Parts = Cases[gaugeParts,
			{_, _, TwoMetricVertex[lIndex1_, lIndex3_, lIndex2_, lIndex4_]} /;
			OrderedQ[{lIndex1, lIndex2, lIndex3, lIndex4}]
		];
		g14g23Parts = Cases[gaugeParts,
			{_, _, TwoMetricVertex[lIndex1_, lIndex4_, lIndex2_, lIndex3_]} /;
			OrderedQ[{lIndex1, lIndex2, lIndex3, lIndex4}]
		];

		Utils`AssertWithMessage[
			Length[g12g34Parts] + Length[g13g24Parts] +
			Length[g14g23Parts] === Length[gaugeParts],
			"CXXDiagrams`CombineQuadrupleVectorPartsInGroup[]: Incompatible \
indices encountered:" <> ToString[gaugeParts]
		];

		parts1And2 = Flatten[Cases[g13g24Parts,
			{part2Scalar_, #[[2]], #[[3]][[{1, 3, 2, 4}]]} :>
			{{#[[1]], part2Scalar}, #[[2]], #[[3]]}] & /@ g12g34Parts, 1];

		vvvvParts = Flatten[Cases[g14g23Parts,
			{part3Scalar_, #[[2]], #[[3]][[{1, 4, 2, 3}]]} :>
			{{#[[1, 1]], #[[1, 2]], part3Scalar}, #[[2]],
				QuadrupleVectorVertex @@ #[[3]]}] & /@
			parts1And2, 1];

		{g12g34Parts, g13g24Parts, g14g23Parts, parts1And2} = Complement[#,
			vvvvParts, SameTest -> (#1[[2]] == #2[[2]] &)] & /@
			{g12g34Parts, g13g24Parts, g14g23Parts, parts1And2};

		parts1And3 = Flatten[Cases[g13g24Parts,
			{part3Scalar_, #[[2]], #[[3]][[{1, 4, 2, 3}]]} :>
			{{#[[1]], part3Scalar}, #[[2]], #[[3]]}] & /@ g14g23Parts, 1];

		parts2And3 = Flatten[Cases[g13g24Parts,
			{part3Scalar_, #[[2]], #[[3]][[{1, 4, 2, 3}]]} :>
			{{#[[1]], part3Scalar}, #[[2]], #[[3]]}] & /@ g13g24Parts, 1];

		Join[vvvvParts,
			{{#[[1, 1]], #[[1, 2]], 0}, #[[2]],
				QuadrupleVectorVertex @@ #[[3]]} & /@ parts1And2,
			{{#[[1, 1]], 0, #[[1, 2]]}, #[[2]],
				QuadrupleVectorVertex @@ #[[3]]} & /@ parts1And3,
			{{0, #[[1, 1]], #[[1, 2]]}, #[[2]],
				QuadrupleVectorVertex @@ #[[3]]} & /@ parts2And3
		]
	]

(** \brief Given a list of gauge (colour and Lorentz) structures
 * combine all `TwoMetricVertex[]` parts to `QuadrupleVectorVertex[]`
 * parts. Pad with zeros if necessary.
 * \param gaugeParts a list of gauge structures of the format
 * `{{scalar1, colourStructure1, LorentzStructure1}, ...}`
 * \returns the list of gauge structures where all `TwoMetricVertex[]`
 * parts are combined to `QuadrupleVectorVertex[]` parts.
 * \note There may be no two gauge parts with the same gauge structures
 * in the given list.
 **)
CombineQuadrupleVectorParts[gaugeParts_List] :=
	Module[{ggParts, ggGroups, allParts},
		ggParts = Cases[gaugeParts,
			{_, Except[ZeroColouredVertex], _TwoMetricVertex}];
		ggGroups = GatherBy[ggParts, Sort[List @@ #[[3]]] &];

		allParts = Join[
			Cases[gaugeParts, Except[{_, _, _TwoMetricVertex}]],
			Flatten[CombineQuadrupleVectorPartsInGroup /@ ggGroups, 1]
		];

		If[Length[allParts] === 0 && Length[gaugeParts] =!= 0,
			allParts = Cases[gaugeParts,
				{_, ZeroColouredVertex, gg_TwoMetricVertex} :>
				{{0, 0, 0}, ZeroColouredVertex,
					QuadrupleVectorVertex @@ gg}][[1;;1]]
		];
		allParts
	]

(** \brief Given a list of gauge (colour and Lorentz) structure
 * combine all parts with gauge structures that are internal to
 * CXXDiagrams to such that are exported by CXXDiagrams. These are:
 * - LeftChiralVertex
 * - RightChiralVertex
 * - TwoMetricVertex
 * Pad with zeros if necessary.
 * \param gaugeParts a list of gauge structures of the format
 * `{{scalar1, colourStructure1, LorentzStructure1}, ...}`
 * \returns the resulting gauge structure where all gauge structures
 * that are internal to CXXDiagrams are combined to such that are
 * exported by CXXDiagrams.
 * \note There may be no two gauge parts with the same gauge structures
 * in the given list.
 * \note The resulting gauge structure has to decompose as
 * `{scalar1, ...} * ColourStructure * LorentzStructure`
 * otherwise a fatal error is issued.
 **)
FullGaugeStructureFromParts[gaugeParts_List] :=
	Module[{gaugeStructures},
		gaugeStructures = CombineChiralParts[gaugeParts];
		gaugeStructures = CombineQuadrupleVectorParts[gaugeStructures];

		Utils`AssertWithMessage[Length[gaugeStructures] === 1,
			"CXXDiagrams`FullGaugeStructureFromParts[]: Vertex with \
composite gauge structure encountered:" <> ToString[gaugeStructures]];

		gaugeStructures[[1]]
	]

(** \brief Given a vertex return its gauge (colour and Lorentz)
 * structure as determined by `FullGaugeStructureFromParts[]`
 * when applied to the gauge parts of the vertex.
 * \param vertex the given vertex as a list of fields
 * \returns the gauge (colour and Lorentz) structure as determined
 * by `FullGaugeStructureFromParts[]` when applied to the gauge parts
 * of the vertex.
 **)
GaugeStructureOfVertex[vertex_] :=
	Module[{lorentzParts, gaugeStructures},
		(* This is arguably a bug in SARAH causing e.g {hh, hh, VZ} to be mapped to {0, 0} *)
		lorentzParts = If[MatchQ[vertex, {_, {_, 0}}],
			Utils`AssertWithMessage[vertex[[2,1]] === 0,
				"CXXDiagrams`SortedVertex[]: ill-formed vertex: " <> ToString[vertex]];
			{{0, MomentumDifferenceVertex[vertex[[1,1]], vertex[[1,2]]]}},
			LabelLorentzPart /@ Rest[vertex]];

		gaugeStructures = GaugeStructureOfVertexLorentzPart /@ lorentzParts;
		FullGaugeStructureFromParts[gaugeStructures]
	]

(** \brief Creates c++ code that makes functions available that
 * numerically evaluate any of the given vertices.
 * \param vertices a list of vertices
 * \param MaximumVerticesLimit An integer option that specify an upper
 * limit of vertices that shall go into a single block of code.
 * \returns a list `{{prototypes1, definitions1}, ...}` containing the
 * corresponding c++ code where no sublist contains more than
 * `MaximumVerticesLimit` number of vertices.
 **)
CreateVertices::errLostVertices =
"Some vertices lost after splitting of cxxVertices into multiple lists.";
CreateVertices::errMaximumVerticesLimit =
"CXXDiagrams`.`MaximumVerticesLimit should be positive integer number.";
CreateVertices[
   vertices:{{__}...},
   OptionsPattern[{MaximumVerticesLimit -> 500}]] :=
Module[{cxxVertices, vertexPartition},
   cxxVertices = CreateVertex /@ DeleteDuplicates[vertices];

   (* Mathematica 7 does not support the `UpTo[n]` notation *)
   vertexPartition = Partition[cxxVertices, OptionValue[MaximumVerticesLimit]];
   (* Partition splits cxxVertices into lists of length MaximumVerticesLimit.
      If the length of cxxVertices is not a multiple of MaximumVerticesLimit,
      some vertices will be discarded! (This is why Complement is used) *)
   If[# =!= {}, AppendTo[vertexPartition,#] ] &@ Complement[cxxVertices, Sequence @@ vertexPartition];

   Utils`AssertOrQuit[Sort[cxxVertices] === Sort[Join@@vertexPartition],CreateVertices::errLostVertices];

   Map[StringJoin[Riffle[#, "\n\n"]] &, Transpose /@ vertexPartition, {2}]
] /; Utils`AssertOrQuit[And[IntegerQ@OptionValue@MaximumVerticesLimit, OptionValue@MaximumVerticesLimit>0],CreateVertices::errMaximumVerticesLimit];
Utils`MakeUnknownInputDefinition@CreateVertices;
(** \brief Creates c++ code that makes a function available that
 * numerically evaluates the given vertex.
 * \param fields a vertex given as a list of fields
 * \returns a list {prototypes, definitions} containing the
 * corresponding c++ code.
 **)
CreateVertex[fields_List] :=
  Module[{fieldSequence},
		fieldSequence = StringJoin @ Riffle[
			CXXNameOfField[#, prefixNamespace -> "fields"] & /@ fields, ", "];

		{
		"template<> struct VertexImpl<" <> fieldSequence <> ">" <> "\n" <>
		"{\n" <> TextFormatting`IndentText[
			"static " <> SymbolName[VertexTypeForFields[fields]] <>
				" evaluate(const std::array<int, " <>
				ToString[Total[NumberOfFieldIndices /@ fields]] <>
			">& indices, const context_base& context);"] <> "\n" <>
		"};"
		,
		SymbolName[VertexTypeForFields[fields]] <>
			" VertexImpl<" <> fieldSequence <> ">::evaluate(\n" <>
				TextFormatting`IndentText["const std::array<int, " <>
					ToString[Total[NumberOfFieldIndices /@ fields]] <> ">& indices, " <>
				"const context_base& context)"] <> "\n" <>
		"{\n" <> TextFormatting`IndentText[
				VertexFunctionBodyForFields[fields] <> "\n"] <>
		"}"
		}
  ]

(** \brief Returns the Lorentz structure of a given vertex.
 * param fields a vertex given as a list of fields
 * \returns the Lorentz structure of a given vertex.
 **)
VertexTypeForFields[fields_List] :=
	AtomHead[GaugeStructureOfVertex[SortedVertex[fields]][[3]]]

(** \brief Creates the c++ code containing the actual numerical
 * evaluation of a given vertex.
 * \param fields a vertex given as a list of fields
 * \returns the c++ code containing the actual numerical
 * evaluation of a given vertex.
 **)
VertexFunctionBodyForFields[fields_List] :=
	Module[{sortedVertex, sortedFields, sortedIndexedFields, indexedFields,
			indexFields, gaugeStructure,
			fieldsOrdering, sortedFieldsOrdering, inverseFOrdering,
			fOrderingWRTSortedF, expr, exprL, exprR, expr1, expr2, expr3,
			vertexRules, incomingScalar, outgoingScalar, incomingGhost,
			lIndex1, lIndex2, lIndex3, lIndex4},
		sortedVertex = SortedVertex[fields, ApplyGUTNormalization -> True];
		sortedIndexedFields = sortedVertex[[1]];
		sortedFields = Vertices`StripFieldIndices /@ sortedIndexedFields;

		indexFields = {SARAH`Cp @@ sortedFields -> SARAH`Cp @@ sortedIndexedFields};

		(* Mathematica 7 does not know about permutations... :'-( *)
		fieldsOrdering = Ordering[fields];
		sortedFieldsOrdering = Ordering[sortedFields];

		inverseFOrdering = Ordering[fieldsOrdering];
		fOrderingWRTSortedF = sortedFieldsOrdering[[inverseFOrdering]];

		indexedFields = sortedIndexedFields[[fOrderingWRTSortedF]];

		gaugeStructure = GaugeStructureOfVertex[sortedVertex];
		Switch[gaugeStructure[[3]],
			ScalarVertex | InverseMetricVertex,
			vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> gaugeStructure[[1]]};

			expr = Vertices`SortCp[SARAH`Cp @@ fields] /. indexFields /. vertexRules;

			DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
			Parameters`CreateLocalConstRefs[expr] <> "\n" <>
			"const " <> GetComplexScalarCType[] <> " result = " <>
			Parameters`ExpressionToString[expr] <> ";\n\n" <>
			"return {result};",

			ChiralVertex,
			vertexRules = {
				(SARAH`Cp @@ sortedIndexedFields)[SARAH`PL] -> gaugeStructure[[1,1]],
				(SARAH`Cp @@ sortedIndexedFields)[SARAH`PR] -> gaugeStructure[[1,2]]
			};

			exprL = Vertices`SortCp[(SARAH`Cp @@ fields)[SARAH`PL]] /.
				indexFields /. vertexRules;
			exprR = Vertices`SortCp[(SARAH`Cp @@ fields)[SARAH`PR]] /.
				indexFields /. vertexRules;

			DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
			Parameters`CreateLocalConstRefs[{exprL, exprR}] <> "\n" <>
			"const " <> GetComplexScalarCType[] <> " left = " <>
				Parameters`ExpressionToString[exprL] <> ";\n\n" <>
			"const " <> GetComplexScalarCType[] <> " right = " <>
			Parameters`ExpressionToString[exprR] <> ";\n\n" <>
			"return {left, right};",

			_MomentumVertex,
			incomingGhost = gaugeStructure[[3,1]];
			vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> gaugeStructure[[1]]};

			expr = Vertices`SortCp[SARAH`Cp @@ fields] /. indexFields /. vertexRules;

			DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
			Parameters`CreateLocalConstRefs[expr] <> "\n" <>
			"const " <> GetComplexScalarCType[] <> " result = " <>
			Parameters`ExpressionToString[expr] <> ";\n\n" <>
			"return {result, " <>
				ToString@Utils`MathIndexToCPP[Position[indexedFields, incomingGhost, {1}][[1,1]]] <>
			"};",

			_TripleVectorVertex,
			{lIndex1, lIndex2, lIndex3} = LorentzIndexOfField /@
				sortedIndexedFields;

      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> gaugeStructure[[1]]};

			expr = Vertices`SortCp[SARAH`Cp @@ fields] /. indexFields /. vertexRules;

			DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return {result, " <>
				If[Or @@ ((And @@
						(RotateLeft[{lIndex1, lIndex2, lIndex3}, #] ===
							LorentzIndexOfField /@ gaugeStructure[[3]])) & /@
						{0, 1, 2}),
					"TripleVectorVertex::even_permutation{}",
					"TripleVectorVertex::odd_permutation{}"] <>
			"};",

			_QuadrupleVectorVertex,
			{lIndex1, lIndex2, lIndex3, lIndex4} = LorentzIndexOfField /@
				sortedIndexedFields;

      vertexRules = {
				(SARAH`Cp @@ sortedIndexedFields)[
					SARAH`g[lIndex1, lIndex2] * SARAH`g[lIndex3, lIndex4]] ->
        gaugeStructure[[1, QuadrupleVectorVertexPartForLorentzIndices[
					lIndex1, lIndex2, lIndex3, lIndex4]]],
				(SARAH`Cp @@ sortedIndexedFields)[
					SARAH`g[lIndex1, lIndex3] * SARAH`g[lIndex2, lIndex4]] ->
        gaugeStructure[[1, QuadrupleVectorVertexPartForLorentzIndices[
					lIndex1, lIndex3, lIndex2, lIndex4]]],
				(SARAH`Cp @@ sortedIndexedFields)[
					SARAH`g[lIndex1, lIndex4] * SARAH`g[lIndex2, lIndex3]] ->
        gaugeStructure[[1, QuadrupleVectorVertexPartForLorentzIndices[
					lIndex1, lIndex4, lIndex2, lIndex3]]]
			};

      expr1 = Vertices`SortCp[(SARAH`Cp @@ fields)[
					SARAH`g[lIndex1, lIndex2] * SARAH`g[lIndex3, lIndex4]]] /.
						indexFields /. vertexRules;
      expr2 = Vertices`SortCp[(SARAH`Cp @@ fields)[
					SARAH`g[lIndex1, lIndex3] * SARAH`g[lIndex2, lIndex4]]] /.
						indexFields /. vertexRules;
      expr3 = Vertices`SortCp[(SARAH`Cp @@ fields)[
					SARAH`g[lIndex1, lIndex4] * SARAH`g[lIndex2, lIndex3]]] /. indexFields
						/. vertexRules;

			DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[{expr1, expr2, expr3}] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " part1 = " <>
      Parameters`ExpressionToString[expr1] <> ";\n\n" <>
      "const " <> GetComplexScalarCType[] <> " part2 = " <>
      Parameters`ExpressionToString[expr2] <> ";\n\n" <>
      "const " <> GetComplexScalarCType[] <> " part3 = " <>
      Parameters`ExpressionToString[expr3] <> ";\n\n" <>
      "return {part1, part2, part3};",

      _MomentumDifferenceVertex,
      {incomingScalar, outgoingScalar} = List @@ gaugeStructure[[3]];
      vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> gaugeStructure[[1]]};

      expr = Vertices`SortCp[SARAH`Cp @@ fields] /. indexFields /. vertexRules;

      "int minuend_index = " <>
        ToString@Utils`MathIndexToCPP[Position[indexedFields, incomingScalar, {1}][[1,1]]] <> ";\n" <>
      "int subtrahend_index = " <>
        ToString@Utils`MathIndexToCPP[Position[indexedFields, outgoingScalar, {1}][[1,1]]] <> ";\n\n" <>
      DeclareIndices[StripUnbrokenGaugeIndices /@ indexedFields, "indices"] <>
      Parameters`CreateLocalConstRefs[expr] <> "\n" <>
      "const " <> GetComplexScalarCType[] <> " result = " <>
      Parameters`ExpressionToString[expr] <> ";\n\n" <>
      "return {result, minuend_index, subtrahend_index};",

      _,
      Print["Unrecognized gauge structure: " <> ToString[gaugeStructure[[3]]]];
      Quit[1]
    ]
  ]

(** \brief Given a sequence of Lorentz indices determine which part of
 * a `QuadrupleVectorVertex` is encoded by
 * `g[lIndex1, lIndex2] * g[lIndex3, lIndex4]`
 * \param lIndex1 the first Lorentz index
 * \param lIndex2 the second Lorentz index
 * \param lIndex3 the third Lorentz index
 * \param lIndex4 the fourth Lorentz index
 * \returns the corresponding part of the `QuadrupleVectorVertex`
 **)
QuadrupleVectorVertexPartForLorentzIndices[
		lIndex1_, lIndex2_, lIndex3_, lIndex4_] :=
	Module[{properlyOrderedindices},
		properlyOrderedindices = Flatten[List @@@ (List @@
			(SARAH`g[lIndex1, lIndex2] * SARAH`g[lIndex3, lIndex4]))];

		Switch[Ordering[properlyOrderedindices],
			{1, 2, 3, 4}, 1,
			{1, 3, 2, 4}, 2,
			{1, 3, 4, 2}, 3,
			_, Print["CXXDiagrams`QuadrupleVectorVertexPartForLorentzIndices[]: \
Cannot properly order Lorentz indices: " <>
ToString[{lIndex1, lIndex2, lIndex3, lIndex4}]]; Quit[1]]
	]

(** \brief a helper function that declares the local indices used
 * by `VertexFunctionBodyForFields[]`
 **)
DeclareIndices[indexedFields_List, arrayName_String] :=
    Module[{p, total = 0, fieldIndexList, decl = ""},
           DeclareIndex[idx_, num_Integer, an_String] := (
               "const int " <> CConversion`ToValidCSymbolString[idx] <>
               " = " <> an <> "[" <> ToString[num] <> "];\n");
           For[p = 1, p <= Length[indexedFields], p++,
               fieldIndexList = Vertices`FieldIndexList[indexedFields[[p]]];
               decl = decl <> StringJoin[DeclareIndex[#, total++, arrayName]& /@ fieldIndexList];
              ];
           Assert[total == Total[Length[Vertices`FieldIndexList[#]]& /@ indexedFields]];
           decl
          ]

(** \brief a helper function that returns a c++ type capable of storing
 * a complex scalar.
 **)
GetComplexScalarCType[] :=
    CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]]

(** \brief Creates c++ code that makes functions available that
 * return tree-level masses of given fields.
 * \returns the corresponding c++ code as a string.
 **)
CreateMassFunctions[] :=
  Module[{massiveFields,
          ghostMappings = SelfEnergies`ReplaceGhosts[FlexibleSUSY`FSEigenstates]},
    massiveFields = TreeMasses`GetParticles[];

    StringJoin @ Riffle[
      Module[{fieldInfo = TreeMasses`FieldInfo[#], numberOfIndices},
             numberOfIndices = Length @ fieldInfo[[5]];

             "template<> inline\n" <>
             "double context_base::mass_impl<" <>
               CXXNameOfField[#, prefixNamespace -> "fields"] <>
             ">(const std::array<int, " <> ToString @ numberOfIndices <>
             ">& indices) const\n" <>
             "{ return model.get_M" <> CXXNameOfField[# /. ghostMappings] <>
             If[TreeMasses`GetDimension[#] === 1, "()", "(indices[0])"] <> "; }"
            ] & /@ massiveFields, "\n\n"]
        ];

(** \brief Creates c++ code that makes functions available that
 * return physical masses of given fields.
 * \returns the corresponding c++ code as a string.
 **)
CreatePhysicalMassFunctions[fieldsNamespace_:""] :=
  Module[{massiveFields,
          ghostMappings = SelfEnergies`ReplaceGhosts[FlexibleSUSY`FSEigenstates]},
    massiveFields = TreeMasses`GetParticles[];

    StringJoin @ Riffle[
      Module[{fieldInfo = TreeMasses`FieldInfo[#], numberOfIndices},
             numberOfIndices = Length @ fieldInfo[[5]];

             "template<> inline\n" <>
             "double context_base::physical_mass_impl<" <>
               CXXNameOfField[#, prefixNamespace -> "fields"] <>
             ">(const std::array<int, " <> ToString @ numberOfIndices <>
             ">& indices) const\n" <>
             "{ return model.get_physical().M" <> CXXNameOfField[# /. ghostMappings] <>
             If[TreeMasses`GetDimension[#] === 1, "", "[indices[0]]"] <> "; }"
            ] & /@ massiveFields, "\n\n"]
        ];

(** \brief Creates the c++ code for a function that returns the
 * numerical value of the electrical charge of the electron.
 * \returns the corresponding c++ code as a string.
 **)
CreateUnitCharge[] :=
  Module[{electron,photon,vertex,vertexBody,
          numberOfElectronIndices,numberOfPhotonIndices},
         electron = AtomHead @ TreeMasses`GetSMElectronLepton[];
         photon = SARAH`Photon;
         vertex = {SARAH`bar[electron], electron, photon};
         vertexBody = VertexFunctionBodyForFields[vertex];
         numberOfElectronIndices = NumberOfFieldIndices[electron];
         numberOfPhotonIndices = NumberOfFieldIndices[photon];

         "ChiralVertex unit_charge(const context_base& context)\n" <>
         "{\n" <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString @ numberOfElectronIndices <> "> electron_indices = {" <>
              If[TreeMasses`GetDimension[electron] =!= 1,
                 " " <> ToString @ (TreeMasses`FieldInfo[electron][[2]]-1) <> (* Electron has the lowest index *)
                 If[numberOfElectronIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfElectronIndices-1}],
                    ""] <> " ",
                 If[numberOfElectronIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfElectronIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString @ numberOfPhotonIndices <> "> photon_indices = {" <>
               If[TreeMasses`GetDimension[photon] =!= 1,
                 " " <> ToString @ (TreeMasses`FieldInfo[photon][[2]]-1) <>
                 If[numberOfPhotonIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfPhotonIndices-1}],
                    ""] <> " ",
                 If[numberOfPhotonIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfPhotonIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString[
              Total[NumberOfFieldIndices /@ {photon,electron,electron}]] <>
            "> indices = " <>
              "concatenate(photon_indices, electron_indices, electron_indices);\n\n") <>

         TextFormatting`IndentText @ vertexBody <> "\n" <>
         "}"
  ]

(** \brief Return the number of indices a field has as would be
 * determined by inspecting the result of ``TreeMasses`FieldInfo[]``.
 * \param field the given field
 * \returns the number of indices a field has as would be
 * determined by inspecting the result of ``TreeMasses`FieldInfo[]``.
 **)
NumberOfFieldIndices[field_] := Length @ TreeMasses`FieldInfo[field][[5]]

(** \brief Return the index bounds a field has as would be
 * determined by inspecting the result of ``TreeMasses`FieldInfo[]``.
 * The result is of the format `{{min1, ...}, {max1, ...}}`.
 * \param field the given field
 * \returns the number of indices a field has as would be
 * determined by inspecting the result of ``TreeMasses`FieldInfo[]``.
 **)
IndexBoundsForField[field_] :=
  Module[{fieldInfo = TreeMasses`FieldInfo[field]},
    If[NumberOfFieldIndices[field] === 0,
       Return[{{},{}}]];
    If[Length @ Cases[fieldInfo[[5]],{SARAH`generation,_}] === 0,
       Transpose[{1,#[[2]]} & /@ fieldInfo[[5]]],
       Transpose[Prepend[
         {1,#[[2]]} & /@ DeleteCases[fieldInfo[[5]],{SARAH`generation,_}],
         {fieldInfo[[2]],fieldInfo[[3]]}]]]]

(** \brief If the SARAH vertices are not yet available, force SARAH to
 * reread its vertex files.
 **)
LoadVerticesIfNecessary[] :=
   If[Head[SARAH`VertexList3] === Symbol || Length[SARAH`VertexList3] === 0,
        SA`CurrentStates = FlexibleSUSY`FSEigenstates;
        SARAH`InitVertexCalculation[FlexibleSUSY`FSEigenstates, False];
        SARAH`partDefinition = ParticleDefinitions[FlexibleSUSY`FSEigenstates];
        SARAH`Particles[SARAH`Current] = SARAH`Particles[FlexibleSUSY`FSEigenstates];
        SARAH`ReadVertexList[FlexibleSUSY`FSEigenstates, False, False, True];
        SARAH`MakeCouplingLists;
   ]

(** \brief Remove any Lorentz and colour indices of a given field
 * \param p the given field
 * \returns the field with any Lorentz and colour indices removed.
 **)
StripUnbrokenGaugeIndices[p_] := StripColourIndices[StripLorentzIndices[p]]

(** \brief Remove any Lorentz indices of a given field
 * \param p the given field
 * \returns the field with any Lorentz indices removed.
 **)
StripLorentzIndices[p_Symbol] := p
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]]
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]]
StripLorentzIndices[p_] :=
	Module[{remainingIndices},
		remainingIndices = Select[p[[1]], (!Vertices`SarahLorentzIndexQ[#] &)];
		If[Length[remainingIndices] === 0, Head[p],
			Head[p][remainingIndices]]
	]

(** \brief Remove any colour indices of a given field
 * \param p the given field
 * \returns the field with any colour indices removed.
 **)
StripColourIndices[p_Symbol] := p
StripColourIndices[SARAH`bar[p_]] := SARAH`bar[StripColourIndices[p]]
StripColourIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripColourIndices[p]]
StripColourIndices[p_] :=
	Module[{remainingIndices},
		remainingIndices = Select[p[[1]], (!Vertices`SarahColorIndexQ[#] &)];
		If[Length[remainingIndices] === 0, Head[p],
			Head[p][remainingIndices]]
	]

ColorFactorForDiagram[topology_, diagram_] :=
   ColourFactorForIndexedDiagramFromGraph[
      CXXDiagrams`IndexDiagramFromGraph[diagram, topology], topology
   ];

(* TODO: generalize the Extraction *)
ExtractColourFactor[colourfactor_ * SARAH`Lam[ctIndex1_, ctIndex2_, ctIndex3_] /; NumericQ[colourfactor]] := 2*colourfactor;
ExtractColourFactor[SARAH`Lam[ctIndex1_, ctIndex2_, ctIndex3_]] := 2;
ExtractColourFactor[colourfactor_ * SARAH`Delta[ctIndex1_, ctIndex2_] /; NumericQ[colourfactor]] := colourfactor;
ExtractColourFactor[SARAH`Delta[ctIndex1_, ctIndex2_]] := 1;
ExtractColourFactor[colourfactor_ /; NumericQ[colourfactor]] := colourfactor;
ExtractColourFactor[args___] :=
   (Print["Error: ExtractColourFactor cannot convert argument ", args]; Quit[1]);

End[];
EndPackage[];
