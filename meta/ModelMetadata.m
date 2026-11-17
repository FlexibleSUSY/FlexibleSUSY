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

(* Machine-readable description of a generated model, written to
   models/<M>/<M>_metadata.json.  It is read by the Python interface
   (python/flexiblesusy) to map input parameters to SLHA blocks and to
   interpret the SLHA output, so that nothing model-specific is
   hard-coded there. *)

BeginPackage["ModelMetadata`", {"CConversion`", "Parameters`", "TreeMasses`", "Utils`"}];

$formatVersion::usage = "Version of the format of the metadata file.  Increase it
when the meaning of an existing field changes.";

TypeToMetadata::usage = "TypeToMetadata[type] converts a CConversion type into an
association <|\"kind\" -> \"real\"|\"complex\"|\"integer\", \"shape\" -> {dims...}|>.";

SLHABlockToMetadata::usage = "SLHABlockToMetadata[block] converts the SLHA block
specification of a parameter, block or {block, index}, into an association
<|\"block\" -> name, \"index\" -> index or Null|>.  Returns Null if the parameter
is not read from or written to SLHA.";

InputParameterToMetadata::usage = "InputParameterToMetadata[{par, block, type}]
describes one input parameter.";

ParticleToMetadata::usage = "ParticleToMetadata[particle] describes one mass
eigenstate: name, multiplicity, PDG codes and whether it is a SM particle.";

OutputParameterToMetadata::usage = "OutputParameterToMetadata[{par, block}]
describes a parameter written to the SLHA output.";

CreateModelMetadata::usage = "CreateModelMetadata[modelName, version,
inputParameters, particles, lesHouchesParameters] returns the metadata as an
association.";

WriteModelMetadata::usage = "WriteModelMetadata[fileName, metadata] writes the
metadata as JSON.";

Begin["`Private`"];

$formatVersion = 1;

ScalarKind[CConversion`realScalarCType]    := "real";
ScalarKind[CConversion`complexScalarCType] := "complex";
ScalarKind[CConversion`integerScalarCType] := "integer";
Utils`MakeUnknownInputDefinition[ScalarKind];

TypeToMetadata[CConversion`ScalarType[t_]] :=
    <|"kind" -> ScalarKind[t], "shape" -> {}|>;

TypeToMetadata[(CConversion`ArrayType | CConversion`VectorType)[t_, n_Integer]] :=
    <|"kind" -> ScalarKind[t], "shape" -> {n}|>;

TypeToMetadata[CConversion`MatrixType[t_, m_Integer, n_Integer]] :=
    <|"kind" -> ScalarKind[t], "shape" -> {m, n}|>;

TypeToMetadata[CConversion`TensorType[t_, dims__Integer]] :=
    <|"kind" -> ScalarKind[t], "shape" -> {dims}|>;

Utils`MakeUnknownInputDefinition[TypeToMetadata];

SLHABlockToMetadata[None | Null | {}] := Null;

SLHABlockToMetadata[{block_, idx_Integer}] :=
    <|"block" -> ToString[block], "index" -> idx|>;

SLHABlockToMetadata[block_ /; Head[block] =!= List] :=
    <|"block" -> ToString[block], "index" -> Null|>;

Utils`MakeUnknownInputDefinition[SLHABlockToMetadata];

InputParameterToMetadata[{par_, block_, type_}] :=
    Join[
        <|"name" -> CConversion`ToValidCSymbolString[par]|>,
        TypeToMetadata[type],
        <|"slha" -> SLHABlockToMetadata[block]|>
    ];

Utils`MakeUnknownInputDefinition[InputParameterToMetadata];

ParticleToMetadata[particle_] :=
    Module[{multiplicity = TreeMasses`GetDimension[particle],
            pdg = Parameters`GetPDGCodesForParticle[particle]},
           Utils`AssertOrQuit[IntegerQ[multiplicity] && multiplicity > 0,
                              ParticleToMetadata::multiplicity, particle, multiplicity];
           Utils`AssertOrQuit[VectorQ[pdg, IntegerQ],
                              ParticleToMetadata::pdg, particle, pdg];
           <|"name" -> CConversion`ToValidCSymbolString[particle],
             "multiplicity" -> multiplicity,
             "pdg" -> pdg,
             "sm" -> TrueQ[TreeMasses`IsSMParticle[particle]]|>
          ];

ParticleToMetadata::multiplicity = "Particle `1` has an invalid multiplicity `2`.";
ParticleToMetadata::pdg = "Particle `1` has invalid PDG codes `2`.";

(* Factor between the value in the SLHA output and the parameter: the gauge
   couplings in Block gauge are written without GUT normalization (cf.
   WriteOut`ReadSLHAOutputBlock). *)
SLHAFactor[par_, {block_, _Integer}] /;
    ToUpperCase[ToString[block]] === "GAUGE" && Parameters`IsGaugeCoupling[par] :=
    N[Parameters`GetGUTNormalization[par]];

SLHAFactor[_, _] := 1;

OutputParameterToMetadata[{par_, block_}] :=
    Join[
        <|"name" -> CConversion`ToValidCSymbolString[par]|>,
        TypeToMetadata[Parameters`GetType[par]],
        <|"slha" -> SLHABlockToMetadata[block],
          "slha_factor" -> SLHAFactor[par, block]|>
    ];

Utils`MakeUnknownInputDefinition[OutputParameterToMetadata];

(* parameters of FSLesHouchesList that are written to the SLHA output:
   mixing matrices (output parameters) and model parameters *)
SelectOutputParameters[lesHouchesParameters_List, pars_List] :=
    Cases[lesHouchesParameters,
          {par_ /; MemberQ[pars, par], block_ /; block =!= None, ___} :> {par, block}];

CreateModelMetadata[modelName_String, version_String, inputParameters_List,
                    particles_List, lesHouchesParameters_List] :=
    <|"format_version" -> $formatVersion,
      "model" -> modelName,
      "flexiblesusy_version" -> version,
      "input_parameters" -> (InputParameterToMetadata /@ inputParameters),
      "particles" -> (ParticleToMetadata /@ particles),
      "mixing_matrices" -> (OutputParameterToMetadata /@
          SelectOutputParameters[lesHouchesParameters, Parameters`GetOutputParameters[]]),
      "model_parameters" -> (OutputParameterToMetadata /@
          SelectOutputParameters[lesHouchesParameters, Parameters`GetModelParameters[]])
    |>;

Utils`MakeUnknownInputDefinition[CreateModelMetadata];

WriteModelMetadata[fileName_String, metadata_Association] :=
    Module[{result = Export[fileName, metadata, "RawJSON", "Compact" -> 2]},
           Utils`AssertOrQuit[result === fileName, WriteModelMetadata::failed, fileName];
          ];

WriteModelMetadata::failed = "Could not write the model metadata to `1`.";

Utils`MakeUnknownInputDefinition[WriteModelMetadata];

End[];

EndPackage[];
