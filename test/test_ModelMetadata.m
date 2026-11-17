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

Needs["TestSuite`", "TestSuite.m"];
Needs["ModelMetadata`", "ModelMetadata.m"];

real    = CConversion`realScalarCType;
complex = CConversion`complexScalarCType;
integer = CConversion`integerScalarCType;

Print["testing TypeToMetadata[] ..."];

TestEquality[TypeToMetadata[CConversion`ScalarType[real]],
             <|"kind" -> "real", "shape" -> {}|>];
TestEquality[TypeToMetadata[CConversion`ScalarType[integer]],
             <|"kind" -> "integer", "shape" -> {}|>];
TestEquality[TypeToMetadata[CConversion`VectorType[complex, 3]],
             <|"kind" -> "complex", "shape" -> {3}|>];
TestEquality[TypeToMetadata[CConversion`ArrayType[real, 2]],
             <|"kind" -> "real", "shape" -> {2}|>];
TestEquality[TypeToMetadata[CConversion`MatrixType[real, 3, 2]],
             <|"kind" -> "real", "shape" -> {3, 2}|>];
TestEquality[TypeToMetadata[CConversion`TensorType[complex, 2, 3, 4]],
             <|"kind" -> "complex", "shape" -> {2, 3, 4}|>];

Print["testing SLHABlockToMetadata[] ..."];

TestEquality[SLHABlockToMetadata[{"MINPAR", 3}], <|"block" -> "MINPAR", "index" -> 3|>];
TestEquality[SLHABlockToMetadata[{MSOFTIN, 1}], <|"block" -> "MSOFTIN", "index" -> 1|>];
TestEquality[SLHABlockToMetadata[MSQ2IN], <|"block" -> "MSQ2IN", "index" -> Null|>];
TestEquality[SLHABlockToMetadata[None], Null];

Print["testing InputParameterToMetadata[] ..."];

TestEquality[InputParameterToMetadata[{TanBeta, {"MINPAR", 3}, CConversion`ScalarType[real]}],
             <|"name" -> "TanBeta", "kind" -> "real", "shape" -> {},
               "slha" -> <|"block" -> "MINPAR", "index" -> 3|>|>];
TestEquality[InputParameterToMetadata[{mq2Input, MSQ2IN, CConversion`MatrixType[real, 3, 3]}],
             <|"name" -> "mq2Input", "kind" -> "real", "shape" -> {3, 3},
               "slha" -> <|"block" -> "MSQ2IN", "index" -> Null|>|>];

Print["testing JSON export ..."];

json = ExportString[
    <|"input_parameters" ->
        {InputParameterToMetadata[{TanBeta, {"MINPAR", 3}, CConversion`ScalarType[real]}],
         InputParameterToMetadata[{x, None, CConversion`ScalarType[real]}]}|>,
    "RawJSON"];

TestEquality[ImportString[json, "RawJSON"],
             <|"input_parameters" ->
                 {<|"name" -> "TanBeta", "kind" -> "real", "shape" -> {},
                    "slha" -> <|"block" -> "MINPAR", "index" -> 3|>|>,
                  <|"name" -> "x", "kind" -> "real", "shape" -> {}, "slha" -> Null|>}|>];

PrintTestSummary[];
