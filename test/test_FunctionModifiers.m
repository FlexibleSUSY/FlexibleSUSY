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
Needs["FunctionModifiers`", "FunctionModifiers.m"];

TestEquality[
   MakeOverride["double get_ewsb_eq_hh_1() const;"],
   "double get_ewsb_eq_hh_1() const override;"
];
TestEquality[
   MakeOverride["double get_ewsb_eq_hh_1() const;\n"],
   "double get_ewsb_eq_hh_1() const override;\n"
];


TestEquality[
   MakeOverride["double get_ewsb_eq_hh_1();"],
   "double get_ewsb_eq_hh_1() override;"
];
TestEquality[
   MakeOverride["double get_ewsb_eq_hh_1();\n"],
   "double get_ewsb_eq_hh_1() override;\n"
];

TestEquality[
   MakeOverride["void set_v(double v_) { v = v_; }"],
   "void set_v(double v_) override { v = v_; }"
];

TestEquality[
   MakeOverride[
      "std::complex<double> CpbarFdFdAhPR(int gI1, int gI2) const;\n" <>
      "std::complex<double> CpbarFdFdAhPL(int gI1, int gI2) const;\n" <>
      "std::complex<double> CpbarFeFeAhPR(int gI1, int gI2) const;"
   ],
   "std::complex<double> CpbarFdFdAhPR(int gI1, int gI2) const override;\n" <>
   "std::complex<double> CpbarFdFdAhPL(int gI1, int gI2) const override;\n" <>
   "std::complex<double> CpbarFeFeAhPR(int gI1, int gI2) const override;"
];

PrintTestSummary[];
