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

t1l    = << "meta/MSSM/tquark_1loop_strong.m";
t1lqcd = << "meta/MSSM/tquark_1loop_qcd.m" /. GSY -> GS;
t2l    = << "meta/MSSM/tquark_2loop_strong.m";
t2lqcd = << "meta/MSSM/tquark_2loop_qcd.m" /. GSY -> GS;

colorCA = 3; colorCF = 4/3; Tf = 1/2; GS = g3;
MGl = mgl; MT = mt; SX = 2 mt Xt;
s2t = SX/(mmst1 - mmst2);

fin[0, args__] := fin[args];

loopFunctions = {
    Hmine[mm1_,mm2_] :> 2 * PolyLog[2, 1-mm1/mm2] + 1/2 * Log[mm1/mm2]^2,
    fin[mm1_,mm2_]   :> 1/2 * ( - (mm1 + mm2) * ( 7 + Zeta[2] )  
                                + 6 * (mm1 * Log[mm1/mmu] + mm2 * Log[mm2/mmu])
                                - 2 * (mm1 * Log[mm1/mmu]^2 + mm2 * Log[mm2/mmu]^2 )
                                +1/2 * (mm1 + mm2) * Log[mm1/mm2]^2 + (mm1-mm2)*Hmine[mm1,mm2] )
};

sqcd = t2l //. loopFunctions;

sqcdexp = Normal[Series[sqcd /. {
      mmsb1 -> mmsusy + x*dsb1, mmsb2 -> mmsusy + x*dsb2,
      mmst1 -> mmsusy + x*dst1, mmst2 -> mmsusy + x*dst2,
      mmgl -> mmsusy + x*dgl
      }, {x, 0, 0}]];

sqcdexp1 = Collect[sqcdexp /. zt2 -> Zeta[2], {Log[__]}, Together];

MakePoint[MS_, Xtt_, Q_] :=
    {
        g3 -> Sqrt[4 Pi 0.1184],
        mmsusy -> MS^2,
        mt -> 173.34,
        Xt -> Xtt MS,
        mmu -> Q^2,
        mmt -> mt^2,
        mgl -> MS,
        log2 -> Log[2],
        zt2 -> Zeta[2],
        zt3 -> Zeta[3]
    }

MakePoint[MG_, Mst1_, Mst2_, MS_, xt_, Q_] :=
    {
        g3 -> Sqrt[4 Pi 0.1184],
        mmsusy -> MS^2,
        mt -> 173.34,
        mmt -> mt^2,
        Xt -> xt,
        mmu -> Q^2,
        mgl -> MG, mmgl -> mgl^2,
        mst1 -> Mst1, mst2 -> Mst2,
        mmst1 -> mst1^2, mmst2 -> mst2^2,
        log2 -> Log[2],
        zt2 -> Zeta[2],
        zt3 -> Zeta[3]
    }

Print["# universal SUSY masses: "];
Print[N[sqcdexp1 / (4 Pi)^4 //. #]]& /@ {
    MakePoint[1000, 0, 1000],
    MakePoint[1000, 0, 1100],
    MakePoint[1000, Sqrt[6], 1000],
    MakePoint[1000, Sqrt[6], 1100],
    MakePoint[1000, -Sqrt[6], 1000],
    MakePoint[1000, -Sqrt[6], 1100]
};

Print[];
Print["# non-universal SUSY masses: "];
Print[N[{ t1lqcd / (4 Pi)^2, (t1l - t1lqcd) / (4 Pi)^2, t2lqcd / (4 Pi)^4, (t2l - t2lqcd) / (4 Pi)^4 } //. loopFunctions //. #]]& /@ {
    MakePoint[1000, 1200, 1300, 1400,               0, 1000],
    MakePoint[1000, 1200, 1300, 1400,               0, 1100],
    MakePoint[1000, 1200, 1300, 1400,  Sqrt[6] * 1300, 1000],
    MakePoint[1000, 1200, 1300, 1400,  Sqrt[6] * 1300, 1100],
    MakePoint[1000, 1200, 1300, 1400, -Sqrt[6] * 1300, 1000],
    MakePoint[1000, 1200, 1300, 1400, -Sqrt[6] * 1300, 1100]
};

Print[];
PrintTestSummary[];
