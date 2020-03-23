
#include "gm2calc/MSSMNoFV_onshell.hpp"
#include "gm2_slha_io.hpp"
#include "error.hpp"

#include <sstream>

char const * const slha_input =
"Block GM2CalcInput\n"
"     1     0.00781653       # alpha(MZ)\n"
"     2     0.00729735       # alpha in Thompson limit\n"
"Block SMINPUTS\n"
"     3     1.17600000E-01   # alpha_s(MZ) SM MSbar\n"
"     4     9.11876000E+01   # MZ(pole)\n"
"     5     4.20000000E+00   # mb(mb) SM MSbar\n"
"     6     1.73300000E+02   # mtop(pole)\n"
"     7     1.77700000E+00   # mtau(pole)\n"
"     8     0.00000000E+00   # mnu3(pole)\n"
"    11     5.02667588E-04   # melectron(pole)\n"
"    12     0.00000000E+00   # mnu1(pole)\n"
"    13     0.1039357        # mmuon(pole)\n"
"    14     0.00000000E+00   # mnu2(pole)\n"
"    21     4.76052706E-03   # md\n"
"    22     2.40534062E-03   # mu\n"
"    23     1.04230487E-01   # ms\n"
"    24     1.27183378E+00   # mc\n"
"Block MASS\n"
"        24     8.04040000E+01   # VWm\n"
"   1000012     3.50280760E+02   # SveL\n"
"   1000014     3.50272662E+02   # SvmL\n"
"   1000016     3.49143898E+02   # SvtL\n"
"   1000024     3.87129550E+02   # Cha(1)\n"
"   1000037     6.38851559E+02   # Cha(2)\n"
"        25     1.14752354E+02   # hh(1)\n"
"        35     7.07297611E+02   # hh(2)\n"
"        37     7.11829537E+02   # Hpm(2)\n"
"        36     7.07024952E+02   # Ah(2)\n"
"   1000001     9.99590838E+02   # Sd(1)\n"
"   2000001     1.04394524E+03   # Sd(2)\n"
"   1000003     9.99587075E+02   # Ss(1)\n"
"   2000003     1.04394342E+03   # Ss(2)\n"
"   1000005     9.56519959E+02   # Sb(1)\n"
"   2000005     9.96560664E+02   # Sb(2)\n"
"   1000011     2.29105585E+02   # Se(1)\n"
"   2000011     3.59250950E+02   # Se(2)\n"
"   1000013     2.29055791E+02   # Sm(1)\n"
"   2000013     3.59259629E+02   # Sm(2)\n"
"   1000015     2.22125984E+02   # Stau(1)\n"
"   2000015     3.60399301E+02   # Stau(2)\n"
"   1000002     1.00291908E+03   # Su(1)\n"
"   2000002     1.04107273E+03   # Su(2)\n"
"   1000004     1.00291298E+03   # Sc(1)\n"
"   2000004     1.04107343E+03   # Sc(2)\n"
"   1000006     7.96708530E+02   # St(1)\n"
"   2000006     1.00222044E+03   # St(2)\n"
"   1000022     2.09850638E+02   # Chi(1)\n"
"   1000023     3.87142196E+02   # Chi(2)\n"
"   1000025    -6.24072971E+02   # Chi(3)\n"
"   1000035     6.38603776E+02   # Chi(4)\n"
"Block HMIX Q= 8.66360379E+02\n"
"     1     6.18568499E+02   # Mu\n"
"     2     10               # vu/vd\n"
"     4     5.18944643E+05   # Sqr(MAh(1))\n"
"Block Au Q= 8.66360379E+02\n"
"  1  1    -1.12397387E+03   # TYu(0,0)/Yu(0,0)\n"
"  2  2    -1.12396894E+03   # TYu(1,1)/Yu(1,1)\n"
"  3  3    -8.70714986E+02   # TYu(2,2)/Yu(2,2)\n"
"Block Ad Q= 8.66360379E+02\n"
"  1  1    -1.37205166E+03   # TYd(0,0)/Yd(0,0)\n"
"  2  2    -1.37204709E+03   # TYd(1,1)/Yd(1,1)\n"
"  3  3    -1.28330100E+03   # TYd(2,2)/Yd(2,2)\n"
"Block Ae Q= 8.66360379E+02\n"
"  1  1    -2.93735295E+02   # TYe(0,0)/Ye(0,0)\n"
"  2  2    -2.93720212E+02   # TYe(1,1)/Ye(1,1)\n"
"  3  3    -2.92154796E+02   # TYe(2,2)/Ye(2,2)\n"
"Block MSOFT Q= 8.66360379E+02\n"
"     1     2.15154384E+02   # MassB\n"
"     2     3.90626858E+02   # MassWB\n"
"     3     1.10300877E+03   # MassG\n"
"    31     3.51653258E+02   # SignedAbsSqrt(ml2(0,0))\n"
"    32     3.51646284E+02   # SignedAbsSqrt(ml2(1,1))\n"
"    33     3.50674223E+02   # SignedAbsSqrt(ml2(2,2))\n"
"    34     2.21215037E+02   # SignedAbsSqrt(me2(0,0))\n"
"    35     2.21192429E+02   # SignedAbsSqrt(me2(1,1))\n"
"    36     2.18022142E+02   # SignedAbsSqrt(me2(2,2))\n"
"    41     1.00711403E+03   # SignedAbsSqrt(mq2(0,0))\n"
"    42     1.00711149E+03   # SignedAbsSqrt(mq2(1,1))\n"
"    43     9.29083096E+02   # SignedAbsSqrt(mq2(2,2))\n"
"    44     9.69369660E+02   # SignedAbsSqrt(mu2(0,0))\n"
"    45     9.69366965E+02   # SignedAbsSqrt(mu2(1,1))\n"
"    46     7.99712943E+02   # SignedAbsSqrt(mu2(2,2))\n"
"    47     9.64756473E+02   # SignedAbsSqrt(md2(0,0))\n"
"    48     9.64753818E+02   # SignedAbsSqrt(md2(1,1))\n"
"    49     9.60016201E+02   # SignedAbsSqrt(md2(2,2))\n"
;

using namespace flexiblesusy;

gm2calc::MSSMNoFV_onshell read(char const * const slha_input)
{
   gm2calc::MSSMNoFV_onshell osmodel;
   gm2calc::GM2_slha_io slha_io;
   std::istringstream sstr(slha_input);

   try {
      slha_io.read_from_stream(sstr);
      slha_io.fill_slha(osmodel);
   } catch (const Error& error) {
      BOOST_FAIL(error.what());
   }

   return osmodel;
}

gm2calc::MSSMNoFV_onshell setup()
{
   gm2calc::MSSMNoFV_onshell osmodel(read(slha_input));

   BOOST_REQUIRE_NO_THROW(osmodel.convert_to_onshell());

   return osmodel;
}
