#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_tsil

#include <boost/test/unit_test.hpp>
#include "tsil_cpp.h"
#include <cmath>

BOOST_AUTO_TEST_CASE(test_individual)
{
   TSIL_REAL x = 1;
   TSIL_REAL y = 2;
   TSIL_REAL z = 3;
   TSIL_REAL u = 4;
   TSIL_REAL v = 5;
   TSIL_COMPLEXCPP s = 10;
   TSIL_REAL qq = 1;

   const auto A    = TSIL_A_(x, qq);
   const auto Ap   = TSIL_Ap_(x, qq);
   const auto Aeps = TSIL_Aeps_(x, qq);
   const auto B    = TSIL_B_(x, y, s, qq);
   const auto Bp   = TSIL_Bp_(x, y, s, qq);
   const auto dBds = TSIL_dBds_(x, y, s, qq);
   const auto Beps = TSIL_Beps_(x, y, s, qq);
   const auto I2   = TSIL_I2_(x, y, z, qq);
   const auto I2p  = TSIL_I2p_(x, y, z, qq);
   const auto I2p2 = TSIL_I2p2_(x, y, z, qq);
   const auto I2pp = TSIL_I2pp_(x, y, z, qq);
   const auto I2p3 = TSIL_I2p3_(x, y, z, qq);

   TSIL_COMPLEXCPP M;
   int is_analytic = TSIL_Manalytic_(x, x, y, y, 0, s, &M);
   BOOST_CHECK(is_analytic != 0);

   TSIL_COMPLEXCPP S;
   is_analytic = TSIL_Sanalytic_(x, y, y, x, qq, &S);
   BOOST_CHECK(is_analytic != 0);

   TSIL_COMPLEXCPP T;
   is_analytic = TSIL_Tanalytic_(x, 0, y, s, qq, &T);
   BOOST_CHECK(is_analytic != 0);

   TSIL_COMPLEXCPP Tbar;
   is_analytic = TSIL_Tbaranalytic_(0, x, y, s, qq, &Tbar);
   BOOST_CHECK(is_analytic != 0);

   TSIL_COMPLEXCPP U;
   is_analytic = TSIL_Uanalytic_(x, 0, 0, y, s, qq, &U);
   BOOST_CHECK(is_analytic != 0);

   TSIL_COMPLEXCPP V;
   is_analytic = TSIL_Vanalytic_(x, y, 0, y, s, qq, &V);
   BOOST_CHECK(is_analytic != 0);

   BOOST_CHECK(std::isfinite(std::real(A   )));
   BOOST_CHECK(std::isfinite(std::real(Ap  )));
   BOOST_CHECK(std::isfinite(std::real(Aeps)));
   BOOST_CHECK(std::isfinite(std::real(B   )));
   BOOST_CHECK(std::isfinite(std::real(Bp  )));
   BOOST_CHECK(std::isfinite(std::real(dBds)));
   BOOST_CHECK(std::isfinite(std::real(Beps)));
   BOOST_CHECK(std::isfinite(std::real(I2  )));
   BOOST_CHECK(std::isfinite(std::real(I2p )));
   BOOST_CHECK(std::isfinite(std::real(I2p2)));
   BOOST_CHECK(std::isfinite(std::real(I2pp)));
   BOOST_CHECK(std::isfinite(std::real(I2p3)));
   BOOST_CHECK(std::isfinite(std::real(M   )));
   BOOST_CHECK(std::isfinite(std::real(S   )));
   BOOST_CHECK(std::isfinite(std::real(T   )));
   BOOST_CHECK(std::isfinite(std::real(Tbar)));
   BOOST_CHECK(std::isfinite(std::real(U   )));
   BOOST_CHECK(std::isfinite(std::real(V   )));

   BOOST_CHECK(std::isfinite(std::imag(A   )));
   BOOST_CHECK(std::isfinite(std::imag(Ap  )));
   BOOST_CHECK(std::isfinite(std::imag(Aeps)));
   BOOST_CHECK(std::isfinite(std::imag(B   )));
   BOOST_CHECK(std::isfinite(std::imag(Bp  )));
   BOOST_CHECK(std::isfinite(std::imag(dBds)));
   BOOST_CHECK(std::isfinite(std::imag(Beps)));
   BOOST_CHECK(std::isfinite(std::imag(I2  )));
   BOOST_CHECK(std::isfinite(std::imag(I2p )));
   BOOST_CHECK(std::isfinite(std::imag(I2p2)));
   BOOST_CHECK(std::isfinite(std::imag(I2pp)));
   BOOST_CHECK(std::isfinite(std::imag(I2p3)));
   BOOST_CHECK(std::isfinite(std::imag(M   )));
   BOOST_CHECK(std::isfinite(std::imag(S   )));
   BOOST_CHECK(std::isfinite(std::imag(T   )));
   BOOST_CHECK(std::isfinite(std::imag(Tbar)));
   BOOST_CHECK(std::isfinite(std::imag(U   )));
   BOOST_CHECK(std::isfinite(std::imag(V   )));

   BOOST_TEST_MESSAGE("A    = " << A   );
   BOOST_TEST_MESSAGE("Ap   = " << Ap  );
   BOOST_TEST_MESSAGE("Aeps = " << Aeps);
   BOOST_TEST_MESSAGE("B    = " << B   );
   BOOST_TEST_MESSAGE("Bp   = " << Bp  );
   BOOST_TEST_MESSAGE("dBds = " << dBds);
   BOOST_TEST_MESSAGE("Beps = " << Beps);
   BOOST_TEST_MESSAGE("I2   = " << I2  );
   BOOST_TEST_MESSAGE("I2p  = " << I2p );
   BOOST_TEST_MESSAGE("I2p2 = " << I2p2);
   BOOST_TEST_MESSAGE("I2pp = " << I2pp);
   BOOST_TEST_MESSAGE("I2p3 = " << I2p3);
   BOOST_TEST_MESSAGE("M    = " << M   );
   BOOST_TEST_MESSAGE("S    = " << S   );
   BOOST_TEST_MESSAGE("T    = " << T   );
   BOOST_TEST_MESSAGE("Tbar = " << Tbar);
   BOOST_TEST_MESSAGE("U    = " << U   );
   BOOST_TEST_MESSAGE("V    = " << V   );
}

BOOST_AUTO_TEST_CASE(test_numeric)
{
   TSIL_REAL x = 1;
   TSIL_REAL y = 2;
   TSIL_REAL z = 3;
   TSIL_REAL u = 4;
   TSIL_REAL v = 5;
   TSIL_COMPLEXCPP s = 6;
   TSIL_REAL qq = 7;

   TSIL_DATA data{};
   TSIL_SetParameters_(&data, x, y, z, u, v, qq);
   TSIL_Evaluate_(&data, std::real(s));

   const auto Bxz   = TSIL_GetFunction_(&data, "Bxz"  );
   const auto Byu   = TSIL_GetFunction_(&data, "Byu"  );
   const auto Svyz  = TSIL_GetFunction_(&data, "Svyz" );
   const auto Suxv  = TSIL_GetFunction_(&data, "Suxv" );
   const auto Sxuv  = TSIL_GetFunction_(&data, "Sxuv" );
   const auto Suvx  = TSIL_GetFunction_(&data, "Suvx" );
   const auto Svux  = TSIL_GetFunction_(&data, "Svux" );
   const auto Sxvu  = TSIL_GetFunction_(&data, "Sxvu" );
   const auto Svxu  = TSIL_GetFunction_(&data, "Svxu" );
   const auto Tuxv  = TSIL_GetFunction_(&data, "Tuxv" );
   const auto Tyzv  = TSIL_GetFunction_(&data, "Tyzv" );
   const auto Txuv  = TSIL_GetFunction_(&data, "Txuv" );
   const auto Tzyv  = TSIL_GetFunction_(&data, "Tzyv" );
   const auto Tvxu  = TSIL_GetFunction_(&data, "Tvxu" );
   const auto Tvyz  = TSIL_GetFunction_(&data, "Tvyz" );
   const auto Uuyxv = TSIL_GetFunction_(&data, "Uuyxv");
   const auto Uzxyv = TSIL_GetFunction_(&data, "Uzxyv");
   const auto Uzxvy = TSIL_GetFunction_(&data, "Uzxvy");
   const auto Uxzuv = TSIL_GetFunction_(&data, "Uxzuv");
   const auto Uyuzv = TSIL_GetFunction_(&data, "Uyuzv");
   const auto Vzxyv = TSIL_GetFunction_(&data, "Vzxyv");
   const auto Vuyxv = TSIL_GetFunction_(&data, "Vuyxv");
   const auto Vxzuv = TSIL_GetFunction_(&data, "Vxzuv");
   const auto Vyuzv = TSIL_GetFunction_(&data, "Vyuzv");
   const auto M     = TSIL_GetFunction_(&data, "M"    );

   BOOST_CHECK_CLOSE_FRACTION(std::real(Bxz  ),   2.1828815497944780, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Byu  ),   1.3097038687359490, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Svyz ), -31.2543090480954149, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Suxv ), -31.5461792884989656, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Sxuv ), -31.5461792884989656, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Suvx ), -31.5461792884989656, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Svux ), -31.5461792884989656, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Sxvu ), -31.5461792884989656, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Svxu ), -31.5461792884989656, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Tuxv ),   0.8763732955387802, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Tyzv ),   0.5123629327539176, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Txuv ),   0.4662594413917679, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Tzyv ),   0.6656794217040922, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Tvxu ),   1.0054064846763813, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Tvyz ),   0.9432272340849203, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Uuyxv),   0.9863678693113318, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Uzxyv),   1.0286440457448356, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Uzxvy),   1.0286440457448356, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Uxzuv),   0.2725852029409196, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Uyuzv),   0.2952750228795177, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Vzxyv),   0.5589744541751330, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Vuyxv),   0.1762338999120509, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Vxzuv),   0.1985617199983827, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(Vyuzv),   0.0588643927348281, 1e-14);
   BOOST_CHECK_CLOSE_FRACTION(std::real(M    ),   0.4877200674519486, 1e-14);
}
