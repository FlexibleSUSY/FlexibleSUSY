#include <algorithm>
#include <math.h>
#include <iostream>

template <typename T>
typename std::enable_if<!std::is_unsigned<T>::value, bool>::type
is_zero(T a, T prec = std::numeric_limits<T>::epsilon()) noexcept
{
   return std::abs(a) <= prec;
}

void sort3(double& s1, double& s12, double& s2, double& m0, double& m1, double& m2)
{
    if (s1 > s12)
    {
        std::swap(s1, s12);
        std::swap(m0, m2);
    }
    if (s12 > s2)
    {
        std::swap(s12, s2);
        std::swap(m0, m1);
    }
    if (s1 > s12)
    {
        std::swap(s1, s12);
        std::swap(m0, m2);
    }
}

double ukC0(double sm1, double sm2, double sm3, double sm4, double sm5, double sm6) {
   sort3(sm1, sm2, sm3, sm4, sm5, sm6);
   if( is_zero(sm1) && (sm2>0) && (sm3>0) && (!is_zero(sm2-sm3)) ) {
      if(is_zero(sm4-sm5)) {
         return -((sm4 - sm6)*(36*std::pow(sm4 - sm6,4) + 9*sm2*std::pow(sm4 - sm6,2)*(sm4 + 5*sm6) +
         9*sm3*std::pow(sm4 - sm6,2)*(sm4 + 5*sm6) +
         4*sm2*sm3*(std::pow(sm4,2) + 19*sm4*sm6 + 10*std::pow(sm6,2))) -
         6*sm6*(6*std::pow(sm4 - sm6,4) + 3*sm2*std::pow(sm4 - sm6,2)*(2*sm4 + sm6) +
         3*sm3*std::pow(sm4 - sm6,2)*(2*sm4 + sm6) +
         2*sm2*sm3*(3*std::pow(sm4,2) + 6*sm4*sm6 + std::pow(sm6,2)))*std::log(sm4/sm6))/
         (36.*std::pow(sm4 - sm6,6));
      }
      if(is_zero(sm4-sm6)) {
         return (-((sm4 - sm5)*(36*sm4*std::pow(sm4 - sm5,4) +
          9*sm2*sm4*std::pow(sm4 - sm5,2)*(sm4 + 5*sm5) +
          3*sm3*std::pow(sm4 - sm5,2)*(std::pow(sm4,2) - 5*sm4*sm5 - 2*std::pow(sm5,2)) +
          sm2*sm3*(std::pow(sm4,3) - 11*std::pow(sm4,2)*sm5 - 47*sm4*std::pow(sm5,2) -
          3*std::pow(sm5,3)))) +
          6*sm4*sm5*(3*std::pow(sm4 - sm5,2)*(2*std::pow(sm4 - sm5,2) - sm3*sm5) +
          sm2*(6*std::pow(sm4,3) - 3*sm4*(2*sm3 + 3*sm4)*sm5 - 4*sm3*std::pow(sm5,2) +
          3*std::pow(sm5,3)))*std::log(sm4/sm5))/(36.*sm4*std::pow(sm4 - sm5,6));
      }
      if(is_zero(sm5-sm6)) {
         return (9*std::pow(sm4 - sm5,3)*sm5*(4*std::pow(sm4 - sm5,2) + sm3*(5*sm4 + sm5)) -
          sm2*(sm4 - sm5)*(3*std::pow(sm4 - sm5,2)*
          (2*std::pow(sm4,2) + 5*sm4*sm5 - std::pow(sm5,2)) +
          sm3*(3*std::pow(sm4,3) + 47*std::pow(sm4,2)*sm5 + 11*sm4*std::pow(sm5,2) -
          std::pow(sm5,3))) + 6*sm4*sm5*
          (3*sm2*sm4*std::pow(sm4 - sm5,2) - 6*std::pow(sm4 - sm5,4) -
          3*sm3*std::pow(sm4 - sm5,2)*(sm4 + 2*sm5) + sm2*sm3*sm4*(4*sm4 + 6*sm5))*
          std::log(sm4/sm5))/(36.*std::pow(sm4 - sm5,6)*sm5);
      }
      if(is_zero(sm4-sm5) && is_zero(sm5-sm6)) {
         return -(2*sm2*sm3 + 15*(sm2 + sm3)*sm4 + 180*std::pow(sm4,2))/(360.*std::pow(sm4,3));
      }
      return -((sm4*std::log(sm4/sm6))/((sm4 - sm5)*(sm4 - sm6))) +
        sm3*((sm4*sm5 - 2*sm4*sm6 + sm5*sm6)/
        (2.*std::pow(sm4 - sm6,2)*(sm4*sm5 - std::pow(sm5,2) - sm4*sm6 + sm5*sm6)) -
        (std::pow(sm5,2)*std::log(sm4/sm5))/(2.*std::pow(-sm4 + sm5,2)*std::pow(sm5 - sm6,2)) -
        (sm6*(-2*sm4*sm5 + sm4*sm6 + std::pow(sm6,2))*std::log(sm4/sm6))/
        (2.*std::pow(sm4 - sm6,3)*std::pow(-sm5 + sm6,2))) -
        (sm5*std::log(sm5/sm6))/((-sm4 + sm5)*(sm5 - sm6)) +
        sm2*(-(sm4*sm5 + sm4*sm6 - 2*sm5*sm6)/
        (2.*std::pow(sm5 - sm6,2)*(std::pow(sm4,2) - sm4*sm5 - sm4*sm6 + sm5*sm6)) +
        (std::pow(sm4,2)*std::log(sm4/sm6))/(2.*std::pow(sm4 - sm5,2)*std::pow(sm4 - sm6,2)) -
        (sm5*(std::pow(sm5,2) - 2*sm4*sm6 + sm5*sm6)*std::log(sm5/sm6))/
        (2.*std::pow(-sm4 + sm5,2)*std::pow(sm5 - sm6,3)) +
        sm3*(-(2*std::pow(sm4,3)*std::pow(sm5,2) + 2*std::pow(sm4,2)*std::pow(sm5,3) +
        5*std::pow(sm4,3)*sm5*sm6 - 22*std::pow(sm4,2)*std::pow(sm5,2)*sm6 +
        5*sm4*std::pow(sm5,3)*sm6 - std::pow(sm4,3)*std::pow(sm6,2) +
        7*std::pow(sm4,2)*sm5*std::pow(sm6,2) + 7*sm4*std::pow(sm5,2)*std::pow(sm6,2) -
        std::pow(sm5,3)*std::pow(sm6,2) - 5*std::pow(sm4,2)*std::pow(sm6,3) +
        6*sm4*sm5*std::pow(sm6,3) - 5*std::pow(sm5,2)*std::pow(sm6,3))/
        (6.*std::pow(sm4 - sm5,2)*std::pow(sm4 - sm6,3)*std::pow(sm5 - sm6,3)) +
        (std::pow(sm4,2)*(2*std::pow(sm4,2) + sm4*sm6 - 3*sm5*sm6)*std::log(sm4/sm6))/
        (3.*std::pow(sm4 - sm5,3)*std::pow(sm4 - sm6,4)) +
        (std::pow(sm5,2)*(2*std::pow(sm5,2) - 3*sm4*sm6 + sm5*sm6)*std::log(sm5/sm6))/
        (3.*std::pow(-sm4 + sm5,3)*std::pow(sm5 - sm6,4))));
   }
   if( (sm1<0.0) && (sm2>0) && (sm3>0) && (is_zero(sm1+sm3)) && (!is_zero(sm2-sm3)) ) {
      if(is_zero(sm4-sm5)) {
         return -(2*std::pow(sm4 - sm6,3)*(6*sm4*std::pow(sm4 - sm6,2) +
            sm3*(std::pow(sm4,2) + 10*sm4*sm6 + std::pow(sm6,2))) +
            sm2*(sm4 - sm6)*(3*sm4*std::pow(sm4 - sm6,2)*(sm4 + 5*sm6) +
            sm3*(sm4 + sm6)*(std::pow(sm4,2) + 28*sm4*sm6 + std::pow(sm6,2))) -
            6*sm4*sm6*(2*std::pow(sm4 - sm6,4) + 2*sm3*std::pow(sm4 - sm6,2)*(sm4 + sm6) +
            sm2*std::pow(sm4 - sm6,2)*(2*sm4 + sm6) +
            2*sm2*sm3*(std::pow(sm4,2) + 3*sm4*sm6 + std::pow(sm6,2)))*std::log(sm4/sm6))/
            (12.*sm4*std::pow(sm4 - sm6,6));
      }
      if(is_zero(sm4-sm6)) {
         return (2*std::pow(sm4 - sm5,3)*(-6*sm4*std::pow(sm4 - sm5,2) +
            sm3*(std::pow(sm4,2) + 10*sm4*sm5 + std::pow(sm5,2))) +
            sm2*(sm4 - sm5)*(-3*sm4*std::pow(sm4 - sm5,2)*(sm4 + 5*sm5) +
            sm3*(sm4 + sm5)*(std::pow(sm4,2) + 28*sm4*sm5 + std::pow(sm5,2))) +
            6*sm4*sm5*(2*std::pow(sm4 - sm5,4) - 2*sm3*std::pow(sm4 - sm5,2)*(sm4 + sm5) +
            sm2*std::pow(sm4 - sm5,2)*(2*sm4 + sm5) -
            2*sm2*sm3*(std::pow(sm4,2) + 3*sm4*sm5 + std::pow(sm5,2)))*std::log(sm4/sm5))/
            (12.*sm4*std::pow(sm4 - sm5,6));
      }
      if(is_zero(sm5-sm6)) {
         return ((sm4 - sm5)*(12*std::pow(sm4 - sm5,2)*sm5 +
            sm2*(-2*std::pow(sm4,2) - 5*sm4*sm5 + std::pow(sm5,2))) +
            6*sm4*(sm2*sm4 - 2*std::pow(sm4 - sm5,2))*sm5*std::log(sm4/sm5))/
            (12.*std::pow(sm4 - sm5,4)*sm5);
      }
      if(is_zero(sm4-sm5) && is_zero(sm5-sm6)) {
         return -(sm2 + 12*sm4)/(24.*std::pow(sm4,2));
      }
      return -((sm4*std::log(sm4/sm6))/((sm4 - sm5)*(sm4 - sm6))) +
        sm3*(-(std::pow(sm4,2)*sm5 + sm4*std::pow(sm5,2) + std::pow(sm4,2)*sm6 -
        6*sm4*sm5*sm6 + std::pow(sm5,2)*sm6 + sm4*std::pow(sm6,2) + sm5*std::pow(sm6,2))/
        (2.*std::pow(sm4 - sm5,2)*std::pow(sm4 - sm6,2)*(sm5 - sm6)) -
        (sm5*(std::pow(sm5,2) - sm4*sm6)*std::log(sm4/sm5))/
        (std::pow(-sm4 + sm5,3)*std::pow(sm5 - sm6,2)) +
        (sm6*(sm4*sm5 - std::pow(sm6,2))*std::log(sm4/sm6))/
        (std::pow(sm4 - sm6,3)*std::pow(-sm5 + sm6,2))) -
        (sm5*std::log(sm5/sm6))/((-sm4 + sm5)*(sm5 - sm6)) +
        sm2*(-(sm4*sm5 + sm4*sm6 - 2*sm5*sm6)/
        (2.*std::pow(sm5 - sm6,2)*(std::pow(sm4,2) - sm4*sm5 - sm4*sm6 + sm5*sm6)) +
        ((2*std::pow(sm4,2)*sm5 - 4*sm4*std::pow(sm5,2) + 3*std::pow(sm5,3) -
        2*std::pow(sm4,2)*sm6 + 2*sm4*sm5*sm6 - std::pow(sm5,2)*sm6)*std::log(sm4/sm5))/
        (2.*std::pow(sm4 - sm5,2)*std::pow(sm5 - sm6,3)) +
        ((2*std::pow(sm4,2)*sm5 - 2*std::pow(sm4,2)*sm6 - 6*sm4*sm5*sm6 +
        4*sm4*std::pow(sm6,2) + 3*sm5*std::pow(sm6,2) - std::pow(sm6,3))*std::log(sm4/sm6))/
        (2.*std::pow(sm4 - sm6,2)*std::pow(-sm5 + sm6,3)) +
        std::log(sm5/sm6)/std::pow(sm5 - sm6,2) +
        sm3*(-(std::pow(sm4,4)*std::pow(sm5,2) - 5*std::pow(sm4,3)*std::pow(sm5,3) -
        2*std::pow(sm4,2)*std::pow(sm5,4) + 10*std::pow(sm4,4)*sm5*sm6 -
        19*std::pow(sm4,3)*std::pow(sm5,2)*sm6 +
        38*std::pow(sm4,2)*std::pow(sm5,3)*sm6 - 5*sm4*std::pow(sm5,4)*sm6 +
        std::pow(sm4,4)*std::pow(sm6,2) - 19*std::pow(sm4,3)*sm5*std::pow(sm6,2) -
        19*sm4*std::pow(sm5,3)*std::pow(sm6,2) + std::pow(sm5,4)*std::pow(sm6,2) -
        5*std::pow(sm4,3)*std::pow(sm6,3) + 38*std::pow(sm4,2)*sm5*std::pow(sm6,3) -
        19*sm4*std::pow(sm5,2)*std::pow(sm6,3) + 10*std::pow(sm5,3)*std::pow(sm6,3) -
        2*std::pow(sm4,2)*std::pow(sm6,4) - 5*sm4*sm5*std::pow(sm6,4) +
        std::pow(sm5,2)*std::pow(sm6,4))/
        (6.*std::pow(sm5 - sm6,3)*
        std::pow(std::pow(sm4,2) - sm4*sm5 - sm4*sm6 + sm5*sm6,3)) -
        (std::pow(sm4,2)*(sm5 - sm6)*(std::pow(sm4,2) - sm5*sm6)*std::log(sm4/sm6))/
        (std::pow(sm4 - sm5,4)*std::pow(sm4 - sm6,4)) +
        (sm5*(std::pow(sm5,4) + std::pow(sm4,2)*sm5*sm6 - 4*sm4*std::pow(sm5,2)*sm6 +
        std::pow(sm5,3)*sm6 + std::pow(sm4,2)*std::pow(sm6,2))*std::log(sm5/sm6))/
        (std::pow(-sm4 + sm5,4)*std::pow(sm5 - sm6,4))));
   }
   std::cout << "error\n";
   return 0.0;
}
