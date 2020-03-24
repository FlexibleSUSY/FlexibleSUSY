// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "config.h"
#include "stopwatch.hpp"

namespace flexiblesusy {

/// returns run-time of a function in seconds
template <class F>
double time_in_seconds(F&& f)
{
   flexiblesusy::Stopwatch sw;
   sw.start();
   f();
   sw.stop();
   return sw.get_time_in_seconds();
}

/// returns run-time of a function, applied to each element in the container, in seconds
template <class F, class Container>
double time_in_seconds(F&&f, Container&& container)
{
   return time_in_seconds([&f, &container] {
      for (const auto& v : container) {
         volatile auto _ = f(v);
      }
   });
}

} // namespace flexiblesusy

#ifdef ENABLE_RANDOM

#include <random>

namespace flexiblesusy {

/// generates vector of uniformly distributed random values
template <class T>
std::vector<T> generate_random_data(unsigned n, T start, T stop)
{
   std::minstd_rand gen;
   std::uniform_real_distribution<T> dist(start, stop);

   std::vector<T> v(n);
   std::generate(begin(v), end(v),
                 [&dist,&gen](){ return dist(gen); });

   return v;
}

} // namespace flexiblesusy

#endif

#endif
