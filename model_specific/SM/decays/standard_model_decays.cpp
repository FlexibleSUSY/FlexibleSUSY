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


/**
 * @file standard_model_decays.cpp
 * @brief implementation of particle decays in the SM
 *
 * Contains the definition of SM decays class methods,
 * which calculate the decay rates for particles in the model.
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#include <cmath>
#include <iomanip>

#include <boost/range/algorithm.hpp>
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_LIMIT_LIST_SIZE 50
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/list.hpp>
#include <boost/hana/size.hpp>
#include <gsl/gsl_monte_miser.h>

#include "standard_model_decays.hpp"
#include "standard_model.hpp"
#include "decays/one_loop_decay_diagrams.hpp"
#include "concatenate.hpp"
#include "decays/decay_functions.hpp"
#include "config.h"
#include "thread_pool.hpp"
#include "wrappers.hpp"
#include "cxx_qft/standard_model_particle_aliases.hpp"
#include "cxx_qft/standard_model_vertices.hpp"
#include "error.hpp"
#include "cxx_qft/standard_model_qft.hpp"

namespace flexiblesusy {

#define CLASSNAME Standard_model_decays
#define PHYSICAL(parameter) model.get_physical().parameter

using namespace standard_model_cxx_diagrams;
using namespace standard_model_cxx_diagrams::fields;
namespace info = standard_model_info;
using namespace std::complex_literals;

const standard_model_decay_table& CLASSNAME::get_decay_table() const
{
   return decay_table;
}

const FlexibleDecay_problems& CLASSNAME::get_problems() const
{
   return problems;
}

void CLASSNAME::clear_problems()
{
   problems.clear();
}

void CLASSNAME::clear()
{
   clear_problems();
   decay_table.clear();
}

void CLASSNAME::calculate_decays()
{
try {
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 1u));

   tp.run_task([this] () { calculate_hh_decays(); });

#else
   calculate_hh_decays();

#endif
}
catch (std::exception& e) {
   problems.add_error(e.what());
}
}

double CLASSNAME::get_alphas(standard_model_cxx_diagrams::context_base const& context) const
{
   return Sqr(context.model.get_g3())/(4.*Pi);
}

double CLASSNAME::get_alpha(standard_model_cxx_diagrams::context_base const& context) const
{
   return Sqr(unit_charge(context))/(4.*Pi);
}

template <int N>
constexpr auto _to_array = [](auto&& ...x) {
   return std::array<int, N>{std::forward<decltype(x)>(x)...};
};

template <typename Field>
int fieldPDG(const typename field_indices<Field>::type indx) {
   return boost::hana::unpack(Field::pdgids, _to_array<boost::hana::size(Field::pdgids)>).at(indx.size() > 0 ? indx.at(0) : 0);
}

/* 1-loop BSM amplitudes
 *
 * Notes:
 *    - Since FS removes a factor of i from every vertex we reintroduce it
 *      calls to one loop amplitudes
 */

// hh -> {VG, VG}
template<>
Decay_amplitude_SVV CLASSNAME::calculate_amplitude_1l<hh, VG, VG>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh >::type const& idx_1,
   typename cxx_diagrams::field_indices<fields::VG >::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::VG >::type const& idx_3) const{

   // amplitude type
   Decay_amplitude_SVV result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_vector_1 = context.physical_mass<VG>(idx_2);
   result.m_vector_2 = context.physical_mass<VG>(idx_3);

   // set the initial value of an amplitude to 0
   result.form_factor_g = std::complex<double>(0., 0.);
   result.form_factor_11 = std::complex<double>(0., 0.);
   result.form_factor_12 = std::complex<double>(0., 0.);
   result.form_factor_21 = std::complex<double>(0., 0.);
   result.form_factor_22 = std::complex<double>(0., 0.);
   result.form_factor_eps = std::complex<double>(0., 0.);

   const double ren_scale {result.m_decay};

   // ----------------- 1-loop contributions to the amplitude -----------------

   // topology T1
   // internal particles in the diagram: Fd, bar[Fd], Fd
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, hh>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VG>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, VG>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {2.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fd>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Fu, bar[Fu], Fu
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, hh>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VG>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, VG>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {2.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fu>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fd], Fd, bar[Fd]
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, hh>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VG>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, VG>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {2.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fd>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fu], Fu, bar[Fu]
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, hh>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VG>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, VG>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {2.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fu>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
               }
            }
         }
      }
   }

   return result;
}

// hh -> {VP, VP}
template<>
Decay_amplitude_SVV CLASSNAME::calculate_amplitude_1l<hh, VP, VP>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh >::type const& idx_1,
   typename cxx_diagrams::field_indices<fields::VP >::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::VP >::type const& idx_3) const{

   // amplitude type
   Decay_amplitude_SVV result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_vector_1 = context.physical_mass<VP>(idx_2);
   result.m_vector_2 = context.physical_mass<VP>(idx_3);

   // set the initial value of an amplitude to 0
   result.form_factor_g = std::complex<double>(0., 0.);
   result.form_factor_11 = std::complex<double>(0., 0.);
   result.form_factor_12 = std::complex<double>(0., 0.);
   result.form_factor_21 = std::complex<double>(0., 0.);
   result.form_factor_22 = std::complex<double>(0., 0.);
   result.form_factor_eps = std::complex<double>(0., 0.);

   // FormCalc's Finite variable
   static constexpr double Finite {1.};

   const double ren_scale {result.m_decay};

   // ----------------- 1-loop contributions to the amplitude -----------------

   // topology T1
   // internal particles in the diagram: Fd, bar[Fd], Fd
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, hh>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VP>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fd>(vertexId2::template indices_of_field<1>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<1>(indexId1) == vertexId1::template indices_of_field<0>(indexId1) &&
                     vertexId2::template indices_of_field<1>(indexId2) == vertexId2::template indices_of_field<0>(indexId2) &&
                     vertexId3::template indices_of_field<1>(indexId3) == vertexId3::template indices_of_field<0>(indexId3)) {
                     const auto correction_S = 1. + get_alphas(context)/Pi * delta_hAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     const auto correction_P = 1. + get_alphas(context)/Pi * delta_AhAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                     temp_result.form_factor_eps = correction_P * temp_result.form_factor_eps;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Fe, bar[Fe], Fe
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fe, hh>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, VP>;
      using vertexId3 = Vertex<typename bar<Fe>::type, Fe, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fe>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fe>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Fu, bar[Fu], Fu
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, hh>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VP>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fu>(vertexId2::template indices_of_field<1>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<1>(indexId1) == vertexId1::template indices_of_field<0>(indexId1) &&
                     vertexId2::template indices_of_field<1>(indexId2) == vertexId2::template indices_of_field<0>(indexId2) &&
                     vertexId3::template indices_of_field<1>(indexId3) == vertexId3::template indices_of_field<0>(indexId3)) {
                     const auto correction_S = 1. + get_alphas(context)/Pi * delta_hAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     const auto correction_P = 1. + get_alphas(context)/Pi * delta_AhAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                     temp_result.form_factor_eps = correction_P * temp_result.form_factor_eps;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[gWp], gWp, bar[gWp]
   {
      using vertexId1 = Vertex<typename bar<gWp>::type, gWp, hh>;
      using vertexId2 = Vertex<typename bar<gWp>::type, gWp, VP>;
      using vertexId3 = Vertex<typename bar<gWp>::type, gWp, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<gWp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<gWp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<gWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[gWpC], gWpC, bar[gWpC]
   {
      using vertexId1 = Vertex<typename bar<gWpC>::type, gWpC, hh>;
      using vertexId2 = Vertex<typename bar<gWpC>::type, gWpC, VP>;
      using vertexId3 = Vertex<typename bar<gWpC>::type, gWpC, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<gWpC>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<gWpC>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<gWpC>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[Hp], Hp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<Hp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g2n2_SSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(0, 1), 1.0i*vertexId3Val.
                        value(1, 0),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[Hp], VWp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g4n4_SSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], VWp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g5n5_SVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g7n7_SVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::odd_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g6n6_VSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(0, 1), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], Hp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<1>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g8n8_VSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[VWp], Hp
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g9n9_VVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[VWp], VWp
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g10n10_VVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(cxx_diagrams::
                        TripleVectorVertex::odd_permutation {}), 1.0i*vertexId3Val.value(
                        cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale,
                     Finite);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fd], Fd, bar[Fd]
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, hh>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VP>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fd>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<0>(indexId1) == vertexId1::template indices_of_field<1>(indexId1) &&
                     vertexId2::template indices_of_field<0>(indexId2) == vertexId2::template indices_of_field<1>(indexId2) &&
                     vertexId3::template indices_of_field<0>(indexId3) == vertexId3::template indices_of_field<1>(indexId3)) {
                     const auto correction_S = 1. + get_alphas(context)/Pi * delta_hAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     const auto correction_P = 1. + get_alphas(context)/Pi * delta_AhAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                     temp_result.form_factor_eps = correction_P * temp_result.form_factor_eps;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fe], Fe, bar[Fe]
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fe, hh>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, VP>;
      using vertexId3 = Vertex<typename bar<Fe>::type, Fe, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fe>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fu], Fu, bar[Fu]
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, hh>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VP>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fu>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<0>(indexId1) == vertexId1::template indices_of_field<1>(indexId1) &&
                     vertexId2::template indices_of_field<0>(indexId2) == vertexId2::template indices_of_field<1>(indexId2) &&
                     vertexId3::template indices_of_field<0>(indexId3) == vertexId3::template indices_of_field<1>(indexId3)) {
                     const auto correction_S = 1. + get_alphas(context)/Pi * delta_hAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     const auto correction_P = 1. + get_alphas(context)/Pi * delta_AhAA_2loopQCD_for_quark_loop(result.m_decay, mInternal1, ren_scale);
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                     temp_result.form_factor_eps = correction_P * temp_result.form_factor_eps;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: gWp, bar[gWp], gWp
   {
      using vertexId1 = Vertex<typename bar<gWp>::type, gWp, hh>;
      using vertexId2 = Vertex<typename bar<gWp>::type, gWp, VP>;
      using vertexId3 = Vertex<typename bar<gWp>::type, gWp, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<gWp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<gWp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<gWp>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: gWpC, bar[gWpC], gWpC
   {
      using vertexId1 = Vertex<typename bar<gWpC>::type, gWpC, hh>;
      using vertexId2 = Vertex<typename bar<gWpC>::type, gWpC, VP>;
      using vertexId3 = Vertex<typename bar<gWpC>::type, gWpC, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<gWpC>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<gWpC>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<gWpC>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<Hp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g2n2_SSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val.
                        value(0, 1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], Hp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g4n4_SSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[VWp], Hp
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g5n5_SVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(0, 1), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[VWp], VWp
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g7n7_SVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[Hp], Hp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<1>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g6n6_VSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[Hp], VWp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<1>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<1>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g8n8_VSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::odd_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], VWp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g9n9_VVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g10n10_VVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(cxx_diagrams::
                        TripleVectorVertex::even_permutation {}), 1.0i*vertexId3Val.value(
                        cxx_diagrams::TripleVectorVertex::odd_permutation {}),
                     ren_scale,
                     Finite);
               }
            }
         }
      }
   }

   // topology T4
   // internal particles in the diagram: Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP, VP>;

      static constexpr double symmetryFac {2.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<3>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
               const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t4g1n15_SS(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T4
   // internal particles in the diagram: VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VP, VWp>;

      static constexpr double symmetryFac {2.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<3>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
               const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t4g2n16_VV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value3(), 1.0i*vertexId2Val.value1
                     (), 1.0i*vertexId2Val.value2(),
                  ren_scale,
                  Finite);
            }
         }
      }
   }

   // topology T6
   // internal particles in the diagram: Hp, conj[VWp]
   {
      using vertexId1 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId2 = Vertex<hh, Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<3>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t6g1n19_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T6
   // internal particles in the diagram: conj[Hp], VWp
   {
      using vertexId1 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId2 = Vertex<hh, typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<2>(indexId1);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<3>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t6g1n19_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T7
   // internal particles in the diagram: Hp, conj[VWp]
   {
      using vertexId1 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId2 = Vertex<hh, Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<3>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t7g1n20_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T7
   // internal particles in the diagram: conj[Hp], VWp
   {
      using vertexId1 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId2 = Vertex<hh, typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<2>(indexId1);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<3>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t7g1n20_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   return result;
}

// hh -> {VP, VZ}
template<>
Decay_amplitude_SVV CLASSNAME::calculate_amplitude_1l<hh, VP, VZ>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<fields::VP>::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::VZ>::type const& idx_3) const {

   // amplitude type
   Decay_amplitude_SVV result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_vector_1 = 0.;
   result.m_vector_2 = qedqcd.displayPoleMZ();

   // set the initial value of an amplitude to 0
   result.form_factor_g = std::complex<double>(0., 0.);
   result.form_factor_11 = std::complex<double>(0., 0.);
   result.form_factor_12 = std::complex<double>(0., 0.);
   result.form_factor_21 = std::complex<double>(0., 0.);
   result.form_factor_22 = std::complex<double>(0., 0.);
   result.form_factor_eps = std::complex<double>(0., 0.);

   // FormCalc's Finite variable
   static constexpr double Finite {1.};

   const double ren_scale {result.m_decay};

   // ----------------- 1-loop contributions to the amplitude -----------------

   // topology T1
   // internal particles in the diagram: Fd, bar[Fd], Fd
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, hh>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VZ>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fd>(vertexId2::template indices_of_field<1>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<1>(indexId1) == vertexId1::template indices_of_field<0>(indexId1) &&
                     vertexId2::template indices_of_field<1>(indexId2) == vertexId2::template indices_of_field<0>(indexId2) &&
                     vertexId3::template indices_of_field<1>(indexId3) == vertexId3::template indices_of_field<0>(indexId3)
                     && result.m_decay/mInternal1 < 0.8 && result.m_vector_2/result.m_decay < 0.75
                     ) {
                     const double correction_S = 1 - get_alphas(context)/Pi;
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Fe, bar[Fe], Fe
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fe, hh>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, VZ>;
      using vertexId3 = Vertex<typename bar<Fe>::type, Fe, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fe>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fe>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Fu, bar[Fu], Fu
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, hh>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VZ>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Fu>(vertexId2::template indices_of_field<1>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<1>(indexId1) == vertexId1::template indices_of_field<0>(indexId1) &&
                     vertexId2::template indices_of_field<1>(indexId2) == vertexId2::template indices_of_field<0>(indexId2) &&
                     vertexId3::template indices_of_field<1>(indexId3) == vertexId3::template indices_of_field<0>(indexId3)
                     && result.m_decay/mInternal1 < 0.8 && result.m_vector_2/result.m_decay < 0.75
                     ) {
                     const double correction_S = 1 - get_alphas(context)/Pi;
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[gWp], gWp, bar[gWp]
   {
      using vertexId1 = Vertex<typename bar<gWp>::type, gWp, hh>;
      using vertexId2 = Vertex<typename bar<gWp>::type, gWp, VZ>;
      using vertexId3 = Vertex<typename bar<gWp>::type, gWp, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<gWp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<gWp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<gWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[gWpC], gWpC, bar[gWpC]
   {
      using vertexId1 = Vertex<typename bar<gWpC>::type, gWpC, hh>;
      using vertexId2 = Vertex<typename bar<gWpC>::type, gWpC, VZ>;
      using vertexId3 = Vertex<typename bar<gWpC>::type, gWpC, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<gWpC>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<gWpC>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<gWpC>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[Hp], Hp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<Hp>::type, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g2n2_SSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(0, 1), 1.0i*vertexId3Val.
                        value(1, 0),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], Hp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g4n4_SSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], VWp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g5n5_SVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g7n7_SVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VZ>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g6n6_VSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(0, 1), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], Hp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VZ>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<1>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g8n8_VSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], VWp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g9n9_VVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[VWp], VWp
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g10n10_VVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(cxx_diagrams::
                        TripleVectorVertex::odd_permutation {}), 1.0i*vertexId3Val.value(
                        cxx_diagrams::TripleVectorVertex::odd_permutation {}),
                     ren_scale,
                     Finite);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fd], Fd, bar[Fd]
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, hh>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VZ>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fd>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<0>(indexId1) == vertexId1::template indices_of_field<1>(indexId1) &&
                     vertexId2::template indices_of_field<0>(indexId2) == vertexId2::template indices_of_field<1>(indexId2) &&
                     vertexId3::template indices_of_field<0>(indexId3) == vertexId3::template indices_of_field<1>(indexId3)
                     && result.m_decay/mInternal1 < 0.8 && result.m_vector_2/result.m_decay < 0.75
                     ) {
                     const double correction_S = 1 - get_alphas(context)/Pi;
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fe], Fe, bar[Fe]
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fe, hh>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, VZ>;
      using vertexId3 = Vertex<typename bar<Fe>::type, Fe, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fe>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fu], Fu, bar[Fu]
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, hh>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VZ>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fu>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  auto temp_result = symmetryFac * colorFac * calculate_diagram_SVV_t1g1n1_FFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale);
                  if (static_cast<int>(flexibledecay_settings.get(FlexibleDecay_settings::include_higher_order_corrections)) > 0 &&
                     vertexId1::template indices_of_field<0>(indexId1) == vertexId1::template indices_of_field<1>(indexId1) &&
                     vertexId2::template indices_of_field<0>(indexId2) == vertexId2::template indices_of_field<1>(indexId2) &&
                     vertexId3::template indices_of_field<0>(indexId3) == vertexId3::template indices_of_field<1>(indexId3)
                     && result.m_decay/mInternal1 < 0.8 && result.m_vector_2/result.m_decay < 0.75
                     ) {
                     const double correction_S = 1 - get_alphas(context)/Pi;
                     temp_result.form_factor_g   = correction_S * temp_result.form_factor_g;
                     temp_result.form_factor_11  = correction_S * temp_result.form_factor_11;
                     temp_result.form_factor_12  = correction_S * temp_result.form_factor_12;
                     temp_result.form_factor_21  = correction_S * temp_result.form_factor_21;
                     temp_result.form_factor_22  = correction_S * temp_result.form_factor_22;
                  }
                  result += temp_result;
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: gWp, bar[gWp], gWp
   {
      using vertexId1 = Vertex<typename bar<gWp>::type, gWp, hh>;
      using vertexId2 = Vertex<typename bar<gWp>::type, gWp, VZ>;
      using vertexId3 = Vertex<typename bar<gWp>::type, gWp, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<gWp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<gWp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<gWp>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: gWpC, bar[gWpC], gWpC
   {
      using vertexId1 = Vertex<typename bar<gWpC>::type, gWpC, hh>;
      using vertexId2 = Vertex<typename bar<gWpC>::type, gWpC, VZ>;
      using vertexId3 = Vertex<typename bar<gWpC>::type, gWpC, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<2>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<gWpC>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<gWpC>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<gWpC>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g3n3_UUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<Hp>::type, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g2n2_SSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val.
                        value(0, 1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[Hp], VWp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g4n4_SSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[VWp], Hp
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g5n5_SVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(0, 1), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: Hp, conj[VWp], VWp
   {
      using vertexId1 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g7n7_SVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::odd_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[Hp], Hp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VZ>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<1>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g6n6_VSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val
                        .value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[Hp], VWp
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VWp, VZ>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<1>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<2>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g8n8_VSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(0, 1), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.
                        value(cxx_diagrams::TripleVectorVertex::odd_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VWp, conj[VWp], Hp
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g9n9_VVS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<2>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t1g10n10_VVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(cxx_diagrams::
                        TripleVectorVertex::even_permutation {}), 1.0i*vertexId3Val.value(
                        cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale,
                     Finite);
               }
            }
         }
      }
   }

   // topology T3
   // internal particles in the diagram: VZ, Hp
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<3>(indexId2))
               continue;

            if (vertexId2::template indices_of_field<0>(indexId2) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
               const double mInternal2 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t3g1n13_VS(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T3
   // internal particles in the diagram: VZ, VWp
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<3>(indexId2))
               continue;

            if (vertexId2::template indices_of_field<0>(indexId2) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
               const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t3g2n14_VV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value2(), 1.0i*vertexId2Val.value3
                     (), 1.0i*vertexId2Val.value1(),
                  ren_scale,
                  Finite);
            }
         }
      }
   }

   // topology T4
   // internal particles in the diagram: Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP, VZ>;

      static constexpr double symmetryFac {2.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<3>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
               const double mInternal2 {context.mass<Hp>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t4g1n15_SS(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T4
   // internal particles in the diagram: VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp, VZ>;

      static constexpr double symmetryFac {2.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<3>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};
               const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t4g2n16_VV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value2(), 1.0i*vertexId2Val.value1
                     (), 1.0i*vertexId2Val.value3(),
                  ren_scale,
                  Finite);
            }
         }
      }
   }

   // topology T6
   // internal particles in the diagram: Hp, conj[VWp]
   {
      using vertexId1 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId2 = Vertex<hh, Hp, typename conj<VWp>::type, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<3>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<2>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t6g1n19_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T6
   // internal particles in the diagram: conj[Hp], VWp
   {
      using vertexId1 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId2 = Vertex<hh, typename conj<Hp>::type, VWp, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<2>(indexId1);
            const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<3>(indexId2);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t6g1n19_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T7
   // internal particles in the diagram: Hp, conj[VWp]
   {
      using vertexId1 = Vertex<typename conj<Hp>::type, VWp, VZ>;
      using vertexId2 = Vertex<hh, Hp, typename conj<VWp>::type, VP>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<3>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<2>(indexId1);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<VWp>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t7g1n20_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T7
   // internal particles in the diagram: conj[Hp], VWp
   {
      using vertexId1 = Vertex<Hp, typename conj<VWp>::type, VZ>;
      using vertexId2 = Vertex<hh, typename conj<Hp>::type, VP, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {

            // skip indices that don't match external indices
            const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
            const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<2>(indexId1);

            if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
               continue;

            // connect internal particles in vertices
            if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
               continue;

            if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<3>(indexId2))
               continue;

            auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
            auto const vertexId2Val = vertexId2::evaluate(indexId2, context);

            if (!vertexId1Val.isZero() && !vertexId2Val.isZero()) {
               // internal masses
               const double mInternal1 {context.mass<Hp>(vertexId1::template indices_of_field<0>(indexId1))};
               const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<1>(indexId1))};

               result += symmetryFac * colorFac * calculate_diagram_SVV_t7g1n20_SV(
                  result.m_decay, result.m_vector_1, result.m_vector_2,
                  mInternal1, mInternal2,
                  1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(),
                  ren_scale);
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, Fd, bar[Fd]
   {
      using vertexId1 = Vertex<typename bar<Fd>::type, Fd, VZ>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VP>;
      using vertexId3 = Vertex<hh, VZ, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Fd>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fd>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g6n36_VFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.right
                        (), -1.0i*vertexId2Val.left(), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, Fe, bar[Fe]
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fe, VZ>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, VP>;
      using vertexId3 = Vertex<hh, VZ, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Fe>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g6n36_VFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.right
                        (), -1.0i*vertexId2Val.left(), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, Fu, bar[Fu]
   {
      using vertexId1 = Vertex<typename bar<Fu>::type, Fu, VZ>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VP>;
      using vertexId3 = Vertex<hh, VZ, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<1>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<1>(indexId1) != vertexId2::template indices_of_field<0>(indexId2))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Fu>(vertexId1::template indices_of_field<1>(indexId1))};
                  const double mInternal3 {context.mass<typename bar<Fu>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g6n36_VFF(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.right
                        (), -1.0i*vertexId2Val.left(), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, bar[gWp], gWp
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<typename bar<gWp>::type, gWp, VZ>;
      using vertexId3 = Vertex<typename bar<gWp>::type, gWp, VP>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<gWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<gWp>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g8n38_VUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, bar[gWpC], gWpC
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<typename bar<gWpC>::type, gWpC, VZ>;
      using vertexId3 = Vertex<typename bar<gWpC>::type, gWpC, VP>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<gWpC>::type>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<gWpC>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g8n38_VUU(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), -1.0i*vertexId2Val.value(1), -1.0i*vertexId3Val.
                        value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, Hp, conj[Hp]
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<Hp, typename conj<Hp>::type, VP>;
      using vertexId3 = Vertex<Hp, typename conj<Hp>::type, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g7n37_VSS(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(1, 0), 1.0i*vertexId3Val.
                        value(0, 1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, conj[Hp], VWp
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<Hp, typename conj<VWp>::type, VP>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VWp, VZ>;

      constexpr double symmetryFac {1.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g9n39_VSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, Hp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<typename conj<Hp>::type, VP, VWp>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VZ>;

      constexpr double symmetryFac {1.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g9n39_VSV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(), 1.0i*vertexId3Val.value()
                        ,
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T9
   // internal particles in the diagram: VZ, VWp, conj[VWp]
   {
      using vertexId1 = Vertex<hh, VZ, VZ>;
      using vertexId2 = Vertex<typename conj<VWp>::type, VP, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VWp, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
               const auto externalFieldIndicesIn2 = vertexId2::template indices_of_field<1>(indexId2);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<1>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SVV_t9g10n40_VVV(
                     result.m_decay, result.m_vector_1, result.m_vector_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.value(), 1.0i*vertexId2Val.value(cxx_diagrams::
                        TripleVectorVertex::even_permutation {}), 1.0i*vertexId3Val.value(
                        cxx_diagrams::TripleVectorVertex::even_permutation {}),
                     ren_scale,
                     Finite);
               }
            }
         }
      }
   }

   return result;
}

// hh -> {VZ, VZ}
template<>
Decay_amplitude_SVV CLASSNAME::calculate_amplitude_tree<hh, VZ, VZ>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<fields::VZ>::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::VZ>::type const& idx_3) const {

   // amplitude type
   Decay_amplitude_SVV result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_vector_1 = qedqcd.displayPoleMZ();
   result.m_vector_2 = qedqcd.displayPoleMZ();

   // set the initial value of an amplitude to 0
   result.form_factor_g = std::complex<double>(0., 0.);
   result.form_factor_11 = std::complex<double>(0., 0.);
   result.form_factor_12 = std::complex<double>(0., 0.);
   result.form_factor_21 = std::complex<double>(0., 0.);
   result.form_factor_22 = std::complex<double>(0., 0.);
   result.form_factor_eps = std::complex<double>(0., 0.);
   // @todo correct prefactors
   // tree-level amplitude
   const auto indices = concatenate(idx_1, idx_2, idx_3);
   const auto vertex =  Vertex<hh, VZ, VZ>::evaluate(indices, context);

   result.form_factor_g += vertex.value();

   return result;
}

// hh -> {VWp, conj[VWp]}
template<>
Decay_amplitude_SVV CLASSNAME::calculate_amplitude_tree<hh, VWp, typename conj<VWp>::type>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<fields::VWp>::type const& idx_2,
   typename cxx_diagrams::field_indices<typename fields::conj<fields::VWp>::type>::type const& idx_3) const {

   // amplitude type
   Decay_amplitude_SVV result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_vector_1 = context.physical_mass<VWp>(idx_2);
   result.m_vector_2 = context.physical_mass<typename conj<VWp>::type>(idx_3);

   // set the initial value of an amplitude to 0
   result.form_factor_g = std::complex<double>(0., 0.);
   result.form_factor_11 = std::complex<double>(0., 0.);
   result.form_factor_12 = std::complex<double>(0., 0.);
   result.form_factor_21 = std::complex<double>(0., 0.);
   result.form_factor_22 = std::complex<double>(0., 0.);
   result.form_factor_eps = std::complex<double>(0., 0.);
   // @todo correct prefactors
   // tree-level amplitude
   const auto indices = concatenate(idx_1, idx_2, idx_3);
   const auto vertex =  Vertex<hh, typename conj<VWp>::type, VWp>::evaluate(indices, context);

   result.form_factor_g += vertex.value();

   return result;
}

// hh -> {bar[Fv], Fv}
template<>
Decay_amplitude_SFF CLASSNAME::calculate_amplitude_1l<hh, typename bar<Fv>::type, Fv>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<typename fields::bar<fields::Fv>::type>::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::Fv>::type const& idx_3) const {

   // amplitude type
   Decay_amplitude_SFF result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_fermion_1 = context.physical_mass<typename bar<Fv>::type>(idx_2);
   result.m_fermion_2 = context.physical_mass<Fv>(idx_3);

   // set the initial value of an amplitude to 0
   result.form_factor_left = std::complex<double>(0., 0.);
   result.form_factor_right = std::complex<double>(0., 0.);

   // FormCalc's Finite variable
   static constexpr double Finite {1.};

   const double ren_scale {result.m_decay};

   // ----------------- 1-loop contributions to the amplitude -----------------

   // topology T1
   // internal particles in the diagram: bar[Fe], Fe, Hp
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fv, typename conj<Hp>::type>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, hh>;
      using vertexId3 = Vertex<typename bar<Fv>::type, Fe, Hp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<0>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fe>::type>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g1n1_FFS(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.right
                        (),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: bar[Fe], Fe, VWp
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fv, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, hh>;
      using vertexId3 = Vertex<typename bar<Fv>::type, Fe, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId3::template indices_of_field<0>(indexId3);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};
                  const double mInternal2 {context.mass<typename bar<Fe>::type>(vertexId2::template indices_of_field<0>(indexId2))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g3n3_FFV(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.left(), -1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left()
                        , 1.0i*vertexId2Val.right(), -1.0i*vertexId3Val.right(), -1.0i*vertexId3Val.
                        left(),
                     ren_scale,
                     Finite);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], Hp, Fe
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fv, typename conj<Hp>::type>;
      using vertexId2 = Vertex<typename bar<Fv>::type, Fe, Hp>;
      using vertexId3 = Vertex<hh, Hp, typename conj<Hp>::type>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<0>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId2::template indices_of_field<2>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g2n2_SSF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left(),
                        1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[Hp], VWp, Fe
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fv, typename conj<Hp>::type>;
      using vertexId2 = Vertex<typename bar<Fv>::type, Fe, VWp>;
      using vertexId3 = Vertex<hh, Hp, typename conj<VWp>::type>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<0>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<Hp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g4n4_SVF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     1.0i*vertexId1Val.left(), 1.0i*vertexId1Val.right(), -1.0i*vertexId2Val.right()
                        , -1.0i*vertexId2Val.left(), 1.0i*vertexId3Val.value(0, 1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], Hp, Fe
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fv, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename bar<Fv>::type, Fe, Hp>;
      using vertexId3 = Vertex<hh, typename conj<Hp>::type, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<0>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId2::template indices_of_field<2>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g5n5_VSF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.left(), -1.0i*vertexId1Val.right(), 1.0i*vertexId2Val.left()
                        , 1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.value(0, 1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: conj[VWp], VWp, Fe
   {
      using vertexId1 = Vertex<typename bar<Fe>::type, Fv, typename conj<VWp>::type>;
      using vertexId2 = Vertex<typename bar<Fv>::type, Fe, VWp>;
      using vertexId3 = Vertex<hh, typename conj<VWp>::type, VWp>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<0>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<typename conj<VWp>::type>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g6n6_VVF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.right
                        (), -1.0i*vertexId2Val.left(), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T1
   // internal particles in the diagram: VZ, VZ, Fv
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId3 = Vertex<hh, VZ, VZ>;

      static constexpr double symmetryFac {1.000000000000000};

      static constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<0>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId2::template indices_of_field<0>(indexId2);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId1::template indices_of_field<0>(indexId1) != vertexId2::template indices_of_field<1>(indexId2))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<VZ>(vertexId2::template indices_of_field<2>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fv>::type>(vertexId1::template indices_of_field<0>(indexId1))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t1g6n6_VVF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.right
                        (), -1.0i*vertexId2Val.left(), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: Fd, bar[Fd], VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<typename bar<Fd>::type, Fd, VZ>;
      using vertexId3 = Vertex<typename bar<Fd>::type, Fd, hh>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Fd>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fd>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g6n18_VFF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.left(
                        ), -1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.
                        right(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: Fe, bar[Fe], VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<typename bar<Fe>::type, Fe, VZ>;
      using vertexId3 = Vertex<typename bar<Fe>::type, Fe, hh>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Fe>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fe>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g6n18_VFF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.left(
                        ), -1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.
                        right(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: Fu, bar[Fu], VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<typename bar<Fu>::type, Fu, VZ>;
      using vertexId3 = Vertex<typename bar<Fu>::type, Fu, hh>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {3.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId3::template indices_of_field<2>(indexId3);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId2::template indices_of_field<2>(indexId2))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Fu>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<Fu>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g6n18_VFF(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), -1.0i*vertexId2Val.left(
                        ), -1.0i*vertexId2Val.right(), 1.0i*vertexId3Val.left(), 1.0i*vertexId3Val.
                        right(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: bar[gWp], gWp, VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<typename bar<gWp>::type, gWp, hh>;
      using vertexId3 = Vertex<typename bar<gWp>::type, gWp, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<gWp>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<gWp>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g8n20_VUU(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), 1.0i*vertexId2Val.value(
                        ), -1.0i*vertexId3Val.value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: bar[gWpC], gWpC, VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<typename bar<gWpC>::type, gWpC, hh>;
      using vertexId3 = Vertex<typename bar<gWpC>::type, gWpC, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<2>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<0>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<gWpC>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<typename bar<gWpC>::type>(vertexId2::template indices_of_field<0>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g8n20_VUU(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), 1.0i*vertexId2Val.value(
                        ), -1.0i*vertexId3Val.value(1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: Hp, conj[Hp], VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<hh, Hp, typename conj<Hp>::type>;
      using vertexId3 = Vertex<Hp, typename conj<Hp>::type, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<2>(indexId2))};
                  const double mInternal3 {context.mass<Hp>(vertexId2::template indices_of_field<1>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g7n19_VSS(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), 1.0i*vertexId2Val.value(
                        ), 1.0i*vertexId3Val.value(0, 1),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: conj[Hp], VWp, VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<hh, Hp, typename conj<VWp>::type>;
      using vertexId3 = Vertex<typename conj<Hp>::type, VWp, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<Hp>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g9n21_VSV(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), 1.0i*vertexId2Val.value(
                        0, 1), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: Hp, conj[VWp], VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<hh, typename conj<Hp>::type, VWp>;
      using vertexId3 = Vertex<Hp, typename conj<VWp>::type, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<Hp>::type>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g9n21_VSV(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), 1.0i*vertexId2Val.value(
                        0, 1), 1.0i*vertexId3Val.value(),
                     ren_scale);
               }
            }
         }
      }
   }

   // topology T10
   // internal particles in the diagram: VWp, conj[VWp], VZ
   {
      using vertexId1 = Vertex<typename bar<Fv>::type, Fv, VZ>;
      using vertexId2 = Vertex<hh, typename conj<VWp>::type, VWp>;
      using vertexId3 = Vertex<typename conj<VWp>::type, VWp, VZ>;

      constexpr double symmetryFac {2.000000000000000};

      constexpr double colorFac {1.000000000000000};

      // loops over vertices' indices
      for (const auto& indexId1: index_range<vertexId1>()) {
         for (const auto& indexId2: index_range<vertexId2>()) {
            for (const auto& indexId3: index_range<vertexId3>()) {

               // skip indices that don't match external indices
               const auto externalFieldIndicesIn1 = vertexId2::template indices_of_field<0>(indexId2);
               const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
               const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<0>(indexId1);

               if (externalFieldIndicesIn1 != idx_1 || externalFieldIndicesIn2 != idx_2 || externalFieldIndicesIn3 != idx_3)
                  continue;

               // connect internal particles in vertices
               if (vertexId1::template indices_of_field<2>(indexId1) != vertexId3::template indices_of_field<2>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<1>(indexId2) != vertexId3::template indices_of_field<1>(indexId3))
                  continue;

               if (vertexId2::template indices_of_field<2>(indexId2) != vertexId3::template indices_of_field<0>(indexId3))
                  continue;

               auto const vertexId1Val = vertexId1::evaluate(indexId1, context);
               auto const vertexId2Val = vertexId2::evaluate(indexId2, context);
               auto const vertexId3Val = vertexId3::evaluate(indexId3, context);

               if (!vertexId1Val.isZero() && !vertexId2Val.isZero() && !vertexId3Val.isZero()) {
                  // internal masses
                  const double mInternal1 {context.mass<VZ>(vertexId1::template indices_of_field<2>(indexId1))};
                  const double mInternal2 {context.mass<typename conj<VWp>::type>(vertexId2::template indices_of_field<1>(indexId2))};
                  const double mInternal3 {context.mass<VWp>(vertexId2::template indices_of_field<2>(indexId2))};

                  result += symmetryFac * colorFac * calculate_diagram_SFF_t5g10n22_VVV(
                     result.m_decay, result.m_fermion_1, result.m_fermion_2,
                     mInternal1, mInternal2, mInternal3,
                     -1.0i*vertexId1Val.right(), -1.0i*vertexId1Val.left(), 1.0i*vertexId2Val.value(
                        ), 1.0i*vertexId3Val.value(cxx_diagrams::TripleVectorVertex::
                        even_permutation {}),
                     ren_scale);
               }
            }
         }
      }
   }

   return result;
}

// hh -> {bar[Fd], Fd}
template<>
Decay_amplitude_SFF CLASSNAME::calculate_amplitude_tree<hh, typename bar<Fd>::type, Fd>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<typename fields::bar<fields::Fd>::type>::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::Fd>::type const& idx_3) const {

   // amplitude type
   Decay_amplitude_SFF result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_fermion_1 = context.physical_mass<typename bar<Fd>::type>(idx_2);
   result.m_fermion_2 = context.physical_mass<Fd>(idx_3);

   // set the initial value of an amplitude to 0
   result.form_factor_left = std::complex<double>(0., 0.);
   result.form_factor_right = std::complex<double>(0., 0.);
   // @todo correct prefactors
   // tree-level amplitude
   const auto indices = concatenate(idx_3, idx_2, idx_1);
   const auto vertex =  Vertex<typename bar<Fd>::type, Fd, hh>::evaluate(indices, context);

   result.form_factor_left += vertex.left();
   result.form_factor_right += vertex.right();

   return result;
}

// hh -> {bar[Fu], Fu}
template<>
Decay_amplitude_SFF CLASSNAME::calculate_amplitude_tree<hh, typename bar<Fu>::type, Fu>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<typename fields::bar<fields::Fu>::type>::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::Fu>::type const& idx_3) const {

   // amplitude type
   Decay_amplitude_SFF result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_fermion_1 = context.physical_mass<typename bar<Fu>::type>(idx_2);
   result.m_fermion_2 = context.physical_mass<Fu>(idx_3);

   // set the initial value of an amplitude to 0
   result.form_factor_left = std::complex<double>(0., 0.);
   result.form_factor_right = std::complex<double>(0., 0.);
   // @todo correct prefactors
   // tree-level amplitude
   const auto indices = concatenate(idx_3, idx_2, idx_1);
   const auto vertex =  Vertex<typename bar<Fu>::type, Fu, hh>::evaluate(indices, context);

   result.form_factor_left += vertex.left();
   result.form_factor_right += vertex.right();

   return result;
}

// hh -> {bar[Fe], Fe}
template<>
Decay_amplitude_SFF CLASSNAME::calculate_amplitude_tree<hh, typename bar<Fe>::type, Fe>(
   const context_base& context,
   typename cxx_diagrams::field_indices<fields::hh>::type const& idx_1,
   typename cxx_diagrams::field_indices<typename fields::bar<fields::Fe>::type>::type const& idx_2,
   typename cxx_diagrams::field_indices<fields::Fe>::type const& idx_3) const {

      // amplitude type
   Decay_amplitude_SFF result;

   // external particles' masses
   result.m_decay = context.physical_mass<hh>(idx_1);
   result.m_fermion_1 = context.physical_mass<typename bar<Fe>::type>(idx_2);
   result.m_fermion_2 = context.physical_mass<Fe>(idx_3);

   // @todo correct prefactors
   // tree-level amplitude
   const auto indices = concatenate(idx_3, idx_2, idx_1);
   const auto vertex =  Vertex<typename bar<Fe>::type, Fe, hh>::evaluate(indices, context);

   result.form_factor_left += vertex.left();
   result.form_factor_right += vertex.right();

   return result;
}


// -------- specializations for decays needing higher order SM corrections --------

// List of potential Z boson decay products excluding pure SM decays
typedef boost::mpl::list<
> BSMForZdecay;

// List of potential W boson decay products excluding pure SM decays
typedef boost::mpl::list<
> BSMForWdecay;

template <typename DecayProds, typename V>
bool check_3body_Vff_decay(const context_base& context, double  mHOS, const typename field_indices<V>::type& idx)  {
   bool found_problem = false;
   const double mVOS {context.physical_mass<V>(idx)};
   boost::mpl::for_each<DecayProds>(
      [mHOS, mVOS, &idx, &context, &found_problem](auto arg) {
         using T = decltype(arg);
         using Field1 = typename boost::mpl::at<T, boost::mpl::int_<0>>::type;
         using Field2 = typename boost::mpl::at<T, boost::mpl::int_<1>>::type;
         using vertexId1 = Vertex<V, Field1, Field2>;
         for (const auto& indexId1: index_range<vertexId1>()) {
            const auto externalFieldIndicesIn1 = vertexId1::template indices_of_field<0>(indexId1);
            if (externalFieldIndicesIn1 != idx) {
               continue;
            }
            const auto externalFieldIndicesIn2 = vertexId1::template indices_of_field<1>(indexId1);
            const auto externalFieldIndicesIn3 = vertexId1::template indices_of_field<2>(indexId1);

            const double mInternal2 {context.physical_mass<Field1>(externalFieldIndicesIn2)};
            const double mInternal3 {context.physical_mass<Field2>(externalFieldIndicesIn3)};

            if (mHOS-mVOS > mInternal2 + mInternal3
                && !Vertex<V, Field1, Field2>::evaluate(indexId1, context).isZero()) {
               found_problem = true;
            }
         }
      }
   );
   return found_problem;
}

#include "decays/specializations/H/decay_H_to_ZZ.inc"
#include "decays/specializations/H/decay_H_to_WpWm.inc"
#include "decays/specializations/H/decay_H_to_gg.inc"
#include "decays/specializations/H/decay_H_to_AA.inc"
#include "decays/specializations/H/decay_H_to_AZ.inc"
#include "decays/specializations/H/decay_H_to_ubaru.inc"
#include "decays/specializations/H/decay_H_to_dbard.inc"
#include "decays/specializations/H/decay_H_to_lbarl.inc"

void CLASSNAME::calculate_hh_decays()
{

   if (run_to_decay_particle_scale) {
      auto decay_mass = PHYSICAL(Mhh);
      if (decay_mass > qedqcd.displayPoleMZ()) {
         try {
            model.run_to(decay_mass);
         }
         catch (const NonPerturbativeRunningError& e) {
            WARNING(std::string(e.what()) + " λ(μ=" + std::to_string(decay_mass) + " GeV) = " + std::to_string(model.get_Lambdax()));
         }
         model.solve_ewsb_tree_level();
         model.calculate_DRbar_masses();
      }
   }

   context_base context {model};

   auto& decays = decay_table.get_hh_decays();

   if (context.physical_mass<hh>(std::array<int, 0>{}) > context.physical_mass<
      VG>(std::array<int, 0> {}) + context.physical_mass<VG>(std::array<int, 0>
      {})) {
      decays.set_decay(partial_width_hh_to_VGVG(model), {standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VG), standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VG)}, create_process_string<hh, VG,
         VG>({},{},{}));
   }

   if (context.physical_mass<hh>(std::array<int, 0>{}) > context.physical_mass<
      VP>(std::array<int, 0> {}) + context.physical_mass<VP>(std::array<int, 0>
      {})) {
      decays.set_decay(partial_width_hh_to_VPVP(model), {standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VP), standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VP)}, create_process_string<hh, VP,
         VP>({},{},{}));
   }

   if (context.physical_mass<hh>(std::array<int, 0>{}) > context.physical_mass<
      VP>(std::array<int, 0> {}) + context.physical_mass<VZ>(std::array<int, 0>
      {})) {
      decays.set_decay(partial_width_hh_to_VPVZ(model), {standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VP), standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VZ)}, create_process_string<hh, VP,
         VZ>({},{},{}));
   }
      decays.set_decay(partial_width_hh_to_VZVZ(model), {standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VZ), standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VZ)}, create_process_string<hh, VZ,
         VZ>({},{},{}));
      decays.set_decay(partial_width_hh_to_conjVWpVWp(model), {-
         standard_model_info::get_pdg_code_for_particle(standard_model_info::VWp), standard_model_info::
         get_pdg_code_for_particle(standard_model_info::VWp)}, create_process_string<hh,
         typename conj<VWp>::type, VWp>({},{},{}));

   for (int gO1 = 0; gO1 < 3; ++gO1) {
      for (int gO2 = 0; gO2 < 3; ++gO2) {
         if (context.physical_mass<hh>(std::array<int, 0>{}) < context.
            physical_mass<typename bar<Fv>::type>(std::array<int, 1> {gO1}) +
            context.physical_mass<Fv>(std::array<int, 1> {gO2})) {
            continue;
         }
         decays.set_decay(partial_width_hh_to_barFvFv(model, gO1, gO2
            ), {-standard_model_info::get_pdg_code_for_particle(standard_model_info::Fv, gO1), standard_model_info
            ::get_pdg_code_for_particle(standard_model_info::Fv, gO2)},
            create_process_string<hh, typename bar<Fv>::type, Fv>({},{gO1},{gO2
            }));

      }
   }

   for (int gO1 = 0; gO1 < 3; ++gO1) {
      for (int gO2 = 0; gO2 < 3; ++gO2) {
         if (context.physical_mass<hh>(std::array<int, 0>{}) < context.
            physical_mass<typename bar<Fd>::type>(std::array<int, 1> {gO1}) +
            context.physical_mass<Fd>(std::array<int, 1> {gO2})) {
            continue;
         }
         decays.set_decay(partial_width_hh_to_barFdFd(model, gO1, gO2
            ), {-standard_model_info::get_pdg_code_for_particle(standard_model_info::Fd, gO1), standard_model_info
            ::get_pdg_code_for_particle(standard_model_info::Fd, gO2)},
            create_process_string<hh, typename bar<Fd>::type, Fd>({},{gO1},{gO2
            }));

      }
   }

   for (int gO1 = 0; gO1 < 3; ++gO1) {
      for (int gO2 = 0; gO2 < 3; ++gO2) {
         decays.set_decay(partial_width_hh_to_barFuFu(model, gO1, gO2
            ), {-standard_model_info::get_pdg_code_for_particle(standard_model_info::Fu, gO1), standard_model_info
            ::get_pdg_code_for_particle(standard_model_info::Fu, gO2)},
            create_process_string<hh, typename bar<Fu>::type, Fu>({},{gO1},{gO2
            }));

      }
   }

   for (int gO1 = 0; gO1 < 3; ++gO1) {
      for (int gO2 = 0; gO2 < 3; ++gO2) {
         if (context.physical_mass<hh>(std::array<int, 0>{}) < context.
            physical_mass<typename bar<Fe>::type>(std::array<int, 1> {gO1}) +
            context.physical_mass<Fe>(std::array<int, 1> {gO2})) {
            continue;
         }
         decays.set_decay(partial_width_hh_to_barFeFe(model, gO1, gO2
            ), {-standard_model_info::get_pdg_code_for_particle(standard_model_info::Fe, gO1), standard_model_info
            ::get_pdg_code_for_particle(standard_model_info::Fe, gO2)},
            create_process_string<hh, typename bar<Fe>::type, Fe>({},{gO1},{gO2
            }));

      }
   }

   if (flexibledecay_settings.get(FlexibleDecay_settings::call_higgstools) ||
      flexibledecay_settings.get(FlexibleDecay_settings::call_lilith) ||
      flexibledecay_settings.get(FlexibleDecay_settings::print_effc_block))
   {
      auto found = std::find_if(std::begin(neutral_higgs_effc), std::end(
         neutral_higgs_effc), [](NeutralHiggsEffectiveCouplings const& effC) {
         return effC.particle == field_as_string<hh>({});});
      found->width = decays.get_total_width();
      found->mass = context.physical_mass<hh>({});
      const auto _indices = concatenate(typename cxx_diagrams::field_indices<hh>::type {}, typename cxx_diagrams::field_indices<hh>::type {}, typename cxx_diagrams::field_indices<hh>::type {}, typename cxx_diagrams::field_indices<hh>::type {});
      const auto h4vertex = Vertex<hh, hh, hh, hh>::evaluate(_indices, context);
      found->lam = std::real(h4vertex.value());
   }
}

double CLASSNAME::partial_width_hh_to_VGVG(Standard_model const& model)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<VG>::type out_1_indices {};
   const typename field_indices<VG>::type out_2_indices {};

   return get_partial_width<hh, VG, VG>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_VPVP(Standard_model const& model)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<VP>::type out_1_indices {};
   const typename field_indices<VP>::type out_2_indices {};

   return get_partial_width<hh, VP, VP>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_VPVZ(Standard_model const& model)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<VP>::type out_1_indices {};
   const typename field_indices<VZ>::type out_2_indices {};

   return get_partial_width<hh, VP, VZ>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_VZVZ(Standard_model const& model)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<VZ>::type out_1_indices {};
   const typename field_indices<VZ>::type out_2_indices {};

   return get_partial_width<hh, VZ, VZ>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_conjVWpVWp(Standard_model const& model)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<VWp>::type out_1_indices {};
   const typename field_indices<typename conj<VWp>::type>::type out_2_indices {};

   return get_partial_width<hh, VWp, typename conj<VWp>::type>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_barFvFv(Standard_model const& model, int gO1, int gO2)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<typename bar<Fv>::type>::type out_1_indices {{gO1}};
   const typename field_indices<Fv>::type out_2_indices {{gO2}};

   return get_partial_width<hh, typename bar<Fv>::type, Fv>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_barFdFd(Standard_model const& model, int gO1, int gO2)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<typename bar<Fd>::type>::type out_1_indices {{gO1}};
   const typename field_indices<Fd>::type out_2_indices {{gO2}};

   return get_partial_width<hh, typename bar<Fd>::type, Fd>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_barFuFu(Standard_model const& model, int gO1, int gO2)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<typename bar<Fu>::type>::type out_1_indices {{gO1}};
   const typename field_indices<Fu>::type out_2_indices {{gO2}};

   return get_partial_width<hh, typename bar<Fu>::type, Fu>(context, in_indices, out_1_indices, out_2_indices);
}

double CLASSNAME::partial_width_hh_to_barFeFe(Standard_model const& model, int gO1, int gO2)
{
   context_base context {model};
   const typename field_indices<hh>::type in_indices {};
   const typename field_indices<typename bar<Fe>::type>::type out_1_indices {{gO1}};
   const typename field_indices<Fe>::type out_2_indices {{gO2}};

   return get_partial_width<hh, typename bar<Fe>::type, Fe>(context, in_indices, out_1_indices, out_2_indices);
}


} // namespace flexiblesusy
