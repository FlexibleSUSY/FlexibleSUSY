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

#ifndef CXXQFT_FIELDS_H
#define CXXQFT_FIELDS_H

#include <array> /**< Only for field_indices and isSMField*/

#include <boost/array.hpp>
#include <boost/version.hpp>

#include <boost/mpl/erase.hpp>
#include <boost/mpl/joint_view.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>

#include <boost/fusion/adapted/boost_array.hpp>
#include <boost/fusion/adapted/mpl.hpp>

#include <boost/hana/config.hpp>
#include <boost/hana/equal.hpp>
#include <boost/hana/fwd/transform.hpp>

#if BOOST_VERSION >= 105800
#include <boost/fusion/include/move.hpp>
#else
#include <boost/fusion/include/copy.hpp>
#endif

namespace flexiblesusy::cxx_diagrams {
   /** @brief Declare a type that can hold the field indices for any given field.
    * @todo Field should have nested numberOfFieldIndices.
    * @returns `::type` evaluates to `std::array<int,int N>`.
    */
   template <class Field>
   struct field_indices {
      /** @b using keyword here is type alias
       * (see https://en.cppreference.com/w/cpp/language/type_alias).
       * It allows to use @code{.cpp}typename field_indices<Field>::type@endcode
       * to get @code{.cpp}std::array<int, Field::numberOfFieldIndices>@endcode
       * type in the code.
       */
      using type = std::array<int, Field::numberOfFieldIndices>;
   };

   namespace detail {
      /** @brief Declare a metafunction that gives the int value for Field
       * during compilation time or `int_<int N>` typename of this index.
       * @returns `::value` gives compile time integer value.
       * @returns `::type` gives `int_<int N>` typename.
       */
      template <class Field>
      struct number_of_field_indices {
         static constexpr int value =
            std::tuple_size<typename field_indices<Field>::type>::value;
         using type = boost::mpl::int_<value>;
      };
      /*
       * @param FieldSequence a Forward Sequence, e.g vector<type1,typ2,...>
       * (see https://www.boost.org/doc/libs/1_69_0/libs/mpl/doc/refmanual/forward-sequence.html).
       * @returns `::type` which accomodates `mpl::int_<N>` with `int N` being the
       * sum of @b number of indices of fields inside FieldSequence.
       * @returns `::value` which represents `N` where `N` being the `int`
       * sum of @b number of indices of fields inside FieldSequence
       */
      template <class FieldSequence>
      struct total_number_of_field_indices {
         using type = typename boost::mpl::fold<
            FieldSequence,
            boost::mpl::int_<0>,
            boost::mpl::plus<boost::mpl::_1,
                             number_of_field_indices<boost::mpl::_2>
                            >
            >::type;
         static constexpr int value = type::value;
      };
   } // namespace detail

   /* @note Is defined by compiler only if the number of generation for a field
    * is not 1.
    */
   template <class Field>
   std::enable_if_t<Field::numberOfGenerations != 1, bool>
   isSMField(const typename field_indices<Field>::type& indices) {
      boost::array<bool, Field::numberOfGenerations> sm_flags;

      #if BOOST_VERSION >= 105800
      boost::fusion::move(typename Field::sm_flags(), sm_flags); /**< Making copy of the information from Field to sm_flags */
      #else
      boost::fusion::copy(typename Field::sm_flags(), sm_flags);
      #endif

      return sm_flags[indices.front()]; /**< front() accesses the first element. Either all are SM indices or none. */
   }

   template <class Field>
   std::enable_if_t<Field::numberOfGenerations == 1, bool>
   isSMField(const typename field_indices<Field>::type&) {
      return boost::mpl::at_c<typename Field::sm_flags, 0>::type::value;
   }

   namespace fields
   {
   enum class ParticleType {
      scalar,
      fermion,
      vector,
      ghost
   };

   template<typename Field>
   struct is_massless {
      static constexpr bool value = Field::massless;
   };
   template<typename Field>
   constexpr bool is_massless_v = is_massless<Field>::value;

   enum class ParticleColorRep {
      singlet,
      triplet,
      anti_triplet,
      sextet,
      octet
   };

   template<typename Field>
   struct is_singlet {
      static constexpr bool value =
         Field::colorRep == ParticleColorRep::singlet;
   };
   template<typename Field>
   constexpr bool is_singlet_v = is_singlet<Field>::value;

   template<typename Field>
   struct is_triplet {
      static constexpr bool value = Field::colorRep == ParticleColorRep::triplet;
   };
   template<typename Field>
   constexpr bool is_triplet_v = is_triplet<Field>::value;

   template<typename Field>
   struct is_anti_triplet {
      static constexpr bool value =
         Field::colorRep == ParticleColorRep::anti_triplet;
   };
   template<typename Field>
   constexpr bool is_anti_triplet_v = is_anti_triplet<Field>::value;

   template<typename Field>
   struct is_octet {
      static constexpr bool value = Field::colorRep == ParticleColorRep::octet;
   };
   template<typename Field>
   constexpr bool is_octet_v = is_octet<Field>::value;

   template<typename Field>
   constexpr std::enable_if_t<
      is_triplet<Field>::value, ParticleColorRep
      >
   color_conj() {
      return ParticleColorRep::anti_triplet;
   }
   template<typename Field>
   constexpr std::enable_if_t<
      is_anti_triplet<Field>::value, ParticleColorRep
      >
   color_conj() {
      return ParticleColorRep::triplet;
   }
   template<typename Field>
   constexpr std::enable_if_t<
      !is_triplet<Field>::value && !is_anti_triplet<Field>::value, ParticleColorRep
      >
   color_conj() {
      return Field::colorRep;
   }

   template<class Field>
   struct bar {
      using index_bounds = typename Field::index_bounds;
      using sm_flags = typename Field::sm_flags;
      using lorentz_conjugate = Field;
      using type = bar<Field>;

      static constexpr int numberOfGenerations = Field::numberOfGenerations;
      static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
      static constexpr double electricCharge = -Field::electricCharge;
      static constexpr auto particleType = Field::particleType;
      static constexpr auto colorRep = color_conj<Field>();
      static constexpr auto massless = Field::massless;
      static constexpr auto pdgids = boost::hana::transform(Field::pdgids, [](int x) {return -x;});
   };

   template<class Field>
   struct conj {
      using index_bounds = typename Field::index_bounds;
      using sm_flags = typename Field::sm_flags;
      using lorentz_conjugate = Field;
      using type = conj<Field>;

      static constexpr int numberOfGenerations = Field::numberOfGenerations;
      static constexpr int numberOfFieldIndices = Field::numberOfFieldIndices;
      static constexpr double electricCharge = -Field::electricCharge;
      static constexpr auto particleType = Field::particleType;
      static constexpr auto colorRep = color_conj<Field>();
      static constexpr auto massless = Field::massless;
      static constexpr auto pdgids = boost::hana::transform(Field::pdgids, [](int x) {return -x;});
   };

   // Double Lorentz conjugation
   template <class Field>
   struct bar<bar<Field>> {
      using type = Field;
   };
   template <class Field>
   struct conj<conj<Field>> {
      using type = Field;
   };

   // Remove Lorentz conjugation
   template <class Field>
   struct remove_lorentz_conjugation {
      using type = Field;
   };

   template <class Field>
   struct remove_lorentz_conjugation<bar<Field>> {
      using type = Field;
   };

   template <class Field>
   struct remove_lorentz_conjugation<conj<Field>> {
      using type = Field;
   };

   template <class Field>
   struct is_scalar : public std::false_type {};

   template <class Field>
   struct is_scalar<bar<Field> > : public is_scalar<Field> {};

   template <class Field>
   struct is_scalar<conj<Field> > : public is_scalar<Field> {};

   template <class Field>
   struct is_fermion : public std::false_type {};

   template <class Field>
   struct is_fermion<bar<Field> > : public is_fermion<Field> {};

   template <class Field>
   struct is_fermion<conj<Field> > : public is_fermion<Field> {};

   template <class Field>
   struct is_vector : public std::false_type {};

   template <class Field>
   struct is_vector<bar<Field> > : public is_vector<Field> {};

   template <class Field>
   struct is_vector<conj<Field> > : public is_vector<Field> {};

   template <class Field>
   struct is_ghost : public std::false_type {};

   template <class Field>
   struct is_ghost<bar<Field> > : public is_ghost<Field> {};

   template <class Field>
   struct is_ghost<conj<Field> > : public is_ghost<Field> {};

   } // namespace generic_fields

   using fields::bar;
   using fields::conj;
   using fields::remove_lorentz_conjugation;

} // namespace flexiblesusy::cxx_diagrams

#endif
