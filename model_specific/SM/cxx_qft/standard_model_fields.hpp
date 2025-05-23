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
 * @file cxx_qft/SM_fields.hpp
 *
 * This file was generated with FlexibleSUSY 2.7.1 and SARAH 4.14.5 .
 */

#ifndef STANDARDMODEL_CXXQFT_FIELDS_H
#define STANDARDMODEL_CXXQFT_FIELDS_H

#include <boost/hana/tuple.hpp>

#include "cxx_qft/fields.hpp"

namespace flexiblesusy {
namespace standard_model_cxx_diagrams::fields {

using cxx_diagrams::fields::ParticleType;
using cxx_diagrams::fields::ParticleColorRep;
using cxx_diagrams::fields::conj;
using cxx_diagrams::fields::bar;

struct VG {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::octet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(21);
   using lorentz_conjugate = VG;
};

struct gG {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::octet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = typename bar<gG>::type;
};

struct Hp {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 1;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = typename conj<Hp>::type;
};

struct Fv {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(12,14,16);
   using lorentz_conjugate = typename bar<Fv>::type;
};

struct Ah {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = Ah;
};

struct hh {
   static constexpr auto particleType = ParticleType::scalar;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(25);
   using lorentz_conjugate = hh;
};

struct VP {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(22);
   using lorentz_conjugate = VP;
};

struct VZ {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(23);
   using lorentz_conjugate = VZ;
};

struct gP {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = true;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = typename bar<gP>::type;
};

struct gZ {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 0;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = typename bar<gZ>::type;
};

struct gWp {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 1;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = typename bar<gWp>::type;
};

struct gWpC {
   static constexpr auto particleType = ParticleType::ghost;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = -1;
   static constexpr auto pdgids = boost::hana::make_tuple(0);
   using lorentz_conjugate = typename bar<gWpC>::type;
};

struct Fd {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -0.3333333333333333;
   static constexpr auto pdgids = boost::hana::make_tuple(1,3,5);
   using lorentz_conjugate = typename bar<Fd>::type;
};

struct Fu {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::triplet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = 0.6666666666666666;
   static constexpr auto pdgids = boost::hana::make_tuple(2,4,6);
   using lorentz_conjugate = typename bar<Fu>::type;
};

struct Fe {
   static constexpr auto particleType = ParticleType::fermion;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int, 0>,
     boost::mpl::vector_c<int, 3>
   >;
   static constexpr int numberOfGenerations = 3;
   using sm_flags = boost::mpl::vector_c<bool, true, true, true>;
   static constexpr int numberOfFieldIndices = 1;
   static constexpr double electricCharge = -1;
   static constexpr auto pdgids = boost::hana::make_tuple(11,13,15);
   using lorentz_conjugate = typename bar<Fe>::type;
};

struct VWp {
   static constexpr auto particleType = ParticleType::vector;
   static constexpr auto colorRep = ParticleColorRep::singlet;
   static constexpr auto massless = false;
   using index_bounds = boost::mpl::pair<
     boost::mpl::vector_c<int>,
     boost::mpl::vector_c<int>
   >;
   static constexpr int numberOfGenerations = 1;
   using sm_flags = boost::mpl::vector_c<bool, true>;
   static constexpr int numberOfFieldIndices = 0;
   static constexpr double electricCharge = 1;
   static constexpr auto pdgids = boost::hana::make_tuple(24);
   using lorentz_conjugate = typename conj<VWp>::type;
};

// Named fields
using Electron = Fe;

using scalars = boost::mpl::vector<Hp, Ah, hh>;
using fermions = boost::mpl::vector<Fv, Fd, Fu, Fe>;
using vectors = boost::mpl::vector<VG, VP, VZ, VWp>;
using ghosts = boost::mpl::vector<gG, gP, gZ, gWp, gWpC>;

} // namespace standard_model_cxx_diagrams::fields

namespace cxx_diagrams::fields{

// Fields that are their own Lorentz conjugates.
template<> struct conj<standard_model_cxx_diagrams::fields::VG> { using type = standard_model_cxx_diagrams::fields::VG; };
template<> struct conj<standard_model_cxx_diagrams::fields::Ah> { using type = standard_model_cxx_diagrams::fields::Ah; };
template<> struct conj<standard_model_cxx_diagrams::fields::hh> { using type = standard_model_cxx_diagrams::fields::hh; };
template<> struct conj<standard_model_cxx_diagrams::fields::VP> { using type = standard_model_cxx_diagrams::fields::VP; };
template<> struct conj<standard_model_cxx_diagrams::fields::VZ> { using type = standard_model_cxx_diagrams::fields::VZ; };

template<>
struct is_vector<flexiblesusy::standard_model_cxx_diagrams::fields::VG > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::standard_model_cxx_diagrams::fields::gG > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::standard_model_cxx_diagrams::fields::Hp > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::standard_model_cxx_diagrams::fields::Fv > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::standard_model_cxx_diagrams::fields::Ah > : public std::true_type {};
template<>
struct is_scalar<flexiblesusy::standard_model_cxx_diagrams::fields::hh > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::standard_model_cxx_diagrams::fields::VP > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::standard_model_cxx_diagrams::fields::VZ > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::standard_model_cxx_diagrams::fields::gP > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::standard_model_cxx_diagrams::fields::gZ > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::standard_model_cxx_diagrams::fields::gWp > : public std::true_type {};
template<>
struct is_ghost<flexiblesusy::standard_model_cxx_diagrams::fields::gWpC > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::standard_model_cxx_diagrams::fields::Fd > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::standard_model_cxx_diagrams::fields::Fu > : public std::true_type {};
template<>
struct is_fermion<flexiblesusy::standard_model_cxx_diagrams::fields::Fe > : public std::true_type {};
template<>
struct is_vector<flexiblesusy::standard_model_cxx_diagrams::fields::VWp > : public std::true_type {};


} // namespace cxx_diagrams::fields

} // namespace flexiblesusy

#endif
