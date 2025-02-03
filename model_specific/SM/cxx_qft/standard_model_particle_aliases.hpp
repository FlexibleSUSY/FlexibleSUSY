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
 * @file standard_model_particle_aliases.hpp
 *
 * This file was generated with FlexibleSUSY @FlexibleSUSYVersion@ and SARAH @SARAHVersion@ .
 */

#include "cxx_qft/standard_model_fields.hpp"

namespace flexiblesusy::standard_model_cxx_diagrams::fields {

using Higgs = hh;
using WpBoson = VWp;
using WmBoson = typename conj<VWp>::type;
using Photon = VP;
using ZBoson = VZ;
using Gluon = VG;
using ChargedLepton = Fe;
using Neutrino = Fv;
using UpTypeQuark = Fu;
using DownTypeQuark = Fd;

}
