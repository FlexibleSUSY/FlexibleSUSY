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

#ifndef SoftsusyMSSM_TWO_SCALE_H
#define SoftsusyMSSM_TWO_SCALE_H

#include "SoftsusyMSSM.hpp"
#include "model.hpp"
#include "softsusy.h"

namespace flexiblesusy {

class Two_scale;

template<>
class SoftsusyMSSM<Two_scale>: public Model, public softsusy::MssmSoftsusy {
public:
   SoftsusyMSSM();
   virtual ~SoftsusyMSSM();

   virtual void calculate_spectrum();
   virtual void clear_problems() {}
   virtual std::string name() const { return "SoftsusyMSSM"; }
   virtual void run_to(double, double eps = -1.0);
   virtual void print(std::ostream& s) const { s << static_cast<softsusy::MssmSoftsusy>(*this); }
   virtual void set_precision(double p) { precision = p; }

   void calculate_DRbar_masses() { calcDrBarPars(); }
   void set_scale(double scale) { setMu(scale); }
   double get_scale() const { return displayMu(); }
   SoftsusyMSSM calc_beta() const { return beta2(); }
   void setSugraBcs(double m0, double m12, double a0) { standardSugra(m0, m12, a0); }

private:
   double precision;
   SoftsusyMSSM(const softsusy::SoftParsMssm&);
};

}

#endif
