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

#ifndef SINGLE_SCALE_CONSTRAINT_H
#define SINGLE_SCALE_CONSTRAINT_H

#include "logger.hpp"
#include <string>

namespace flexiblesusy {

class Model;

class Single_scale_constraint {
public:
   virtual ~Single_scale_constraint() = default;
   virtual void apply() = 0;                    ///< apply constraint
   virtual double get_scale() const = 0;        ///< get scale where to apply
   virtual std::string name() const { return "unknown constraint"; } ///< name of constraint
   virtual void set_model(Model*) = 0; ///< set model where to apply the constraint
};

} // namespace flexiblesusy

#endif
