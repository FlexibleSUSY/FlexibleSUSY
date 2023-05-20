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

#include "global_thread_pool.hpp"
#include "thread_pool.hpp"
#include "logger.hpp"

namespace flexiblesusy {

namespace {
struct Global_thread_pool_msg {
   Global_thread_pool_msg() {
      VERBOSE_MSG("Creating global thread pool ...");
   }
};
} // anonymous namespace

Thread_pool& global_thread_pool()
{
   static Global_thread_pool_msg msg;
   static Thread_pool tp;
   return tp;
}

} // namespace flexiblesusy
