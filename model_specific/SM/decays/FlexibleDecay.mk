#  ====================================================================
#  This file is part of FlexibleSUSY.
#
#  FlexibleSUSY is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  FlexibleSUSY is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FlexibleSUSY.  If not, see
#  <http://www.gnu.org/licenses/>.
#  ====================================================================

LIBSM_SRC += \
		$(DIR)/decays/SM_decay_table.cpp \
		$(DIR)/decays/SM_decays.cpp

LIBSM_HDR += \
		$(DIR)/decays/SM_decay_table.hpp \
		$(DIR)/decays/SM_decays.hpp \
		$(DIR)/decays/SM_decay_amplitudes.hpp

EXESM_SRC += \
		$(DIR)/run_decays_SM.cpp
