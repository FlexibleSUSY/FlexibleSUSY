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

#ifndef SLHA_FORMAT_H
#define SLHA_FORMAT_H

namespace flexiblesusy {

/// SLHA line formatter for the MASS block entries
extern const char * const mass_formatter;
/// SLHA line formatter for the mixing matrix entries ;
extern const char * const mixing_matrix_formatter;
/// SLHA line formatter for vector entries
extern const char * const vector_formatter;
/// SLHA number formatter
extern const char * const number_formatter;
/// SLHA line formatter for entries with three indices
extern const char * const tensor_formatter;
/// SLHA scale formatter
extern const char * const scale_formatter;
/// SLHA line formatter for the one-element entries ;
extern const char * const single_element_formatter;
/// SLHA line formatter for the SPINFO block entries
extern const char * const spinfo_formatter;

namespace {
   /// maximum line length in SLHA output
   constexpr unsigned SLHA_MAX_LINE_LENGTH = 200;
} // namespace

#define FORMAT_MASS(pdg, mass, name)                                           \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int pdg_ = (pdg);                                                  \
      const double mass_ = (mass);                                             \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, mass_formatter, pdg_, mass_,    \
                    name_.c_str());                                            \
      return std::string(buf);                                                 \
   }()

#define FORMAT_MIXING_MATRIX(i, k, entry, name)                                \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int i_ = (i);                                                      \
      const int k_ = (k);                                                      \
      const double entry_ = (entry);                                           \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, mixing_matrix_formatter, i_,    \
                    k_, entry_, name_.c_str());                                \
      return std::string(buf);                                                 \
   }()

#define FORMAT_VECTOR(i, entry, name)                                          \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int i_ = (i);                                                      \
      const double entry_ = (entry);                                           \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, vector_formatter, i_, entry_,   \
                    name_.c_str());                                            \
      return std::string(buf);                                                 \
   }()

#define FORMAT_ELEMENT(pdg, value, name)                                       \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int pdg_ = (pdg);                                                  \
      const double value_ = (value);                                           \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, single_element_formatter, pdg_, \
                    value_, name_.c_str());                                    \
      return std::string(buf);                                                 \
   }()

#define FORMAT_SCALE(n)                                                        \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const double n_ = (n);                                                   \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, scale_formatter, n_);           \
      return std::string(buf);                                                 \
   }()

#define FORMAT_NUMBER(n, str)                                                  \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const double n_ = (n);                                                   \
      const std::string str_ = (str);                                          \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, number_formatter, n_,           \
                    str_.c_str());                                             \
      return std::string(buf);                                                 \
   }()

#define FORMAT_SPINFO(n, str)                                                  \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int n_ = (n);                                                      \
      const std::string str_ = (str);                                          \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, spinfo_formatter, n_,           \
                    str_.c_str());                                             \
      return std::string(buf);                                                 \
   }()

#define FORMAT_RANK_THREE_TENSOR(i, j, k, entry, name)                         \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int i_ = (i);                                                      \
      const int j_ = (j);                                                      \
      const int k_ = (k);                                                      \
      const double entry_ = (entry);                                           \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, tensor_formatter, i_, j_, k_,   \
                    entry_, name_.c_str());                                    \
      return std::string(buf);                                                 \
   }()

} // namespace flexiblesusy

#endif
