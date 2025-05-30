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

#include <string>

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
/// SLHA line formatter for the OBSINFO block entries
extern const char * const obsinfo_formatter;
/// FLHA line formatter for FWCOEF, IMFCOEF block entries
extern const char * const wilson_formatter;
/// SLHA line formatter for the DECAY block
extern const char * const format_total_width;
/// SLHA line formatter for the EFFECTIVECOUPLINGS block
extern const char * const format_effectivecouplings;

namespace {
   /// maximum line length in SLHA output
   constexpr unsigned SLHA_MAX_LINE_LENGTH = 200;
} // namespace

template <typename Container>
std::string format_decay(double br, const Container& pids, const std::string& name)
{
   const int nda = pids.size();
   char buf[SLHA_MAX_LINE_LENGTH] = { 0 };
   int written = 0;
   written += std::snprintf(buf + written, SLHA_MAX_LINE_LENGTH - written,
                            "   %16.8E  %2d  ", br, nda);
   for (int i = 0; i < nda; i++) {
      written += std::snprintf(buf + written, SLHA_MAX_LINE_LENGTH - written,
                               " %9d", pids[i]);
   }
   written += std::snprintf(buf + written, SLHA_MAX_LINE_LENGTH - written,
                            "  # %s\n", name.c_str());
   return std::string(buf);
}

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

#define FORMAT_OBSINFO(i, j, str)                                              \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int i_ = (i);                                                      \
      const int j_ = (j);                                                      \
      const std::string str_ = (str);                                          \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, obsinfo_formatter, i_, j_,      \
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

#define FORMAT_WILSON_COEFFICIENTS(f, m, x, y, ph, entry, name)                \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const std::string f_ = (f);                                              \
      const std::string m_ = (m);                                              \
      const int x_ = (x);                                                      \
      const int y_ = (y);                                                      \
      const int ph_ = (ph);                                                    \
      const double entry_ = (entry);                                           \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, wilson_formatter, f_.c_str(),   \
                    m_.c_str(), x_, y_, ph_, entry_, name_.c_str());           \
      return std::string(buf);                                                 \
   }()

#define FORMAT_TOTAL_WIDTH(pdg, width, name)                                   \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int pdg_ = (pdg);                                                  \
      const double width_ = (width);                                           \
      const std::string name_ = (name);                                        \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, format_total_width,             \
                    pdg_, width_, name_.c_str());                              \
      return std::string(buf);                                                 \
   }()

#define FORMAT_EFFECTIVECOUPLINGS(pdg1, pdg2, pdg3, width, comment)            \
   [&] {                                                                       \
      char buf[SLHA_MAX_LINE_LENGTH];                                          \
      const int pdg1_ = (pdg1);                                                \
      const int pdg2_ = (pdg2);                                                \
      const int pdg3_ = (pdg3);                                                \
      const double width_ = (width);                                           \
      const std::string comment_ = (comment);                                  \
      std::snprintf(buf, SLHA_MAX_LINE_LENGTH, format_effectivecouplings,      \
                    pdg1_, pdg2_, pdg3_, width_, comment_.c_str());            \
      return std::string(buf);                                                 \
   }()

} // namespace flexiblesusy

#endif

