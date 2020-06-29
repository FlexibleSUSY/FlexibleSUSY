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

#ifndef SLHA_IO_H
#define SLHA_IO_H

#include "slha_format.hpp"

#include <complex>
#include <functional>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace SLHAea {
   class Coll;
   class Line;
} // namespace SLHAea

namespace flexiblesusy {

   class Spectrum_generator_settings;
   class Physical_input;
   struct PMNS_parameters;

/**
 * @class SLHA_io
 * @brief Handles reading and writing of SLHA files
 *
 * Reading: There are two ways to read block entries from SLHA files:
 * a) using the read_block() function with a %SLHA_io::Tuple_processor
 * or b) using the read_entry() function for each entry.  Note, that
 * a) is much faster than b) (more than 1000 times) because b) needs
 * to search for the block each time read_entry() is called.
 *
 * Example how to use a tuple processor (fast!):
 * \code{.cpp}
void process_tuple(double* array, int key, double value) {
   array[key] = value;
}

void read_file() {
   double array[1000];

   SLHA_io reader;
   reader.read_from_file("file.slha");

   SLHA_io::Tuple_processor processor = [&array] (int key, double value) {
      return process_tuple(array, key, value);
   };

   reader.read_block("MyBlock", processor);
}
 * \endcode
 *
 * Example how to use a for loop (slow!):
 * \code{.cpp}
void read_file() {
   double array[1000];

   SLHA_io reader;
   reader.read_from_file("file.slha");

   for (int i = 0; i < 1000; i++) {
      array[i] = reader.read_entry("MyBlock", i);
   }
}
 * \endcode
 */
class SLHA_io {
public:
   using Tuple_processor = std::function<void(int, double)>;
   enum Position { front, back };
   struct Modsel {
      bool quark_flavour_violated{false};   ///< MODSEL[6]
      bool lepton_flavour_violated{false};  ///< MODSEL[6]
      double parameter_output_scale{0.};    ///< MODSEL[12]
      void clear() { *this = Modsel(); }
   };

   struct CKM_wolfenstein {
      double lambdaW{0.}, aCkm{0.}, rhobar{0.}, etabar{0.};
      void clear() { *this = CKM_wolfenstein(); }
   };

   SLHA_io();
   SLHA_io(const SLHA_io&);
   SLHA_io(SLHA_io&&) noexcept;
   ~SLHA_io();

   SLHA_io& operator=(const SLHA_io&);
   SLHA_io& operator=(SLHA_io&&) noexcept;

   void clear();

   // reading functions
   bool block_exists(const std::string&) const;
   void fill(softsusy::QedQcd&) const;
   void fill(Spectrum_generator_settings&) const;
   void fill(Physical_input&) const;
   const Modsel& get_modsel() const { return modsel; }
   const SLHAea::Coll& get_data() const;
   void read_from_file(const std::string&);
   void read_from_source(const std::string&);
   void read_from_stream(std::istream&);
   double read_block(const std::string&, const Tuple_processor&) const;
   template <class Derived>
   double read_block(const std::string&, Eigen::PlainObjectBase<Derived>&) const;
   double read_block(const std::string&, double&) const;
   double read_entry(const std::string&, int) const;
   double read_scale(const std::string&) const;

   // writing functions
   void set_data(const SLHAea::Coll&);
   void set_block(const std::ostringstream&, Position position = back);
   void set_block(const std::string&, Position position = back);
   void set_blocks(const std::vector<std::string>&, Position position = back);
   void set_block(const std::string&, double, const std::string&, double scale = 0.);
   template <class Derived>
   void set_block(const std::string&, const Eigen::MatrixBase<Derived>&, const std::string&, double scale = 0.);
   template <class Derived>
   void set_block_imag(const std::string&, const Eigen::MatrixBase<Derived>&, const std::string&, double scale = 0.);
   void set_modsel(const Modsel&);
   void set_physical_input(const Physical_input&);
   void set_settings(const Spectrum_generator_settings&);
   void set_sminputs(const softsusy::QedQcd&);
   void write_to_file(const std::string&) const;
   void write_to_stream() const;
   void write_to_stream(std::ostream&) const;

private:
   std::unique_ptr<SLHAea::Coll> data; ///< SHLA data
   Modsel modsel{};            ///< data from block MODSEL

   static std::string block_head(const std::string& name, double scale);
   static bool read_scale(const SLHAea::Line& line, double& scale);

   void read_modsel();
   double read_matrix(const std::string&, double*, int, int) const;
   double read_matrix(const std::string&, std::complex<double>*, int, int) const;
   double read_vector(const std::string&, double*, int) const;
   double read_vector(const std::string&, std::complex<double>*, int) const;

   void set_vector(const std::string&, const double*, const std::string&, double, int);
   void set_vector(const std::string&, const std::complex<double>*, const std::string&, double, int);
   void set_matrix(const std::string&, const double*, const std::string&, double, int, int);
   void set_matrix(const std::string&, const std::complex<double>*, const std::string&, double, int, int);

   void set_vector_imag(const std::string&, const double*, const std::string&, double, int);
   void set_vector_imag(const std::string&, const std::complex<double>*, const std::string&, double, int);
   void set_matrix_imag(const std::string&, const double*, const std::string&, double, int, int);
   void set_matrix_imag(const std::string&, const std::complex<double>*, const std::string&, double, int, int);
};

/**
 * Fills a matrix or vector from a SLHA block
 *
 * @param block_name block name
 * @param dense matrix or vector to be filled
 *
 * @return scale (or 0 if no scale is defined)
 */
template <class Derived>
double SLHA_io::read_block(const std::string& block_name, Eigen::PlainObjectBase<Derived>& dense) const
{
   return dense.cols() == 1
      ? read_vector(block_name, dense.data(), dense.rows())
      : read_matrix(block_name, dense.data(), dense.rows(), dense.cols());
}

/**
 * Writes real part of a matrix or vector to SLHA object
 *
 * @param name bloch name
 * @param dense matrix ox vector
 * @param symbol symbol name
 * @param scale renormalization scale
 */
template<class Derived>
void SLHA_io::set_block(const std::string& name,
                        const Eigen::MatrixBase<Derived>& dense,
                        const std::string& symbol, double scale)
{
   dense.cols() == 1
      ? set_vector(name, dense.eval().data(), symbol, scale, dense.rows())
      : set_matrix(name, dense.eval().data(), symbol, scale, dense.rows(), dense.cols());
}

/**
 * Writes imaginary part of a matrix or vector to SLHA object
 *
 * @param name bloch name
 * @param dense matrix ox vector
 * @param symbol symbol name
 * @param scale renormalization scale
 */
template<class Derived>
void SLHA_io::set_block_imag(const std::string& name,
                             const Eigen::MatrixBase<Derived>& dense,
                             const std::string& symbol, double scale)
{
   dense.cols() == 1
      ? set_vector_imag(name, dense.eval().data(), symbol, scale, dense.rows())
      : set_matrix_imag(name, dense.eval().data(), symbol, scale, dense.rows(), dense.cols());
}

} // namespace flexiblesusy

#endif
