
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_slha_io

#include <boost/test/unit_test.hpp>

#include "linalg2.hpp"
#include "mixings.hpp"
#include "slha_io.hpp"
#include "slhaea.h"
#include "stopwatch.hpp"
#include "wrappers.hpp"
#include <Eigen/Core>
#include <complex>
#include <string>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_read_entry )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block A\n"
      "   1   1.2      # comment 1\n"
      "   2   1.3E+2   # comment 2\n"
      "   2   2.3E+2   # comment 3";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   BOOST_CHECK_EQUAL(reader.read_entry("A",1), 1.2);
   BOOST_CHECK_EQUAL(reader.read_entry("A",2), 230.);
}

BOOST_AUTO_TEST_CASE( test_read_block )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block Matrix\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::MatrixXd matrix(Eigen::MatrixXd::Zero(2,2));
   reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(matrix(0,0), 1.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 2.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 3.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 4.0);
}

BOOST_AUTO_TEST_CASE( test_read_block_dense )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block Matrix\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::Matrix<double,2,2> matrix(Eigen::Matrix<double,2,2>::Zero());
   reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(matrix(0,0), 1.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 2.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 3.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 4.0);
}

BOOST_AUTO_TEST_CASE( test_read_entry_doubled )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   std::string str =\
      "Block DOUBLE\n"
      "   1  1\n";
   block.str(str);
   coll.push_back(block);

   str =
      "Block DOUBLE\n"
      "   1  2\n";
   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   double entry = reader.read_entry("DOUBLE", 1);

   BOOST_CHECK_EQUAL(entry, 2.0);
}

BOOST_AUTO_TEST_CASE( test_read_block_doubled )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   std::string str =
      "Block Matrix\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";
   block.str(str);
   coll.push_back(block);

   str =
      "Block Matrix\n"
      "   1  1  5.0      # element 1,1\n"
      "   1  2  6.0      # element 1,2\n"
      "#  2  1  7.0      # element 2,1\n"
      "   2  2  8.0      # element 2,2\n";
   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::Matrix<double,2,2> matrix(Eigen::Matrix<double,2,2>::Zero());
   reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(matrix(0,0), 5.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 6.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 3.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 8.0);
}

BOOST_AUTO_TEST_CASE( test_read_scale )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   std::string str =
      "Block Matrix Q= 1234.56\n"
      "Block Matrix Q= 2000\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   str =
      "Block Matrix Q= 2000\n"
      "   1  1  11.0      # element 1,1\n"
      "   1  2  12.0      # element 1,2\n"
      "   2  1  13.0      # element 2,1\n"
      "   2  2  14.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   const double scale = reader.read_scale("Matrix");

   BOOST_CHECK_EQUAL(scale, 2000.);
}

BOOST_AUTO_TEST_CASE( test_read_scale_from_block )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   std::string str =
      "Block Matrix Q= 1234.56\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   str =
      "Block Matrix Q= 2000\n"
      "   1  1  11.0      # element 1,1\n"
      "   1  2  12.0      # element 1,2\n"
      "   2  1  13.0      # element 2,1\n"
      "   2  2  14.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::MatrixXd matrix(Eigen::MatrixXd::Zero(2,2));
   const double scale = reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(scale, 2000.0);
   BOOST_CHECK_EQUAL(matrix(0,0), 11.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 12.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 13.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 14.0);
}

BOOST_AUTO_TEST_CASE( test_write_read_matrix )
{
   const double scale = 100.0;
   const char* block_name = "Matrix";

   Eigen::MatrixXd matrix(2,2);
   matrix << 1, 2, 3, 4;

   SLHA_io slha_io;
   slha_io.set_block(block_name, matrix, "M", scale);

   Eigen::MatrixXd new_matrix(2,2);
   const double new_scale = slha_io.read_block(block_name, new_matrix);

   BOOST_CHECK_EQUAL(scale, new_scale);
   BOOST_CHECK_EQUAL(matrix(0,0), new_matrix(0,0));
   BOOST_CHECK_EQUAL(matrix(0,1), new_matrix(0,1));
   BOOST_CHECK_EQUAL(matrix(1,0), new_matrix(1,0));
   BOOST_CHECK_EQUAL(matrix(1,1), new_matrix(1,1));
}

// non-square matrices: the element (i,k) must be written as entry "i k" and
// read back to (i,k)
BOOST_AUTO_TEST_CASE(test_write_read_matrix_non_square)
{
   using C = std::complex<double>;
   Eigen::Matrix<double, 3, 2> m32;
   m32 << 11, 12, 21, 22, 31, 32;
   Eigen::Matrix<C, 2, 3> m23;
   m23 << C(11, -11), C(12, -12), C(13, -13), C(21, -21), C(22, -22),
      C(23, -23);

   SLHA_io slha_io;
   slha_io.set_block("M32", m32, "M32");
   slha_io.set_block("M23", m23, "M23");
   slha_io.set_block_imag("IMM23", m23, "M23");

   // entry "i k" must hold 10*i + k (or its negative for IMM23)
   const auto check_entries = [&slha_io](const std::string& name, int rows,
                                         int cols, double sign) {
      const auto& block = slha_io.get_data().at(name);
      for (int i = 1; i <= rows; i++) {
         for (int k = 1; k <= cols; k++) {
            const auto& line = block.at({std::to_string(i), std::to_string(k)});
            BOOST_CHECK_EQUAL(std::stod(line.at(2)), sign * (10 * i + k));
         }
      }
   };

   check_entries("M32", 3, 2, 1.);
   check_entries("M23", 2, 3, 1.);
   check_entries("IMM23", 2, 3, -1.);

   Eigen::Matrix<double, 3, 2> m32_read = Eigen::Matrix<double, 3, 2>::Zero();
   slha_io.read_block("M32", m32_read);
   BOOST_CHECK(m32_read == m32);

   Eigen::Matrix<double, 2, 3> m23_re = Eigen::Matrix<double, 2, 3>::Zero();
   Eigen::Matrix<double, 2, 3> m23_im = Eigen::Matrix<double, 2, 3>::Zero();
   slha_io.read_block("M23", m23_re);
   slha_io.read_block("IMM23", m23_im);
   BOOST_CHECK(m23_re == m23.real());
   BOOST_CHECK(m23_im == m23.imag());
}

BOOST_AUTO_TEST_CASE( test_write_read_matrix_complex )
{
   const double scale = 100.0;
   const char* block_name_real = "Matrix";
   const char* block_name_imag = "ImMatrix";
   const std::complex<double> i(0.0, 1.0);

   Eigen::MatrixXcd matrix(2,2);
   matrix << std::complex<double>(1.0,1.0),
             std::complex<double>(2.0,2.0),
             std::complex<double>(3.0,3.0),
             std::complex<double>(4.0,4.0);

   SLHA_io slha_io;
   slha_io.set_block(block_name_real, matrix, "M", scale);
   slha_io.set_block(block_name_imag, matrix, "M", scale);

   Eigen::MatrixXd new_matrix_real(2,2), new_matrix_imag(2,2);
   const double new_scale_real = slha_io.read_block(block_name_real, new_matrix_real);
   const double new_scale_imag = slha_io.read_block(block_name_imag, new_matrix_imag);

   Eigen::MatrixXcd new_matrix = new_matrix_real + new_matrix_imag*i;

   BOOST_CHECK_EQUAL(scale, new_scale_real);
   BOOST_CHECK_EQUAL(scale, new_scale_imag);
   BOOST_CHECK_EQUAL(matrix(0,0), new_matrix(0,0));
   BOOST_CHECK_EQUAL(matrix(0,1), new_matrix(0,1));
   BOOST_CHECK_EQUAL(matrix(1,0), new_matrix(1,0));
   BOOST_CHECK_EQUAL(matrix(1,1), new_matrix(1,1));
}

BOOST_AUTO_TEST_CASE( test_write_read_vector )
{
   const double scale = 100.0;
   const char* block_name = "Vector";

   Eigen::VectorXd vector(2);
   vector << 1, 2;

   SLHA_io slha_io;
   slha_io.set_block(block_name, vector, "V", scale);

   Eigen::VectorXd new_vector(2);
   const double new_scale = slha_io.read_block(block_name, new_vector);

   BOOST_CHECK_EQUAL(scale, new_scale);
   BOOST_CHECK_EQUAL(vector(0), new_vector(0));
   BOOST_CHECK_EQUAL(vector(1), new_vector(1));
}

BOOST_AUTO_TEST_CASE( test_write_read_vector_complex )
{
   const double scale = 100.0;
   const char* block_name_real = "Vector";
   const char* block_name_imag = "ImVector";
   const std::complex<double> i(0.0, 1.0);

   Eigen::VectorXcd vector(2);
   vector << std::complex<double>(1.0,1.0),
             std::complex<double>(2.0,2.0);

   SLHA_io slha_io;
   slha_io.set_block(block_name_real, vector, "M", scale);
   slha_io.set_block(block_name_imag, vector, "M", scale);

   Eigen::VectorXd new_vector_real(2), new_vector_imag(2);
   const double new_scale_real = slha_io.read_block(block_name_real, new_vector_real);
   const double new_scale_imag = slha_io.read_block(block_name_imag, new_vector_imag);

   Eigen::VectorXcd new_vector = new_vector_real + new_vector_imag*i;

   BOOST_CHECK_EQUAL(scale, new_scale_real);
   BOOST_CHECK_EQUAL(scale, new_scale_imag);
   BOOST_CHECK_EQUAL(vector(0), new_vector(0));
   BOOST_CHECK_EQUAL(vector(1), new_vector(1));
}

/**
 * Creates a SLHAea block with name `TestBlock' with
 *  `number_of_entries' entries of the form
 *
 *  Block TestBlock
 *     0  0 + offset
 *     1  1 + offset
 *     2  2 + offset
 *     3  3 + offset
 *
 * @param number_of_entries number of block entries
 * @param offset offset
 *
 * @return SLHAea block with block name `TestBlock'
 */
SLHAea::Block create_block(int number_of_entries, int offset, double scale = 0.)
{
   SLHAea::Block block;
   std::string str = "Block TestBlock";
   if (scale != 0.)
      str += " Q= " + std::to_string(scale);
   str += '\n';

   for (int i = 0; i < number_of_entries; i++) {
      const std::string key(std::to_string(i));
      const std::string num(std::to_string(i + offset));
      str += "   " + key + "  " + num + "\n";
   }

   block.str(str);

   return block;
}

void process_tuple(double* array, int key, double value)
{
   array[key] = value;
}

void check_array(double* array, int number_of_entries, int offset = 0)
{
   for (int i = 0; i < number_of_entries; i++)
      BOOST_CHECK_EQUAL(array[i], i + offset);
}

void reset_array(double* array, int number_of_entries, double value = 0.)
{
   for (int i = 0; i < number_of_entries; i++)
      array[i] = value;
}

BOOST_AUTO_TEST_CASE( test_processor_doubled )
{
   using namespace std::placeholders;
   const int number_of_entries = 10;
   SLHA_io reader;
   SLHAea::Coll coll;

   SLHAea::Block block = create_block(number_of_entries, 0, 1234.56);
   coll.push_back(block);

   block = create_block(number_of_entries, 10, 2345.67);
   coll.push_back(block);

   reader.set_data(coll);

   double array[number_of_entries];
   reset_array(array, number_of_entries, 0.);

   SLHA_io::Tuple_processor processor
      = std::bind(&process_tuple, array, _1, _2);

   const double scale = reader.read_block("TestBlock", processor);

   BOOST_CHECK_EQUAL(scale, 2345.67);
   check_array(array, number_of_entries, 10);
}

/**
 * This test compares the speed of reading an array from a long SLHA
 * block using a) a tuple processor and b) a for loop.  The tuple
 * processor is more complicated to use, but is more than 1000 times
 * faster than the for loop.
 */
BOOST_AUTO_TEST_CASE( test_processor_vs_loop )
{
   using namespace std::placeholders;
   Stopwatch timer;
   double processor_time = 0., loop_time = 0.;
   const int number_of_entries = 1000;
   SLHA_io reader;
   SLHAea::Coll coll;

   SLHAea::Block block = create_block(number_of_entries, 0);

   coll.push_back(block);
   reader.set_data(coll);

   double array[number_of_entries];
   reset_array(array, number_of_entries);

   SLHA_io::Tuple_processor processor
      = std::bind(&process_tuple, array, _1, _2);

   // measure time to read the array using a tuple processor
   timer.start();
   reader.read_block("TestBlock", processor);
   timer.stop();
   processor_time = timer.get_time_in_seconds();

   check_array(array, number_of_entries);
   reset_array(array, number_of_entries);

   // measure time to read the array using a for loop
   timer.start();
   for (int i = 0; i < number_of_entries; i++) {
      array[i] = reader.read_entry("TestBlock", i);
   }
   timer.stop();
   loop_time = timer.get_time_in_seconds();

   check_array(array, number_of_entries);

   BOOST_TEST_MESSAGE("time using the tuple processor: " << processor_time << " s");
   BOOST_TEST_MESSAGE("time using the for loop: " << loop_time << " s");

   BOOST_CHECK_LT(10 * processor_time, loop_time);
}

BOOST_AUTO_TEST_CASE( test_slha_mixing_matrix_convention )
{
   Eigen::Matrix<double, 2, 2> mass_matrix;
   Eigen::Matrix<std::complex<double>, 2, 2> Z;
   Eigen::Array<double, 2, 1> eigenvalues;

   mass_matrix(0,0) = 0;
   mass_matrix(0,1) = 1;
   mass_matrix(1,0) = 1;
   mass_matrix(1,1) = 0;

   // diagonalize with SARAH convention
   fs_diagonalize_symmetric(mass_matrix, eigenvalues, Z);

   BOOST_CHECK_GT(eigenvalues(0), 0.);
   BOOST_CHECK_GT(eigenvalues(1), 0.);

   BOOST_CHECK_GT(Z.imag().cwiseAbs().maxCoeff(), 0.);

   BOOST_CHECK((Z.row(0).imag().cwiseAbs().maxCoeff() > 0. &&
                Z.row(1).imag().cwiseAbs().maxCoeff() == 0.) ||
               (Z.row(1).imag().cwiseAbs().maxCoeff() > 0. &&
                Z.row(0).imag().cwiseAbs().maxCoeff() == 0.));

   // convert to SLHA convention
   convert_symmetric_fermion_mixings_to_slha(eigenvalues, Z);

   BOOST_CHECK(eigenvalues(0) < 0. || eigenvalues(1) < 0.);
   BOOST_CHECK_EQUAL(Z.imag().cwiseAbs().maxCoeff(), 0.);

   // reconstruct mass matrix
   Eigen::Matrix<std::complex<double>, 2, 2> reconstructed_mass_matrix;
   reconstructed_mass_matrix = Z.transpose() * eigenvalues.matrix().asDiagonal() * Z;

   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(0,0)), mass_matrix(0,0), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(0,1)), mass_matrix(0,1), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(1,0)), mass_matrix(1,0), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(1,1)), mass_matrix(1,1), 1.0e-15);

   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(0,0)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(0,1)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(1,0)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(1,1)), 0.);

   // convert to HK convention
   convert_symmetric_fermion_mixings_to_hk(eigenvalues, Z);

   BOOST_CHECK(eigenvalues(0) > 0. && eigenvalues(1) > 0.);
   BOOST_CHECK_GT(Z.imag().cwiseAbs().maxCoeff(), 0.);

   // reconstruct mass matrix
   reconstructed_mass_matrix = Z.transpose() * eigenvalues.matrix().asDiagonal() * Z;

   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(0,0)), mass_matrix(0,0), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(0,1)), mass_matrix(0,1), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(1,0)), mass_matrix(1,0), 1.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(Re(reconstructed_mass_matrix(1,1)), mass_matrix(1,1), 1.0e-15);

   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(0,0)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(0,1)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(1,0)), 0.);
   BOOST_CHECK_EQUAL(Im(reconstructed_mass_matrix(1,1)), 0.);
}

// A real row with imaginary parts from rounding errors must not be
// multiplied by i: before, such a row became purely imaginary (printed as
// zero in the SLHA block NMIX) and the sign of its mass was flipped.
BOOST_AUTO_TEST_CASE(test_slha_mixing_matrix_convention_rounding_noise)
{
   Eigen::Matrix<double, 3, 3> mass_matrix;
   // clang-format off
   mass_matrix << 300.,  0., -40.,
                    0., 600.,  50.,
                  -40.,  50., -200.;
   // clang-format on

   Eigen::Array<double, 3, 1> eigenvalues;
   Eigen::Matrix<std::complex<double>, 3, 3> Z;
   fs_diagonalize_symmetric(mass_matrix, eigenvalues, Z);

   // add rounding noise to the imaginary parts of the real rows
   for (int i = 0; i < 3; i++) {
      if (Z.row(i).imag().cwiseAbs().maxCoeff() == 0.) {
         Z(i, 2) += std::complex<double>(0., 1e-15);
      }
   }

   convert_symmetric_fermion_mixings_to_slha(eigenvalues, Z);

   // the mixing matrix is real up to the noise and diagonalizes the mass matrix
   BOOST_CHECK_SMALL(Z.imag().cwiseAbs().maxCoeff(), 1e-14);

   const Eigen::Matrix<double, 3, 3> reconstructed =
      Z.real().transpose() * eigenvalues.matrix().asDiagonal() * Z.real();

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_SMALL(std::abs(reconstructed(i, k) - mass_matrix(i, k)),
                           1e-10);
      }
   }

   // exactly one negative mass (det < 0)
   BOOST_CHECK_EQUAL((eigenvalues < 0.).count(), 1);
}

// A complex symmetric mass matrix (CP violation) can't be diagonalized by a
// real mixing matrix: the masses must stay positive and the mixing matrix
// unchanged (SLHA2 convention for CP violation), and has_imaginary_part()
// must report it, so that its imaginary part is written to the SLHA output.
BOOST_AUTO_TEST_CASE(test_slha_mixing_matrix_convention_cp_violation)
{
   using C = std::complex<double>;
   Eigen::Matrix<C, 3, 3> mass_matrix;
   // clang-format off
   mass_matrix << C(300., 50.),  C(0.),  C(-40.),
                  C(0.),         C(600.), C(50., -20.),
                  C(-40.),       C(50., -20.), C(-200., 30.);
   // clang-format on

   Eigen::Array<double, 3, 1> eigenvalues;
   Eigen::Matrix<C, 3, 3> Z;
   fs_diagonalize_symmetric(mass_matrix, eigenvalues, Z);

   const Eigen::Array<double, 3, 1> eigenvalues_hk(eigenvalues);
   const Eigen::Matrix<C, 3, 3> Z_hk(Z);

   convert_symmetric_fermion_mixings_to_slha(eigenvalues, Z);

   BOOST_CHECK((eigenvalues == eigenvalues_hk).all());
   BOOST_CHECK(Z == Z_hk);
   BOOST_CHECK((eigenvalues > 0.).all());
   BOOST_CHECK(has_imaginary_part(Z));

   // the mixing matrix still diagonalizes the mass matrix: M = Z^T m Z
   const Eigen::Matrix<C, 3, 3> reconstructed =
      Z.transpose() * eigenvalues.matrix().cast<C>().asDiagonal() * Z;

   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_SMALL(std::abs(reconstructed(i, k) - mass_matrix(i, k)),
                           1e-10);
      }
   }

   // 1x1 case: a complex phase that is not a multiple of Pi/2
   double m1 = 100.;
   Eigen::Matrix<C, 1, 1> z1;
   z1(0, 0) = std::polar(1., 1.);
   const auto z1_hk = z1;
   convert_symmetric_fermion_mixings_to_slha(m1, z1);
   BOOST_CHECK_EQUAL(m1, 100.);
   BOOST_CHECK(z1 == z1_hk);
   BOOST_CHECK(has_imaginary_part(z1));

   // real matrices and rounding noise are not reported
   BOOST_CHECK(!has_imaginary_part(Eigen::Matrix<double, 2, 2>::Identity()));
   Eigen::Matrix<C, 2, 2> noisy = Eigen::Matrix<C, 2, 2>::Identity();
   noisy(0, 1) = C(0., 1e-15);
   BOOST_CHECK(!has_imaginary_part(noisy));
}

// The rows of a complex scalar mixing matrix have arbitrary phases.  For a
// real mass matrix, convert_scalar_mixings_to_slha() must make them real
// (otherwise the real part in the SLHA output isn't the mixing matrix), and
// for a complex hermitian mass matrix the matrix must still diagonalize it.
BOOST_AUTO_TEST_CASE(test_slha_scalar_mixing_matrix_convention)
{
   using C = std::complex<double>;

   const auto check_diagonalizes = [](const Eigen::Matrix<C, 3, 3>& mass_matrix,
                                      const Eigen::Array<double, 3, 1>& m2,
                                      const Eigen::Matrix<C, 3, 3>& Z) {
      // M = Z^dagger m2 Z
      const Eigen::Matrix<C, 3, 3> reconstructed =
         Z.adjoint() * m2.matrix().cast<C>().asDiagonal() * Z;
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            BOOST_CHECK_SMALL(std::abs(reconstructed(i, k) - mass_matrix(i, k)),
                              1e-8);
         }
      }
   };

   // real mass matrix, diagonalized with random phases of the rows
   Eigen::Matrix<C, 3, 3> mass_matrix;
   // clang-format off
   mass_matrix << C(9.e4),  C(-3.e3), C(0.),
                  C(-3.e3), C(4.e5),  C(2.e4),
                  C(0.),    C(2.e4),  C(2.5e5);
   // clang-format on

   Eigen::Array<double, 3, 1> m2;
   Eigen::Matrix<C, 3, 3> Z;
   fs_diagonalize_hermitian(mass_matrix, m2, Z);
   for (int i = 0; i < 3; i++) {
      Z.row(i) *= std::polar(1., 0.7 + 1.3 * i);
   }
   BOOST_REQUIRE(has_imaginary_part(Z));

   convert_scalar_mixings_to_slha(Z);

   BOOST_CHECK(!has_imaginary_part(Z));
   check_diagonalizes(mass_matrix, m2, Z);
   for (int i = 0; i < 3; i++) {
      Eigen::Index k = 0;
      Z.row(i).cwiseAbs().maxCoeff(&k);
      BOOST_CHECK_GT(std::real(Z(i, k)), 0.);
   }

   // complex hermitian mass matrix (CP violation): Z stays complex
   mass_matrix(0, 1) = C(-3.e3, 2.e3);
   mass_matrix(1, 0) = std::conj(mass_matrix(0, 1));
   fs_diagonalize_hermitian(mass_matrix, m2, Z);
   convert_scalar_mixings_to_slha(Z);
   BOOST_CHECK(has_imaginary_part(Z));
   check_diagonalizes(mass_matrix, m2, Z);

   // a real matrix is left unchanged
   Eigen::Matrix<double, 2, 2> Zr;
   Zr << 0.6, -0.8, -0.8, -0.6;
   const auto Zr_orig = Zr;
   convert_scalar_mixings_to_slha(Zr);
   BOOST_CHECK(Zr == Zr_orig);
}

// The sfermion mixing matrix in the super-CKM basis must diagonalize the
// sfermion mass matrix written with the soft-breaking parameters in SLHA
// convention (m_Q^2 -> vL m_Q^2 vL^dagger, m_d^2 -> vR^* m_d^2 vR^T,
// T_d -> vR^* T_d vL^dagger).
BOOST_AUTO_TEST_CASE(test_slha_sfermion_mixing_super_ckm_basis)
{
   using C = std::complex<double>;
   using M3 = Eigen::Matrix<C, 3, 3>;
   using M6 = Eigen::Matrix<C, 6, 6>;

   // hermitian soft masses and a trilinear coupling with flavour mixing
   M3 mq2, md2, Td;
   // clang-format off
   mq2 << C(1.0e6),        C(2.0e4, 1.0e4), C(-3.0e3, 5.0e2),
          C(2.0e4, -1.0e4), C(1.1e6),       C(4.0e4, -2.0e3),
          C(-3.0e3, -5.0e2), C(4.0e4, 2.0e3), C(0.9e6);
   md2 << C(8.0e5),        C(-1.0e4, 3.0e3), C(2.0e3),
          C(-1.0e4, -3.0e3), C(8.5e5),       C(5.0e3, 1.0e3),
          C(2.0e3),         C(5.0e3, -1.0e3), C(7.0e5);
   Td  << C(10., 1.), C(2., -3.), C(-5.),
          C(1., 2.),  C(30.),     C(4., 1.),
          C(-2.),     C(3., 3.),  C(500., -20.);
   // clang-format on

   // rotations of the left- and right-handed fermions (as from fs_svd)
   M3 vL, vR;
   Eigen::Array<double, 3, 1> w;
   fs_diagonalize_hermitian(M3(mq2 * 1e-6 + M3::Identity()), w, vL);
   fs_diagonalize_hermitian(M3(md2.conjugate() * 1e-6 + 2. * M3::Identity()), w,
                            vR);

   // sfermion mass matrix in the interaction basis (f_L, f_R)
   const double v = 20.;
   M6 mass_matrix;
   mass_matrix << mq2, v * Td.adjoint(), v * Td, md2;

   Eigen::Array<double, 6, 1> m2;
   M6 Z;
   fs_diagonalize_hermitian(mass_matrix, m2, Z);

   convert_sfermion_mixings_to_slha(Z, vL, vR);

   // mass matrix with the parameters in SLHA convention
   const M3 mq2_slha = vL * mq2 * vL.adjoint();
   const M3 md2_slha = vR.conjugate() * md2 * vR.transpose();
   const M3 Td_slha = vR.conjugate() * Td * vL.adjoint();
   M6 mass_matrix_slha;
   mass_matrix_slha << mq2_slha, v * Td_slha.adjoint(), v * Td_slha, md2_slha;

   const M6 diag = Z * mass_matrix_slha * Z.adjoint();
   for (int i = 0; i < 6; i++) {
      for (int k = 0; k < 6; k++) {
         BOOST_CHECK_SMALL(std::abs(diag(i, k) - (i == k ? C(m2(i)) : C(0.))),
                           1e-6 * m2.maxCoeff());
      }
   }

   // left-handed sfermions only (sneutrinos)
   M3 ZV;
   Eigen::Array<double, 3, 1> mv2;
   fs_diagonalize_hermitian(mq2, mv2, ZV);
   convert_sfermion_mixings_to_slha(ZV, vL);
   const M3 diag_v = ZV * mq2_slha * ZV.adjoint();
   for (int i = 0; i < 3; i++) {
      for (int k = 0; k < 3; k++) {
         BOOST_CHECK_SMALL(
            std::abs(diag_v(i, k) - (i == k ? C(mv2(i)) : C(0.))),
            1e-6 * mv2.maxCoeff());
      }
   }

   // a real mixing matrix with real rotations stays real
   Eigen::Matrix<double, 2, 2> Zr;
   Zr << 0.6, -0.8, 0.8, 0.6;
   Eigen::Matrix<double, 1, 1> one;
   one << -1.;
   convert_sfermion_mixings_to_slha(Zr, one, one);
   BOOST_CHECK_EQUAL(Zr(0, 0), -0.6);
   BOOST_CHECK_EQUAL(Zr(1, 1), -0.6);
}

template<int N>
void convert_symmetric_fermion_mixings_to_slha_forloop(
   Eigen::Array<double, N, 1>& m,
   Eigen::Matrix<std::complex<double>, N, N>& z)
{
   for (int i = 0; i < N; i++) {
      // check if i'th row contains non-zero imaginary parts
      if (!is_zero(z.row(i).imag().cwiseAbs().maxCoeff())) {
         z.row(i) *= std::complex<double>(0.0,1.0);
         m(i) *= -1;
#ifdef ENABLE_DEBUG
         if (!is_zero(z.row(i).imag().cwiseAbs().maxCoeff())) {
            WARNING("Row " << i << " of the following fermion mixing matrix"
                    " contains entries which have non-zero real and imaginary"
                    " parts:\nZ = " << z);
         }
#endif
      }
   }
}

template<int N>
void convert_symmetric_fermion_mixings_to_slha_rediagonalization(
   Eigen::Array<double, N, 1>& m,
   Eigen::Matrix<std::complex<double>, N, N>& z)
{
   Eigen::Matrix<std::complex<double>, N, N> y =
      z.transpose() * m.matrix().asDiagonal() * z;
#ifdef ENABLE_DEBUG
   if (!is_zero(y.imag().cwiseAbs().maxCoeff())) {
      WARNING("The following symmetric fermion mass matrix contains entries"
          " with non-zero imaginary parts:\nY = " << y);
   }
#endif
   Eigen::Matrix<double, N, N> real_z;
   fs_diagonalize_hermitian(y.real().eval(), m, real_z);
   z = real_z.template cast<std::complex<double> >();
}

#define MEASURE(type,iterations)                                   \
   do {                                                            \
      Stopwatch stopwatch;                                         \
      double time = 0.;                                            \
      Eigen::Matrix<double, 4, 4> mass_matrix;                     \
      Eigen::Matrix<std::complex<double>, 4, 4> z;                 \
      Eigen::Array<double, 4, 1> m;                                \
      mass_matrix << 0, 1, 0, 1,                                   \
                     1, 0, 0, 0,                                   \
                     0, 0, 1, 0,                                   \
                     1, 0, 0, 0;                                   \
      fs_diagonalize_symmetric(mass_matrix, m, z);                 \
      for (int i = 0; i < iterations; i++) {                       \
         stopwatch.start();                                        \
         convert_symmetric_fermion_mixings_to_slha_##type(m,z);    \
         stopwatch.stop();                                         \
         time += stopwatch.get_time_in_seconds();                  \
      }                                                            \
      BOOST_TEST_MESSAGE("conversion via " #type ": " << time << " s"); \
   } while (0)

BOOST_AUTO_TEST_CASE( test_slha_mixing_matrix_conversion_speed )
{
   const int number_of_iterations = 1000000;

   MEASURE(forloop          , number_of_iterations);
   MEASURE(rediagonalization, number_of_iterations);
}
