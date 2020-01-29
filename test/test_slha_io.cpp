
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_slha_io

#include <boost/test/unit_test.hpp>

#include "slha_io.hpp"
#include "linalg2.hpp"
#include "stopwatch.hpp"
#include <Eigen/Core>
#include <boost/lexical_cast.hpp>

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

   const std::string str = "Block Matrix Q= 1234.56\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::MatrixXd matrix(Eigen::MatrixXd::Zero(2,2));
   const double scale = reader.read_scale("Matrix");

   BOOST_CHECK_EQUAL(scale, 1234.56);
}

BOOST_AUTO_TEST_CASE( test_read_scale_from_block )
{
   SLHAea::Coll coll;
   SLHAea::Block block;

   const std::string str = "Block Matrix Q= 1234.56\n"
      "   1  1  1.0      # element 1,1\n"
      "   1  2  2.0      # element 1,2\n"
      "   2  1  3.0      # element 2,1\n"
      "   2  2  4.0      # element 2,2\n";

   block.str(str);
   coll.push_back(block);

   SLHA_io reader;
   reader.set_data(coll);

   Eigen::MatrixXd matrix(Eigen::MatrixXd::Zero(2,2));
   const double scale = reader.read_block("Matrix", matrix);

   BOOST_CHECK_EQUAL(scale, 1234.56);
   BOOST_CHECK_EQUAL(matrix(0,0), 1.0);
   BOOST_CHECK_EQUAL(matrix(0,1), 2.0);
   BOOST_CHECK_EQUAL(matrix(1,0), 3.0);
   BOOST_CHECK_EQUAL(matrix(1,1), 4.0);
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
      str += " Q= " + boost::lexical_cast<std::string>(scale);
   str += '\n';

   for (int i = 0; i < number_of_entries; i++) {
      const std::string key(boost::lexical_cast<std::string>(i));
      const std::string num(boost::lexical_cast<std::string>(i + offset));
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
   const int number_of_entries = 10000;
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

   BOOST_CHECK_LT(100. * processor_time, loop_time);
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
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(eigenvalues, Z);

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
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(eigenvalues, Z);

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
   const int number_of_iterations = 10000000;

   MEASURE(forloop          , number_of_iterations);
   MEASURE(rediagonalization, number_of_iterations);
}
