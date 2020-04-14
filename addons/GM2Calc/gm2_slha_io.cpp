// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "gm2_slha_io.hpp"

#include "gm2calc/MSSMNoFV_onshell.hpp"

#include "gm2_config_options.hpp"
#include "gm2_log.hpp"
#include "gm2_numerics.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include <boost/lexical_cast.hpp>
#include <Eigen/Core>

namespace gm2calc {

namespace {

   struct HMIX_data {
      double mu{0.0};
      double tanb{0.0};
      double v{0.0};
      double mA2{0.0};
   };

   struct GM2CalcInput_data {
      double alpha_MZ{0.0};
      double alpha_thompson{0.0};
   };

   void process_gm2calcconfig_tuple(Config_options& /*config_options*/, int /*key*/, double /*value*/);
   void process_gm2calcinput_tuple(GM2CalcInput_data& /*data*/, int /*key*/, double /*value*/);
   void process_gm2calcinput_tuple(MSSMNoFV_onshell& /*model*/, int /*key*/, double /*value*/);
   void process_sminputs_tuple(MSSMNoFV_onshell& /*model*/, int /*key*/, double /*value*/);
   void process_hmix_tuple(HMIX_data& /*data*/, int /*key*/, double /*value*/);
   void process_mass_tuple(MSSMNoFV_onshell_physical& /*physical*/, int /*key*/, double /*value*/);
   void process_msoft_tuple(MSSMNoFV_onshell& /*model*/, int /*key*/, double /*value*/);

} // anonymous namespace

/**
 * @brief reads from source
 *
 * If source is "-", then read_from_stream() is called.  Otherwise,
 * read_from_file() is called.
 *
 * @param source string that specifies the source
 */
void GM2_slha_io::read_from_source(const std::string& source)
{
   if (source == "-") {
      read_from_stream(std::cin);
   } else {
      read_from_file(source);
   }
}

/**
 * @brief opens SLHA input file and reads the content
 * @param file_name SLHA input file name
 */
void GM2_slha_io::read_from_file(const std::string& file_name)
{
   std::ifstream ifs(file_name);
   if (ifs.good()) {
      data.clear();
      data.read(ifs);
   } else {
      throw EReadError("cannot read input file: \"" + file_name + "\"");
   }
}

/**
 * @brief reads SLHA data from a stream
 * @param istr input stream
 */
void GM2_slha_io::read_from_stream(std::istream& istr)
{
   data.read(istr);
}

/**
 * Reads scale definition from an SLHA block.
 *
 * @param block block
 *
 * @return scale (or 0 if no scale is defined)
 */
double GM2_slha_io::read_scale(const SLHAea::Block& block)
{
   double scale = 0.0;

   for (const auto& line : block) {
      // read scale from block definition
      if (line.is_block_def() && line.size() > 3 && line[2] == "Q=") {
         scale = convert_to<double>(line[3]);
      }
   }

   return scale;
}

/**
 * Reads scale definition from SLHA block.
 *
 * @param block_name block name
 *
 * @return scale (or 0 if no scale is defined)
 */
double GM2_slha_io::read_scale(const std::string& block_name) const
{
   double scale = 0.;
   auto block = SLHAea::Coll::find(data.cbegin(), data.cend(), block_name);

   while (block != data.cend()) {
      scale = GM2_slha_io::read_scale(*block);
      ++block;
      block = SLHAea::Coll::find(block, data.cend(), block_name);
   }

   return scale;
}

/**
 * Returns true if the block scale after Q= matches \a scale, false
 * otherwise.  If scale == 0, the functions returns true.
 *
 * @param block SLHA block
 * @param scale scale
 * @param eps absolute tolerance to treat two scales being the same
 */
bool GM2_slha_io::is_at_scale(const SLHAea::Block& block, double scale, double eps)
{
   if (is_zero(scale, std::numeric_limits<double>::epsilon())) {
      return true;
   }

   const auto block_scale = GM2_slha_io::read_scale(block);

   return is_equal(scale, block_scale, eps);
}

/**
 * Applies processor to each (key, value) pair of a SLHA block.
 * Non-data lines are ignored.
 *
 * @param block the block
 * @param processor tuple processor to be applied
 */
void GM2_slha_io::read_block(const SLHAea::Block& block, const Tuple_processor& processor)
{
   for (const auto& line : block) {
      if (line.is_data_line() && line.size() >= 2) {
         const auto key = convert_to<int>(line[0]);
         const auto value = convert_to<double>(line[1]);
         processor(key, value);
      }
   }
}

/**
 * Applies processor to each (key, value) pair of a SLHA block.
 * Non-data lines are ignored.
 *
 * @param block_name block name
 * @param processor tuple processor to be applied
 * @param scale (or 0 if scale should be ignored)
 */
void GM2_slha_io::read_block(const std::string& block_name,
                             const Tuple_processor& processor,
                             double scale) const
{
   auto block = SLHAea::Coll::find(data.cbegin(), data.cend(), block_name);

   while (block != data.cend()) {
      if (is_at_scale(*block, scale)) {
         GM2_slha_io::read_block(*block, processor);
      }

      ++block;
      block = SLHAea::Coll::find(block, data.cend(), block_name);
   }
}

void GM2_slha_io::write_to_file(const std::string& file_name)
{
   std::ofstream ofs(file_name);
   write_to_stream(ofs);
}

void GM2_slha_io::write_to_stream(std::ostream& ostr)
{
   if (ostr.good()) {
      ostr << data;
   } else {
      ERROR("cannot write SLHA file");
   }
}

/**
 * Fills a block entry with a value.  If the block or the entry do not
 * exist, the block / entry is created.
 *
 * @param block_name block name
 * @param entry number of the entry
 * @param value value
 * @param description comment
 */
void GM2_slha_io::fill_block_entry(const std::string& block_name,
                                   unsigned entry, double value,
                                   const std::string& description)
{
   auto block = SLHAea::Coll::find(data.cbegin(), data.cend(), block_name);

   if (block == data.cend()) {
      SLHAea::Block block;
      block.str("Block " + block_name);
      data.push_back(block);
   }

   data[block_name][entry] = (FORMAT_ELEMENT(entry, value, description)).str();
}

/**
 * Fills a block entry with a string.  If the block or the entry do
 * not exist, the block / entry is created.
 *
 * @param block_name block name
 * @param entry number of the entry
 * @param description comment
 */
void GM2_slha_io::fill_block_entry(const std::string& block_name,
                                   unsigned entry,
                                   const std::string& description)
{
   auto block = SLHAea::Coll::find(data.cbegin(), data.cend(), block_name);

   if (block == data.cend()) {
      SLHAea::Block block;
      block.str("Block " + block_name);
      data.push_front(block);
   }

   data[block_name][entry] = (FORMAT_SPINFO(entry, description)).str();
}

void GM2_slha_io::fill_from_msoft(MSSMNoFV_onshell& model) const
{
   GM2_slha_io::Tuple_processor processor = [&model] (int key, double value) {
      return process_msoft_tuple(model, key, value);
   };

   read_block("MSOFT", processor, model.get_scale());
}

void GM2_slha_io::fill_from_A(MSSMNoFV_onshell& model) const
{
   const double scale = model.get_scale();

   {
      Eigen::Matrix<double,3,3> Ae(Eigen::Matrix<double,3,3>::Zero());
      read_block("AE", Ae, scale);
      model.set_Ae(Ae);
   }
   {
      Eigen::Matrix<double,3,3> Au(Eigen::Matrix<double,3,3>::Zero());
      read_block("AU", Au, scale);
      model.set_Au(Au);
   }
   {
      Eigen::Matrix<double,3,3> Ad(Eigen::Matrix<double,3,3>::Zero());
      read_block("AD", Ad, scale);
      model.set_Ad(Ad);
   }
}

void GM2_slha_io::fill_from_hmix(MSSMNoFV_onshell& model) const
{
   HMIX_data hmix;

   GM2_slha_io::Tuple_processor processor = [&hmix] (int key, double value) {
      return process_hmix_tuple(hmix, key, value);
   };

   read_block("HMIX", processor, model.get_scale());

   const double tanb = hmix.tanb;
   const double scb = tanb / (1 + tanb*tanb); // sin(beta)*cos(beta)

   model.set_Mu(hmix.mu);
   model.set_TB(hmix.tanb);
   model.set_BMu(hmix.mA2 * scb);
}

void GM2_slha_io::fill_scale(MSSMNoFV_onshell& model) const
{
   const double eps = std::numeric_limits<double>::epsilon();
   const double scale = read_scale("HMIX");

   if (is_zero(scale, eps)) {
      throw EInvalidInput("Could not determine renormalization scale"
                          " from HMIX block");
   }

   model.set_scale(scale);
}

void GM2_slha_io::fill_from_sminputs(MSSMNoFV_onshell& model) const
{
   GM2_slha_io::Tuple_processor processor = [&model] (int key, double value) {
      return process_sminputs_tuple(model, key, value);
   };

   read_block("SMINPUTS", processor);
}

void GM2_slha_io::fill_from_mass(MSSMNoFV_onshell_physical& physical) const
{
   GM2_slha_io::Tuple_processor processor = [&physical] (int key, double value) {
      return process_mass_tuple(physical, key, value);
   };

   read_block("MASS", processor);
   read_block("NMIX", physical.ZN);
   read_block("SMUMIX", physical.ZM);

   physical.convert_to_hk();
}

void GM2_slha_io::fill_alpha_from_gm2calcinput(MSSMNoFV_onshell& model) const
{
   GM2CalcInput_data data;

   GM2_slha_io::Tuple_processor processor = [&data] (int key, double value) {
      return process_gm2calcinput_tuple(data, key, value);
   };

   read_block("GM2CalcInput", processor);

   if (data.alpha_MZ > std::numeric_limits<double>::epsilon()) {
      model.set_alpha_MZ(data.alpha_MZ);
   }

   if (data.alpha_thompson > std::numeric_limits<double>::epsilon()) {
      model.set_alpha_thompson(data.alpha_thompson);
   }
}

/**
 * Reads the GM2CalcInput block and fills the model parameter class.
 *
 * This function assumes that MW(pole) and MZ(pole) are non-zero.
 */
void GM2_slha_io::fill_from_gm2calcinput(MSSMNoFV_onshell& model) const
{
   GM2_slha_io::Tuple_processor processor = [&model] (int key, double value) {
      return process_gm2calcinput_tuple(model, key, value);
   };

   read_block("GM2CalcInput", processor);
}

/**
 * Reads model parameters in GM2Calc format from GM2CalcInput and
 * SMINPUTS blocks
 *
 * @param model model
 */
void GM2_slha_io::fill_gm2calc(MSSMNoFV_onshell& model) const
{
   fill_from_sminputs(model);
   fill_from_gm2calcinput(model);
}

/**
 * Reads model parameters in SLHA format (from SLHA and GM2CalcInput
 * input blocks)
 *
 * @param model model
 */
void GM2_slha_io::fill_slha(MSSMNoFV_onshell& model) const
{
   // read all pole masses (including MW) from SMINPUTS
   fill_from_sminputs(model);
   fill_from_mass(model.get_physical());
   fill_scale(model);
   fill_from_hmix(model);
   fill_from_A(model);
   fill_from_msoft(model);
   fill_alpha_from_gm2calcinput(model);
}

/**
 * Reads configuration from GM2CalcConfig block
 *
 * @param config_options configuration settings
 */
void GM2_slha_io::fill(Config_options& config_options) const
{
   GM2_slha_io::Tuple_processor processor = [&config_options] (int key, double value) {
      return process_gm2calcconfig_tuple(config_options, key, value);
   };

   read_block("GM2CalcConfig", processor);
}

namespace {

bool is_integer(double value)
{
   double intpart;
   return std::modf(value, &intpart) == 0.0;
}

template <class Source>
std::string to_string(Source arg)
{
   return boost::lexical_cast<std::string>(arg);
}

void read_bool(double value, bool& result, const char* error_msg)
{
   if (value == 0.0 || value == 1.0) {
      result = (value != 0.0);
   } else {
      throw EInvalidInput(std::string(error_msg) + ": " +
                          gm2calc::to_string(value) +
                          " (allowed values: 0 or 1)");
   }
}

template <typename T>
void read_integer(double value, T& result, T min, T max, const char* error_msg)
{
   if (is_integer(value) && value >= min && value <= max) {
      result = static_cast<T>(value);
   } else {
      throw EInvalidInput(
         std::string(error_msg) + ": " + gm2calc::to_string(value) +
         " (allowed integer values: " + gm2calc::to_string(min) + ",...," +
         gm2calc::to_string(max) + ")");
   }
}

void read_double_non_zero(double value, double& result)
{
   const double eps = std::numeric_limits<double>::epsilon();

   if (!is_zero(value, eps)) {
      result = value;
   }
}

void process_gm2calcconfig_tuple(
   Config_options& config_options, int key, double value)
{
   switch (key) {
   case 0:
      read_integer(value,
                   config_options.output_format,
                   static_cast<Config_options::E_output_format>(0),
                   static_cast<Config_options::E_output_format>(
                      Config_options::NUMBER_OF_OUTPUT_FORMATS - 1),
                   "unsupported output format in GM2CalcConfig[0]");
      break;
   case 1:
      read_integer(value, config_options.loop_order, 0u, 2u,
                   "unsupported loop order in GM2CalcConfig[1]");
      break;
   case 2:
      read_bool(value, config_options.tanb_resummation,
                "unsupported tan(beta) resummation flag value in GM2CalcConfig[2]");
      break;
   case 3:
      read_bool(value, config_options.force_output,
                "unsupported force output flag value in GM2CalcConfig[3]");
      break;
   case 4:
      read_bool(value, config_options.verbose_output,
                "unsupported verbose output flag value in GM2CalcConfig[4]");
      break;
   case 5:
      read_bool(value, config_options.calculate_uncertainty,
                "unsupported uncertainty flag value in GM2CalcConfig[5]");
      break;
   default:
      WARNING("Unrecognized entry in block GM2CalcConfig: " << key);
      break;
   }
}

void process_gm2calcinput_tuple(
   MSSMNoFV_onshell& model, int key, double value)
{
   switch (key) {
   case  0: model.set_scale(value);          break;
   case  1: model.set_alpha_MZ(value);       break;
   case  2: model.set_alpha_thompson(value); break;
   case  3: model.set_TB(    value);         break;
   case  4: model.set_Mu(    value);         break;
   case  5: model.set_MassB( value);         break;
   case  6: model.set_MassWB(value);         break;
   case  7: model.set_MassG( value);         break;
   case  8: model.set_MA0(   value);         break;
   case  9: model.set_ml2(0, 0, signed_sqr(value)); break;
   case 10: model.set_ml2(1, 1, signed_sqr(value)); break;
   case 11: model.set_ml2(2, 2, signed_sqr(value)); break;
   case 12: model.set_me2(0, 0, signed_sqr(value)); break;
   case 13: model.set_me2(1, 1, signed_sqr(value)); break;
   case 14: model.set_me2(2, 2, signed_sqr(value)); break;
   case 15: model.set_mq2(0, 0, signed_sqr(value)); break;
   case 16: model.set_mq2(1, 1, signed_sqr(value)); break;
   case 17: model.set_mq2(2, 2, signed_sqr(value)); break;
   case 18: model.set_mu2(0, 0, signed_sqr(value)); break;
   case 19: model.set_mu2(1, 1, signed_sqr(value)); break;
   case 20: model.set_mu2(2, 2, signed_sqr(value)); break;
   case 21: model.set_md2(0, 0, signed_sqr(value)); break;
   case 22: model.set_md2(1, 1, signed_sqr(value)); break;
   case 23: model.set_md2(2, 2, signed_sqr(value)); break;
   case 24: model.set_Ae( 0, 0, value); break;
   case 25: model.set_Ae( 1, 1, value); break;
   case 26: model.set_Ae( 2, 2, value); break;
   case 27: model.set_Ad( 0, 0, value); break;
   case 28: model.set_Ad( 1, 1, value); break;
   case 29: model.set_Ad( 2, 2, value); break;
   case 30: model.set_Au( 0, 0, value); break;
   case 31: model.set_Au( 1, 1, value); break;
   case 32: model.set_Au( 2, 2, value); break;
   default:
      WARNING("Unrecognized entry in block GM2CalcInput: " << key);
      break;
   }
}

void process_gm2calcinput_tuple(
   GM2CalcInput_data& data, int key, double value)
{
   switch (key) {
   case  0: /* scale */                  break;
   case  1: data.alpha_MZ = value;       break;
   case  2: data.alpha_thompson = value; break;
   default:
      break;
   }
}

void process_sminputs_tuple(
   MSSMNoFV_onshell& model, int key, double value)
{
   const double Pi = 3.14159265358979323846;
   MSSMNoFV_onshell_physical& physical = model.get_physical();

   switch (key) {
   case  1: /* alpha_em(MZ) */      break;
   case  2: /* G_F */               break;
   case  3: model.set_g3(std::sqrt(4*Pi*value)); break;
   case  4: physical.MVZ = value;   break;
   case  5: physical.MFb = value;   break;
   case  6: physical.MFt = value;   break;
   case  7: physical.MFtau = value; break;
   case  8: physical.MFvt = value;  break;
   case  9: physical.MVWm = value;  break;
   case 11: physical.MFe = value;   break;
   case 12: physical.MFve = value;  break;
   case 13: physical.MFm = value;   break;
   case 14: physical.MFvm = value;  break;
   case 21: physical.MFd = value;   break;
   case 23: physical.MFs = value;   break;
   case 22: physical.MFu = value;   break;
   case 24: physical.MFc = value;   break;
   default:
      WARNING("Unrecognized entry in block SMINPUTS: " << key);
      break;
   }
}

void process_hmix_tuple(
   HMIX_data& data, int key, double value)
{
   switch (key) {
   case 1: data.mu = value  ; break;
   case 2: data.tanb = value; break;
   case 3: data.v = value   ; break;
   case 4: data.mA2 = value ; break;
   default:
      WARNING("Unrecognized entry in block HMIX: " << key);
      break;
   }
}

void process_msoft_tuple(
   MSSMNoFV_onshell& model, int key, double value)
{
   switch (key) {
   case 21: model.set_mHd2(value)                 ; break;
   case 22: model.set_mHu2(value)                 ; break;
   case 31: model.set_ml2(0, 0, signed_sqr(value)); break;
   case 32: model.set_ml2(1, 1, signed_sqr(value)); break;
   case 33: model.set_ml2(2, 2, signed_sqr(value)); break;
   case 34: model.set_me2(0, 0, signed_sqr(value)); break;
   case 35: model.set_me2(1, 1, signed_sqr(value)); break;
   case 36: model.set_me2(2, 2, signed_sqr(value)); break;
   case 41: model.set_mq2(0, 0, signed_sqr(value)); break;
   case 42: model.set_mq2(1, 1, signed_sqr(value)); break;
   case 43: model.set_mq2(2, 2, signed_sqr(value)); break;
   case 44: model.set_mu2(0, 0, signed_sqr(value)); break;
   case 45: model.set_mu2(1, 1, signed_sqr(value)); break;
   case 46: model.set_mu2(2, 2, signed_sqr(value)); break;
   case 47: model.set_md2(0, 0, signed_sqr(value)); break;
   case 48: model.set_md2(1, 1, signed_sqr(value)); break;
   case 49: model.set_md2(2, 2, signed_sqr(value)); break;
   case  1: model.set_MassB(               value) ; break;
   case  2: model.set_MassWB(              value) ; break;
   case  3: model.set_MassG(               value) ; break;
   default:
      WARNING("Unrecognized entry in block MSOFT: " << key);
      break;
   }
}

void process_mass_tuple(
   MSSMNoFV_onshell_physical& physical, int key, double value)
{
   switch (key) {
   case 1000012: physical.MSveL = value;    break;
   case 1000014: physical.MSvmL = value;    break;
   case 1000016: physical.MSvtL = value;    break;
   case 1000001: physical.MSd(0) = value;   break;
   case 2000001: physical.MSd(1) = value;   break;
   case 1000002: physical.MSu(0) = value;   break;
   case 2000002: physical.MSu(1) = value;   break;
   case 1000011: physical.MSe(0) = value;   break;
   case 2000011: physical.MSe(1) = value;   break;
   case 1000013: physical.MSm(0) = value;   break;
   case 2000013: physical.MSm(1) = value;   break;
   case 1000015: physical.MStau(0) = value; break;
   case 2000015: physical.MStau(1) = value; break;
   case 1000003: physical.MSs(0) = value;   break;
   case 2000003: physical.MSs(1) = value;   break;
   case 1000004: physical.MSc(0) = value;   break;
   case 2000004: physical.MSc(1) = value;   break;
   case 1000005: physical.MSb(0) = value;   break;
   case 2000005: physical.MSb(1) = value;   break;
   case 1000006: physical.MSt(0) = value;   break;
   case 2000006: physical.MSt(1) = value;   break;
   case 24     : read_double_non_zero(value, physical.MVWm); break;
   case 25     : physical.Mhh(0) = value;   break;
   case 35     : physical.Mhh(1) = value;   break;
   case 36     : physical.MAh(1) = value;   break;
   case 37     : physical.MHpm(1) = value;  break;
   case 1000021: physical.MGlu = value;     break;
   case 1000022: physical.MChi(0) = value;  break;
   case 1000023: physical.MChi(1) = value;  break;
   case 1000025: physical.MChi(2) = value;  break;
   case 1000035: physical.MChi(3) = value;  break;
   case 1000024: physical.MCha(0) = value;  break;
   case 1000037: physical.MCha(1) = value;  break;
   default:
      break;
   }
}

} // anonymous namespace

} // namespace gm2calc
