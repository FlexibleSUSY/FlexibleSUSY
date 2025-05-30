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

#include "slha_io.hpp"
#include "error.hpp"
#include "ew_input.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "physical_input.hpp"
#include "pmns.hpp"
#include "slhaea.h"
#include "spectrum_generator_settings.hpp"
#include "observables/l_to_l_conversion/settings.hpp"
#include "decays/flexibledecay_settings.hpp"
#include "string_conversion.hpp"
#include "string_format.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace flexiblesusy {

namespace {

int round(double a) noexcept
{
   return static_cast<int>(a >= 0. ? a + 0.5 : a - 0.5);
}

/**
 * fill Modsel struct from given key - value pair
 *
 * @param modsel MODSEL data
 * @param key SLHA key in MODSEL
 * @param value value corresponding to key
 */
void process_modsel_tuple(SLHA_io::Modsel& modsel, int key, double value)
{
   switch (key) {
   case 1: // SUSY breaking model (defined in FlexibleSUSY model file)
   case 3: // SUSY model (defined in SARAH model file)
   case 4: // R-parity violation (defined in SARAH model file)
   case 5: // CP-parity violation (defined in SARAH model file)
   case 11:
   case 21:
      WARNING("Key " << key << " in Block MODSEL currently not supported");
      break;
   case 6: // Flavour violation (defined in SARAH model file)
   {
      const int ivalue = flexiblesusy::round(value);

      if (ivalue < 0 || ivalue > 3) {
         WARNING("Value " << ivalue << " in MODSEL block entry 6 out of range");
      } else {
         const auto uvalue = static_cast<unsigned>(ivalue);
         modsel.quark_flavour_violated = ((uvalue & 0x1u) != 0);
         modsel.lepton_flavour_violated = ((uvalue & 0x2u) != 0);
      }
   }
      break;
   case 12:
      modsel.parameter_output_scale = value;
      break;
   default:
      WARNING("Unrecognized entry in block MODSEL: " << key);
      break;
   }
}

/**
 * fill qedqcd from given key - value pair
 *
 * @param qedqcd low-energy data set
 * @param key SLHA key in SMINPUTS
 * @param value value corresponding to key
 */
void process_sminputs_tuple(softsusy::QedQcd& qedqcd, int key, double value)
{
   switch (key) {
   case 1:
      qedqcd.setAlpha(softsusy::ALPHA, 1.0 / value);
      qedqcd.setAlphaEmInput(1.0 / value);
      break;
   case 2:
      qedqcd.setFermiConstant(value);
      break;
   case 3:
      qedqcd.setAlpha(softsusy::ALPHAS, value);
      qedqcd.setAlphaSInput(value);
      break;
   case 4:
      qedqcd.setPoleMZ(value);
      qedqcd.set_scale(value);
      break;
   case 5:
      qedqcd.setMass(softsusy::mBottom, value);
      qedqcd.setMbMb(value);
      break;
   case 6:
      qedqcd.setPoleMt(value);
      break;
   case 7:
      qedqcd.setMass(softsusy::mTau, value);
      qedqcd.setPoleMtau(value);
      break;
   case 8:
      qedqcd.setNeutrinoPoleMass(3, value);
      break;
   case 9:
      qedqcd.setPoleMW(value);
      break;
   case 11:
      qedqcd.setMass(softsusy::mElectron, value);
      qedqcd.setPoleMel(value);
      break;
   case 12:
      qedqcd.setNeutrinoPoleMass(1, value);
      break;
   case 13:
      qedqcd.setMass(softsusy::mMuon, value);
      qedqcd.setPoleMmuon(value);
      break;
   case 14:
      qedqcd.setNeutrinoPoleMass(2, value);
      break;
   case 21:
      qedqcd.setMass(softsusy::mDown, value);
      qedqcd.setMd2GeV(value);
      break;
   case 22:
      qedqcd.setMass(softsusy::mUp, value);
      qedqcd.setMu2GeV(value);
      break;
   case 23:
      qedqcd.setMass(softsusy::mStrange, value);
      qedqcd.setMs2GeV(value);
      break;
   case 24:
      qedqcd.setMass(softsusy::mCharm, value);
      qedqcd.setMcMc(value);
      break;
   default:
      WARNING("Unrecognized entry in block SMINPUTS: " << key);
      break;
   }
}

void process_flexiblesusy_tuple(Spectrum_generator_settings& settings,
                                int key, double value)
{
   if (0 <= key && key < static_cast<int>(Spectrum_generator_settings::NUMBER_OF_OPTIONS)) {
      settings.set(static_cast<Spectrum_generator_settings::Settings>(key), value);
   } else {
      WARNING("Unrecognized entry in block FlexibleSUSY: " << key);
   }
}

void process_ltolconversion_tuple(LToLConversion_settings& settings,
                                int key, double value)
{
   if (0 <= key && key < static_cast<int>(LToLConversion_settings::NUMBER_OF_OPTIONS)) {
      settings.set(static_cast<LToLConversion_settings::Settings>(key), value);
   } else {
      WARNING("Unrecognized entry in block LToLConversion: " << key);
   }
}

void process_flexibledecay_tuple(FlexibleDecay_settings& settings,
                                int key, double value)
{
   if (0 <= key && key < static_cast<int>(FlexibleDecay_settings::NUMBER_OF_OPTIONS)) {
      settings.set(static_cast<FlexibleDecay_settings::Settings>(key), value);
   } else {
      WARNING("Unrecognized entry in block FlexibleDecay: " << key);
   }
}

void process_flexiblesusyinput_tuple(
   Physical_input& input,
   int key, double value)
{
   if (0 <= key && key < static_cast<int>(Physical_input::NUMBER_OF_INPUT_PARAMETERS)) {
      input.set(static_cast<Physical_input::Input>(key), value);
   } else {
      WARNING("Unrecognized entry in block FlexibleSUSYInput: " << key);
   }
}

/**
 * fill CKM_wolfenstein from given key - value pair
 *
 * @param ckm_wolfenstein Wolfenstein parameters
 * @param key SLHA key in SMINPUTS
 * @param value value corresponding to key
 */
void process_vckmin_tuple(SLHA_io::CKM_wolfenstein& ckm_wolfenstein, int key, double value)
{
   switch (key) {
   case 1:
      ckm_wolfenstein.lambdaW = value;
      break;
   case 2:
      ckm_wolfenstein.aCkm = value;
      break;
   case 3:
      ckm_wolfenstein.rhobar = value;
      break;
   case 4:
      ckm_wolfenstein.etabar = value;
      break;
   default:
      WARNING("Unrecognized entry in block VCKMIN: " << key);
      break;
   }
}

/**
 * fill PMNS_parameters from given key - value pair
 *
 * @param pmns_parameters PMNS matrix parameters
 * @param key SLHA key in SMINPUTS
 * @param value value corresponding to key
 */
void process_upmnsin_tuple(PMNS_parameters& pmns_parameters, int key, double value)
{
   switch (key) {
   case 1:
      pmns_parameters.theta_12 = value;
      break;
   case 2:
      pmns_parameters.theta_23 = value;
      break;
   case 3:
      pmns_parameters.theta_13 = value;
      break;
   case 4:
      pmns_parameters.delta = value;
      break;
   case 5:
      pmns_parameters.alpha_1 = value;
      break;
   case 6:
      pmns_parameters.alpha_2 = value;
      break;
   default:
      WARNING("Unrecognized entry in block UPMNSIN: " << key);
      break;
   }
}

int column_major_index(int r, int c, int /* rows */ , int cols)
{
   return c*cols + r;
}

template <typename T>
struct real_prefix {
   constexpr static const char* const value = "Re(";
};

template <>
struct real_prefix<double> {
   constexpr static const char* const value = "";
};

template <typename T>
struct real_suffix {
   constexpr static const char* const value = ")";
};

template <>
struct real_suffix<double> {
   constexpr static const char* const value = "";
};

} // anonymous namespace

namespace detail {

bool read_scale(const SLHAea::Line& line, double& scale)
{
   if (line.is_block_def() && line.size() > 3 && line[2] == "Q=") {
      scale = to_double(line[3].c_str());
      return true;
   }
   return false;
}

template <typename T>
double read_matrix_(const SLHAea::Coll& data, const std::string& block_name, T* a, int rows, int cols)
{
   auto block = SLHAea::Coll::find(data.cbegin(), data.cend(), block_name);

   double scale = 0.;

   while (block != data.cend()) {
      for (const auto& line: *block) {
         detail::read_scale(line, scale);

         if (line.is_data_line() && line.size() >= 3) {
            const int i = to_int(line[0].c_str()) - 1;
            const int k = to_int(line[1].c_str()) - 1;
            if (0 <= i && i < rows && 0 <= k && k < cols) {
               a[k*cols + i] = to_double(line[2].c_str());
            }
         }
      }

      ++block;
      block = SLHAea::Coll::find(block, data.cend(), block_name);
   }

   return scale;
}

template <typename T>
double read_vector_(const SLHAea::Coll& data, const std::string& block_name, T* a, int len)
{
   auto block = SLHAea::Coll::find(data.cbegin(), data.cend(), block_name);

   double scale = 0.;

   while (block != data.cend()) {
      for (const auto& line: *block) {
         detail::read_scale(line, scale);

         if (line.is_data_line() && line.size() >= 2) {
            const int i = to_int(line[0].c_str()) - 1;
            if (0 <= i && i < len) {
               a[i] = to_double(line[1].c_str());
            }
         }
      }

      ++block;
      block = SLHAea::Coll::find(block, data.cend(), block_name);
   }

   return scale;
}


template <typename T>
std::string format_vector(const std::string& name, const T* a, const std::string& symbol, int rows)
{
   constexpr const char* const prefix = real_prefix<T>::value;
   constexpr const char* const suffix = real_suffix<T>::value;

   std::ostringstream ss;
   ss << name;

   for (int i = 1; i <= rows; ++i) {
      ss << FORMAT_VECTOR(i, std::real(a[i-1]), (prefix + symbol + "(" + flexiblesusy::to_string(i) + ")" + suffix));
   }

   return ss.str();
}


template <typename T>
std::string format_matrix(const std::string& name, const T* a, const std::string& symbol, int rows, int cols)
{
   constexpr const char* const prefix = real_prefix<T>::value;
   constexpr const char* const suffix = real_suffix<T>::value;

   std::ostringstream ss;
   ss << name;

   for (int i = 1; i <= rows; ++i) {
      for (int k = 1; k <= cols; ++k) {
         const int idx = column_major_index(i-1, k-1, rows, cols);
         ss << FORMAT_MIXING_MATRIX(i, k, std::real(a[idx]),
               (prefix + symbol + "(" + flexiblesusy::to_string(i) + "," + flexiblesusy::to_string(k) + ")" + suffix));
      }
   }

   return ss.str();
}


template <typename T>
std::string format_vector_imag(const std::string& name, const T* a, const std::string& symbol, int rows)
{
   std::ostringstream ss;
   ss << name;

   for (int i = 1; i <= rows; ++i) {
      ss << FORMAT_VECTOR(i, std::imag(a[i-1]), ("Im(" + symbol + "(" + flexiblesusy::to_string(i) + "))"));
   }

   return ss.str();
}


template <typename T>
std::string format_matrix_imag(const std::string& name, const T* a, const std::string& symbol, int rows, int cols)
{
   std::ostringstream ss;
   ss << name;

   for (int i = 1; i <= rows; ++i) {
      for (int k = 1; k <= cols; ++k) {
         const int idx = column_major_index(i-1, k-1, rows, cols);
         ss << FORMAT_MIXING_MATRIX(i, k, std::imag(a[idx]),
               ("Im(" + symbol + "(" + flexiblesusy::to_string(i) + "," + flexiblesusy::to_string(k) + "))"));
      }
   }

   return ss.str();
}


} // namespace detail


SLHA_io::SLHA_io()
   : data(std::make_unique<SLHAea::Coll>())
{
}


SLHA_io::SLHA_io(const SLHA_io& other)
   : data(std::make_unique<SLHAea::Coll>(*other.data))
   , modsel(other.modsel)
{
}


SLHA_io::SLHA_io(SLHA_io&& other) noexcept
   : data(std::move(other.data))
   , modsel(std::move(other.modsel))
{
}


SLHA_io::~SLHA_io() = default;


SLHA_io& SLHA_io::operator=(const SLHA_io& other)
{
   return *this = SLHA_io(other);
}


SLHA_io& SLHA_io::operator=(SLHA_io&& other) noexcept
{
   data = std::move(other.data);
   modsel = std::move(other.modsel);
   return *this;
}


void SLHA_io::clear()
{
   data->clear();
   modsel.clear();
}


const SLHAea::Coll& SLHA_io::get_data() const
{
   return *data;
}


void SLHA_io::set_data(const SLHAea::Coll& data_)
{
   data.reset(new SLHAea::Coll(data_));
}


std::string SLHA_io::block_head(const std::string& name, double scale)
{
   const double eps = std::numeric_limits<double>::epsilon();

   std::string result("Block " + name);

   if (!is_zero(scale, eps)) {
      result += " Q= " + FORMAT_SCALE(scale);
   }

   result += '\n';

   return result;
}

bool SLHA_io::block_exists(const std::string& block_name) const
{
   return data->find(block_name) != data->cend();
}

/**
 * @brief reads from source
 *
 * If source is "-", then read_from_stream() is called.  Otherwise,
 * read_from_file() is called.
 *
 * @param source string that specifies the source
 */
void SLHA_io::read_from_source(const std::string& source)
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
void SLHA_io::read_from_file(const std::string& file_name)
{
   std::ifstream ifs(file_name);
   if (ifs.good()) {
      read_from_stream(ifs);
   } else {
      throw ReadError(R"(cannot read SLHA file: ")" + file_name + R"(")");
   }
}

/**
 * @brief clears stored data and reads SLHA data from a stream
 * @param istr input stream
 */
void SLHA_io::read_from_stream(std::istream& istr)
{
   data->clear();
   data->read(istr);
   read_modsel();
}

/**
 * Reads the scale from the line, if the line is a block head and
 * contains a scale definition.  Otherwise, the scale is not read.
 *
 * @param line the line
 * @param scale the scale to write the value to
 *
 * @return true if scale has been read; false otherwise
 */
bool SLHA_io::read_scale(const SLHAea::Line& line, double& scale)
{
   return detail::read_scale(line, scale);
}

void SLHA_io::read_modsel()
{
   Tuple_processor modsel_processor = [this] (int key, double value) {
      return process_modsel_tuple(modsel, key, value);
   };

   read_block("MODSEL", modsel_processor);
}

void SLHA_io::fill(softsusy::QedQcd& qedqcd) const
{
   CKM_wolfenstein ckm_wolfenstein;
   PMNS_parameters pmns_parameters;

   Tuple_processor sminputs_processor = [&qedqcd] (int key, double value) {
      return process_sminputs_tuple(qedqcd, key, value);
   };

   read_block("SMINPUTS", sminputs_processor);

   if (modsel.quark_flavour_violated) {
      Tuple_processor vckmin_processor = [&ckm_wolfenstein] (int key, double value) {
         return process_vckmin_tuple(ckm_wolfenstein, key, value);
      };

      read_block("VCKMIN", vckmin_processor);
   }

   if (modsel.lepton_flavour_violated) {
      Tuple_processor upmnsin_processor = [&pmns_parameters] (int key, double value) {
         return process_upmnsin_tuple(pmns_parameters, key, value);
      };

      read_block("UPMNSIN", upmnsin_processor);
   }

   // fill CKM parameters in qedqcd
   CKM_parameters ckm_parameters;
   ckm_parameters.set_from_wolfenstein(
      ckm_wolfenstein.lambdaW,
      ckm_wolfenstein.aCkm,
      ckm_wolfenstein.rhobar,
      ckm_wolfenstein.etabar);
   qedqcd.setCKM(ckm_parameters);

   // fill PMNS parameters in qedqcd
   qedqcd.setPMNS(pmns_parameters);
}

/**
 * Fill struct of extra physical input parameters from SLHA object
 * (FlexibleSUSYInput block)
 *
 * @param input struct of physical input parameters
 */
void SLHA_io::fill(Physical_input& input) const
{
   Tuple_processor processor = [&input] (int key, double value) {
      return process_flexiblesusyinput_tuple(input, key, value);
   };

   read_block("FlexibleSUSYInput", processor);
}

/**
 * Fill struct of decay settings from SLHA object
 * (LToLConversion block)
 *
 * @param settings struct of decay settings
 */
void SLHA_io::fill(LToLConversion_settings& settings) const
{
   Tuple_processor processor = [&settings] (int key, double value) {
      return process_ltolconversion_tuple(settings, key, value);
   };

   read_block("LToLConversion", processor);
}

/**
 * Fill struct of decay settings from SLHA object
 * (FlexibleDecay block)
 *
 * @param settings struct of decay settings
 */
void SLHA_io::fill(FlexibleDecay_settings& settings) const
{
   Tuple_processor flexibledecay_processor = [&settings] (int key, double value) {
      return process_flexibledecay_tuple(settings, key, value);
   };

   read_block("FlexibleDecay", flexibledecay_processor);
}

/**
 * Fill struct of spectrum generator settings from SLHA object
 * (FlexibleSUSY block)
 *
 * @param settings struct of spectrum generator settings
 */
void SLHA_io::fill(Spectrum_generator_settings& settings) const
{
   Tuple_processor flexiblesusy_processor = [&settings] (int key, double value) {
      return process_flexiblesusy_tuple(settings, key, value);
   };

   read_block("FlexibleSUSY", flexiblesusy_processor);
}

/**
 * Applies processor to each (key, value) pair of a SLHA block.
 * Non-data lines are ignored.
 *
 * @param block_name block name
 * @param processor tuple processor to be applied
 *
 * @return scale (or 0 if no scale is defined)
 */
double SLHA_io::read_block(const std::string& block_name, const Tuple_processor& processor) const
{
   auto block = SLHAea::Coll::find(data->cbegin(), data->cend(), block_name);
   double scale = 0.;

   while (block != data->cend()) {
      for (const auto& line: *block) {
         read_scale(line, scale);

         if (line.is_data_line() && line.size() >= 2) {
            const auto key = to_int(line[0].c_str());
            const auto value = to_double(line[1].c_str());
            processor(key, value);
         }
      }

      ++block;
      block = SLHAea::Coll::find(block, data->cend(), block_name);
   }

   return scale;
}

/**
 * Fills an entry from a SLHA block
 *
 * @param block_name block name
 * @param entry entry to be filled
 *
 * @return scale (or 0 if no scale is defined)
 */
double SLHA_io::read_block(const std::string& block_name, double& entry) const
{
   auto block = SLHAea::Coll::find(data->cbegin(), data->cend(), block_name);
   double scale = 0.;

   while (block != data->cend()) {
      for (const auto& line: *block) {
         read_scale(line, scale);

         if (line.is_data_line()) {
            entry = to_double(line[0].c_str());
         }
      }

      ++block;
      block = SLHAea::Coll::find(block, data->cend(), block_name);
   }

   return scale;
}

double SLHA_io::read_entry(const std::string& block_name, int key) const
{
   auto block = SLHAea::Coll::find(data->cbegin(), data->cend(), block_name);
   double entry = 0.;
   const SLHAea::Block::key_type keys(1, flexiblesusy::to_string(key));

   while (block != data->cend()) {
      auto line = block->find(keys);

      while (line != block->end()) {
         if (line->is_data_line() && line->size() > 1) {
            entry = to_double(line->at(1).c_str());
         }

         ++line;
         line = block->find(line, block->end(), keys);
      }

      ++block;
      block = SLHAea::Coll::find(block, data->cend(), block_name);
   }

   return entry;
}

/**
 * Reads scale definition from SLHA block.
 *
 * @param block_name block name
 *
 * @return scale (or 0 if no scale is defined)
 */
double SLHA_io::read_scale(const std::string& block_name) const
{
   double scale = 0.;
   auto block = SLHAea::Coll::find(data->cbegin(), data->cend(), block_name);

   while (block != data->cend()) {
      for (const auto& line: *block) {
         read_scale(line, scale);
      }

      ++block;
      block = SLHAea::Coll::find(block, data->cend(), block_name);
   }

   return scale;
}

void SLHA_io::set_block(const std::ostringstream& lines, Position position)
{
   set_block(lines.str(), position);
}

void SLHA_io::set_block(const std::string& lines, Position position)
{
   SLHAea::Block block;
   block.str(lines);

   data->erase(block.name());

   if (position == front) {
      data->push_front(std::move(block));
   } else {
      data->push_back(std::move(block));
   }
}

void SLHA_io::set_blocks(const std::vector<std::string>& blocks, Position position)
{
   for (const auto& block: blocks) {
      set_block(block, position);
   }
}

/**
 * This function treats a given scalar as 1x1 matrix.  Such a case is
 * not defined in the SLHA standard, but we still handle it to avoid
 * problems.
 */
void SLHA_io::set_block(const std::string& name, double value,
                        const std::string& symbol, double scale)
{
   std::ostringstream ss;
   ss << block_head(name, scale);
   ss << FORMAT_MIXING_MATRIX(1, 1, value, symbol);

   set_block(ss);
}

void SLHA_io::set_modsel(const Modsel& modsel_)
{
   modsel = modsel_;
   const auto qfv = static_cast<unsigned>(modsel.quark_flavour_violated);
   const auto lfv = 2u * static_cast<unsigned>(modsel.lepton_flavour_violated);

   std::ostringstream ss;
   ss << block_head("MODSEL", 0.0);
   ss << FORMAT_ELEMENT(6 , qfv | lfv, "quark/lepton flavour violation");
   ss << FORMAT_ELEMENT(12, modsel.parameter_output_scale, "running parameter output scale (GeV)");

   set_block(ss);
}

void SLHA_io::set_physical_input(const Physical_input& input)
{
   const auto& names = flexiblesusy::Physical_input::get_names();

   std::ostringstream ss;
   ss << block_head("FlexibleSUSYInput", 0.0);

   for (std::size_t i = 0; i < names.size(); i++) {
      ss << FORMAT_ELEMENT(i, input.get(static_cast<Physical_input::Input>(i)),
                           names[i]);
   }

   set_block(ss);
}

void SLHA_io::set_settings(const Spectrum_generator_settings& settings)
{
   std::ostringstream ss;
   ss << block_head("FlexibleSUSY", 0.0);

   for (int i = 0; i < Spectrum_generator_settings::NUMBER_OF_OPTIONS; i++) {
      ss << FORMAT_ELEMENT(i, settings.get(static_cast<Spectrum_generator_settings::Settings>(i)),
                           settings.get_description(static_cast<Spectrum_generator_settings::Settings>(i)));
   }

   set_block(ss);
}

void SLHA_io::set_LToLConversion_settings(const LToLConversion_settings& settings)
{
   std::ostringstream ss;
   ss << block_head("LToLConversion", 0.0);

   for (int i = 0; i < LToLConversion_settings::NUMBER_OF_OPTIONS; i++) {
      ss << FORMAT_ELEMENT(i, settings.get(static_cast<LToLConversion_settings::Settings>(i)),
                           settings.get_description(static_cast<LToLConversion_settings::Settings>(i)));
   }

   set_block(ss);
}

void SLHA_io::set_FlexibleDecay_settings(const FlexibleDecay_settings& settings)
{
   std::ostringstream ss;
   ss << block_head("FlexibleDecay", 0.0);

   for (int i = 0; i < FlexibleDecay_settings::NUMBER_OF_OPTIONS; i++) {
      ss << FORMAT_ELEMENT(i, settings.get(static_cast<FlexibleDecay_settings::Settings>(i)),
                           settings.get_description(static_cast<FlexibleDecay_settings::Settings>(i)));
   }

   set_block(ss);
}

void SLHA_io::set_sminputs(const softsusy::QedQcd& qedqcd)
{
   std::ostringstream ss;

   ss << block_head("SMINPUTS", 0.0);
   ss << FORMAT_ELEMENT( 1, 1./qedqcd.displayAlphaEmInput()  , "alpha_em^(-1)(MZ) SM(5) MSbar");
   ss << FORMAT_ELEMENT( 2, qedqcd.displayFermiConstant()    , "G_Fermi");
   ss << FORMAT_ELEMENT( 3, qedqcd.displayAlphaSInput()      , "alpha_s(MZ) SM(5) MSbar");
   ss << FORMAT_ELEMENT( 4, qedqcd.displayPoleMZ()           , "MZ(pole)");
   ss << FORMAT_ELEMENT( 5, qedqcd.displayMbMb()             , "mb(mb) SM(5) MSbar");
   ss << FORMAT_ELEMENT( 6, qedqcd.displayPoleMt()           , "Mtop(pole)");
   ss << FORMAT_ELEMENT( 7, qedqcd.displayPoleMtau()         , "Mtau(pole)");
   ss << FORMAT_ELEMENT( 8, qedqcd.displayNeutrinoPoleMass(3), "Mv3(pole)");
   ss << FORMAT_ELEMENT( 9, qedqcd.displayPoleMW()           , "MW(pole)");
   ss << FORMAT_ELEMENT(11, qedqcd.displayPoleMel()          , "Melectron(pole)");
   ss << FORMAT_ELEMENT(12, qedqcd.displayNeutrinoPoleMass(1), "Mv1(pole)");
   ss << FORMAT_ELEMENT(13, qedqcd.displayPoleMmuon()        , "Mmuon(pole)");
   ss << FORMAT_ELEMENT(14, qedqcd.displayNeutrinoPoleMass(2), "Mv2(pole)");
   ss << FORMAT_ELEMENT(21, qedqcd.displayMd2GeV()           , "md(2GeV)");
   ss << FORMAT_ELEMENT(22, qedqcd.displayMu2GeV()           , "mu(2GeV)");
   ss << FORMAT_ELEMENT(23, qedqcd.displayMs2GeV()           , "ms(2GeV)");
   ss << FORMAT_ELEMENT(24, qedqcd.displayMcMc()             , "mc(mc) SM(4) MSbar");

   set_block(ss);
}

void SLHA_io::set_unitarity_infinite_s(
   Spectrum_generator_settings const& spectrum_generator_settings, UnitarityInfiniteS const& unitarity)
{
   if (spectrum_generator_settings.get(Spectrum_generator_settings::calculate_observables)) {
      std::ostringstream block;
      block << "Block FlexibleSUSYUnitarity Q= " << FORMAT_SCALE(unitarity.renScale) << '\n'
            << FORMAT_ELEMENT(0, unitarity.allowed, "Tree-level unitarity limits fulfilled or not")
            << FORMAT_ELEMENT(1, unitarity.maxAbsReEigenval, "max(|re(eigenvalues(a0))|)");
      set_block(block);
   }
}

void SLHA_io::write_to_file(const std::string& file_name) const
{
   std::ofstream ofs(file_name);
   write_to_stream(ofs);
}

void SLHA_io::write_to_stream(std::ostream& ostr) const
{
   if (ostr.good()) {
      ostr << *data;
   } else {
      ERROR("cannot write SLHA file");
   }
}


void SLHA_io::write_to_stream() const
{
   write_to_stream(std::cerr);
}


double SLHA_io::read_vector(const std::string& block_name, double* a, int len) const
{
   return detail::read_vector_(*data, block_name, a, len);
}


double SLHA_io::read_vector(const std::string& block_name, std::complex<double>* a, int len) const
{
   return detail::read_vector_(*data, block_name, a, len);
}


double SLHA_io::read_matrix(const std::string& block_name, double* a, int rows, int cols) const
{
   return detail::read_matrix_(*data, block_name, a, rows, cols);
}


double SLHA_io::read_matrix(const std::string& block_name, std::complex<double>* a, int rows, int cols) const
{
   return detail::read_matrix_(*data, block_name, a, rows, cols);
}


void SLHA_io::set_vector(const std::string& name, const double* a, const std::string& symbol, double scale, int rows)
{
   set_block(detail::format_vector(block_head(name, scale), a, symbol, rows));
}


void SLHA_io::set_vector(const std::string& name, const std::complex<double>* a, const std::string& symbol, double scale, int rows)
{
   set_block(detail::format_vector(block_head(name, scale), a, symbol, rows));
}


void SLHA_io::set_matrix(const std::string& name, const double* a, const std::string& symbol, double scale, int rows, int cols)
{
   set_block(detail::format_matrix(block_head(name, scale), a, symbol, rows, cols));
}


void SLHA_io::set_matrix(const std::string& name, const std::complex<double>* a, const std::string& symbol, double scale, int rows, int cols)
{
   set_block(detail::format_matrix(block_head(name, scale), a, symbol, rows, cols));
}


void SLHA_io::set_vector_imag(const std::string& name, const double* a, const std::string& symbol, double scale, int rows)
{
   set_block(detail::format_vector_imag(block_head(name, scale), a, symbol, rows));
}


void SLHA_io::set_vector_imag(const std::string& name, const std::complex<double>* a, const std::string& symbol, double scale, int rows)
{
   set_block(detail::format_vector_imag(block_head(name, scale), a, symbol, rows));
}


void SLHA_io::set_matrix_imag(const std::string& name, const double* a, const std::string& symbol, double scale, int rows, int cols)
{
   set_block(detail::format_matrix_imag(block_head(name, scale), a, symbol, rows, cols));
}


void SLHA_io::set_matrix_imag(const std::string& name, const std::complex<double>* a, const std::string& symbol, double scale, int rows, int cols)
{
   set_block(detail::format_matrix_imag(block_head(name, scale), a, symbol, rows, cols));
}

void SLHA_io::set_hs_or_lilith(std::string const& block_name, const std::size_t ndof, const double chi2, const double chi2SMmin, const double mhSM, const double pval)
{
   std::ostringstream ss;

   ss << block_head(block_name, 0.0);
   ss << FORMAT_ELEMENT(1, ndof, "number of degrees of freedom");
   ss << FORMAT_ELEMENT(2, chi2, "𝜒²");
   ss << FORMAT_ELEMENT(3, chi2SMmin, "SM 𝜒² for mh = " + std::to_string(mhSM) + " GeV");
   // SLHA doesn't print nicelly numbers with 3 digit exponent
   ss << FORMAT_ELEMENT(4, pval > 1e-100 ? pval : 0., "p-value");

   set_block(ss);
}

void SLHA_io::set_higgsbounds(std::vector<std::tuple<int, double, double, std::string>> const& v)
{
   std::ostringstream ss;

   ss << block_head("HIGGSBOUNDS", 0.0);
   for (auto const& el : v) {
      ss << FORMAT_MIXING_MATRIX(std::get<0>(el), 1, std::get<1>(el), std::get<3>(el));
      ss << FORMAT_MIXING_MATRIX(std::get<0>(el), 2, std::get<2>(el), "expRatio");
   }

   set_block(ss);
}

void SLHA_io::set_effectivecouplings_block(const std::vector<std::tuple<int, int, int, double, std::string>>& effCouplings)
{
   std::ostringstream decay;
   decay << "Block EFFHIGGSCOUPLINGS\n";

   for (auto const& effC : effCouplings) {
      decay << FORMAT_EFFECTIVECOUPLINGS(std::get<0>(effC),  std::get<1>(effC), std::get<2>(effC), std::get<3>(effC), std::get<4>(effC));
   }

   set_block(decay);
}

#define DECAY_FERMION_RE(PDG1, PDG2, CHANNEL) (ss << (!effC.CHANNEL.first.empty() ? FORMAT_EFFECTIVECOUPLINGS(effC.pdgid, PDG1,  PDG2, std::real(effC.CHANNEL.second), effC.CHANNEL.first + "/SM with mhSM = m" + effC.particle) : ""))
#define DECAY_VBOSON(PDG1, PDG2, CHANNEL) (ss << (!effC.CHANNEL.first.empty() ? FORMAT_EFFECTIVECOUPLINGS(effC.pdgid, PDG1,  PDG2, effC.CHANNEL.second, effC.CHANNEL.first + "/SM with mhSM = m" + effC.particle) : ""))

void SLHA_io::set_normalized_effectivecouplings_block(const EffectiveCoupling_list& effCouplings) {
   std::ostringstream ss;
   ss << "Block NORMALIZEDEFFHIGGSCOUPLINGS\n";
   for (auto const& effC : effCouplings) {
      ss << FORMAT_EFFECTIVECOUPLINGS(effC.pdgid, 0,  0, effC.width_sm, "SM Higgs width for mhSM = m" + effC.particle);
      if (effC.CP != -1) {
         DECAY_FERMION_RE(-1, 1, uu);
         DECAY_FERMION_RE(-2, 2, dd);
         DECAY_FERMION_RE(-3, 3, ss);
         DECAY_FERMION_RE(-4, 4, cc);
         DECAY_FERMION_RE(-5, 5, bb);
         DECAY_FERMION_RE(-6, 6, tt);
         DECAY_FERMION_RE(-11, 11, ee);
         DECAY_FERMION_RE(-13, 13, mumu);
         DECAY_FERMION_RE(-15, 15, tautau);

         DECAY_VBOSON(-24, 24, WW);
         DECAY_VBOSON(23, 23, ZZ);
      }

      DECAY_VBOSON(21, 21, gg);
      DECAY_VBOSON(22, 22, gamgam);
      DECAY_VBOSON(23, 22, Zgam);
   }

   set_block(ss);
}

#define DECAY_FERMION_IM(PDG1, PDG2, CHANNEL) (ss << (!effC.CHANNEL.first.empty() ? FORMAT_EFFECTIVECOUPLINGS(effC.pdgid, PDG1,  PDG2, std::imag(effC.CHANNEL.second), effC.CHANNEL.first + "/SM with mhSM = m" + effC.particle) : ""))

void SLHA_io::set_imnormalized_effectivecouplings_block(const EffectiveCoupling_list& effCouplings) {
   std::ostringstream ss;
   ss << "Block IMNORMALIZEDEFFHIGGSCOUPLINGS\n";
   for (auto const& effC : effCouplings) {
      if (effC.CP == 1) continue;
      DECAY_FERMION_IM(-1, 1, uu);
      DECAY_FERMION_IM(-2, 2, dd);
      DECAY_FERMION_IM(-3, 3, ss);
      DECAY_FERMION_IM(-4, 4, cc);
      DECAY_FERMION_IM(-5, 5, bb);
      DECAY_FERMION_IM(-6, 6, tt);
      DECAY_FERMION_IM(-11, 11, ee);
      DECAY_FERMION_IM(-13, 13, mumu);
      DECAY_FERMION_IM(-15, 15, tautau);
   }

   set_block(ss);
}

} // namespace flexiblesusy
