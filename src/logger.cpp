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

#include "logger.hpp"
#include "config.h"
#include "error.hpp"

namespace flexiblesusy {

namespace {

struct Print_green {
   Print_green() {
      std::cerr << "\033[0;32m";
   }
   ~Print_green() {
      std::cerr << "\033[0m";
   }
};

struct Print_blue {
   Print_blue() {
      std::cerr << "\033[0;34m";
   }
   ~Print_blue() {
      std::cerr << "\033[0m";
   }
};

struct Print_red {
   Print_red() {
      std::cerr << "\033[0;31m";
   }
   ~Print_red() {
      std::cerr << "\033[0m";
   }
};

struct Print_red_bold {
   Print_red_bold() {
      std::cerr << "\033[1;31m";
   }
   ~Print_red_bold() {
      std::cerr << "\033[0m";
   }
};

struct Print_red_background {
   Print_red_background() {
      std::cerr << "\033[41;1;37m";
   }
   ~Print_red_background() {
      std::cerr << "\033[0m";
   }
};

struct Append_endl {
   ~Append_endl() {
      std::cerr << std::endl;;
   }
};

} // anonymous namespace

#ifdef ENABLE_SILENT

void print_verbose(std::function<void()>&&, const char*, int) {}
void print_debug(std::function<void()>&&, const char*, int) {}
void print_info(std::function<void()>&&, const char*, int) {}
void print_warning(std::function<void()>&&, const char*, int) {}
void print_error(std::function<void()>&&, const char*, int) {}
void print_fatal(std::function<void()>&&, const char*, int) {}

#else

#ifdef ENABLE_VERBOSE
void print_verbose(std::function<void()>&& f, const char* /* filename */, int /* line */ )
{
   const auto en = Append_endl();
#ifdef ENABLE_COLORS
   const auto cp = Print_green();
#endif
   f();
}
#else
void print_verbose(std::function<void()>&&, const char*, int) {}
#endif

#ifdef ENABLE_DEBUG
void print_debug(std::function<void()>&& f, const char* /* filename */, int /* line */)
{
   const auto en = Append_endl();
#ifdef ENABLE_COLORS
   const auto cp = Print_blue();
#endif
   f();
}
#else
void print_debug(std::function<void()>&&, const char*, int) {}
#endif

void print_info(std::function<void()>&& f, const char* /* filename */, int /* line */)
{
   const auto en = Append_endl();
   f();
}

void print_warning(std::function<void()>&& f, const char* /* filename */, int /* line */)
{
   const auto en = Append_endl();
#ifdef ENABLE_COLORS
   const auto cp = Print_red();
#endif
   f();
}

void print_error(std::function<void()>&& f, const char* /* filename */, int /* line */)
{
   const auto en = Append_endl();
#ifdef ENABLE_COLORS
   const auto cp = Print_red_bold();
#endif
   f();
}

void print_fatal(std::function<void()>&& f, const char* filename, int line)
{
   const auto en = Append_endl();
#ifdef ENABLE_COLORS
   const auto cp = Print_red_background();
#endif

   std::cerr << "Fatal: (file: " << filename << ", line: " << line << ") ";
   f();
   throw FatalError();
}

#endif

} // namespace flexiblesusy
