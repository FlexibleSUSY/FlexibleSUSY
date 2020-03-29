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

#ifndef GM2_ERROR_HPP
#define GM2_ERROR_HPP

#include <stdexcept>
#include <string>

namespace gm2calc {

class Error : public std::runtime_error {
public:
   explicit Error(const char* msg) : std::runtime_error(msg) {}
   explicit Error(const std::string& msg) : std::runtime_error(msg) {}
};

/**
 * @class ESetupError
 * @brief Spectrum generator was not setup correctly
 */
class ESetupError : public Error {
public:
   explicit ESetupError(const char* msg) : Error(msg) {}
   explicit ESetupError(const std::string& msg) : Error(msg) {}
};

class EInvalidInput : public Error {
public:
   explicit EInvalidInput(const char* msg) : Error(msg) {}
   explicit EInvalidInput(const std::string& msg) : Error(msg) {}
};

class EPhysicalProblem : public Error {
public:
   explicit EPhysicalProblem(const char* msg) : Error(msg) {}
   explicit EPhysicalProblem(const std::string& msg) : Error(msg) {}
};

class EReadError : public Error {
public:
   explicit EReadError(const char* msg) : Error(msg) {}
   explicit EReadError(const std::string& msg) : Error(msg) {}
};

} // namespace gm2calc

#endif
