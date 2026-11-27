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

/**
 * @file lilith_interface.hpp
 *
 * @brief interface to the Python class lilith.Lilith
 *
 * Replaces the C API of Lilith (lilith/c-api/lilith.c), whose functions
 * have the same names.  Unlike that API, these functions
 *
 * - throw a flexiblesusy::Error with Python's message instead of
 *   ending the process with exit(1) or returning -1,
 * - release the Python objects they create,
 * - return -2 log(L) in double precision, not as a float, and
 * - initialize Lilith and read its experimental database only once per
 *   database (see initialize_lilith).
 *
 * The Python interpreter must be initialized before they are called.
 */

#ifndef DECAYS_LILITH_INTERFACE_H
#define DECAYS_LILITH_INTERFACE_H

#include "config.h"

#ifdef ENABLE_LILITH

#include <Python.h>

namespace flexiblesusy
{

/**
 * Returns an object of the class lilith.Lilith that has read the
 * experimental results listed in the file @a experimental_input (the
 * default list of Lilith if it is empty).
 *
 * The object is created at the first call for a database and returned
 * again at later calls, so the database is read only once.  The returned
 * reference is borrowed: it belongs to the interface and is never
 * released, so the caller must not release it either.
 *
 * @throw SetupError if Lilith can't be imported or the database can't be
 * read
 */
PyObject* initialize_lilith(const char* experimental_input);

/**
 * Passes the XML input @a XMLinputstring to @a lilithcalc.  Replaces the
 * input of the previous call.
 *
 * @return @a lilithcalc (borrowed)
 * @throw Error if Lilith can't read the input
 */
PyObject* lilith_readuserinput(PyObject* lilithcalc,
                               const char* XMLinputstring);

/**
 * Returns -2 log(L) for the input of the last lilith_readuserinput call.
 *
 * @throw Error if Lilith raises an error
 */
double lilith_computelikelihood(PyObject* lilithcalc);

/**
 * Returns the number of degrees of freedom of the experimental results.
 *
 * @throw Error if Lilith raises an error
 */
double lilith_exp_ndf(PyObject* lilithcalc);

} // namespace flexiblesusy

#endif // ENABLE_LILITH

#endif
