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

#include "config.h"

#ifdef ENABLE_LILITH

#include "decays/lilith_interface.hpp"

#include "error.hpp"

#include <map>
#include <string>

#ifndef FLEXIBLESUSY_LILITH_DIR
#error "FLEXIBLESUSY_LILITH_DIR (Lilith's root directory) is not defined"
#endif

namespace flexiblesusy
{

namespace
{

/// owns a new reference to a Python object
class Py_ref
{
public:
   explicit Py_ref(PyObject* obj = nullptr) noexcept : obj(obj) {}
   ~Py_ref() { Py_XDECREF(obj); }
   Py_ref(const Py_ref&) = delete;
   Py_ref(Py_ref&&) = delete;
   Py_ref& operator=(const Py_ref&) = delete;
   Py_ref& operator=(Py_ref&&) = delete;

   PyObject* get() const noexcept { return obj; }
   /// gives up the ownership
   PyObject* release() noexcept
   {
      PyObject* tmp = obj;
      obj = nullptr;
      return tmp;
   }
   explicit operator bool() const noexcept { return obj != nullptr; }

private:
   PyObject* obj{nullptr};
};

/// returns str(obj), or an empty string if that fails
std::string to_string(PyObject* obj)
{
   if (obj == nullptr) {
      return {};
   }
   const Py_ref str(PyObject_Str(obj));
   const char* chars = str ? PyUnicode_AsUTF8(str.get()) : nullptr;
   if (chars == nullptr) {
      PyErr_Clear();
      return {};
   }
   return chars;
}

/// returns "type: message" of the current Python exception and clears it
std::string fetch_python_error()
{
   PyObject* type = nullptr;
   PyObject* value = nullptr;
   PyObject* traceback = nullptr;
   PyErr_Fetch(&type, &value, &traceback);
   PyErr_NormalizeException(&type, &value, &traceback);
   const Py_ref type_ref(type);
   const Py_ref value_ref(value);
   const Py_ref traceback_ref(traceback);

   std::string name;
   if (type != nullptr) {
      const Py_ref type_name(PyObject_GetAttrString(type, "__name__"));
      name = to_string(type_name.get());
      PyErr_Clear();
   }
   const std::string message = to_string(value);

   if (name.empty()) {
      return message.empty() ? "unknown error" : message;
   }
   return message.empty() ? name : name + ": " + message;
}

[[noreturn]] void throw_python_error(const std::string& what)
{
   throw Error("Lilith: " + what + ": " + fetch_python_error());
}

[[noreturn]] void throw_python_setup_error(const std::string& what)
{
   throw SetupError("Lilith: " + what + ": " + fetch_python_error());
}

/// returns the tuple (item,), taking over the reference to item, or
/// nullptr on error
PyObject* make_args(PyObject* item)
{
   if (item == nullptr) {
      return nullptr;
   }
   PyObject* args = PyTuple_New(1);
   if (args == nullptr) {
      Py_DECREF(item);
      return nullptr;
   }
   PyTuple_SET_ITEM(args, 0, item); // steals the reference
   return args;
}

/// calls obj.name(*args), returns a new reference or nullptr on error
PyObject* call_method(PyObject* obj, const char* name, PyObject* args)
{
   const Py_ref method(PyObject_GetAttrString(obj, name));
   if (!method) {
      return nullptr;
   }
   return PyObject_CallObject(method.get(), args);
}

/// returns obj.name as a double
double get_double_attr(PyObject* obj, const char* name)
{
   const Py_ref attr(PyObject_GetAttrString(obj, name));
   if (!attr) {
      throw_python_error(std::string("cannot get ") + name);
   }
   const double value = PyFloat_AsDouble(attr.get());
   if (value == -1. && PyErr_Occurred()) {
      throw_python_error(std::string("cannot convert ") + name +
                         " to a number");
   }
   return value;
}

/// imports the module lilith from Lilith's root directory
PyObject* import_lilith()
{
   PyObject* sys_path = PySys_GetObject("path"); // borrowed
   const Py_ref lilith_dir(PyUnicode_FromString(FLEXIBLESUSY_LILITH_DIR));
   if (sys_path == nullptr || !lilith_dir) {
      throw_python_setup_error("cannot get sys.path");
   }
   const int contained = PySequence_Contains(sys_path, lilith_dir.get());
   if (contained < 0 ||
       (contained == 0 && PyList_Append(sys_path, lilith_dir.get()) < 0)) {
      throw_python_setup_error("cannot add " FLEXIBLESUSY_LILITH_DIR
                               " to sys.path");
   }

   PyObject* module = PyImport_ImportModule("lilith");
   if (module == nullptr) {
      throw_python_setup_error(
         "cannot import lilith from " FLEXIBLESUSY_LILITH_DIR);
   }
   return module;
}

/// returns a new lilith.Lilith object that has read the database
PyObject* create_lilith(const std::string& experimental_input)
{
   const Py_ref module(import_lilith());
   const Py_ref lilith_class(PyObject_GetAttrString(module.get(), "Lilith"));
   if (!lilith_class) {
      throw_python_setup_error("cannot find the class Lilith");
   }

   // Lilith(verbose=False)
   const Py_ref args(make_args(PyBool_FromLong(0)));
   Py_ref lilithcalc(args ? PyObject_CallObject(lilith_class.get(), args.get())
                          : nullptr);
   if (!lilithcalc) {
      throw_python_setup_error("cannot create a Lilith object");
   }

   // Lilith's default list is used if no list was given
   const Py_ref read_args(
      experimental_input.empty()
         ? PyTuple_New(0)
         : make_args(PyUnicode_FromString(experimental_input.c_str())));
   const Py_ref result(
      read_args ? call_method(lilithcalc.get(), "readexpinput", read_args.get())
                : nullptr);
   if (!result) {
      throw_python_setup_error("cannot read the experimental results \"" +
                               experimental_input + "\"");
   }

   return lilithcalc.release();
}

} // anonymous namespace

PyObject* initialize_lilith(const char* experimental_input)
{
   // The objects are never released: they would have to be released
   // before Python is finalized, and Python stays initialized for the
   // lifetime of the process (see call_lilith).
   static std::map<std::string, PyObject*> lilith_objects;

   const std::string input(experimental_input != nullptr ? experimental_input
                                                         : "");

   const auto it = lilith_objects.find(input);
   if (it != lilith_objects.end()) {
      return it->second;
   }

   PyObject* lilithcalc = create_lilith(input);
   lilith_objects.emplace(input, lilithcalc);

   return lilithcalc;
}

PyObject* lilith_readuserinput(PyObject* lilithcalc, const char* XMLinputstring)
{
   const Py_ref args(make_args(PyUnicode_FromString(XMLinputstring)));
   const Py_ref result(
      args ? call_method(lilithcalc, "readuserinput", args.get()) : nullptr);
   if (!result) {
      throw_python_error("cannot read the input");
   }
   return lilithcalc;
}

double lilith_computelikelihood(PyObject* lilithcalc)
{
   const Py_ref result(call_method(lilithcalc, "computelikelihood", nullptr));
   if (!result) {
      throw_python_error("cannot compute the likelihood");
   }
   return get_double_attr(lilithcalc, "l");
}

double lilith_exp_ndf(PyObject* lilithcalc)
{
   return get_double_attr(lilithcalc, "exp_ndf");
}

} // namespace flexiblesusy

#endif // ENABLE_LILITH
