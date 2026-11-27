#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lilith_interface

// Python.h must be included before the standard headers
#include "decays/lilith_interface.hpp"

#include <boost/test/unit_test.hpp>

#include "decays/decay.hpp"
#include "decays/experimental_constraints.hpp"
#include "error.hpp"
#include "physical_input.hpp"

#include <cmath>
#include <complex>
#include <optional>
#include <string>
#include <utility>

using namespace flexiblesusy;

namespace
{

/// Python stays initialized, as in call_lilith
struct Python_fixture {
   Python_fixture() { Py_Initialize(); }
};

BOOST_GLOBAL_FIXTURE(Python_fixture);

/// SM-like Higgs boson, with the bottom coupling scaled by kb
EffectiveCoupling_list make_higgs(double kb)
{
   NeutralHiggsEffectiveCouplings h;
   h.particle = "hh(1)";
   h.mass = 125.09;
   h.width = 4.1e-3;
   h.CP = CP_property::even;
   h.gamgam.second = 1.;
   h.Zgam.second = 1.;
   h.gg.second = 1.;
   h.ZZ.second = 1.;
   h.WW.second = 1.;
   h.tt.second = 1.;
   h.cc.second = 1.;
   h.bb.second = kb;
   h.tautau.second = 1.;
   h.mumu.second = 1.;

   EffectiveCoupling_list list;
   list.push_back(std::move(h));
   return list;
}

SignalResult lilith(double kb)
{
   const auto res = call_lilith(make_higgs(kb), Physical_input(), "");
   BOOST_REQUIRE(res.has_value());
   return res.value_or(SignalResult{});
}

std::string user_input(double kb)
{
   return R"(<?xml version="1.0"?>
<lilithinput>
<reducedcouplings>
<mass>125.09</mass>
<C to="tt">1.0</C>
<C to="cc">1.0</C>
<C to="bb">)" +
          std::to_string(kb) + R"(</C>
<C to="tautau">1.0</C>
<C to="mumu">1.0</C>
<C to="ZZ">1.0</C>
<C to="WW">1.0</C>
<precision>BEST-QCD</precision>
</reducedcouplings>
</lilithinput>)";
}

/// -2 log(L) from a new lilith.Lilith object in Python, without the interface
double likelihood_from_new_object(const std::string& input)
{
   PyObject* globals = PyDict_New();
   PyObject* builtins = PyImport_ImportModule("builtins");
   PyDict_SetItemString(globals, "__builtins__", builtins);
   PyObject* py_input = PyUnicode_FromString(input.c_str());
   PyDict_SetItemString(globals, "user_input", py_input);
   PyObject* result = PyRun_String("import lilith\n"
                                   "calc = lilith.Lilith(False)\n"
                                   "calc.readexpinput()\n"
                                   "calc.readuserinput(user_input)\n"
                                   "calc.computelikelihood()\n"
                                   "l = calc.l\n",
                                   Py_file_input, globals, globals);
   BOOST_REQUIRE(result != nullptr);
   const double l = PyFloat_AsDouble(PyDict_GetItemString(globals, "l"));
   Py_XDECREF(result);
   Py_XDECREF(py_input);
   Py_XDECREF(builtins);
   Py_XDECREF(globals);
   return l;
}

double likelihood_from_interface(const std::string& input)
{
   PyObject* calc = initialize_lilith("");
   lilith_readuserinput(calc, input.c_str());
   return lilith_computelikelihood(calc);
}

} // anonymous namespace

BOOST_AUTO_TEST_CASE(test_sm_like_higgs)
{
   const auto res = lilith(1.);

   BOOST_CHECK_GT(res.ndof, 0U);
   BOOST_CHECK_GT(res.chi2BSM, 0.);
   // the SM reference has ZZ = 1.00001 and no tt, cc, ... couplings
   BOOST_CHECK_CLOSE_FRACTION(res.chi2BSM, res.chi2SM, 1e-3);
}

BOOST_AUTO_TEST_CASE(test_repeated_calls_give_the_same_result)
{
   const auto first = lilith(1.);
   const auto other = lilith(0.8);
   const auto again = lilith(1.);

   // a different point gives a different result
   BOOST_CHECK_GT(std::abs(other.chi2BSM - first.chi2BSM), 1.);
   BOOST_CHECK_EQUAL(other.chi2SM, first.chi2SM);
   BOOST_CHECK_EQUAL(again.chi2BSM, first.chi2BSM);
   BOOST_CHECK_EQUAL(again.chi2SM, first.chi2SM);
   BOOST_CHECK_EQUAL(again.ndof, first.ndof);
}

BOOST_AUTO_TEST_CASE(test_reused_object_agrees_with_new_object)
{
   // the cached object has read other inputs before
   likelihood_from_interface(user_input(0.7));
   const double l_interface = likelihood_from_interface(user_input(0.9));
   const double l_new = likelihood_from_new_object(user_input(0.9));

   // exact, as the likelihood is no longer returned as a float
   BOOST_CHECK_EQUAL(l_interface, l_new);
}

BOOST_AUTO_TEST_CASE(test_invalid_database_throws)
{
   BOOST_CHECK_THROW(initialize_lilith("/nonexistent/lilith.list"), SetupError);
   BOOST_CHECK_THROW(
      call_lilith(make_higgs(1.), Physical_input(), "/nonexistent/lilith.list"),
      SetupError);

   // a valid database still works afterwards
   BOOST_CHECK_EQUAL(lilith(1.).chi2BSM, lilith(1.).chi2BSM);
   BOOST_CHECK(!PyErr_Occurred());
}

BOOST_AUTO_TEST_CASE(test_invalid_input_throws)
{
   const double before = lilith(1.).chi2BSM;

   PyObject* calc = initialize_lilith("");
   BOOST_CHECK_THROW(lilith_readuserinput(calc, "<lilithinput><broken"), Error);
   BOOST_CHECK(!PyErr_Occurred());

   // the next point isn't affected
   BOOST_CHECK_EQUAL(lilith(1.).chi2BSM, before);
}
