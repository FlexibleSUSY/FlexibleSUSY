
/** \file rge.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief RGE objects consisting of energy scale and parameters 
   and loops (order in perturbation theory) and thresholds
   (accuracy of calculation, handles decoupling ususally).
*/

#ifndef RGE_H
#define RGE_H

#include "def.h"
#include "utils.h"
#include "linalg.h"
#include "numerics.h"
#include "rk_legacy.hpp"

namespace softsusy {

/// Describes a set of parameters and RGEs in general. 
class RGE
{
private:
  double mu; ///< Renormalisation scale
  int numPars; ///< Number of parameters
  int loops; ///< To what order does the RG evolution run
  int thresholds; ///< Threshold flag: not always used

public:
  /// Default constructor fills data with zeroes
  RGE(): mu(0.0), numPars(0), loops(0), thresholds(0) {};
  virtual ~RGE() {};
  
  /// Sets renormalisation scale to e
  void setMu(double e) { mu = e; };
  /// Default is to just set mu, but if you've got a derived class, 
  /// you may like to call several setMu's.
  virtual void setmu(double f) { setMu(f); };
  /// Set number of parameters in RGE object
  void setPars(int i) { numPars = i; };
  /// Set number of loops used
  void setLoops(int l) { loops = l; };
  /// Set level of threshold approximation used. It's typically used in the
  /// beta functions or boundary conditions.
  void setThresholds(int t) { thresholds = t; };

  /// Return renomalisation scale
  double displayMu() const { return mu; };
  /// Return number of loops
  int displayLoops() const { return loops; };
  /// Return number of parameters
  int howMany() const { return numPars; };
  /// Return level of threshold approximation
  int displayThresholds() const { return thresholds; };  
  /// Displays all RGE parameters in a double vector. Obligatory function for
  /// derived classes.
  virtual const DoubleVector display() const = 0;
  /// Sets all parameters in an RGE object to the contents of s. It's an
  /// obligatory function for a derived class.
  virtual void set(const DoubleVector & s) = 0;
  /// Beta functions - all RGE objects must have them.
  virtual DoubleVector beta() const = 0;

  /// Runs parameters from scale "from" to "to, eps is fractional accuracy
  virtual int run(double from, double to, double eps= -1.0);
  /// Runs parameters to scale "to", eps is accuracy
  virtual int runto(double to, double eps = -1.0);
  
  /// Does the actual calling of Runge Kutta: default precision is TOLERANCE.
  /// Returns >0 if there's a problem with the running.
  /// IO parameters: x1=start renormalisation scale, x2=end renormalisation
  /// scale, v=boundary condition on parameters to be run on input, but gets
  /// changed to calculated final values on output, function to
  /// calculate derivatives are encoded in *derivs, eps is the numerical
  /// accuracy 
  int callRK(double x1, double x2, DoubleVector & v, 
	     DoubleVector (*derivs)(double, const DoubleVector &), 
	     double eps = -1.0);

};  

/// Runge-Kutte user defined routine: given log renorm scale x and the
/// dependent variables y of an RGE, will calculate the derivitives dydx.
DoubleVector allDerivs(double, const DoubleVector &);

} // namespace softsusy

#endif
