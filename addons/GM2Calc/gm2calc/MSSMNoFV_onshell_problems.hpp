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

#ifndef GM2_MSSMNoFV_ONSHELL_PROBLEMS_HPP
#define GM2_MSSMNoFV_ONSHELL_PROBLEMS_HPP

#include <iosfwd>
#include <string>
#include <vector>

namespace gm2calc {

/**
 * @class MSSMNoFV_onshell_problems
 * @brief contains problem and warning flags
 *
 * If a problem has occurred, the physical particle spectrum cannot be
 * trusted (for example tachyons are present).  A warning means an
 * imprecision has occurred and care must be taken when interpreting
 * the particle spectrum.
 */
class MSSMNoFV_onshell_problems {
public:
   struct Convergence_problem {
      void clear() { *this = Convergence_problem(); }
      double precision{0.0};  ///< achieved accuracy (in GeV)
      unsigned iterations{0}; ///< used number of iterations
   };

   void clear();              ///< delete all problems and warnings
   void clear_problems();     ///< delete all problems
   void clear_warnings();     ///< delete all warnings
   void flag_no_convergence_Mu_MassB_MassWB(double, unsigned);
   void flag_no_convergence_me2(double, unsigned);
   void flag_tachyon(const std::string&);
   void unflag_no_convergence_Mu_MassB_MassWB();
   void unflag_no_convergence_me2();
   bool no_Mu_MassB_MassWB_convergence() const;
   bool no_me2_convergence() const;
   Convergence_problem get_Mu_MassB_MassWB_convergence_problem() const;
   Convergence_problem get_me2_convergence_problem() const;
   bool have_tachyon() const; ///< returns true if tachyon exists
   bool have_problem() const; ///< returns true if problem has occurred
   bool have_warning() const; ///< returns true if there is a warning
   std::string get_warnings() const; ///< get warnings as string
   std::string get_problems() const; ///< get problems as string
   void print(std::ostream&) const; ///< print problems and warnings to stream
   void print_problems(std::ostream&) const; ///< print problems to stream
   void print_warnings(std::ostream&) const; ///< print warnings to stream

private:
   bool have_no_convergence_Mu_MassB_MassWB{false};
   bool have_no_convergence_me2{false};
   std::vector<std::string> tachyons;
   Convergence_problem convergence_problem_Mu_MassB_MassWB;
   Convergence_problem convergence_problem_me2;
};

std::ostream& operator<<(std::ostream&, const MSSMNoFV_onshell_problems&);

} // namespace gm2calc

#endif
