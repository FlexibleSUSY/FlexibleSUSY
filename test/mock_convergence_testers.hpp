#ifndef MOCK_CONVERGENCE_TESTERS_H
#define MOCK_CONVERGENCE_TESTERS_H

#include "convergence_tester.hpp"

namespace flexiblesusy {

class Counting_convergence_tester : public Convergence_tester {
public:
   Counting_convergence_tester(int max_iterations_)
      : iteration(0), maximum_iterations(max_iterations_) {}
   virtual ~Counting_convergence_tester() {}
   virtual bool accuracy_goal_reached() {
      return false;
   }
   virtual int max_iterations() const {
      return maximum_iterations;
   }
   virtual void restart() {
      iteration = 0;
   }
private:
   int iteration;
   int maximum_iterations;
};

class Static_convergence_tester : public Convergence_tester {
public:
   Static_convergence_tester(int max_iterations_)
      : iteration(0), maximum_iterations(max_iterations_) {}
   virtual ~Static_convergence_tester() {}
   virtual bool accuracy_goal_reached() {
      return false;
   }
   virtual int max_iterations() const {
      return maximum_iterations;
   }
   virtual void restart() {
      iteration = 0;
   }
private:
   int iteration;
   int maximum_iterations;
};

} // namespace flexiblesusy

#endif
