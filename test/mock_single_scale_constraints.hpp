#ifndef MOCK_SINGLE_SCALE_CONSTRAINTS_H
#define MOCK_SINGLE_SCALE_CONSTRAINTS_H

#include "model.hpp"
#include "single_scale_constraint.hpp"

namespace flexiblesusy {

class Counting_constraint : public Single_scale_constraint {
public:
   Counting_constraint(double scale_)
      : scale(scale_)
      , number_of_apply_calls(0) {}
   virtual ~Counting_constraint() {}
   virtual void apply() { ++number_of_apply_calls; }
   virtual double get_scale() const { return scale; }
   virtual void set_model(Model*) {}

   int get_number_of_apply_calls() const {
      return number_of_apply_calls;
   }

private:
   double scale;
   int number_of_apply_calls;
};

class Static_constraint : public Single_scale_constraint {
public:
   Static_constraint(double scale_)
      : scale(scale_) {}
   virtual ~Static_constraint() {}
   virtual void apply() {}
   virtual double get_scale() const { return scale; }
   virtual void set_model(Model*) {}
private:
   double scale;
};

} // namespace flexiblesusy

#endif
