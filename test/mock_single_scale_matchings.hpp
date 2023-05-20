#ifndef MOCK_SINGLE_SCALE_MATCHINGS_H
#define MOCK_SINGLE_SCALE_MATCHINGS_H

#include "mock_models.hpp"
#include "single_scale_matching.hpp"

namespace flexiblesusy {

class Trivial_matching_condition: public Single_scale_matching {
public:
   Trivial_matching_condition() = default;
   Trivial_matching_condition(double scale_)
      : scale(scale_)
      {}
   virtual ~Trivial_matching_condition() {}
   virtual void match() {
      mHigh->set_parameters(mLow->get_parameters());
   }
   virtual double get_scale() const {
      return scale;
   }
   virtual void set_models(Model* mLow_, Model* mHigh_) {
      mLow = cast_model<Static_model*>(mLow_);
      mHigh = cast_model<Static_model*>(mHigh_);
   }
private:
   Static_model* mLow{nullptr};
   Static_model* mHigh{nullptr};
   double scale{100.};
};

class Counting_matching_condition: public Single_scale_matching {
public:
   Counting_matching_condition(double scale_)
      : scale(scale_)
      , number_of_matches(0)
      , number_of_get_scale(0)
      {}
   virtual ~Counting_matching_condition() {}
   virtual void match() {
      ++number_of_matches;
   }
   virtual double get_scale() const {
      ++number_of_get_scale;
      return scale;
   }
   virtual void set_models(Model*, Model*) {
   }
   int get_number_of_matches() const {
      return number_of_matches;
   }
   int get_number_of_get_scale() const {
      return number_of_get_scale;
   }
private:
   double scale;
   int number_of_matches;
   mutable int number_of_get_scale;
};

} // namespace flexiblesusy

#endif
