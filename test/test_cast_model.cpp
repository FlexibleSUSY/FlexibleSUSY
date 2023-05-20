
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_cast_model

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "model.hpp"

using namespace flexiblesusy;

class MyModel : public Model {
public:
   virtual ~MyModel() {}
   virtual void calculate_spectrum() {}
   virtual void clear_problems() {}
   virtual std::string name() const { return "MyModel"; }
   virtual void print(std::ostream& out = std::cout) const {
      out << "Model: " << name() << '\n';
   }
   virtual void run_to(double, double eps = -1.0) {}
   virtual void set_precision(double) {}
};

BOOST_AUTO_TEST_CASE( test_cast_model )
{
   MyModel model;

   Model* model_ptr = &model;

   MyModel* mymodel_ptr = cast_model<MyModel*>(model_ptr);

   BOOST_CHECK(mymodel_ptr != NULL);
}
