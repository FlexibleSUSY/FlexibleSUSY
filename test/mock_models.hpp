#ifndef MOCK_MODELS_H
#define MOCK_MODELS_H

#include "model.hpp"

#include <Eigen/Core>

namespace flexiblesusy {

class Static_model: public Model {
public:
   Static_model() = default;
   Static_model(const Eigen::VectorXd& pars) : parameters(pars) {}
   virtual ~Static_model() {}
   virtual void calculate_spectrum() {}
   virtual void clear_problems() {}
   virtual std::string name() const { return "Static_model"; }
   virtual void print(std::ostream& out = std::cout) const {
      out << "Model: " << name() << '\n';
   }
   virtual void run_to(double, double) {}
   virtual void set_parameters(const Eigen::VectorXd& v) { parameters = v; }
   virtual Eigen::VectorXd get_parameters() const { return parameters; }
   virtual void set_precision(double) {}
private:
   Eigen::VectorXd parameters{1};
};

class Counting_model: public Model {
public:
   Counting_model() = default;
   virtual ~Counting_model() {}
   virtual void calculate_spectrum() {}
   virtual void clear_problems() {}
   virtual std::string name() const { return "Counting_model"; }
   virtual void print(std::ostream& out = std::cout) const {
      out << "Model: " << name() << '\n';
   }
   virtual void run_to(double, double) { ++number_of_runs; }
   virtual void set_precision(double) {}
   int get_number_of_runs() const {
      return number_of_runs;
   }
private:
   int number_of_runs{0};
};

} // namespace flexiblesusy

#endif
