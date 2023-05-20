#include "lattice_model.hpp"
#include "lattice_solver.hpp"


namespace flexiblesusy {

using namespace std;

Real  Lattice_translator::u(size_t i)
{
    return f->efts[T].units[i];
}

Real& Lattice_translator::y(size_t i)
{
    return f->y(T, m, i);
}

} // namespace flexiblesusy
