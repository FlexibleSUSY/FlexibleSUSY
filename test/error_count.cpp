
#include "error_count.hpp"

namespace flexiblesusy {

int gErrors = 0;

int get_errors() {
   return gErrors;
}

void clear_errors() {
   gErrors = 0;
}

}
