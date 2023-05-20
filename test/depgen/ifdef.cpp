#define SYM1
#ifdef SYM1
#include "base.hpp"
#endif

#undef SYM1
#ifdef SYM1
#include "comment.hpp"
#endif

#ifdef SYM2
#include "does_not_exist.hpp"
#endif
