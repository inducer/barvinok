#include <barvinok/evalue.h>

#if defined(__cplusplus)
extern "C" {
#endif

struct barvinok_options;

evalue *euler_summate(Param_Polyhedron *PP, Polyhedron *TC,
    evalue *e, unsigned nvar, struct barvinok_options *options);

#if defined(__cplusplus)
}
#endif
