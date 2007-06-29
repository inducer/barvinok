#include <barvinok/util.h>
#include <barvinok/options.h>
#include "polysign.h"

static enum order_sign interval_minmax(Polyhedron *I)
{
    int i;
    int min = 1;
    int max = -1;
    assert(I->Dimension == 1);
    POL_ENSURE_VERTICES(I);
    for (i = 0; i < I->NbRays; ++i) {
	if (value_pos_p(I->Ray[i][1]))
	    max = 1;
	else if (value_neg_p(I->Ray[i][1]))
	    min = -1;
	else {
	    if (max < 0)
		max = 0;
	    if (min > 0)
		min = 0;
	}
    }
    if (min > 0)
	return order_gt;
    if (max < 0)
	return order_lt;
    if (min == max)
	return order_eq;
    if (max == 0)
	return order_le;
    if (min == 0)
	return order_ge;
    return order_unknown;
}

/* Returns the sign of the affine function specified by T on the polyhedron D */
enum order_sign PL_polyhedron_affine_sign(Polyhedron *D, Matrix *T,
					    struct barvinok_options *options)
{
    enum order_sign sign;
    Polyhedron *I = Polyhedron_Image(D, T, options->MaxRays);
    if (POL_ISSET(options->MaxRays, POL_INTEGER))
	I = DomainConstraintSimplify(I, options->MaxRays);
    if (emptyQ2(I)) {
	Polyhedron_Free(I);
	I = Polyhedron_Image(D, T, options->MaxRays);
    }
    sign = interval_minmax(I);
    Polyhedron_Free(I);
    return sign;
}

enum lp_result PL_polyhedron_range(Polyhedron *D, Value *obj, Value denom,
				Value *min, Value *max,
				struct barvinok_options *options)
{
    Polyhedron *I;
    int bounded;

    Matrix *T = Matrix_Alloc(2, D->Dimension+1);
    Vector_Copy(obj, T->p[0], D->Dimension+1);
    value_assign(T->p[1][D->Dimension], denom);

    I = Polyhedron_Image(D, T, options->MaxRays);
    Matrix_Free(T);
    POL_ENSURE_VERTICES(I);

    if (emptyQ(I)) {
	Polyhedron_Free(I);
	return lp_empty;
    }
    bounded = line_minmax(I, min, max);

    if (!bounded)
	return lp_unbounded;
    if (value_gt(*min, *max))
	return lp_empty;
    return lp_ok;
}
