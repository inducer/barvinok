#include <assert.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include "bernoulli.h"

#define ALLOC(type) (type*)malloc(sizeof(type))
#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

static struct bernoulli_coef bernoulli_coef;
static struct poly_list faulhaber;

struct bernoulli_coef *bernoulli_coef_compute(int n)
{
    int i, j;
    Value factor, tmp;

    if (n < bernoulli_coef.n)
	return &bernoulli_coef;

    if (n >= bernoulli_coef.size) {
	int size = 3*(n + 5)/2;
	Vector *b;

	b = Vector_Alloc(size);
	if (bernoulli_coef.n) {
	    Vector_Copy(bernoulli_coef.num->p, b->p, bernoulli_coef.n);
	    Vector_Free(bernoulli_coef.num);
	}
	bernoulli_coef.num = b;
	b = Vector_Alloc(size);
	if (bernoulli_coef.n) {
	    Vector_Copy(bernoulli_coef.den->p, b->p, bernoulli_coef.n);
	    Vector_Free(bernoulli_coef.den);
	}
	bernoulli_coef.den = b;
	b = Vector_Alloc(size);
	if (bernoulli_coef.n) {
	    Vector_Copy(bernoulli_coef.lcm->p, b->p, bernoulli_coef.n);
	    Vector_Free(bernoulli_coef.lcm);
	}
	bernoulli_coef.lcm = b;

	bernoulli_coef.size = size;
    }
    value_init(factor);
    value_init(tmp);
    for (i = bernoulli_coef.n; i <= n; ++i) {
	if (i == 0) {
	    value_set_si(bernoulli_coef.num->p[0], 1);
	    value_set_si(bernoulli_coef.den->p[0], 1);
	    value_set_si(bernoulli_coef.lcm->p[0], 1);
	    continue;
	}
	value_set_si(bernoulli_coef.num->p[i], 0);
	value_set_si(factor, -(i+1));
	for (j = i-1; j >= 0; --j) {
	    mpz_mul_ui(factor, factor, j+1);
	    mpz_divexact_ui(factor, factor, i+1-j);
	    value_division(tmp, bernoulli_coef.lcm->p[i-1],
			   bernoulli_coef.den->p[j]);
	    value_multiply(tmp, tmp, bernoulli_coef.num->p[j]);
	    value_multiply(tmp, tmp, factor);
	    value_addto(bernoulli_coef.num->p[i], bernoulli_coef.num->p[i], tmp);
	}
	mpz_mul_ui(bernoulli_coef.den->p[i], bernoulli_coef.lcm->p[i-1], i+1);
	Gcd(bernoulli_coef.num->p[i], bernoulli_coef.den->p[i], &tmp);
	if (value_notone_p(tmp)) {
	    value_division(bernoulli_coef.num->p[i],
			    bernoulli_coef.num->p[i], tmp);
	    value_division(bernoulli_coef.den->p[i],
			    bernoulli_coef.den->p[i], tmp);
	}
	value_lcm(bernoulli_coef.lcm->p[i-1], bernoulli_coef.den->p[i],
		  &bernoulli_coef.lcm->p[i]);
    }
    bernoulli_coef.n = n+1;
    value_clear(factor);
    value_clear(tmp);

    return &bernoulli_coef;
}

struct poly_list *faulhaber_compute(int n)
{
    int i, j;
    Value factor;
    struct bernoulli_coef *bc;

    if (n < faulhaber.n)
	return &faulhaber;

    if (n >= faulhaber.size) {
	int size = 3*(n + 5)/2;
	Vector **poly;

	poly = ALLOCN(Vector *, size);
	for (i = 0; i < faulhaber.n; ++i)
	    poly[i] = faulhaber.poly[i];
	free(faulhaber.poly);
	faulhaber.poly = poly;

	faulhaber.size = size;
    }

    bc = bernoulli_coef_compute(n);

    value_init(factor);
    for (i = faulhaber.n; i <= n; ++i) {
	faulhaber.poly[i] = Vector_Alloc(i+3);
	value_assign(faulhaber.poly[i]->p[i+1], bc->lcm->p[i]);
	mpz_mul_ui(faulhaber.poly[i]->p[i+2], bc->lcm->p[i], i+1);
	value_set_si(factor, 1);
	for (j = 1; j <= i; ++j) {
	    mpz_mul_ui(factor, factor, i+2-j);
	    mpz_divexact_ui(factor, factor, j);
	    value_division(faulhaber.poly[i]->p[i+1-j],
			   bc->lcm->p[i], bc->den->p[j]);
	    value_multiply(faulhaber.poly[i]->p[i+1-j],
			   faulhaber.poly[i]->p[i+1-j], bc->num->p[j]);
	    value_multiply(faulhaber.poly[i]->p[i+1-j],
			   faulhaber.poly[i]->p[i+1-j], factor);
	}
	Vector_Normalize(faulhaber.poly[i]->p, i+3);
    }
    value_clear(factor);
    pl->n = n+1;

    return &faulhaber;
}

/* shift variables in polynomial one down */
static void shift(evalue *e)
{
    int i;
    if (value_notzero_p(e->d))
	return;
    assert(e->x.p->type == polynomial);
    assert(e->x.p->pos > 1);
    e->x.p->pos--;
    for (i = 0; i < e->x.p->size; ++i)
	shift(&e->x.p->arr[i]);
}

static evalue *shifted_copy(evalue *src)
{
    evalue *e = ALLOC(evalue);
    value_init(e->d);
    evalue_copy(e, src);
    shift(e);
    return e;
}

static evalue *power_sums(struct poly_list *faulhaber, evalue *poly,
			  Vector *arg, Value denom, int neg, int alt_neg)
{
    int i;
    evalue *base = affine2evalue(arg->p, denom, arg->Size-1);
    evalue *sum = evalue_zero();

    for (i = 1; i < poly->x.p->size; ++i) {
	evalue *term = evalue_polynomial(faulhaber->poly[i], base);
	evalue *factor = shifted_copy(&poly->x.p->arr[i]);
	emul(factor, term);
	if (alt_neg && (i % 2))
	    evalue_negate(term);
	eadd(term, sum);
	free_evalue_refs(factor);
	free_evalue_refs(term);
	free(factor);
	free(term);
    }
    if (neg)
	evalue_negate(sum);
    free_evalue_refs(base);
    free(base);

    return sum;
}

struct section { Polyhedron *D; evalue *E; };

struct Bernoulli_data {
    unsigned MaxRays;
    struct section *s;
    int size;
    int ns;
    evalue *e;
};

static void Bernoulli_cb(Matrix *M, Value *lower, Value *upper, void *cb_data)
{
    struct Bernoulli_data *data = (struct Bernoulli_data *)cb_data;
    Matrix *M2;
    Polyhedron *T;
    evalue *factor = NULL;
    evalue *linear = NULL;
    int constant = 0;
    Value tmp;
    unsigned dim = M->NbColumns-2;
    Vector *row;

    assert(data->ns < data->size);

    M2 = Matrix_Copy(M);
    T = Constraints2Polyhedron(M2, data->MaxRays);
    Matrix_Free(M2);

    POL_ENSURE_VERTICES(T);
    if (emptyQ(T)) {
	Polyhedron_Free(T);
	return;
    }

    assert(lower != upper);

    row = Vector_Alloc(dim+1);
    value_init(tmp);
    if (value_notzero_p(data->e->d)) {
	factor = data->e;
	constant = 1;
    } else {
	assert(data->e->x.p->type == polynomial);
	if (data->e->x.p->pos > 1) {
	    factor = shifted_copy(data->e);
	    constant = 1;
	} else
	    factor = shifted_copy(&data->e->x.p->arr[0]);
    }
    if (!EVALUE_IS_ZERO(*factor)) {
	value_absolute(tmp, upper[1]);
	/* upper - lower */
	Vector_Combine(lower+2, upper+2, row->p, tmp, lower[1], dim+1);
	value_multiply(tmp, tmp, lower[1]);
	/* upper - lower + 1 */
	value_addto(row->p[dim], row->p[dim], tmp);

	linear = affine2evalue(row->p, tmp, dim);
	emul(factor, linear);
    } else
	linear = evalue_zero();

    if (constant) {
	data->s[data->ns].E = linear;
	data->s[data->ns].D = T;
	++data->ns;
    } else {
	evalue *poly_u = NULL, *poly_l = NULL;
	Polyhedron *D;
	struct poly_list *faulhaber;
	assert(data->e->x.p->type == polynomial);
	assert(data->e->x.p->pos == 1);
	faulhaber = faulhaber_compute(data->e->x.p->size-1);
	/* lower is the constraint
	 *			 b i - l' >= 0		i >= l'/b = l
	 * upper is the constraint
	 *			-a i + u' >= 0		i <= u'/a = u
	 */
	M2 = Matrix_Alloc(3, 2+T->Dimension);
	value_set_si(M2->p[0][0], 1);
	value_set_si(M2->p[1][0], 1);
	value_set_si(M2->p[2][0], 1);
	/* Case 1:
	 *		1 <= l		l' - b >= 0
	 */
	Vector_Oppose(lower+2, M2->p[0]+1, T->Dimension+1);
	value_subtract(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension],
		      lower[1]);
	D = AddConstraints(M2->p_Init, 1, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    evalue *extra;
	    if (!poly_u) {
		Vector_Copy(upper+2, row->p, dim+1);
		value_oppose(tmp, upper[1]);
		value_addto(row->p[dim], row->p[dim], tmp);
		poly_u = power_sums(faulhaber, data->e, row, tmp, 0, 0);
	    }
	    Vector_Oppose(lower+2, row->p, dim+1);
	    extra = power_sums(faulhaber, data->e, row, lower[1], 1, 0);
	    eadd(poly_u, extra);
	    eadd(linear, extra);

	    data->s[data->ns].E = extra;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 2:
	 *		1 <= -u		-u' - a >= 0
	 */
	Vector_Oppose(upper+2, M2->p[0]+1, T->Dimension+1);
	value_addto(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension],
		      upper[1]);
	D = AddConstraints(M2->p_Init, 1, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    evalue *extra;
	    if (!poly_l) {
		Vector_Copy(lower+2, row->p, dim+1);
		value_addto(row->p[dim], row->p[dim], lower[1]);
		poly_l = power_sums(faulhaber, data->e, row, lower[1], 0, 1);
	    }
	    Vector_Oppose(upper+2, row->p, dim+1);
	    value_oppose(tmp, upper[1]);
	    extra = power_sums(faulhaber, data->e, row, tmp, 1, 1);
	    eadd(poly_l, extra);
	    eadd(linear, extra);

	    data->s[data->ns].E = extra;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 3:
	 *		u >= 0		u' >= 0
	 *		-l >= 0		-l' >= 0
	 */
	Vector_Copy(upper+2, M2->p[0]+1, T->Dimension+1);
	Vector_Copy(lower+2, M2->p[1]+1, T->Dimension+1);
	D = AddConstraints(M2->p_Init, 2, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    if (!poly_l) {
		Vector_Copy(lower+2, row->p, dim+1);
		value_addto(row->p[dim], row->p[dim], lower[1]);
		poly_l = power_sums(faulhaber, data->e, row, lower[1], 0, 1);
	    }
	    if (!poly_u) {
		Vector_Copy(upper+2, row->p, dim+1);
		value_oppose(tmp, upper[1]);
		value_addto(row->p[dim], row->p[dim], tmp);
		poly_u = power_sums(faulhaber, data->e, row, tmp, 0, 0);
	    }
	
	    data->s[data->ns].E = ALLOC(evalue);
	    value_init(data->s[data->ns].E->d);
	    evalue_copy(data->s[data->ns].E, poly_u);
	    eadd(poly_l, data->s[data->ns].E);
	    eadd(linear, data->s[data->ns].E);
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 4:
	 *		l < 1		-l' + b - 1 >= 0
	 *		0 < l		l' - 1 >= 0
	 */
	Vector_Copy(lower+2, M2->p[0]+1, T->Dimension+1);
	value_addto(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension], lower[1]);
	value_decrement(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension]);
	Vector_Oppose(lower+2, M2->p[1]+1, T->Dimension+1);
	value_decrement(M2->p[1][1+T->Dimension], M2->p[1][1+T->Dimension]);
	D = AddConstraints(M2->p_Init, 2, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    if (!poly_u) {
		Vector_Copy(upper+2, row->p, dim+1);
		value_oppose(tmp, upper[1]);
		value_addto(row->p[dim], row->p[dim], tmp);
		poly_u = power_sums(faulhaber, data->e, row, tmp, 0, 0);
	    }

	    eadd(linear, poly_u);
	    data->s[data->ns].E = poly_u;
	    poly_u = NULL;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	/* Case 5:
	 * 		-u < 1		u' + a - 1 >= 0
	 *		0 < -u		-u' - 1 >= 0
	 *		l <= 0		-l' >= 0
	 */
	Vector_Copy(upper+2, M2->p[0]+1, T->Dimension+1);
	value_subtract(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension],
			upper[1]);
	value_decrement(M2->p[0][1+T->Dimension], M2->p[0][1+T->Dimension]);
	Vector_Oppose(upper+2, M2->p[1]+1, T->Dimension+1);
	value_decrement(M2->p[1][1+T->Dimension], M2->p[1][1+T->Dimension]);
	Vector_Copy(lower+2, M2->p[2]+1, T->Dimension+1);
	D = AddConstraints(M2->p_Init, 3, T, data->MaxRays);
	if (emptyQ2(D))
	    Polyhedron_Free(D);
	else {
	    if (!poly_l) {
		Vector_Copy(lower+2, row->p, dim+1);
		value_addto(row->p[dim], row->p[dim], lower[1]);
		poly_l = power_sums(faulhaber, data->e, row, lower[1], 0, 1);
	    }

	    eadd(linear, poly_l);
	    data->s[data->ns].E = poly_l;
	    poly_l = NULL;
	    data->s[data->ns].D = D;
	    ++data->ns;
	}

	Matrix_Free(M2);
	Polyhedron_Free(T);
	if (poly_l) {
	    free_evalue_refs(poly_l);
	    free(poly_l);
	}
	if (poly_u) {
	    free_evalue_refs(poly_u);
	    free(poly_u);
	}
	free_evalue_refs(linear);
	free(linear);
    }
    if (factor != data->e) {
	free_evalue_refs(factor);
	free(factor);
    }
    value_clear(tmp);
    Vector_Free(row);
}

evalue *Bernoulli_sum_evalue(evalue *e, unsigned nvar,
			     struct barvinok_options *options)
{
    struct Bernoulli_data data;
    int i, j;
    evalue *sum = evalue_zero();

    if (EVALUE_IS_ZERO(*e))
	return sum;

    if (nvar == 0) {
	eadd(e, sum);
	return sum;
    }

    assert(value_zero_p(e->d));
    assert(e->x.p->type == partition);

    data.size = e->x.p->size * 2 + 16;
    data.s = ALLOCN(struct section, data.size);
    data.MaxRays = options->MaxRays;

    for (i = 0; i < e->x.p->size/2; ++i) {
	Polyhedron *D;
	for (D = EVALUE_DOMAIN(e->x.p->arr[2*i]); D; D = D->next) {
	    unsigned dim = D->Dimension - 1;
	    Polyhedron *next = D->next;
	    evalue tmp;
	    D->next = NULL;

	    value_init(tmp.d);
	    value_set_si(tmp.d, 0);

	    if (value_zero_p(D->Constraint[0][0]) &&
		    value_notzero_p(D->Constraint[0][1])) {
		tmp.x.p = new_enode(partition, 2, dim);
		EVALUE_SET_DOMAIN(tmp.x.p->arr[0], Polyhedron_Project(D, dim));
		evalue_copy(&tmp.x.p->arr[1], &e->x.p->arr[2*i+1]);
		reduce_evalue_in_domain(&tmp.x.p->arr[1], D);
		shift(&tmp.x.p->arr[1]);
	    } else {
		data.ns = 0;
		data.e = &e->x.p->arr[2*i+1];

		for_each_lower_upper_bound(D, Bernoulli_cb, &data);

		if (data.ns == 0)
		    evalue_set_si(&tmp, 0, 1);
		else {
		    tmp.x.p = new_enode(partition, 2*data.ns, dim);
		    for (j = 0; j < data.ns; ++j) {
			EVALUE_SET_DOMAIN(tmp.x.p->arr[2*j], data.s[j].D);
			value_clear(tmp.x.p->arr[2*j+1].d);
			tmp.x.p->arr[2*j+1] = *data.s[j].E;
			free(data.s[j].E);
		    }
		}
	    }

	    if (nvar > 1) {
		evalue *res = Bernoulli_sum_evalue(&tmp, nvar-1, options);
		eadd(res, sum);
		free_evalue_refs(res);
		free(res);
	    } else
		eadd(&tmp, sum);

	    free_evalue_refs(&tmp);
	    D->next = next;;
	}
    }

    free(data.s);

    reduce_evalue(sum);
    return sum;
}

evalue *Bernoulli_sum(Polyhedron *P, Polyhedron *C,
			struct barvinok_options *options)
{
    evalue e;
    evalue *sum;

    value_init(e.d);
    e.x.p = new_enode(partition, 2, P->Dimension);
    EVALUE_SET_DOMAIN(e.x.p->arr[0], Polyhedron_Copy(P));
    evalue_set_si(&e.x.p->arr[1], 1, 1);
    sum = Bernoulli_sum_evalue(&e, P->Dimension - C->Dimension, options);
    free_evalue_refs(&e);
    return sum;
}
