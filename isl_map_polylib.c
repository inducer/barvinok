/*
 * Copyright 2008-2009 Katholieke Universiteit Leuven
 *
 * Use of this software is governed by the GNU GPLv2 license
 *
 * Written by Sven Verdoolaege, K.U.Leuven, Departement
 * Computerwetenschappen, Celestijnenlaan 200A, B-3001 Leuven, Belgium
 */

#include <isl/val.h>
#include <isl/val_gmp.h>
#include <isl/space.h>
#include <isl/set.h>
#include <isl/map.h>
#include <isl/constraint.h>
#include "isl_set_polylib.h"
#include "isl_map_polylib.h"

struct isl_basic_set *isl_basic_set_new_from_polylib(Polyhedron *P,
			struct isl_space *dim)
{
	isl_ctx *ctx;

	if (!dim)
		return NULL;
	ctx = isl_space_get_ctx(dim);
	isl_assert(ctx, isl_space_dim(dim, isl_dim_in) == 0, return NULL);

	return (struct isl_basic_set *)
		isl_basic_map_new_from_polylib(P, dim);
}

/* Return the number of equality constraints in the polyhedron description "P".
 * The equality constraints have a zero in the first column.
 * They also appear before the inequality constraints, but this code
 * does not rely on this order.
 */
static int polyhedron_n_eq(Polyhedron *P)
{
	int i, n = 0;

	for (i = 0; i < P->NbConstraints; ++i)
		if (value_zero_p(P->Constraint[i][0]))
			++n;

	return n;
}

/* Set the row "row" of "dst" to the values in array "src".
 */
static __isl_give isl_mat *set_row(__isl_take isl_mat *dst, int row,
	Value *src)
{
	int i, n;
	isl_ctx *ctx;

	ctx = isl_mat_get_ctx(dst);
	n = isl_mat_cols(dst);
	for (i = 0; i < n; ++i) {
		isl_val *v;

		v = isl_val_int_from_gmp(ctx, src[i]);
		dst = isl_mat_set_element_val(dst, row, i, v);
	}

	return dst;
}

/* Extract the "n_eq" equality constraints from "P", dropping the column
 * that identifies equality constraints.
 */
static __isl_give isl_mat *extract_equalities(isl_ctx *ctx, Polyhedron *P,
	int n_eq)
{
	int i, j;
	isl_mat *eq;

	eq = isl_mat_alloc(ctx, n_eq, P->Dimension + 1);
	for (i = 0, j = 0; i < P->NbConstraints; ++i) {
		if (!value_zero_p(P->Constraint[i][0]))
			continue;
		eq = set_row(eq, j++, P->Constraint[i] + 1);
	}

	return eq;
}

/* Extract the "n_ineq" inequality constraints from "P", dropping the column
 * that identifies equality constraints.
 */
static __isl_give isl_mat *extract_inequalities(isl_ctx *ctx, Polyhedron *P,
	int n_ineq)
{
	int i, j;
	isl_mat *ineq;

	ineq = isl_mat_alloc(ctx, n_ineq, P->Dimension + 1);
	for (i = 0, j = 0; i < P->NbConstraints; ++i) {
		if (value_zero_p(P->Constraint[i][0]))
			continue;
		ineq = set_row(ineq, j++, P->Constraint[i] + 1);
	}

	return ineq;
}

__isl_give isl_basic_map *isl_basic_map_new_from_polylib(Polyhedron *P,
	__isl_take isl_space *space)
{
	isl_ctx *ctx;
	isl_mat *eq, *ineq;
	unsigned n_eq, n_ineq;

	if (!space)
		return NULL;

	ctx = isl_space_get_ctx(space);
	isl_assert(ctx, P, goto error);
	isl_assert(ctx, P->Dimension >= isl_space_dim(space, isl_dim_all),
		    goto error);

	n_eq = polyhedron_n_eq(P);
	n_ineq = P->NbConstraints - n_eq;
	eq = extract_equalities(ctx, P, n_eq);
	ineq = extract_inequalities(ctx, P, n_ineq);

	return isl_basic_map_from_constraint_matrices(space, eq, ineq,
	    isl_dim_in, isl_dim_out, isl_dim_div, isl_dim_param, isl_dim_cst);
error:
	isl_space_free(space);
	return NULL;
}

struct isl_set *isl_set_new_from_polylib(Polyhedron *D, struct isl_space *dim)
{
	isl_ctx *ctx;
	struct isl_set *set = NULL;
	Polyhedron *P;

	if (!dim)
		return NULL;
	ctx = isl_space_get_ctx(dim);
	isl_assert(ctx, isl_space_dim(dim, isl_dim_in) == 0, return NULL);

	set = isl_set_empty(isl_space_copy(dim));
	if (!set)
		goto error;

	for (P = D; P; P = P->next)
		set = isl_set_union_disjoint(set,
		    isl_set_from_basic_set(
		    isl_basic_set_new_from_polylib(P, isl_space_copy(dim))));
	isl_space_free(dim);
	return set;
error:
	isl_space_free(dim);
	return NULL;
}

static isl_stat count_constraints(__isl_take isl_constraint *c, void *user)
{
	int *n = (int *)user;
	(*n)++;
	isl_constraint_free(c);
	return isl_stat_ok;
}

struct isl_poly_copy {
	int n;
	Matrix *M;
};

static isl_stat copy_constraint_to(__isl_take isl_constraint *c, void *user)
{
	int i, j, k;
	enum isl_dim_type types[] = { isl_dim_in, isl_dim_out,
					isl_dim_div, isl_dim_param };
	struct isl_poly_copy *data = (struct isl_poly_copy *)user;
	isl_val *v;

	if (isl_constraint_is_equality(c))
		value_set_si(data->M->p[data->n][0], 0);
	else
		value_set_si(data->M->p[data->n][0], 1);
	k = 1;
	for (i = 0; i < 4; ++i) {
		int n = isl_constraint_dim(c, types[i]);
		for (j = 0; j < n; ++j, ++k) {
			v = isl_constraint_get_coefficient_val(c, types[i], j);
			isl_val_get_num_gmp(v, data->M->p[data->n][k]);
			isl_val_free(v);
		}
	}
	v = isl_constraint_get_constant_val(c);
	isl_val_get_num_gmp(v, data->M->p[data->n][k]);
	isl_val_free(v);
	isl_constraint_free(c);
	data->n++;
	return isl_stat_ok;
}

Polyhedron *isl_basic_map_to_polylib(__isl_keep isl_basic_map *bmap)
{
	Polyhedron *P;
	unsigned nparam;
	unsigned n_in;
	unsigned n_out;
	unsigned max_rays;
	unsigned n_div;
	int n = 0;
	struct isl_poly_copy data;

	if (!bmap)
		return NULL;

	if (isl_basic_map_is_rational(bmap))
		max_rays = POL_NO_DUAL;
	else
		max_rays = POL_NO_DUAL | POL_INTEGER;

	if (isl_basic_map_foreach_constraint(bmap, &count_constraints, &n) < 0)
		return NULL;

	nparam = isl_basic_map_n_param(bmap);
	n_in = isl_basic_map_n_in(bmap);
	n_out = isl_basic_map_n_out(bmap);
	n_div = isl_basic_map_dim(bmap, isl_dim_div);
	data.M = Matrix_Alloc(n, 1 + n_in + n_out + n_div + nparam + 1);
	data.n = 0;
	if (isl_basic_map_foreach_constraint(bmap,
					    &copy_constraint_to, &data) < 0) {
		Matrix_Free(data.M);
		return NULL;
	}
	P = Constraints2Polyhedron(data.M, max_rays);
	Matrix_Free(data.M);

	return P;
}

static isl_stat add_basic_map(__isl_take isl_basic_map *bmap, void *user)
{
	Polyhedron ***next = user;

	**next = isl_basic_map_to_polylib(bmap);
	*next = &(**next)->next;

	isl_basic_map_free(bmap);
	return isl_stat_ok;
}

Polyhedron *isl_map_to_polylib(struct isl_map *map)
{
	Polyhedron *R = NULL;
	Polyhedron **next = &R;

	if (!map)
		return NULL;

	if (isl_map_foreach_basic_map(map, &add_basic_map, &next) < 0)
		goto error;

	return R ? R : Empty_Polyhedron(isl_map_dim(map, isl_dim_all));
error:
	Domain_Free(R);
	return NULL;
}

Polyhedron *isl_basic_set_to_polylib(struct isl_basic_set *bset)
{
	return isl_basic_map_to_polylib((struct isl_basic_map *)bset);
}

Polyhedron *isl_set_to_polylib(struct isl_set *set)
{
	return isl_map_to_polylib((struct isl_map *)set);
}
