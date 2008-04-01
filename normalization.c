#include <assert.h>
#include "normalization.h"

static int is_unit_row(Value *row, int pos, int len)
{
    if (!value_one_p(row[pos]) && !value_mone_p(row[pos]))
	return 0;
    return First_Non_Zero(row+pos+1, len-(pos+1)) == -1;
}

/* Transform the constraints of P to "standard form".
 * In particular, if P is described by
 *		A x + b(p) >= 0
 * then this function returns a matrix H = A U, A = H Q, such
 * that D x' = D Q x >= -b(p), with D a diagonal matrix with
 * positive entries.  The calling function can then construct
 * the standard form H' x' - I s + b(p) = 0, with H' the rows of H
 * that are not positive multiples of unit vectors
 * (since those correspond to D x' >= -b(p)).
 * The number of rows in H' is returned in *rows_p.
 * Optionally returns the matrix that maps the new variables
 * back to the old variables x = U x'.
 * Note that the rows of H (and P) may be reordered by this function.
 */
Matrix *standard_constraints(Polyhedron *P, unsigned nparam, int *rows_p,
			     Matrix **T)
{
    unsigned nvar = P->Dimension - nparam;
    int i, j, d;
    int rows;
    Matrix *M;
    Matrix *H, *U, *Q;

    assert(P->NbEq == 0);

    /* move constraints only involving parameters down
     * and move unit vectors (if there are any) to the right place.
     */
    for (d = 0, j = P->NbConstraints; d < j; ++d) {
	int index;
	index = First_Non_Zero(P->Constraint[d]+1, nvar);
	if (index != -1) {
	    if (index != d &&
		is_unit_row(P->Constraint[d]+1, index, nvar)) {
		Vector_Exchange(P->Constraint[d], P->Constraint[index],
				P->Dimension+2);
		if (index > d)
		    --d;
	    }
	    continue;
	}
	while (d < --j && First_Non_Zero(P->Constraint[j]+1, nvar) == -1)
	    ;
	if (d >= j)
	    break;
	Vector_Exchange(P->Constraint[d], P->Constraint[j], P->Dimension+2);
    }
    M = Matrix_Alloc(d, nvar);
    for (j = 0; j < d; ++j)
	Vector_Copy(P->Constraint[j]+1, M->p[j], nvar);

    neg_left_hermite(M, &H, &Q, &U);
    Matrix_Free(M);
    Matrix_Free(Q);
    if (T)
	*T = U;
    else
	Matrix_Free(U);

    /* Rearrange rows such that top of H is lower diagonal and
     * compute the number of non (multiple of) unit-vector rows.
     */
    rows = H->NbRows-nvar;
    for (i = 0; i < H->NbColumns; ++i) {
	for (j = i; j < H->NbRows; ++j)
	    if (value_notzero_p(H->p[j][i]))
		break;
	if (j != i) {
	    Vector_Exchange(P->Constraint[i], P->Constraint[j], P->Dimension+2);
	    Vector_Exchange(H->p[i], H->p[j], H->NbColumns);
	}
	if (First_Non_Zero(H->p[i], i) != -1)
	    rows++;
    }
    *rows_p = rows;

    return H;
}
