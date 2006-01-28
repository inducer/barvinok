#include <unistd.h>
#include <stdlib.h>
#include <polylib/polylibgmp.h>
#include <barvinok/util.h>
#include <barvinok/barvinok.h>
#include "config.h"
#ifdef HAVE_OMEGA
#include "omega/convert.h"
#endif

/* The input of this example program is a polytope in combined
 * data and parameter space followed by two lines indicating
 * the number of existential variables and parameters respectively.
 * The first lines starts with "E ", followed by a number.
 * The second lines starts with "P ", followed by a number.
 * These two lines are (optionally) followed by the names of the parameters.
 * The polytope is in PolyLib notation.
 */

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    POL_NO_DUAL
#else
#define MAXRAYS  600
#endif

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
#ifdef HAVE_OMEGA
    { "omega",   no_argument,  0,  'o' },
#endif
    { "pip",   no_argument,  0,  'p' },
    { "convert",   no_argument,  0,  'c' },
    { "floor",     no_argument,  0,  'f' },
    { "range-reduction",	no_argument,	0,  'R' },
    { "version",   no_argument,  0,  'V' },
    { 0, 0, 0, 0 }
};
#endif

#ifdef HAVE_OMEGA
#define OMEGA_OPT "o"

Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, char **parms)
{
    varvector varv;
    varvector paramv;
    Relation r = Polyhedron2relation(P, exist, nparam, parms);
    Polyhedron_Free(P);
    return relation2Domain(r, varv, paramv);
}
#else
#define OMEGA_OPT ""
Polyhedron *Omega_simplify(Polyhedron *P, 
			    unsigned exist, unsigned nparam, char **parms)
{
    return P;
}
#endif

int main(int argc, char **argv)
{
    Polyhedron *A;
    Matrix *M;
    char **param_name;
    int exist, nparam, nvar;
    char s[128];
    evalue *EP;
    int c, ind = 0;
    int range = 0;
    int convert = 0;
    int omega = 0;
    int pip = 0;
    int floor = 0;

    while ((c = getopt_long(argc, argv, OMEGA_OPT"pfcRV", options, &ind)) != -1) {
	switch (c) {
	case 'o':
	    omega = 1;
	    break;
	case 'p':
	    pip = 1;
	    break;
	case 'f':
	    floor = 1;
	    break;
	case 'c':
	    convert = 1;
	    break;
	case 'R':
	    range = 1;
	    break;
	case 'V':
	    printf(barvinok_version());
	    exit(0);
	    break;
	}
    }

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, MAXRAYS);
    Matrix_Free(M);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "E %d", &exist)<1))
	fgets(s, 128, stdin);

    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "P %d", &nparam)<1))
	fgets(s, 128, stdin);

    Polyhedron_Print(stdout, P_VALUE_FMT, A);
    printf("exist: %d, nparam: %d\n", exist, nparam);
    param_name = Read_ParamNames(stdin, nparam);
    nvar = A->Dimension - exist - nparam;
    if (omega) {
	A = Omega_simplify(A, exist, nparam, param_name);
	assert(!A->next);
	exist = A->Dimension - nvar - nparam;
    }
    if (pip && exist > 0)
	EP = barvinok_enumerate_pip(A, exist, nparam, MAXRAYS);
    else
	EP = barvinok_enumerate_e(A, exist, nparam, MAXRAYS);
    reduce_evalue(EP);
    evalue_combine(EP);
    if (range)
	evalue_range_reduction(EP);
    print_evalue(stdout, EP, param_name);
    if (floor) {
	fprintf(stderr, "WARNING: floor conversion not supported\n");
	evalue_frac2floor(EP);
	print_evalue(stdout, EP, param_name);
    } else if (convert) {
	evalue_mod2table(EP, nparam);
	print_evalue(stdout, EP, param_name);
    }
    free_evalue_refs(EP);
    free(EP);
    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    return 0;
}
