#include <unistd.h>
#include <sys/times.h>
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
#include <util.h>
#include <barvinok.h>

/* The input of this example program is a polytope in combined
 * data and parameter space followed by two lines indicating
 * the number of existential variables and parameters respectively.
 * The first lines starts with "E ", followed by a number.
 * The second lines starts with "P ", followed by a number.
 * These two lines are (optionally) followed by the names of the parameters.
 * The polytope is in PolyLib notation.
 */

int main()
{
    Polyhedron *A;
    Matrix *M;
    char **param_name;
    int exist, nparam;
    char s[128];
    evalue *EP;

    M = Matrix_Read();
    A = Constraints2Polyhedron(M, 600);
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
    EP = barvinok_enumerate_e(A, exist, nparam, 600);
    print_evalue(stdout, EP, param_name);
    free_evalue_refs(EP);
    Free_ParamNames(param_name, nparam);
    Polyhedron_Free(A);
    return 0;
}
