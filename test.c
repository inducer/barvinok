#include <stdio.h>
#include <polylib/polylibgmp.h>
#include <barvinok.h>

int main()
{
    int i, nbPol, nbMat, func;
    Polyhedron *A, *B, *C;
    char s[128];

    nbPol = nbMat = 0;
    fgets(s, 128, stdin);
    while ((*s=='#') ||
	    ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
	fgets(s, 128, stdin);

    for (i = 0; i < nbPol; ++i) {
	Matrix *M = Matrix_Read();
	A = Constraints2Polyhedron(M, 600);
	Matrix_Free(M);
	fgets(s, 128, stdin);
	while ((*s=='#') || (sscanf(s, "F %d", &func)<1))
	    fgets(s, 128, stdin);

	switch(func) {
	case 1:
	    Polyhedron_Print(stdout, P_VALUE_FMT, A);
	    B = Polyhedron_Polar(A, 600);
	    break;
	}
    }

    return 0;
}
