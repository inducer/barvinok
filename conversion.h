#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

void value2zz(Value v, ZZ& z);
void zz2value(const ZZ& z, Value& v);
void values2zz(Value *p, vec_ZZ& v, int len);
void zz2values(const vec_ZZ& v, Value *p);
void matrix2zz(Matrix *M, mat_ZZ& m, unsigned nr, unsigned nc);
Matrix *rays(Polyhedron *C);
Matrix *rays2(Polyhedron *C);
void randomvector(Polyhedron *P, vec_ZZ& lambda, int nvar);
