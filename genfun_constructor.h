#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
}
#include <barvinok/util.h>
#include <barvinok/genfun.h>
#include "reducer.h"
#include "bfcounter.h"

#ifdef NTL_STD_CXX
using namespace NTL;
#endif

/* base for generating function counting */
struct gf_base {
    np_base *base;
    gen_fun *gf;

    gf_base(np_base *npb, Polyhedron *context) : base(npb) {
	gf = new gen_fun(context);
    }
    virtual ~gf_base() {
    };
    static gf_base *create(Polyhedron *context, unsigned dim, unsigned nparam);

    void start_gf(Polyhedron *P, unsigned MaxRays) {
	base->start(P, MaxRays);
    }
};

struct partial_ireducer : public ireducer, public gf_base {
    partial_ireducer(Polyhedron *context, unsigned dim, unsigned nparam) : 
	    ireducer(dim), gf_base(this, context) {
	lower = nparam;
    }
    ~partial_ireducer() {
    }
    virtual void base(ZZ& c, ZZ& cd, vec_ZZ& num, mat_ZZ& den_f);
};

struct partial_reducer : public reducer, public gf_base {
    vec_ZZ lambda;
    vec_ZZ tmp;

    partial_reducer(Polyhedron *context, unsigned dim, unsigned nparam) : 
	    reducer(dim), gf_base(this, context) {
	lower = nparam;

	tmp.SetLength(dim - nparam);
    }
    virtual void init(Polyhedron *P) {
	randomvector(P, lambda, dim - lower);
    }
    ~partial_reducer() {
    }
    virtual void base(ZZ& c, ZZ& cd, vec_ZZ& num, mat_ZZ& den_f);

    virtual void split(vec_ZZ& num, ZZ& num_s, vec_ZZ& num_p,
		       mat_ZZ& den_f, vec_ZZ& den_s, mat_ZZ& den_r);
};

struct partial_bfcounter : public bfcounter_base, public gf_base {
    partial_bfcounter(Polyhedron *context, unsigned dim, unsigned nparam) : 
	    bfcounter_base(dim), gf_base(this, context) {
	lower = nparam;
    }
    ~partial_bfcounter() {
    }
    virtual void base(mat_ZZ& factors, bfc_vec& v);
};
