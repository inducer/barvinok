#ifndef DECOMPOSER_H
#define DECOMPOSER_H

#include <gmp.h>
#include <barvinok/polylib.h>
#include <barvinok/options.h>

struct signed_cone_consumer {
    virtual void handle(Polyhedron *P, int sign) = 0;
};

struct vertex_decomposer {
    Polyhedron *P;
    unsigned nbV;	// number of vertices
    Param_Vertices *V;	// current vertex
    int vert;		// current vertex index
    signed_cone_consumer& scc;

    vertex_decomposer(Polyhedron *P, unsigned nbV, signed_cone_consumer& scc) :
			P(P), nbV(nbV), scc(scc) {}
    void decompose_at_vertex(Param_Vertices *V, int _i, barvinok_options *options);
};

void barvinok_decompose(Polyhedron *C, signed_cone_consumer& scc,
			barvinok_options *options);

#endif
