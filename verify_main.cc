#include <gmp.h>
#include <NTL/mat_ZZ.h>
extern "C" {
#include <polylib/polylibgmp.h>
#include "ev_operations.h"
}
#include <util.h>
#include <barvinok.h>
#include <barvinok2.h>

#include "verif_ehrhart.h"

#ifdef HAVE_GROWING_CHERNIKOVA
#define MAXRAYS    0
#else
#define MAXRAYS  600
#endif

#include "config.h"
#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "explicit",  no_argument,  0,  'e' },
    { "series",  no_argument,  0,  's' },
    { "verbose",  no_argument,  0,  'v' },
    { 0, 0, 0, 0 }
};
#endif

/* RANGE : normal range for evalutations (-RANGE -> RANGE) */
#define RANGE 50

/* SRANGE : small range for evalutations */
#define SRANGE 15

/* if dimension >= BIDDIM, use SRANGE */
#define BIGDIM 5

/* VSRANGE : very small range for evalutations */
#define VSRANGE 5

/* if dimension >= VBIDDIM, use VSRANGE */
#define VBIGDIM 8

int check_series(Polyhedron *S,Polyhedron *C, gen_fun *gf,
	         int nparam,int pos,Value *z)
{
    int k;
    Value c, tmp;
    value_init(c);
    value_init(tmp);

    if(pos == nparam) {
	/* Computes the coefficient */
	gf->coefficient(&z[S->Dimension-nparam+1], &c);

	/* if c=0 we may be out of context. */
	/* scanning is useless in this case*/
	if(!in_domain(C,&z[S->Dimension-nparam+1])) {
	    /* ok */ ;
	} else {

#ifdef PRINT_ALL_RESULTS
	    printf("EP( ");
	    value_print(stdout,VALUE_FMT,z[S->Dimension-nparam+1]);
	    for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
		printf(", ");
		value_print(stdout,VALUE_FMT,z[k]);
	    }
	    printf(" ) = ");
	    value_print(stdout,VALUE_FMT,c);
	    printf(" ");
#endif
	    /* Manually count the number of points */
	    count_points(1,S,z,&tmp);
#ifdef PRINT_ALL_RESULTS
	    printf(", count = ");
	    value_print(stdout, P_VALUE_FMT, tmp);
	    printf(". ");
#endif

	    if(value_ne(tmp,c)) {
		printf("\n"); 
		fflush(stdout);
		fprintf(stderr,"Error !\n");
		fprintf(stderr,"EP( ");
		value_print(stderr,VALUE_FMT,z[S->Dimension-nparam+1]);
		for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
		  fprintf(stderr,", ");
		  value_print(stderr,VALUE_FMT,z[k]);
		}
		fprintf(stderr," ) should be ");
		value_print(stderr,VALUE_FMT,tmp);
		fprintf(stderr,", while EP eval gives ");
		value_print(stderr,VALUE_FMT,c);
		fprintf(stderr,".\n");
#ifndef DONT_BREAK_ON_ERROR
		value_clear(c); value_clear(tmp);
		return 0;
#endif
	    }
#ifdef PRINT_ALL_RESULTS
	    else
		printf("OK.\n");
#endif
	}
    } else
	for(value_assign(tmp,Min); value_le(tmp,Max); 
				   value_increment(tmp,tmp)) {
#ifndef PRINT_ALL_RESULTS
	  k = VALUE_TO_INT(tmp);
	  if(!pos && !(k%st)) {
	    printf("o");
	    fflush(stdout);
	  }
#endif
	  value_assign(z[pos+S->Dimension-nparam+1],tmp);
	  if(!check_series(S,C,gf,nparam,pos+1,z)) {
	    value_clear(c); value_clear(tmp);
	    return(0);
	  }
	}

    value_clear(c);
    value_clear(tmp);
    return 1;
}

int main(int argc,char *argv[]) {
	
  Matrix *C1, *P1;
  Polyhedron *C, *P, *S;
  Polyhedron *CC, *PP;
  Enumeration *en;
  Value *p, tmp;
  int i,j;
  int m,M;
    int c, ind = 0;
    int series = 0;
    int verbose = 0;
    int function = 0;
    int result = 0;

/******* Read the input *********/
  P1 = Matrix_Read();
  C1 = Matrix_Read();

  if(C1->NbColumns < 2) {
    fprintf(stderr,"Not enough parameters !\n");
    exit(0);
  }
  
  P = Constraints2Polyhedron(P1,MAXRAYS);
  C = Constraints2Polyhedron(C1,MAXRAYS);
  params = Read_ParamNames(stdin, C->Dimension);
  Matrix_Free(C1);
  Matrix_Free(P1);

  /******* Read the options: initialize Min and Max ********/
  if(P->Dimension >= VBIGDIM)
    M = VSRANGE;
  else if(P->Dimension >= BIGDIM)
    M = SRANGE;
  else
    M = RANGE;
  m = -M;

    while ((c = getopt_long(argc, argv, "m:M:r:sve", options, &ind)) != -1) {
	switch (c) {
	case 'e':
	    function = 1;
	    break;
	case 's':
	    series = 1;
	    break;
	case 'v':
	    verbose = 1;
	    break;
	case 'm':
	    m = atoi(optarg);
	    break;
	case 'M':
	    M = atoi(optarg);
	    break;
	case 'r':
	    M = atoi(optarg);
	    m = -M;
	    break;
	}
    }

  if(m > M) {
    fprintf(stderr,"Nothing to do: Min > Max !\n");
    return(0);
  }
  value_init(Min);
  value_init(Max);
  value_set_si(Min,m);
  value_set_si(Max,M);
  value_init(tmp);

  /******* Compute true context *******/
  CC = align_context(C,P->Dimension,MAXRAYS);
  PP = DomainIntersection(P,CC,MAXRAYS);
  Domain_Free(CC);
  C1 = Matrix_Alloc(C->Dimension+1,P->Dimension+1);

  for(i=0;i<C1->NbRows;i++)
    for(j=0;j<C1->NbColumns;j++)
      if(i==j-P->Dimension+C->Dimension)
	value_set_si(C1->p[i][j],1);
      else
	value_set_si(C1->p[i][j],0);
  CC = Polyhedron_Image(PP,C1,MAXRAYS);
  Domain_Free(C);
  C = CC;

    gen_fun *gf = 0;

    /******* Compute EP *********/
    if (!series)
        en = barvinok_enumerate(P,C,MAXRAYS);
    else {
	evalue *EP;
	gf = barvinok_series(P, C, MAXRAYS);
	if (verbose) {
	    gf->print(C->Dimension, params);
	    puts("");
	}
	if (function) {
	    EP = *gf;
	    en =  partition2enumeration(EP);
	}
    }
  
  /******* Initializations for check *********/
  p = (Value *)malloc(sizeof(Value) * (P->Dimension+2));
  for(i=0;i<=P->Dimension;i++) {
    value_init(p[i]);
    value_set_si(p[i],0);
  }
  value_init(p[i]);
  value_set_si(p[i],1);

  /* S = scanning list of polyhedra */
  S = Polyhedron_Scan(P,C,MAXRAYS);

#ifndef PRINT_ALL_RESULTS
  if(C->Dimension > 0) {
    value_substract(tmp,Max,Min);
    if (VALUE_TO_INT(tmp) > 80)
      st = 1+(VALUE_TO_INT(tmp))/80;
    else
      st=1;
    for(i=VALUE_TO_INT(Min);i<=VALUE_TO_INT(Max);i+=st)
      printf(".");
    printf( "\r" );
    fflush(stdout);
  }
#endif

    /******* CHECK NOW *********/
    if(S) {
	if (!series || function) {
	    if (!check_poly(S,C,en,C->Dimension,0,p))
		result = -1;
	} else {
	    if (!check_series(S,C,gf,C->Dimension,0,p))
		result = -1;
	}
    }

    if (result == -1)
	fprintf(stderr,"Check failed !\n");
    
#ifndef PRINT_ALL_RESULTS
  printf( "\n" );
#endif
  
    if (gf)
	delete gf;

  for(i=0;i<=(P->Dimension+1);i++) 
    value_clear(p[i]);
  value_clear(tmp);
  Free_ParamNames(params, C->Dimension);
  return result;
} /* main */



