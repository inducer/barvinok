#include <assert.h>
#include <iostream>
#include <bernstein/bernstein.h>
#include <barvinok/evalue.h>
#include <barvinok/options.h>
#include <barvinok/util.h>
#include <barvinok/bernstein.h>
#include "argp.h"
#include "progname.h"
#include "evalue_convert.h"
#include "evalue_read.h"
#include "verify.h"
#include "range.h"

using std::cout;
using std::cerr;
using std::endl;
using namespace GiNaC;
using namespace bernstein;
using namespace barvinok;

#define ALLOCN(type,n) (type*)malloc((n) * sizeof(type))

#define OPT_VARS  	    (BV_OPT_LAST+1)
#define OPT_SPLIT  	    (BV_OPT_LAST+2)
#define OPT_LOWER  	    (BV_OPT_LAST+3)
#define OPT_METHOD  	    (BV_OPT_LAST+4)

struct argp_option argp_options[] = {
    { "split",		    OPT_SPLIT,	"int" },
    { "variables",	    OPT_VARS,  	"list",	0,
	"comma separated list of variables over which to compute a bound" },
    { "lower",	    	    OPT_LOWER, 	0, 0,	"compute lower bound instead of upper bound"},
    { "optimization-method",	OPT_METHOD,	"bernstein|propagation",
	0, "optimization method to use" },
    { 0 }
};

struct options {
    struct convert_options   convert;
    struct verify_options    verify;
    char* var_list;
    int split;
    int lower;
#define	METHOD_BERNSTEIN	0
#define METHOD_PROPAGATION	1
    int method;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
    struct options *options = (struct options*) state->input;

    switch (key) {
    case ARGP_KEY_INIT:
	state->child_inputs[0] = &options->convert;
	state->child_inputs[1] = &options->verify;
	state->child_inputs[2] = options->verify.barvinok;
	options->var_list = NULL;
	options->split = 0;
	options->lower = 0;
	options->method = METHOD_BERNSTEIN;
	break;
    case OPT_VARS:
	options->var_list = strdup(arg);
	break;
    case OPT_SPLIT:
	options->split = atoi(arg);
	break;
    case OPT_LOWER:
	options->lower = 1;
	break;
    case OPT_METHOD:
	if (!strcmp(arg, "bernstein"))
	    options->method = METHOD_BERNSTEIN;
	else if (!strcmp(arg, "propagation"))
	    options->method = METHOD_PROPAGATION;
	else
	    argp_error(state, "unknown value for --optimization-method option");
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static int check_poly_max(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options);

struct check_poly_max_data : public check_EP_data {
    piecewise_lst	 *pl;

    check_poly_max_data(evalue *EP, piecewise_lst *pl) : pl(pl) {
	this->EP = EP;
	this->cp.check = check_poly_max;
    }
};

static int check_poly_max(const struct check_poly_data *data,
			  int nparam, Value *z,
			  const struct verify_options *options)
{
    int k;
    int ok;
    const check_poly_max_data *max_data;
    max_data = (const check_poly_max_data *)data;
    const char *minmax;
    Value m, n, d;
    value_init(m);
    value_init(n);
    value_init(d);
    int sign;

    if (options->barvinok->bernstein_optimize == BV_BERNSTEIN_MAX) {
	minmax = "max";
	sign = 1;
    } else {
	minmax = "min";
	sign = -1;
    }

    max_data->pl->evaluate(nparam, z, &n, &d);
    if (sign > 0)
	mpz_fdiv_q(m, n, d);
    else
	mpz_cdiv_q(m, n, d);

    if (options->print_all) {
	printf("%s(", minmax);
	value_print(stdout, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    printf(", ");
	    value_print(stdout, VALUE_FMT, z[k]);
	}
	printf(") = ");
	value_print(stdout, VALUE_FMT, n);
	if (value_notone_p(d)) {
	    printf("/");
	    value_print(stdout, VALUE_FMT, d);
	}
	printf(" (");
	value_print(stdout, VALUE_FMT, m);
	printf(")");
    }

    evalue_optimum(max_data, &n, sign);

    if (sign > 0)
	ok = value_ge(m, n);
    else
	ok = value_le(m, n);

    if (options->print_all) {
	printf(", %s(EP) = ", minmax);
	value_print(stdout, VALUE_FMT, n);
	printf(". ");
    }

    if (!ok) {
	printf("\n"); 
	fflush(stdout);
	fprintf(stderr, "Error !\n");
	fprintf(stderr, "%s(", minmax);
	value_print(stderr, VALUE_FMT, z[0]);
	for (k = 1; k < nparam; ++k) {
	    fprintf(stderr,", ");
	    value_print(stderr, VALUE_FMT, z[k]);
	}
	fprintf(stderr, ") should be ");
	if (sign > 0)
	    fprintf(stderr, "greater");
	else
	    fprintf(stderr, "smaller");
	fprintf(stderr, " than or equal to ");
	value_print(stderr, VALUE_FMT, n);
	fprintf(stderr, ", while pl eval gives ");
	value_print(stderr, VALUE_FMT, m);
	fprintf(stderr, ".\n");
	cerr << *max_data->pl << endl;
    } else if (options->print_all)
	printf("OK.\n");

    value_clear(m);
    value_clear(n);
    value_clear(d);

    return ok;
}

static int verify(piecewise_lst *pl, evalue *EP, unsigned nvar, unsigned nparam,
		  struct verify_options *options)
{
    check_poly_max_data data(EP, pl);
    return !check_EP(&data, nvar, nparam, options);
}

static int optimize(evalue *EP, char **all_vars, unsigned nvar, unsigned nparam,
		    struct options *options)
{
    Polyhedron *U;
    piecewise_lst *pl = NULL;
    U = Universe_Polyhedron(nparam);
    int print_solution = 1;
    int result = 0;

    exvector params;
    params = constructParameterVector(all_vars+nvar, nparam);

    if (options->verify.verify) {
	verify_options_set_range(&options->verify, nvar+nparam);
	if (!options->verify.barvinok->verbose)
	    print_solution = 0;
    }

    if (options->lower)
	options->verify.barvinok->bernstein_optimize = BV_BERNSTEIN_MIN;
    else
	options->verify.barvinok->bernstein_optimize = BV_BERNSTEIN_MAX;
    if (options->method == METHOD_BERNSTEIN) {
	pl = evalue_bernstein_coefficients(NULL, EP, U, params,
					   options->verify.barvinok);
	if (options->lower)
	    pl->minimize();
	else
	    pl->maximize();
    } else
	pl = evalue_range_propagation(NULL, EP, params,
				      options->verify.barvinok);
    assert(pl);
    if (print_solution)
	cout << *pl << endl;
    if (options->verify.verify) {
	result = verify(pl, EP, nvar, nparam, &options->verify);
    }
    delete pl;

    Polyhedron_Free(U);

    return result;
}

int main(int argc, char **argv)
{
    evalue *EP;
    char **all_vars = NULL;
    unsigned nvar;
    unsigned nparam;
    struct options options;
    struct barvinok_options *bv_options = barvinok_options_new_with_defaults();
    static struct argp_child argp_children[] = {
	{ &convert_argp,    	0,	"input conversion",	1 },
	{ &verify_argp,    	0,	"verification",		2 },
	{ &barvinok_argp,    	0,	"barvinok options",	3 },
	{ 0 }
    };
    static struct argp argp = { argp_options, parse_opt, 0, 0, argp_children };
    int result = 0;

    options.verify.barvinok = bv_options;
    set_program_name(argv[0]);
    argp_parse(&argp, argc, argv, 0, 0, &options);

    EP = evalue_read_from_file(stdin, options.var_list, &all_vars,
			       &nvar, &nparam, bv_options->MaxRays);
    assert(EP);

    if (options.split)
	evalue_split_periods(EP, options.split, bv_options->MaxRays);

    evalue_convert(EP, &options.convert, bv_options->verbose, nparam, all_vars);

    if (EVALUE_IS_ZERO(*EP))
	print_evalue(stdout, EP, all_vars);
    else
	result = optimize(EP, all_vars, nvar, nparam, &options);

    evalue_free(EP);

    if (options.var_list)
	free(options.var_list);
    Free_ParamNames(all_vars, nvar+nparam);
    barvinok_options_free(bv_options);
    return result;
}