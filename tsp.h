#ifndef _TSP_H
#define _TSP_H

#include <cplex.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <omp.h>

#define VERB_PRINT_FILE_INPUT_NODES 0
#define VERB_PRINT_CLINE_ARGUMENTS 0
#define VERB_PRINT_CPLEX_PROBLEM 0
#define VERB_PRINT_CPLEX_SOLUTION 0
#define VERB_PRINT_SEC_ITERATIONS 1
#define VERB_PRINT_CPLEX_INFOS 0
#define VERB_PRINT_COMP_CUT 1
#define VERB_PRINT_SOL_INFO 1

#define GNUPLOT_OUTPUT 1
#define GNUPLOT_PROGRESS 0
#define GNUPLOT_EXECUTABLE "gnuplot"
#define GNUPLOT_LINECOLOR "#0080cd"

#define STRING_MAX_SIZE 1024
#define SHORT_STRING_MAX_SIZE 32
#define COMMAND_LINE_DEFAULT_SIZE 65536
#define TESTSET_MAX_SIZE 512
#define GRASP_DISTRIBUTION_VAL1 0.33
#define GRASP_DISTRIBUTION_VAL2 0.33
#define GRASP_DISTRIBUTION_VAL3 0.33
#define CHECKS_NO_LIMITS -1

#define INT_ROUND(x) ((int)((x) + 0.5))
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define INLINE_SQUARED_DISTANCE(x1,x2,y1,y2) (((x1)-(x2))*((x1)-(x2))+((y1)-(y2))*((y1)-(y2)))
#define TRUE 1
#define FALSE 0

typedef struct
{
	//input data
	int nnodes;
	double* xcoord;
	double* ycoord;

	// execution parameters
	double time_limit; 
	int random_seed;

	// parameters 
	int model_type;
	char input_file[STRING_MAX_SIZE + 1];	
	int integer_costs;
	char* cplex_model_file;
	char cplex_write_model;

	// model-specific     
	int xstart;
	int ustart;
	int ystart;
	int neighbors;							// number of neighbors to be considered in model9
	double lambda;
} instance;

typedef struct
{
	int index;
	int node;
} node_pair;

typedef struct
{
	int node;
	double dist;
} dist_pair_t;

typedef struct
{
	instance* inst;
	int* succ;
	int* comp;
	int ncomp;
	node_pair* nodes_order;
	FILE* gp_pipe;
	int sec_steps;
	int* tour;
	void (*heuristic_method)(void*);
	int trials;
	double* xstar;						// best sol. available
	double	best_lb;						// best lower bound available
	double latest_cost;
	double execution_time;					// total execution time
	double latest_improvement;
	int checks;
} solution;

typedef struct
{
	solution* sol;
	CPXENVptr env;
	CPXLPptr lp;
} data_wrap;

void	add_sec(solution* sol, CPXENVptr env, CPXLPptr lp);
void	add_sec_cb(solution* sol, int nnodes, CPXCENVptr env, void* cbdata, int wherefrom);
void	build_model(instance* inst, CPXENVptr env, CPXLPptr lp);
void	build_model_0(instance* inst, CPXENVptr env, CPXLPptr lp);
void	build_model_1(instance* inst, CPXENVptr env, CPXLPptr lp);
void	build_model_2(instance* inst, CPXENVptr env, CPXLPptr lp);
void	build_model_3(instance* inst, CPXENVptr env, CPXLPptr lp);

void	build_sol(const double* xstar, solution* sol);
void	build_sol_1(const double* xstar, solution* sol);
void	build_sol_2(const double* xstar, solution* sol);

void	close_instance(instance* inst);
void	close_gp_pipe(solution* sol);
void	close_solution(solution* sol);
int		comparator_nodes(const void* n1, const void* n2);
double	dist(int i, int j, const instance* inst);
int		find_in_distribution(double* v, int size, double x);
void	first_nearest_neighbor(solution* sol);
void	first_nearest_neighbor_wmarkers(solution* sol);
double	get_cost(solution* sol);
char*	get_model_name(int m);
void	get_sol(solution* sol, CPXENVptr env, CPXLPptr lp);
void	get_xstar_format_sol(solution* sol);
void	GRASP(solution* sol);
void	greedy_refinement_2opt(solution* sol);
void	greedy_refinement_3opt(solution* sol);
void	gnuplot_plt(solution* sol);
void	gnuplot_plt_0(solution* sol, char* plot_command, char* style);
void	gnuplot_plt_1(solution* sol, char* plot_command, char* style);
void	gnuplot_plt_2(solution* sol, char* plot_command, char* style);

void	hard_fixing(solution* sol, CPXENVptr env, CPXLPptr lp);
int		hf_set_bounds(solution* sol, CPXENVptr env, CPXLPptr lp, double prob_of_fixing);
void	identify_connected_components(const double* xstar, solution* sol, int nnodes);
void	insertion_method(solution* sol);
void	k_nearest_neighbors(solution* sol);
void	local_branching(solution* sol, CPXENVptr env, CPXLPptr lp);
int		lb_set_constraint(solution* sol, CPXENVptr env, CPXLPptr lp, int k, int replace_constraint);

void	loop_callback(solution* sol, CPXENVptr env, CPXLPptr lp);
void	loop(solution* sol, CPXENVptr env, CPXLPptr lp);
void	loop_gap(solution* sol, CPXENVptr env, CPXLPptr lp);

void	multi_start(solution* sol);
void	multi_start_parallel(solution* sol);
double	obtain_SA_mean_iteration_per_sec(solution* sol, double initial_temperature, double final_temperature, double normalizing_factor);
void	print_cplex_solution(const double* xstar, solution* sol);
void	print_error(const char* err);
void	provide_heuristic_mipstart(solution* sol, CPXENVptr env, CPXLPptr lp);
void	random_solution(solution* sol);

void	refine(solution* sol);
void	setup_gp_pipe(solution* sol);
void	setup_solution(solution* sol);

void	simulated_annealing(solution* sol);
double	average_distance(instance* inst);

void	solve_lp(solution* sol, CPXENVptr env, CPXLPptr lp);
static int CPXPUBLIC subtour_elimination_lazycb(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
static inline void swap(int* a, int* b);
int		TSPopt(solution* sol);

void	VNS(solution* sol);
void	diversificate(solution* sol, int kick_size);

void	greedy_opt_search(solution* sol);

long	xpos(long i, long j, long n);
int		xpos_asymm(int i, int j, int n);
int		upos_mtz(int i, int ustart);
int		ypos_gg(int i, int j, int n, int ystart);
#endif