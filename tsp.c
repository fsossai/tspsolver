#include "tsp.h"
#include <time.h>

char benchmark = 0;

void add_sec(solution* sol, CPXENVptr env, CPXLPptr lp)
{
	char name[STRING_MAX_SIZE + 1];
	sprintf(name, "subtour_elimination");
	char* cname[] = { name };

	sprintf(cname[0], "subtour_elimination");
	int n = sol->inst->nnodes;
	int lastrow;
	double rhs;
	char sense = 'L';
	node_pair* pairs = (node_pair*)malloc(n * sizeof(node_pair));
	int* comp_counter = (int*)calloc(sol->ncomp + 1, sizeof(int));

	for (int i = 0; i < sol->inst->nnodes; i++)
	{
		pairs[i].index = sol->comp[i];
		pairs[i].node = i;
		comp_counter[sol->comp[i]]++;
	}
	qsort(pairs, sol->inst->nnodes, sizeof(node_pair), comparator_nodes);

	// i keeps track of the first node in the k-th conn. component
	for (int k = 1, i = 0; k <= sol->ncomp; k++)
	{
		rhs = comp_counter[k] - 1.0;
		lastrow = CPXgetnumrows(env, lp);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [subtour_elimination]");

		// enumerating all archs in the k-th conn. component
		for (int j1 = 0; j1 < comp_counter[k] - 1; j1++)
		{
			// j2 swipe all nodes in the conn. comp. to get every possibile pair
			for (int j2 = j1 + 1; j2 < comp_counter[k]; j2++)
			{
				if (CPXchgcoef(env, lp, lastrow, xpos(pairs[i + j1].node, pairs[i + j2].node, n), 1.0))
					print_error(" wrong CPXchgcoef [subtour_elimination]");
			}
		}
		i += comp_counter[k];
	}

	free(pairs);
	free(comp_counter);
}

void add_sec_cb(solution * sol, int nnodes, CPXCENVptr env, void* cbdata, int wherefrom)
{
	double rhs;
	char sense = 'L';
	node_pair* pairs = (node_pair*)malloc(nnodes * sizeof(node_pair));
	int* comp_counter = (int*)calloc(sol->ncomp + 1, sizeof(int));
	int* cut_indexes = NULL;
	double* values = NULL;
	int max_nodes_in_comp = 0;

	for (int i = 0; i < nnodes; i++)
	{
		pairs[i].index = sol->comp[i];
		pairs[i].node = i;
		comp_counter[sol->comp[i]]++;		// counting how many nodes for every conn. component
		max_nodes_in_comp = max(max_nodes_in_comp, comp_counter[sol->comp[i]]);
	}

	/* allocate space for every pair of nodes in the maxmimal connected component
	and thus enough for any other connected component.
	the number of pair of nodes = (max_nodes_in_comp CHOOSE 2)
	*/
	int max_archs = (max_nodes_in_comp * (max_nodes_in_comp - 1)) / 2;
	cut_indexes = (int*)malloc(max_archs * sizeof(int));
	values = (double*)malloc(max_archs * sizeof(double));
	for (int i = 0; i < max_archs; i++)
		values[i] = 1.0;

	qsort(pairs, nnodes, sizeof(node_pair), comparator_nodes);

	// i keeps track of the first node in the k-th conn. component
	for (int k = 1, i = 0, nzcnt; k <= sol->ncomp; k++)
	{
		rhs = comp_counter[k] - 1.0;
		nzcnt = 0;	// non-zero variables counter
		// creating a new cut
		// enumerating all archs in the k-th conn. component
		for (int j1 = 0; j1 < comp_counter[k] - 1; j1++)
		{
			// j2 swipe all nodes in the conn. comp. to get every possibile pair
			for (int j2 = j1 + 1; j2 < comp_counter[k]; j2++)
			{
				// an arch has been identified
				cut_indexes[nzcnt++] = xpos(pairs[i + j1].node, pairs[i + j2].node, nnodes);
			}
		}
		// applying cut
		if (CPXcutcallbackadd(env, cbdata, wherefrom, nzcnt, rhs, sense, cut_indexes, values, 0))
			print_error("CPXcutcallbackadd() error");

		i += comp_counter[k];
	}
	free(pairs);
	free(comp_counter);
	free(cut_indexes);
}

void build_model(instance * inst, CPXENVptr env, CPXLPptr lp)
{
	switch (inst->model_type)
	{
	case 0:
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
		build_model_0(inst, env, lp);
		break;
	case 1:
		build_model_1(inst, env, lp);
		break;
	case 2:
		build_model_2(inst, env, lp);
		break;
	case 3:
		build_model_3(inst, env, lp);
		break;
	case 9:
	case 10:
	case 11:
	case 12:
	case 13:
	case 14:
	case 15:
		// Nothing to do here
		break;
	default:
		print_error("unknown model type!");
		break;
	}
}

void build_model_0(instance * inst, CPXENVptr env, CPXLPptr lp)
{
	int nnodes = inst->nnodes;
	int nvars = nnodes * (nnodes - 1) / 2;

	printf("Building basic model ...\n");
	printf(" Adding binary variables ...\n");

	double* obj = (double*)calloc(nvars, sizeof(double));
	double* lb = (double*)calloc(nvars, sizeof(double));
	double* ub = (double*)malloc(nvars * sizeof(double));
	char* xctype = (char*)malloc(nvars * sizeof(char));
	char** names_pp = (char*)calloc(nvars, sizeof(char*));
	char* name_p = (char*)calloc(nvars * SHORT_STRING_MAX_SIZE, sizeof(char));
	char* name_pbackup = name_p;

	int k = 0;
	for (int i = 0; i < nnodes; i++)
	{
		for (int j = i + 1; j < nnodes; j++)
		{
			obj[k] = dist(i, j, inst);
			ub[k] = 1.0;
			xctype[k] = 'B';
			sprintf(name_p, "x(%d,%d)", i + 1, j + 1);
			names_pp[k] = name_p;
			name_p += SHORT_STRING_MAX_SIZE;
			k++;
		}
	}
	if (CPXnewcols(env, lp, nvars, obj, lb, ub, xctype, names_pp))
		print_error(" wrong CPXnewcols on x var.s");

	free(obj);
	free(lb);
	free(ub);
	free(xctype);
	free(names_pp);
	free(name_pbackup);

	char name[SHORT_STRING_MAX_SIZE];
	char* cname[] = { name };

	printf(" Adding degree constraints: %.2lf%%", 0.0);
	for (int h = 0; h < inst->nnodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 2.0;
		char sense = 'E';
		sprintf(name, "degree(%d)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree]");
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h)
				continue;
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst->nnodes), 1.0))
				print_error(" wrong CPXchgcoef [degree]");
		}
		printf("\r Adding degree constraints: %.2lf%%", (double)h / (double)inst->nnodes * 100.0);
	}
	printf("\r Adding degree constraints: %.2lf%%", 100.0);

	if (inst->cplex_write_model)
		CPXwriteprob(env, lp, inst->cplex_model_file, NULL);

	printf("\nModel built\n");
}

void build_model_1(instance * inst, CPXENVptr env, CPXLPptr lp)
{
	int n = inst->nnodes;
	double zero = 0.0;
	double one = 1.0;
	char binary = 'B';
	char integer = 'I';

	char** cname = (char**)calloc(1, sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary var.s x(i,j) 

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			double obj = dist(i, j, inst);
			double ub = (i == j) ? zero : one; //eliminating self-loops
			double lb = zero;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
				print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos_asymm(i, j, n))
				print_error(" wrong position for x var.s");
		}
	}

	// Add degree constraints (archs entering the node h)
	for (int h = 0; h < n; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = 'E';                            // 'E' for equality constraint
		sprintf(cname[0], "degree_in(%d)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree_in]");

		for (int i = 0; i < n; i++)
		{
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(i, h, n), 1.0))
				print_error(" wrong CPXchgcoef [degree_in]");
		}
	}
	// Add degree constraints (archs leaving the node h)
	for (int h = 0; h < n; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		double rhs = 1.0;
		char sense = 'E';                            // 'E' for equality constraint
		sprintf(cname[0], "degree_out(%d)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree_out]");

		for (int i = 0; i < n; i++)
		{
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(h, i, n), 1.0))
				print_error(" wrong CPXchgcoef [degree_out]");
		}
	}

	inst->ustart = CPXgetnumcols(env, lp);
	// Add continuous variables u(1) = 0
	sprintf(cname[0], "u(1)");
	if (CPXnewcols(env, lp, 1, &zero, &zero, &zero, &integer, cname))
		print_error(" wrong CPXnewcols on u var.s");

	// Add continuous variables u(i), i = 2, ..., n
	for (int i = 1; i < n; i++)
	{
		sprintf(cname[0], "u(%d)", i + 1);
		double obj = zero;
		double ub = n - 2.0;
		double lb = zero;
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
			print_error(" wrong CPXnewcols on u var.s");
		if (CPXgetnumcols(env, lp) - 1 != upos_mtz(i, inst->ustart))
			print_error(" wrong position for u var.s");
	}
	// Add subtours elimination constraints: 
	// u(i) - u(j) + n*x(i,j) <= n - 1,  i,j in {2,...,n}, i!=j
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			if (i == j)
				continue;

			int lastrow = CPXgetnumrows(env, lp);
			double rhs = n - 1.0;
			char sense = 'L';                            // Less or equal
			sprintf(cname[0], "subtour(%d,%d)", i + 1, j + 1);

			// Creating constraint
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
				print_error(" wrong CPXnewrows [subtour]");
			// Setting constraint coefficients
			if (CPXchgcoef(env, lp, lastrow, upos_mtz(i, inst->ustart), 1.0))
				print_error(" wrong CPXchgcoef [subtour]");
			if (CPXchgcoef(env, lp, lastrow, upos_mtz(j, inst->ustart), -1.0))
				print_error(" wrong CPXchgcoef [subtour]");
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(i, j, n), n))
				print_error(" wrong CPXchgcoef [subtour]");
		}
	}

	if (inst->cplex_write_model)
		CPXwriteprob(env, lp, inst->cplex_model_file, NULL);

	free(cname[0]);
	free(cname);
}

void build_model_2(instance * inst, CPXENVptr env, CPXLPptr lp)
{
	int n = inst->nnodes;
	int lastrow;
	double zero = 0.0;
	double one = 1.0;
	char binary = 'B';
	char integer = 'I';

	double obj;
	char sense;
	double rhs;
	double ub;
	double lb;

	char** cname = (char**)calloc(1, sizeof(char*));
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary var.s x(i,j) 
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			obj = dist(i, j, inst);
			ub = (i == j) ? zero : one; //eliminating self-loops
			lb = zero;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
				print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos_asymm(i, j, n))
				print_error(" wrong position for x var.s");
		}
	}

	// Add degree constraints (archs entering the node h)
	rhs = 1.0;
	sense = 'E';
	for (int h = 0; h < n; h++)
	{
		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "degree_in(%d)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree_in]");

		for (int i = 0; i < n; i++)
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(i, h, n), 1.0))
				print_error(" wrong CPXchgcoef [degree_in]");
	}

	// Add degree constraints (archs leaving the node h)
	rhs = 1.0;
	sense = 'E';
	for (int h = 0; h < n; h++)
	{
		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "degree_out(%d)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree_out]");

		for (int i = 0; i < n; i++)
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(h, i, n), 1.0))
				print_error(" wrong CPXchgcoef [degree_out]");
	}

	// Add integer variables y(i,j)
	inst->ystart = CPXgetnumcols(env, lp);
	obj = zero;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);
			ub = (i == j || j == 0) ? zero : n - 1.0;
			lb = zero;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
				print_error(" wrong CPXnewcols on y var.s");
			if (CPXgetnumcols(env, lp) - 1 != ypos_gg(i, j, n, inst->ystart))
				print_error(" wrong position for y var.s");
		}
	}

	// Adding special constraints for y(1,j)
	rhs = n - 1.0;
	lastrow = CPXgetnumrows(env, lp);
	sense = 'E';
	sprintf(cname[0], "sink_flow");
	if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
		print_error(" wrong CPXnewrows [sink_flow]");
	for (int j = 1; j < n; j++)
		if (CPXchgcoef(env, lp, lastrow, ypos_gg(0, j, n, inst->ystart), 1.0))
			print_error(" wrong CPXchgcoef [sink_flow]");

	// Adding flow constraints for y(i,j) foreach i,j i!=0, j!=0
	rhs = 1.0;
	sense = 'E';
	for (int h = 1; h < n; h++)
	{
		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "flow(%i)", h + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [flow]");
		for (int i = 0; i < n; i++)
		{
			if (h == i)
				continue;
			if (CPXchgcoef(env, lp, lastrow, ypos_gg(i, h, n, inst->ystart), 1.0))
				print_error(" wrong CPXchgcoef [flow]");
			if (CPXchgcoef(env, lp, lastrow, ypos_gg(h, i, n, inst->ystart), -1.0))
				print_error(" wrong CPXchgcoef [flow]");
		}
	}

	// big-m constraints
	rhs = 0.0;
	sense = 'L';
	for (int j = 1; j < n; j++)
	{
		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "sink_bigm_xy(1,%d)", j + 1);
		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [sink_bigm_xy]");
		if (CPXchgcoef(env, lp, lastrow, xpos_asymm(0, j, n), 1.0 - n) ||
			CPXchgcoef(env, lp, lastrow, ypos_gg(0, j, n, inst->ystart), 1.0))
			print_error(" wrong CPXchgcoef [sink_bigm_xy]");
	}
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			if (i == j)
				continue;

			lastrow = CPXgetnumrows(env, lp);
			sprintf(cname[0], "bigm_xy(%d,%d)", i + 1, j + 1);

			// Creating constraint
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
				print_error(" wrong CPXnewrows [bigm_xy]");
			// Setting constraint coefficients
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(i, j, n), 2.0 - n) ||
				CPXchgcoef(env, lp, lastrow, ypos_gg(i, j, n, inst->ystart), 1.0))
				print_error(" wrong CPXchgcoef [bigm_xy]");
		}
	}

	if (inst->cplex_write_model)
		CPXwriteprob(env, lp, inst->cplex_model_file, NULL);

	free(cname[0]);
	free(cname);
}

void build_model_3(instance * inst, CPXENVptr env, CPXLPptr lp)
{
	int n = inst->nnodes;
	int lastrow;
	char** cname = (char**)calloc(1, sizeof(char*));
	cname[0] = (char*)calloc(100, sizeof(char));

	double zero = 0.0;
	double one = 1.0;
	char binary = 'B';
	char integer = 'I';

	double obj;
	char sense;
	double ub;
	double lb;
	double rhs;

	// Add binary var.s x(i,j) 
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);
			obj = dist(i, j, inst);
			ub = (i == j) ? zero : one; //eliminating self-loops
			lb = zero;

			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname))
				print_error(" wrong CPXnewcols on x var.s");
			if (CPXgetnumcols(env, lp) - 1 != xpos_asymm(i, j, n))
				print_error(" wrong position for x var.s");
		}
	}

	// Add degree constraints (archs entering the node h)
	rhs = 1.0;
	sense = 'E';
	for (int h = 0; h < n; h++)
	{
		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "degree_in(%d)", h + 1);

		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree_in]");

		for (int i = 0; i < n; i++)
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(i, h, n), 1.0))
				print_error(" wrong CPXchgcoef [degree_in]");
	}

	// Add degree constraints (archs leaving the node h)
	rhs = 1.0;
	sense = 'E';
	for (int h = 0; h < n; h++)
	{
		lastrow = CPXgetnumrows(env, lp);
		sprintf(cname[0], "degree_out(%d)", h + 1);

		if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname))
			print_error(" wrong CPXnewrows [degree_out]");

		for (int i = 0; i < n; i++)
			if (CPXchgcoef(env, lp, lastrow, xpos_asymm(h, i, n), 1.0))
				print_error(" wrong CPXchgcoef [degree_out]");
	}

	// Adding second set of variables
	inst->ustart = CPXgetnumcols(env, lp);
	// Add integer variables u(1) = 0
	sprintf(cname[0], "u(1)");
	if (CPXnewcols(env, lp, 1, &zero, &zero, &zero, &integer, cname))
		print_error(" wrong CPXnewcols on u var.s");

	// Add integer variables u(i), i = 2,...,n
	obj = zero;
	ub = n - 2.0;
	lb = zero;
	for (int i = 1; i < n; i++)
	{
		sprintf(cname[0], "u(%d)", i + 1);

		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname))
			print_error(" wrong CPXnewcols on u var.s");
		if (CPXgetnumcols(env, lp) - 1 != upos_mtz(i, inst->ustart))
			print_error(" wrong position for u var.s");
	}

	// Add subtours elimination constraints: 
	// u(i) - u(j) + M*x(i,j) <= M - 1,  i,j in {2,...,n}, i!=j
	int index[3];
	double value[3];
	int izero = 0;
	int M = n;
	rhs = M - 1.0;
	sense = 'L';

	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
		{
			if (i == j)
				continue;

			lastrow = CPXgetnumrows(env, lp);
			sprintf(cname[0], "u_consistency_arc(%d,%d)", i + 1, j + 1);

			index[0] = upos_mtz(i, inst->ustart);
			value[0] = 1.0;
			index[1] = upos_mtz(j, inst->ustart);
			value[1] = -1.0;
			index[2] = xpos_asymm(i, j, n);
			value[2] = M;

			if (CPXaddlazyconstraints(env, lp, 1, 3, &rhs, &sense, &izero, index, value, cname))
				print_error(" wrong CPXaddlazyconstraints [u-consistency]");
		}
	}

	// Add constraints x(i,j) + x(j,i) <= 1 for all i,j in 1,...,n i<j
	int index_2[2];
	double value_2[2];
	izero = 0;
	rhs = 1.0;
	sense = 'L';

	for (int i = 0; i < n; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			lastrow = CPXgetnumrows(env, lp);
			sprintf(cname[0], "two_nodes_loop(%d,%d)", i + 1, j + 1);

			index_2[0] = xpos_asymm(i, j, n);
			value_2[0] = 1.0;
			index_2[1] = xpos_asymm(j, i, n);
			value_2[1] = 1.0;

			if (CPXaddlazyconstraints(env, lp, 1, 2, &rhs, &sense, &izero, index_2, value_2, cname))
				print_error(" wrong CPXaddlazyconstraints [u-consistency]");
		}
	}


	if (inst->cplex_write_model)
		CPXwriteprob(env, lp, inst->cplex_model_file, NULL);

	free(cname[0]);
	free(cname);
}

// deprecabile
void build_sol(const double* xstar, solution * sol)
{
	switch (sol->inst->model_type)
	{
	case 0:
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
		identify_connected_components(xstar, sol, sol->inst->nnodes);
		break;
	case 1:
	case 3:
		build_sol_1(xstar, sol);
		break;
	case 2:
		build_sol_2(xstar, sol);
		break;
	case 9:
	case 10:
	case 11:
	case 12:
	case 13:
	case 14:
	case 15:
		// Nothing to do here
		break;
	default:
		print_error("unknown model type!");
		break;
	}
}

// deprecabile
void build_sol_1(const double* xstar, solution * sol)
{
	int n = sol->inst->nnodes;
	int ustart = sol->inst->ustart;
	if (sol->nodes_order == NULL)
		sol->nodes_order = (node_pair*)calloc(n, sizeof(node_pair));
	node_pair * pairs = sol->nodes_order;
	for (int i = 1; i < n; i++)
	{
		pairs[i].index = INT_ROUND(xstar[upos_mtz(i, ustart)]);
		pairs[i].node = i;
	}
	qsort(&pairs[1], n - 1, sizeof(node_pair), comparator_nodes);
}

// deprecabile
void build_sol_2(const double* xstar, solution * sol)
{
	int n = sol->inst->nnodes;
	int ystart = sol->inst->ystart;
	sol->nodes_order = (node_pair*)calloc(n * n, sizeof(node_pair));
	node_pair* pairs = sol->nodes_order;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			pairs[i * n + j].index = -INT_ROUND(xstar[ypos_gg(i, j, n, ystart)]);
			pairs[i * n + j].node = j;
		}
	}
	qsort(pairs, n * n, sizeof(node_pair), comparator_nodes);
}

void build_tour(solution * sol)
{
	int* tour = sol->tour;
	int* succ = sol->succ;
	int nnodes = sol->inst->nnodes;
	if (tour == NULL)
	{
		tour = (int*)calloc(nnodes, sizeof(int));
		sol->tour = tour;
	}
	tour[0] = 0;
	for (int i = 0; i < nnodes; i++)
		tour[i + 1] = succ[tour[i]];
}

void close_gp_pipe(solution * sol)
{
	if (sol->gp_pipe)
		_pclose(sol->gp_pipe);
	sol->gp_pipe = NULL;
}

void close_instance(instance * inst)
{
	if (inst)
	{
		if (inst->xcoord)
			free(inst->xcoord);
		if (inst->ycoord)
			free(inst->ycoord);
		if (inst->cplex_model_file)
			free(inst->cplex_model_file);
	}
	memset(inst, 0x00, sizeof(instance));
}

void close_solution(solution * sol)
{
	if (sol->comp)
		free(sol->comp);
	if (sol->succ)
		free(sol->succ);
	if (sol->nodes_order)
		free(sol->nodes_order);
	if (sol->tour)
		free(sol->tour);
	if (sol->xstar)
		free(sol->xstar);
	if (sol->inst)
	{
		close_instance(sol->inst);
		free(sol->inst);
	}
	if (sol->gp_pipe)
		_pclose(sol->gp_pipe);
	memset(sol, 0x00, sizeof(solution));
	close_gp_pipe(sol);
}

int comparator_nodes(const void* n1, const void* n2)
{
	return (*(node_pair*)n1).index - (*(node_pair*)n2).index;
}

int comparator_dist_pair(const void* a, const void* b)
{
	if ((*(dist_pair_t*)a).dist < (*(dist_pair_t*)b).dist)
		return -1;
	if ((*(dist_pair_t*)a).dist > (*(dist_pair_t*)b).dist)
		return 1;
	return 0;
}

double dist(int i, int j, const instance * inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if (!inst->integer_costs)
		return sqrt(dx * dx + dy * dy);
	return (double)INT_ROUND(sqrt(dx * dx + dy * dy));
}

int find_in_distribution(double* v, int size, double x)
{
	int a = 0, b = size, middle;
	while (1)
	{
		middle = (a + b) / 2;
		if (x >= v[middle] && x < v[middle + 1])
			return middle;
		if (x < v[middle])
			b = middle;
		else
			a = middle;
		if (a == b + 1)
			return -1;
	}
}

void first_nearest_neighbor(solution * sol)
{
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	long long n_var = (nnodes * (nnodes - 1LL)) / 2LL;

	// 'available' lists all node that have not been inserted yet
	int* available = (int*)calloc(nnodes, sizeof(int));
	for (int i = 0; i < nnodes; i++)
		available[i] = i;

	// tour[0] --> tour[1] --> ... --> tour[nnodes] will be the final tour
	int* tour = (int*)calloc(nnodes + 1, sizeof(int));

	// Storing the tour in the comp/succ format as well
	int* comp = (int*)malloc(nnodes * sizeof(int));
	int* succ = (int*)malloc(nnodes * sizeof(int));
	sol->ncomp = 1;

	// Selecting randomly the first node of the tour
	int selected;
	selected = (int)(((double)rand() / (double)RAND_MAX) * (double)(nnodes - 1));
	tour[nnodes] = tour[0] = available[selected];
	swap(&available[selected], &available[nnodes - 1]);
	comp[0] = 1;

	double min_dist, current_dist, total_cost = 0.0;

	// Creating the tour.
	// From now on, 'nodes - i' can be seen as the number of remaining nodes to selected.
	for (int i = 1, current, nearest; i < nnodes; i++)
	{
		current = tour[i - 1];
		min_dist = DBL_MAX;
		// Finding the nearest point to 'current'
		for (int j = 0; j < nnodes - i; j++)
		{
			// squared distance inlined for performace: d(current,j)
			current_dist = INLINE_SQUARED_DISTANCE(xcoord[current], xcoord[available[j]], ycoord[current], ycoord[available[j]]);
			if (current_dist < min_dist)
			{
				min_dist = current_dist;
				nearest = j;
			}
		}

		tour[i] = available[nearest];
		comp[i] = 1;
		succ[current] = available[nearest];
		swap(&available[nearest], &available[nnodes - i - 1]);
		total_cost += sqrt(min_dist);
	}
	succ[tour[nnodes - 1]] = tour[0];
	total_cost += dist(tour[nnodes - 1], tour[0], sol->inst);

	sol->tour = tour;
	sol->comp = comp;
	sol->succ = succ;
	sol->latest_cost = total_cost;

	free(available);

}

void first_nearest_neighbor_wmarkers(solution * sol)
{
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	long long n_var = (nnodes * (nnodes - 1LL)) / 2LL;

	// 'used' keeps track of nodes that have been already selected
	int* used = (int*)calloc(nnodes, sizeof(int));

	// tour[0] --> tour[1] --> ... --> tour[nnodes] will be the final tour
	int* tour = (int*)calloc(nnodes + 1, sizeof(int));

	// Storing the tour in the comp/succ format as well
	int* comp = (int*)malloc(nnodes * sizeof(int));
	int* succ = (int*)malloc(nnodes * sizeof(int));
	sol->ncomp = 1;

	srand(sol->inst->random_seed);

	// Selecting randomly the first node of the tour
	tour[0] = (int)(((double)rand() / (double)RAND_MAX) * (double)(nnodes - 1));
	used[tour[0]] = 1;
	tour[nnodes] = tour[0];
	comp[0] = 1;

	double min_dist, current_dist, total_cost = 0.0;

	// Creating the tour.
	// From now on, 'nodes - i' can be seen as the number of remaining nodes to selected.
	for (int i = 1, current, nearest; i < nnodes; i++)
	{
		current = tour[i - 1];
		min_dist = DBL_MAX;
		// Finding the nearest point to 'current'
		for (int j = 0; j < nnodes; j++)
		{
			if (!used[j])
			{
				// squared distance inlined for performace: d(current,j)
				current_dist = INLINE_SQUARED_DISTANCE(xcoord[current], xcoord[j], ycoord[current], ycoord[j]);
				if (current_dist < min_dist)
				{
					min_dist = current_dist;
					nearest = j;
				}
			}
		}

		tour[i] = nearest;
		comp[i] = 1;
		used[nearest] = 1;
		succ[current] = nearest;
		total_cost += sqrt(min_dist);
	}
	succ[tour[nnodes - 1]] = tour[0];
	total_cost += dist(tour[nnodes - 1], tour[0], sol->inst);

	sol->tour = tour;
	sol->comp = comp;
	sol->succ = succ;
	sol->latest_cost = total_cost;

	free(used);
}

double get_cost(solution * sol)
{
	double cost = 0.0;
	int current = 0;
	do
	{
		cost += dist(current, sol->succ[current], sol->inst);
	} while (current = sol->succ[current]);
	return cost;
}

double get_cost_tour(solution * sol)
{
	double cost = 0.0;
	int nnodes = sol->inst->nnodes;
	int* tour = sol->tour;
	for (int i = 0; i < nnodes; i++)
	{
		cost += dist(tour[i], tour[i + 1], sol->inst);
	}
	return cost;
}

char* get_model_name(int m)
{
	switch (m)
	{
	case 0:
		return "TOURS";
	case 1:
		return "MTZ";
	case 2:
		return "GG";
	case 3:
		return "MTZ-lazyc";
	case 4:
		return "LOOP";
	case 5:
		return "LOOP_gap";
	case 6:
		return "loop-callback";
	case 7:
		return "hard-fixing";
	case 8:
		return "local-branching";
	case 9:
		return "1-nn";
	case 10:
		return "k-nn";
	case 11:
		return "GRASP";
	case 12:
		return "insertion";
	case 13:
		return "simulated-annealing";
	case 14:
		return "VNS";
	case 15:
		return "greedy-opt-search";
	default:
		return "unknown";
	}
}

void get_sol(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	switch (sol->inst->model_type)
	{
	case 0:
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
	case 7:
	{
		int ncols = CPXgetnumcols(env, lp);
		CPXgetobjval(env, lp, &(sol->latest_cost));
		CPXgetbestobjval(env, lp, &(sol->best_lb));
		sol->xstar = (double*)calloc(ncols, sizeof(double));
		CPXgetx(env, lp, sol->xstar, 0, ncols - 1);
		break;
	}
	case 8:
	case 9:
	case 10:
	case 11:
	case 12:
	case 13:
	case 14:
	case 15:
		// Nothing to do here
		break;
	default:
		print_error("unknown model type!");
		break;
	}

}

void get_xstar_format_sol(solution * sol)
{
	int nnodes = sol->inst->nnodes;
	double* xstar = (double*)calloc((nnodes - 1) * nnodes / 2, sizeof(double));
	int* succ = sol->succ;

	int current = 0;
	for (int i = 0; i < nnodes; i++)
	{
		xstar[xpos(current, succ[current], nnodes)] = 1.0;
		current = succ[current];
	}

	sol->xstar = xstar;
}

void GRASP(solution * sol)
{
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	long long n_var = (nnodes * (nnodes - 1LL)) / 2LL;

	// 'used' keeps track of nodes that have been already selected
	int* used = (int*)calloc(nnodes, sizeof(int));

	// tour[0] --> tour[1] --> ... --> tour[nnodes] will be the final tour
	int* tour = (int*)calloc(nnodes + 1, sizeof(int));

	// Storing the tour in the comp/succ format as well
	int* comp = (int*)malloc(nnodes * sizeof(int));
	int* succ = (int*)malloc(nnodes * sizeof(int));
	sol->ncomp = 1;

	// Creating a cumulative distribution function to randomize which point to choose
	double cumulative[3];
	cumulative[0] = GRASP_DISTRIBUTION_VAL1;
	cumulative[1] = cumulative[0] + GRASP_DISTRIBUTION_VAL2;
	cumulative[2] = cumulative[1] + GRASP_DISTRIBUTION_VAL3;

	// Selecting randomly the first node of the tour
	tour[0] = (int)(((double)rand() / (double)RAND_MAX) * (double)(nnodes - 1));
	used[tour[0]] = 1;
	tour[nnodes] = tour[0];
	comp[0] = 1;

	double current_dist, total_cost = 0.0, uniform;

	// First, second and third smallest distances
	double min_dist_1, min_dist_2, min_dist_3;
	// First, second and third nearest points
	int nearest_1, nearest_2, nearest_3;
	nearest_1 = nearest_2 = nearest_3 = -1; // representing no valid point
	int current, selected, range;
	range = 3; // Means that initially we're gonna consider three neighbors

	// Creating the tour.
	// From now on, 'nodes - i' can be seen as the number of remaining nodes to selected.
	for (int i = 1; i < nnodes; i++)
	{
		current = tour[i - 1];
		min_dist_1 = min_dist_2 = min_dist_3 = DBL_MAX;
		// Finding the three nearest point to 'current'
		for (int j = 0; j < nnodes; j++)
		{
			if (!used[j])
			{
				// squared distance inlined for performace
				current_dist = (xcoord[current] - xcoord[j]) * (xcoord[current] - xcoord[j]) +
					(ycoord[current] - ycoord[j]) * (ycoord[current] - ycoord[j]);
				if (current_dist < min_dist_1)			// found the first nearest candidate
				{
					// cascade updating
					min_dist_3 = min_dist_2;
					min_dist_2 = min_dist_1;
					min_dist_1 = current_dist;
					nearest_3 = nearest_2;
					nearest_2 = nearest_1;
					nearest_1 = j;
				}
				else if (current_dist < min_dist_2)		// found the second nearest candidate
				{
					// cascade updating
					min_dist_3 = min_dist_2;
					min_dist_2 = current_dist;
					nearest_3 = nearest_2;
					nearest_2 = j;
				}
				else if (current_dist < min_dist_3)		// found the third nearest candidate
				{
					min_dist_3 = current_dist;
					nearest_3 = j;
				}
			}
		}

		// Updating range: cannot consider more neighbor than how many we're left with
		range = MIN(3, (nnodes - i));

		// Selecting randomly the next point in the tour according to 'cumulative'
		uniform = ((double)rand() / (double)RAND_MAX) * cumulative[range - 1];
		if (uniform < cumulative[0])
		{
			selected = nearest_1;
			current_dist = min_dist_1;
		}
		else if (uniform < cumulative[1])
		{
			selected = nearest_2;
			current_dist = min_dist_2;
		}
		else
		{
			selected = nearest_3;
			current_dist = min_dist_3;
		}

		tour[i] = selected;
		comp[i] = 1;
		used[selected] = 1;
		succ[current] = selected;
		total_cost += sqrt(current_dist);
	}
	succ[tour[nnodes - 1]] = tour[0];
	total_cost += dist(tour[nnodes - 1], tour[0], sol->inst);

	sol->tour = tour;
	sol->comp = comp;
	sol->succ = succ;
	sol->latest_cost = total_cost;

	free(used);
}

void gnuplot_plt(solution * sol)
{
	setup_gp_pipe(sol);

	fprintf(sol->gp_pipe, "set style line 1 lc rgb '%s' lt 1 lw 2 pt 7 pi -0.5 ps 0.5\n", GNUPLOT_LINECOLOR);
	fprintf(sol->gp_pipe, "set pointintervalbox 3\n");

	char* plot_command = (char*)malloc(COMMAND_LINE_DEFAULT_SIZE + 1);
	switch (sol->inst->model_type)
	{
	case 0:
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
	case 9:
	case 10:
	case 11:
	case 12:
	case 13:
	case 14:
		gnuplot_plt_0(sol, plot_command, "with linespoints ls 1 notitle");
		break;
	case 1:
	case 3:
		gnuplot_plt_1(sol, plot_command, "with linespoints ls 1 notitle");
		break;
	case 2:
		gnuplot_plt_2(sol, plot_command, "with linespoints ls 1 notitle");
		break;
	default:
		break;
	}
	free(plot_command);
}

void gnuplot_plt_0(solution * sol, char* plot_command, char* style)
{
	FILE* gp_pipe = sol->gp_pipe;
	int* comp_sink = (int*)calloc(sol->ncomp + 1, sizeof(int));
	comp_sink[sol->comp[0]] = 1;
	int comp_found = 1;
	// Seeking a 'sink' node for each connected component
	printf("Building plot ...\n");
	fflush(stdout);
	for (int i = 1; i < sol->inst->nnodes; i++)
	{
		if (!comp_sink[sol->comp[i]])
		{
			comp_sink[sol->comp[i]] = i + 1;
			if (++comp_found == sol->ncomp)
				break;
		}
	}

	// Gathering and accumulating data points for each subtour
	for (int tour = 1; tour <= sol->ncomp; tour++)
	{
		int sink, j;
		sink = j = comp_sink[tour] - 1;
		fprintf(gp_pipe, "$tour%i << END\n", tour);
		do
		{
			fprintf(gp_pipe, "%15.3lf %15.3lf\n", sol->inst->xcoord[j], sol->inst->ycoord[j]);
			j = sol->succ[j];
		} while (j != sink);
		fprintf(gp_pipe, "%15.3lf %15.3lf\nEND\n", sol->inst->xcoord[sink], sol->inst->ycoord[sink]);
	}
	// Building the whole gnuplot plot command
	fprintf(gp_pipe, "plot\\\n");
	for (int tour = 1; tour <= sol->ncomp - 1; tour++)
		fprintf(gp_pipe, " '$tour%i' %s,\\\n", tour, style);
	fprintf(gp_pipe, " '$tour%i' %s\n", sol->ncomp, style);
	free(comp_sink);
}

void gnuplot_plt_1(solution * sol, char* plot_command, char* style)
{
	char* current = plot_command;
	sprintf(current, "$tsptour << END\n");
	current += strlen(current);
	int n = sol->inst->nnodes;
	node_pair* sorted_nodes = sol->nodes_order;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	for (int i = 0; i < n; i++)
	{
		sprintf(current, "%15.3lf %15.3lf\n",
			xcoord[sorted_nodes[i].node],
			ycoord[sorted_nodes[i].node]);
		current += strlen(current);
	}
	sprintf(current, "%15.3lf %15.3lf\nEND\n",
		xcoord[sorted_nodes[0].node],
		ycoord[sorted_nodes[0].node]);
	current += strlen(current);
	sprintf(current, "plot '$tsptour' %s\n", style);
	fprintf(sol->gp_pipe, plot_command);
}

void gnuplot_plt_2(solution * sol, char* plot_command, char* style)
{
	char* current = plot_command;
	sprintf(current, "$tsptour << END\n");
	current += strlen(current);
	int n = sol->inst->nnodes;
	node_pair* sorted_archs = sol->nodes_order;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	sprintf(current, "%15.3lf %15.3lf\n", xcoord[0], ycoord[0]);
	current += strlen(current);
	int i = 0;
	do
	{
		sprintf(current, "%15.3lf %15.3lf\n",
			xcoord[sorted_archs[i].node],
			ycoord[sorted_archs[i].node]);
		current += strlen(current);
	} while (sorted_archs[++i].index < 0);
	sprintf(current, "%15.3lf %15.3lf\nEND\n", xcoord[0], ycoord[0]);
	current += strlen(current);
	sprintf(current, "plot '$tsptour' %s\n", style);
	fprintf(sol->gp_pipe, plot_command);
}

void hard_fixing(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	// Variables for counting time
	double iteration_time;
	double start = clock();

	// Getting first heuristic solution
	printf("Getting an intial feasible solution ...\n");
	first_nearest_neighbor(sol);
	greedy_refinement_2opt(sol);
	get_xstar_format_sol(sol);

	// Imposing time limit
	double const num_iterations = 10;
	double const initialization_time = (clock() - start) / CLOCKS_PER_SEC;
	double residual_time = sol->inst->time_limit - initialization_time;
	printf("Initialization time: %.1f s\n", initialization_time);
	double iter_time_limit = (residual_time - initialization_time) / num_iterations;
	CPXsetdblparam(env, CPXPARAM_TimeLimit, iter_time_limit);
	printf("Total time limit %.1f s \n", residual_time);

	int num_fixed_arches, ncols, nnodes;
	ncols = CPXgetnumcols(env, lp);
	nnodes = sol->inst->nnodes;
	double initial_cost, objval, latest_objval, bestbound, fixing_probability;
	initial_cost = objval = bestbound = sol->latest_cost;
	fixing_probability = 0.9;

	printf("Cost of heuristic incumbent solution: %lf\n", sol->latest_cost);

	for (int i = 0; i < num_iterations; i++)
	{
		iteration_time = -clock();

		// Hard Fixing
		num_fixed_arches = hf_set_bounds(sol, env, lp, fixing_probability);

		//Printing initial info about current iteration
		printf("************************************************************\n");
		printf("Iteration %d (limit %.1f s)\n", i + 1, iter_time_limit);
		printf("Hard fixing percentage: %.2f%%,\t actually fixed %.2f%%\n",
			fixing_probability * 100,
			((double)num_fixed_arches / (double)nnodes) * 100);
		printf("Solving problem... \n");

		// Solving the problem
		loop_callback(sol, env, lp);

		// Getting info about solution
		latest_objval = objval;
		CPXgetobjval(env, lp, &objval);
		CPXgetbestobjval(env, lp, &bestbound);

		// Taking time elapsed in current iteration 
		iteration_time += clock();
		iteration_time = iteration_time / CLOCKS_PER_SEC;

		//Printing info about solution of current iteration
		printf("Current solution cost: %.2lf (improved by %.2lf%%)\n",
			objval, (latest_objval - objval) / latest_objval * 100.0);
		printf("Elapsed time: %.1f s \n", iteration_time);

		CPXgetx(env, lp, sol->xstar, 0, ncols - 1);

		// Conditionally updating fixing probability
		if (iteration_time <= iter_time_limit)
			fixing_probability -= 0.1;
		else
			fixing_probability += 0.05;

		residual_time -= iteration_time;
		iter_time_limit = residual_time / (num_iterations - i - 1);
		CPXsetdblparam(env, CPXPARAM_TimeLimit, iter_time_limit);
	}
	sol->latest_improvement = (initial_cost - objval) / initial_cost * 100.0;
	printf("\n Hardfixing improved cost by %.2lf%%\n", sol->latest_improvement);
	sol->latest_cost = objval;
}

int hf_set_bounds(solution * sol, CPXENVptr env, CPXLPptr lp, double prob_of_fixing)
{
	long const nnodes = sol->inst->nnodes;
	char which_bound = 'L';
	double zero = 0.0;
	double one = 1.0;

	// Get solution of previous iteration
	double* xstar = sol->xstar;

	double random_number;
	int num_fixed_arches = 0;
	int index;
	for (int i = 0; i < nnodes; i++)
		for (int j = i + 1; j < nnodes; j++)
		{
			index = xpos(i, j, nnodes);

			// For all variables we set lb = 0
			if (CPXchgbds(env, lp, 1, &index, &which_bound, &zero))
				print_error("Error in hf_set_bounds() !");

			// For some variables we set lb = 1 
			random_number = (double)rand() / (double)RAND_MAX;
			if (random_number < prob_of_fixing && xstar[index] > 0.5)
			{
				if (CPXchgbds(env, lp, 1, &index, &which_bound, &one))
					print_error("Error in hf_set_bounds() !");
				num_fixed_arches++;
			}
		}

	return num_fixed_arches;
}

void local_branching(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	// Variables for counting time
	double iteration_time;
	double start = clock();

	// Getting first heuristic solution
	printf("Getting an intial feasible solution ...\n");
	first_nearest_neighbor(sol);
	greedy_refinement_2opt(sol);
	get_xstar_format_sol(sol);

	// Imposing time limit
	double const initialization_time = (clock() - start) / CLOCKS_PER_SEC;
	double residual_time = sol->inst->time_limit - initialization_time;
	printf("Initialization time: %.1f s\n", initialization_time);
	residual_time -= initialization_time;
	CPXsetdblparam(env, CPXPARAM_TimeLimit, residual_time);
	printf("Total time limit %.1f s \n", residual_time);

	double initial_cost, objval, latest_objval, bestbound;
	initial_cost = objval = bestbound = sol->latest_cost;

	printf("Cost of heuristic incumbent solution: %lf\n", sol->latest_cost);

	int k = 3, iteration = 1;
	while (residual_time > 0)
	{
		iteration_time = -clock();

		// Local Branching
		lb_set_constraint(sol, env, lp, k, iteration != 1);

		//Printing initial info about current iteration
		printf("************************************************************\n");
		printf("Iteration %d (residual %.1f s)\n", iteration, residual_time);
		printf("Number of free variables: %i\n", k);
		printf("Solving problem... \n");

		// Solving the problem
		loop_callback(sol, env, lp);

		// Getting info about solution
		latest_objval = objval;
		CPXgetobjval(env, lp, &objval);
		CPXgetbestobjval(env, lp, &bestbound);

		// Taking time elapsed in current iteration 
		iteration_time += clock();
		iteration_time = iteration_time / CLOCKS_PER_SEC;

		//Printing info about solution of current iteration
		printf("Current solution cost: %.2lf (improved by %.2lf%%)\n",
			objval, (latest_objval - objval) / latest_objval * 100.0);
		printf("Elapsed time: %.1f s \n", iteration_time);

		residual_time -= iteration_time;
		CPXsetdblparam(env, CPXPARAM_TimeLimit, residual_time);

		iteration++; k++;
	}
	double improvement = (initial_cost - objval) / initial_cost * 100.0;
	printf("\n Local branching improved cost by %.2lf%%\n", improvement);
	sol->latest_improvement = improvement;
	sol->latest_cost = objval;
}

int lb_set_constraint(solution * sol, CPXENVptr env, CPXLPptr lp, int k, int replace_constraint)
{
	// Get solution of previous iteration
	int ncols = CPXgetnumcols(env, lp);
	int lastrow = CPXgetnumrows(env, lp) - 1;

	// Deleting constraint related to previous iteration	
	if (replace_constraint)
		CPXdelrows(env, lp, lastrow, lastrow);

	CPXgetx(env, lp, sol->xstar, 0, ncols - 1);
	double* xstar = sol->xstar;

	double zero = 0.0;
	double one = 1.0;
	char name[STRING_MAX_SIZE + 1];
	char* cname[] = { name };

	int nnodes = sol->inst->nnodes;

	double rhs = (double)nnodes - k;
	char sense = 'E';
	sprintf(name, "local_branching");
	CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname);
	lastrow = CPXgetnumrows(env, lp) - 1;

	for (int i = 0; i < nnodes; i++)
		for (int j = i + 1; j < nnodes; j++)
			if (xstar[xpos(i, j, nnodes)] > 0.5)
				CPXchgcoef(env, lp, lastrow, xpos(i, j, nnodes), one);

	return lastrow;
}

void insertion_method(solution * sol)
{
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	long long n_var = (nnodes * (nnodes - 1LL)) / 2LL;

	// 'available' lists all node that have not been inserted yet
	int* available = (int*)calloc(nnodes, sizeof(int));
	for (int i = 0; i < nnodes; i++)
		available[i] = i;

	// Storing the tour in the comp/succ format as well
	int* comp = (int*)malloc(nnodes * sizeof(int));
	int* succ = (int*)malloc(nnodes * sizeof(int));
	sol->ncomp = 1;

	int first, selected;
	// selecting first random point
	selected = (int)(((double)rand() / (double)RAND_MAX) * (double)(nnodes - 1));
	first = available[selected];
	swap(&available[selected], &available[nnodes - 1]);
	// selecting second random point
	selected = (int)(((double)rand() / (double)RAND_MAX) * (double)(nnodes - 2));
	succ[first] = available[selected];
	succ[available[selected]] = first;
	swap(&available[selected], &available[nnodes - 2]);
	comp[0] = comp[1] = 1;

	double current_extra_mileage, min_extra_mileage, total_cost;
	int h, current, fusion_node;
	total_cost = 2.0 * dist(first, succ[first], sol->inst);

	for (int i = 2; i < nnodes; i++)
	{
		// Selecting the next point to be inserted
		selected = (int)(((double)rand() / (double)RAND_MAX) * (double)(nnodes - i - 1));
		h = available[selected];
		swap(&available[selected], &available[nnodes - i - 1]);

		// finding the edge with the minimum extra mileage
		min_extra_mileage = DBL_MAX;
		current = first;
		do
		{
			// inlining the following:
			// current_extra_mileage = d(current, h) + d(h, succ[current]) - d(current, succ[current])
			current_extra_mileage =
				sqrt(INLINE_SQUARED_DISTANCE(xcoord[current], xcoord[h],
					ycoord[current], ycoord[h]))
				+ sqrt(INLINE_SQUARED_DISTANCE(xcoord[h], xcoord[succ[current]],
					ycoord[h], ycoord[succ[current]]))
				- sqrt(INLINE_SQUARED_DISTANCE(xcoord[current], xcoord[succ[current]],
					ycoord[current], ycoord[succ[current]]));
			if (current_extra_mileage < min_extra_mileage)
			{
				min_extra_mileage = current_extra_mileage;
				fusion_node = current;
			}
		} while ((current = succ[current]) != first);

		// finally performing the insertion
		succ[h] = succ[fusion_node];
		succ[fusion_node] = h;
		comp[i] = 1;

		total_cost += min_extra_mileage;
	}

	free(available);
	sol->succ = succ;
	sol->comp = comp;
	sol->latest_cost = total_cost;
}

void k_nearest_neighbors(solution * sol)
{
	int neighbors = sol->inst->neighbors;
	int nnodes = sol->inst->nnodes;
	int n_var = (nnodes * (nnodes - 1)) / 2;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;

	// 'used' keeps track of nodes that have been already selected
	int* used = (int*)calloc(nnodes, sizeof(int));

	// 'cumulative' represents a discrete cumulative distribution function
	double* cumulative = (double*)calloc(nnodes + 1, sizeof(double));

	// tour[0] --> tour[1] --> ... --> tour[nnodes] will be the final tour
	int* tour = (int*)calloc(nnodes + 1, sizeof(int));

	// Storing the tour in the comp/succ format as well
	int* comp = (int*)malloc(nnodes * sizeof(int));
	int* succ = (int*)malloc(nnodes * sizeof(int));
	sol->ncomp = 1;

	// 'pairs' stores pairs (i,d) meaning that the distance from a some node to 'i' is 'd' 
	dist_pair_t * pairs = (dist_pair_t*)calloc(nnodes - 1, sizeof(dist_pair_t));

	double lambda = sol->inst->lambda;	// parameter of the exponential distribution
	double average_selection = 0.0;		// what is the index of the random choise on average
	double total_cost = 0.0;
	printf("  k-nn: using lambda = %.2f\n", lambda);

	// Selecting randomly the first node of the tour
	tour[0] = ((double)rand() / (double)RAND_MAX) * (double)(nnodes - 1);
	used[tour[0]] = 1;
	tour[nnodes] = tour[0];
	comp[0] = 1;

	// Initializing the cumulative distribution according to an
	// exponential distribution function of parameter 'lambda'.
	// The domain is restricted to [0, nnodes-1].
	cumulative[0] = 0.0;
	for (int i = 1; i < nnodes; i++)
		cumulative[i] = cumulative[i - 1] + lambda * exp(-lambda * (i - 1));

	// Creating the tour.
	// From now on, 'nodes - i' can be seen as the number of remaining nodes to selected.
	for (int i = 1, current, selected = 0; i < nnodes; i++)
	{
		current = tour[i - 1];
		// Storing all distances from 'current' node to the j-th
		for (int j = 0, k = 0; j < nnodes; j++)
		{
			if (!used[j])
			{
				pairs[k].dist = INLINE_SQUARED_DISTANCE(xcoord[current], xcoord[j], ycoord[current], ycoord[j]);
				pairs[k].node = j;
				k++;
			}
		}
		// Cannot consider more neighbor than how many we're left with
		neighbors = MIN(neighbors, (nnodes - i));

		// Sorting all pairs according to the value of distance
		qsort(pairs, nnodes - i, sizeof(dist_pair_t), comparator_dist_pair);

		// Finally selecting a random node
		selected = find_in_distribution(
			cumulative, neighbors,
			((double)rand() / ((double)RAND_MAX + 1.0) * cumulative[neighbors]));
		average_selection += selected;
		tour[i] = pairs[selected].node;
		used[tour[i]] = 1;
		succ[current] = tour[i];
		comp[i] = 1;
		total_cost += sqrt(pairs[selected].dist);
		if (!(i % 1000))
			printf("  k-nn: processed %i/%i nodes \n", i, nnodes);
	}
	succ[tour[nnodes - 1]] = tour[0];
	total_cost += dist(tour[nnodes - 1], tour[0], sol->inst);

	average_selection /= nnodes - 1;
	printf("  k-nn: average selection: (%.1f)-th nearest\n", average_selection + 1.0);

	sol->tour = tour;
	sol->comp = comp;
	sol->succ = succ;
	sol->latest_cost = total_cost;

	free(used);
	free(cumulative);
	free(pairs);
}

void identify_connected_components(const double* xstar, solution * sol, int nnodes)
{
	sol->ncomp = 0;
	if (sol->succ == NULL)
	{
		sol->succ = (int*)malloc(nnodes * sizeof(int));
		sol->comp = (int*)malloc(nnodes * sizeof(int));
	}
	for (int i = 0; i < nnodes; i++)
		sol->succ[i] = sol->comp[i] = -1;

	for (int start = 0; start < nnodes; start++)
	{
		if (sol->comp[start] >= 0) // node "start" was already visited, just skip it
			continue;

		// a new component is found
		sol->ncomp++;
		int i = start;
		int done = 0;
		while (!done)  // go and visit the current component
		{
			sol->comp[i] = sol->ncomp;
			done = 1;
			for (int j = 0; j < nnodes; j++)
			{
				if (i != j && xstar[xpos(i, j, nnodes)] > 0.5 && sol->comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before
				{
					sol->succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		sol->succ[i] = start;  // last arc to close the cycle

		// go to the next component...
	}
}

void loop_callback(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	int nthreads = 1;
	data_wrap dw = { .env = env,.lp = lp,.sol = sol };
	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	// let MIP callbacks work on the original model
	CPXsetlazyconstraintcallbackfunc(env, subtour_elimination_lazycb, &dw);
	CPXgetnumcores(env, &nthreads);
	if (!benchmark)
		printf("Using %i threads\n", nthreads);
	CPXsetintparam(env, CPX_PARAM_THREADS, nthreads);

	if (CPXmipopt(env, lp))
		print_error("CPXmipopt() error");
}

void loop(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	int ncols = CPXgetnumcols(env, lp);

	double* xstar = (double*)calloc(ncols, sizeof(double));

	int sec_steps = 0;
	double objval, bestbound;

	while (1)
	{
		CPXmipopt(env, lp);
		CPXgetx(env, lp, xstar, 0, ncols - 1);
		build_sol(xstar, sol);

		// Print main infos
		if (VERB_PRINT_SEC_ITERATIONS && !benchmark)
		{
			CPXgetobjval(env, lp, &objval);
			CPXgetbestobjval(env, lp, &bestbound);
			printf("  Iteration num: %d \n", sec_steps + 1);
			printf("  Number of connected components: %d \n", sol->ncomp);
			printf("  Value of objective function: %.6f \t", objval);
			printf("  Value of best bound: %.6f \n\n", bestbound);
		}

		if (sol->ncomp == 1)
		{
			break;
		}
		else
		{
			if (GNUPLOT_PROGRESS && !benchmark)
				gnuplot_plt(sol);
			add_sec(sol, env, lp);
			sec_steps++;
		}
	}
	sol->sec_steps = sec_steps;
	if (!benchmark)
		printf("Subtour constraints added in %i steps\n", sol->sec_steps);
	free(xstar);
}

void loop_gap(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	int ncols = CPXgetnumcols(env, lp);
	long long nnodes = sol->inst->nnodes;

	double* xstar = (double*)calloc(ncols, sizeof(double));
	int sec_steps = 0;
	double objval, bestbound;

	// Parameters for heuristic approach
	int use_gap_heuristic = 1;
	double mipgap_default = 0.0001; //0.01% tolerance
	double mipgap = ((double)nnodes) / (10000 * 100);

	if (use_gap_heuristic)
	{
		CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, mipgap);
		if (VERB_PRINT_SEC_ITERATIONS && !benchmark)
			printf("** Heuristic activated. Value of gap tolerance: %.4f %%**\n", mipgap * 100);
	}

	while (1)
	{
		CPXmipopt(env, lp);
		CPXgetx(env, lp, xstar, 0, ncols - 1);
		build_sol(xstar, sol);

		// Print main infos
		if (VERB_PRINT_SEC_ITERATIONS && !benchmark)
		{
			CPXgetobjval(env, lp, &objval);
			CPXgetbestobjval(env, lp, &bestbound);
			printf("  Iteration num: %d \n", sec_steps + 1);
			printf("  Number of connected components: %d \n", sol->ncomp);
			printf("  Value of objective function: %.6f \t", objval);
			printf("  Value of best bound: %.6f \n\n", bestbound);
		}

		if (sol->ncomp == 1)
		{
			/* If we enter here, its the first time we find a solution with 1
			connected component using the loose gap tolerance "mip_gap".
			From now, the next time we will find a solution with ncomp = 1
			it will be optimal because the tolerance was re-set to "mipgap_default".
			(So the loop will break entering the usual condition).
			*/
			if (use_gap_heuristic)
			{
				CPXsetdblparam(env, CPXPARAM_MIP_Tolerances_MIPGap, mipgap_default);
				use_gap_heuristic = 0;
				if (VERB_PRINT_SEC_ITERATIONS)
					printf("** Heuristic de-activated. Value of gap tolerance: %.4f  %%**\n", 100 * mipgap_default);
				continue;
			}
			else
				break;
		}
		else
		{
			if (GNUPLOT_PROGRESS && !benchmark)
				gnuplot_plt(sol);
			add_sec(sol, env, lp);
			sec_steps++;
		}
	}
	sol->sec_steps = sec_steps;
	if (!benchmark)
		printf("Subtour constraints added in %i steps\n", sol->sec_steps);
	free(xstar);
}

void print_cplex_solution(const double* xstar, solution * sol)
{
	int n = sol->inst->nnodes;
	switch (sol->inst->model_type)
	{
	case 0:
	case 4:
	case 6:
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				if (xstar[xpos(i, j, n)] > 0.5)
					printf(" x[%4i,%4i]\t= 1\n", i + 1, j + 1);
		break;
	case 1:
	case 3:
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (xstar[xpos_asymm(i, j, n)] > 0.5)
					printf(" x[%4i,%4i]\t= 1\n", i + 1, j + 1);
		for (int i = 0; i < n; i++)
			printf(" u[%4i]\t= %.4f\n", i + 1, xstar[upos_mtz(i, sol->inst->ustart)]);
		break;
	case 2:
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (xstar[xpos_asymm(i, j, n)] > 0.5)
					printf(" x[%4i,%4i]\t= 1\n", i + 1, j + 1);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (xstar[ypos_gg(i, j, n, sol->inst->ystart)] > 0)
					printf(" y[%4i,%4i]\t= %.4lf\n", i + 1, j + 1, xstar[ypos_gg(i, j, n, sol->inst->ystart)]);
		break;
	default:
		print_error("unknown model type!");
		break;
	}
}

void print_error(const char* err)
{
	printf("\n\n ERROR: %s \n\n", err);
	fflush(stdout);
	exit(1);
}

void provide_heuristic_mipstart(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	int nnodes = sol->inst->nnodes;
	int nvars = nnodes * (nnodes - 1) / 2; // it depends on the model type
	int ncols = CPXgetnumcols(env, lp);

	// Providing CPLEX a starting solution
	const int beg[2] = { 0, nnodes };
	int* varindeces = (int*)calloc(nnodes, sizeof(int));
	double* values = (double*)calloc(nnodes, sizeof(double));
	int effortlevel[1] = { CPX_MIPSTART_AUTO };
	for (int i = 0; i < nnodes; i++)
		values[i] = 1.0;
	int counter = 0;
	for (int i = 0; i < nvars; i++)
	{
		if (sol->xstar[i] == 1.0)
			varindeces[counter++] = i;
	}
	int error = CPXaddmipstarts(env, lp, 1, nnodes, beg, varindeces, values, effortlevel, NULL);
	printf("Adding MIP start: %s\n", (error) ? "FAILED" : "SUCCESSFUL");

	free(varindeces);
	free(values);
}

void random_solution(solution * sol)
{
	int nnodes = sol->inst->nnodes;
	int* available = (int*)calloc(nnodes, sizeof(int));
	int* tour = sol->tour;

	if (sol->succ == NULL)
		sol->succ = (int*)malloc(nnodes * sizeof(int));
	if (sol->comp == NULL)
		sol->comp = (int*)malloc(nnodes * sizeof(int));
	if (sol->tour == NULL)
		tour = (int*)calloc(nnodes + 1, sizeof(int));

	for (int i = 0, node; i < nnodes; i++)
		available[i] = i;

	// Creating a random tour in sequential format
	int selected;
	double cost = 0.0;
	selected = (int)((double)rand() / (double)RAND_MAX * (nnodes - 1));
	tour[0] = available[selected];
	swap(&available[selected], &available[nnodes - 1]);

	for (int i = 1, selected; i < nnodes; i++)
	{
		selected = (int)((double)rand() / (double)RAND_MAX * (nnodes - 1 - i));
		tour[i] = available[selected];
		swap(&available[selected], &available[nnodes - 1 - i]);
		cost += dist(tour[i - 1], tour[i], sol->inst);
	}
	tour[nnodes] = tour[0];
	cost += dist(tour[nnodes - 1], tour[nnodes], sol->inst);

	// Converting in succ/comp format
	sol->ncomp = 1;
	for (int i = 0; i < nnodes; i++)
	{
		sol->succ[tour[i]] = tour[i + 1];
		sol->comp[i] = 1;
	}
	sol->latest_cost = cost;
	sol->tour = tour;

	free(available);
}

void greedy_refinement_2opt(solution * sol)
{
	int nnodes = sol->inst->nnodes;
	int* tour = (int*)malloc((nnodes + 1) * sizeof(int));
	double initial_sol_cost = sol->latest_cost;
	double sol_cost = initial_sol_cost;
	long long residual_checks = (long long)sol->checks * (long long)1000000;

	// converting succ/comp format to sequential tour format
	tour[0] = 0;
	for (int i = 0; i < nnodes; i++)
		tour[i + 1] = sol->succ[tour[i]];

	int improvable, swapped = 0;
	double delta;
	do
	{
		improvable = 0;
		for (int i = 0; i < nnodes - 2; i++)
		{
			for (int j = i + 2; j < nnodes; j++)
			{
				if (i == 0 && j == nnodes - 1)
					continue;

				delta =
					-dist(tour[i], tour[i + 1], sol->inst)
					- dist(tour[j], tour[(j + 1) % nnodes], sol->inst)
					+ dist(tour[i], tour[j], sol->inst)
					+ dist(tour[i + 1], tour[(j + 1) % nnodes], sol->inst);
				residual_checks--;
				if (delta < -1e-6)
				{
					sol_cost += delta;
					for (int k = 0; k < (j - i) / 2; k++)
						swap(&tour[i + k + 1], &tour[j - k]);
					improvable = 1;
					printf("\r2-opt refining (greedy): %4i tuples swapped", ++swapped);
				}
				if (residual_checks != CHECKS_NO_LIMITS && residual_checks == 0)
					goto end;
			}
		}
	} while (improvable);
end:;
	printf("\n");
	// converting back to succ/comp format
	sol->ncomp = 1;
	for (int i = 0; i < nnodes; i++)
	{
		sol->succ[tour[i]] = tour[i + 1];
		sol->comp[i] = 1;
	}

	sol->latest_cost = sol_cost;
	free(tour);
}

void greedy_refinement_2opt_test(solution* sol)
{
	int nnodes = sol->inst->nnodes;
	int* tour = (int*)malloc((nnodes + 1) * sizeof(int));
	double initial_sol_cost = sol->latest_cost;
	double sol_cost = initial_sol_cost;
	long long residual_checks = (long long)sol->checks * (long long)1000000;

	// converting succ/comp format to sequential tour format
	tour[0] = 0;
	for (int i = 0; i < nnodes; i++)
		tour[i + 1] = sol->succ[tour[i]];

	int improvable, swapped = 0;
	double delta;
	do
	{
		improvable = 0;
		for (int i = 0; i < nnodes - 2; i++)
		{
			for (int j = i + 2; j < nnodes; j++)
			{
				if (i == 0 && j == nnodes - 1)
					continue;

				delta =
					-dist(tour[i], tour[i + 1], sol->inst)
					- dist(tour[j], tour[(j + 1) % nnodes], sol->inst)
					+ dist(tour[i], tour[j], sol->inst)
					+ dist(tour[i + 1], tour[(j + 1) % nnodes], sol->inst);
				residual_checks--;
				if (delta < -1e-6)
				{
					sol_cost += delta;
					for (int k = 0; k < (j - i) / 2; k++)
						swap(&tour[i + k + 1], &tour[j - k]);
					improvable = 1;
					printf("\r2-opt refining (greedy): %4i tuples swapped", ++swapped);
				}
				if (residual_checks != CHECKS_NO_LIMITS && residual_checks == 0)
					goto end;
			}
		}
	} while (improvable);
end:;
	printf("\n");
	// converting back to succ/comp format
	sol->ncomp = 1;
	for (int i = 0; i < nnodes; i++)
	{
		sol->succ[tour[i]] = tour[i + 1];
		sol->comp[i] = 1;
	}

	sol->latest_cost = sol_cost;
	free(tour);
}

void greedy_refinement_3opt(solution * sol)
{
	// Getting pointers from solution
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;

	int* succ = sol->succ;
	long long residual_checks = (long long)sol->checks * (long long)1000000;

	double delta;
	int improvable, swapped = 0, start, current;
	do
	{
		improvable = 0;
	loop_from_scratch:;
		int node_a = start = rand() % (nnodes - 1);
		//int node_a = start = 0; //D
		int node_b = succ[node_a];
		while (node_b != start && succ[node_b] != start)
		{
			int node_c = succ[node_b];
			int node_d = succ[node_c];
			while (node_d != start)
			{
				int node_e = succ[node_d];
				int node_f = succ[node_e];

				while (node_f != node_a && node_e != start)
				{
					// Computing distances 
					delta =
						-sqrt(INLINE_SQUARED_DISTANCE(xcoord[node_a], xcoord[node_b], ycoord[node_a], ycoord[node_b]))
						- sqrt(INLINE_SQUARED_DISTANCE(xcoord[node_c], xcoord[node_d], ycoord[node_c], ycoord[node_d]))
						- sqrt(INLINE_SQUARED_DISTANCE(xcoord[node_e], xcoord[node_f], ycoord[node_e], ycoord[node_f]))
						+ sqrt(INLINE_SQUARED_DISTANCE(xcoord[node_a], xcoord[node_d], ycoord[node_a], ycoord[node_d]))
						+ sqrt(INLINE_SQUARED_DISTANCE(xcoord[node_e], xcoord[node_b], ycoord[node_e], ycoord[node_b]))
						+ sqrt(INLINE_SQUARED_DISTANCE(xcoord[node_c], xcoord[node_f], ycoord[node_c], ycoord[node_f]));

					residual_checks--;
					if (delta < -1e-6)
					{
						// Crossing edges
						succ[node_a] = node_d;
						succ[node_e] = node_b;
						succ[node_c] = node_f;
						sol->latest_cost += delta;
						printf("\r3-opt refining (greedy): %4i tuples swapped", ++swapped);
						improvable = 1;
						if (residual_checks != CHECKS_NO_LIMITS && residual_checks == 0)
							goto end;
						goto loop_from_scratch;
					}
					if (residual_checks != -1 && residual_checks == 0)
						goto end;
					node_e = succ[node_e];
					node_f = succ[node_e];
				}
				node_c = succ[node_c];
				node_d = succ[node_c];
			}
			node_a = succ[node_a];
			node_b = succ[node_a];
		}
	} while (improvable);
end:;
	printf("\n");
}

void setup_gp_pipe(solution * sol)
{
	if (sol->gp_pipe != NULL)
		return;

	char gp_path[STRING_MAX_SIZE + 1];
	sprintf(gp_path, "%s -p", GNUPLOT_EXECUTABLE);

	sol->gp_pipe = _popen(gp_path, "w");
	fprintf(sol->gp_pipe, "set style line 1 lc rgb '%s' lt 1 lw 2 pt 7 pi -0.5 ps 0.5\n", GNUPLOT_LINECOLOR);
	fprintf(sol->gp_pipe, "set pointintervalbox 3\n");
}

void setup_solution(solution * sol)
{
	memset(sol, 0x00, sizeof(solution));
	sol->inst = (instance*)calloc(1, sizeof(instance));
}

void simulated_annealing(solution * sol)
{
	const int nnodes = sol->inst->nnodes;
	double delta, sol_cost;

	printf("Getting the initial solution ...\n");
	random_solution(sol);
	sol_cost = sol->latest_cost;
	printf("Inital solution cost:\t%lf\n", sol_cost);

	printf("Iteration estimation ...\n");
	const double initial_temperature = 1e+3;
	const double final_temperature = 1e-3;
	const double normalizing_factor = 1.0 / average_distance(sol->inst);
	double temperature = initial_temperature;
	const long long T =
		(long long)(obtain_SA_mean_iteration_per_sec(sol, initial_temperature, final_temperature, normalizing_factor)
			* sol->inst->time_limit * 1.00); //1.05
	//const long long T = (int)1e9;
	const long long print_step = T / 100LL;
	const double cooling = pow(final_temperature / initial_temperature, 1.0 / (double)T);
	printf("Estimated %.1lfM iterations in %.2lf s\n", (double)T / 1e6, sol->inst->time_limit);

	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;
	instance * inst = sol->inst;

	int* tour = sol->tour;
	double remaining_time = sol->inst->time_limit;
	double timer;

	printf("Intial temperature:\t%lf\n", temperature);
	printf("Cooling coefficient:\t%.15lf\n", cooling);
	//timer = -clock();
	int a, b;
	for (long long t = 0; t < T; t++)
	{
		// creating a new candidate solution from the neighborhood

		// Selecting randomly two indexes (a,b) both in [0; nnodes-1] such that a<b
		while ((a = (int)((double)rand() / (double)RAND_MAX * (nnodes - 1))) >=
			(b = (int)((double)rand() / (double)RAND_MAX * (nnodes - 1))) ||
			b - a >= nnodes - 2);

		// calculating the result of applying the 'inverse' operator on the tour cost
		delta =
			-dist(tour[(a - 1 + nnodes) % nnodes], tour[a], inst)
			- dist(tour[b], tour[(b + 1 + nnodes) % nnodes], inst)
			+ dist(tour[(a - 1 + nnodes) % nnodes], tour[b], inst)
			+ dist(tour[a], tour[(b + 1 + nnodes) % nnodes], inst);

		if (delta <= 0 ||
			(((double)rand() / (double)RAND_MAX) < exp(-delta * normalizing_factor / temperature)))
		{
			sol_cost += delta;

			// creating the tour i.e. applying the 'inverse' operator to the tour
			for (int i = 0; i < (b - a + 1) / 2; i++)
				swap(&tour[a + i], &tour[b - i]);
			tour[nnodes] = tour[0];
		}

		temperature *= cooling;
		if (!(t % print_step))
		{
			printf("\rAnnealing ... %3.0f%%", (double)t / (double)T * 100.0);
			//remaining_time -= (timer + clock()) / CLOCKS_PER_SEC;
			//if (remaining_time < 0)
			//	break;
			//timer = -clock();
		}

	}
	printf("\rAnnealing ... 100%%\n");
	printf("Final temperature:\t%.15lf\n", temperature);

	// converting back to succ/comp format
	for (int i = 0; i < nnodes; i++)
		sol->succ[tour[i]] = tour[i + 1];
	printf("Solution improved by %.3f%%\n", (1.0 - sol_cost / sol->latest_cost) * 100.0);
	sol->latest_cost = sol_cost;
}

double average_distance(instance * inst)
{
	int nnodes = inst->nnodes;
	double* xcoord = inst->xcoord;
	double* ycoord = inst->ycoord;

	double distances = 0.0;
	int selected_1, selected_2;
	for (int i = 0; i < nnodes; i++)
	{
		selected_1 = rand() % (nnodes - 1);
		selected_2 = rand() % (nnodes - 1);
		distances += sqrt(INLINE_SQUARED_DISTANCE(xcoord[selected_1], xcoord[selected_2],
			ycoord[selected_1], ycoord[selected_2]));
	}

	return distances / (double)nnodes;
}

void solve_lp(solution * sol, CPXENVptr env, CPXLPptr lp)
{
	switch (sol->inst->model_type)
	{
	case 0:
	case 1:
	case 2:
	case 3:
		if (CPXmipopt(env, lp))
			print_error("CPXmipopt() error");
		break;
	case 4:
		loop(sol, env, lp);
		break;
	case 5:
		loop_gap(sol, env, lp);
		break;
	case 6:
		loop_callback(sol, env, lp);
		break;
	case 7:
		hard_fixing(sol, env, lp);
		break;
	case 8:
		local_branching(sol, env, lp);
		break;
	case 9:
		sol->heuristic_method = first_nearest_neighbor;
		multi_start_parallel(sol);
		refine(sol);
		break;
	case 10:
		sol->heuristic_method = k_nearest_neighbors;
		multi_start_parallel(sol);
		break;
	case 11:
		sol->heuristic_method = GRASP;
		multi_start_parallel(sol);
		break;
	case 12:
		sol->heuristic_method = insertion_method;
		multi_start_parallel(sol);
		break;
	case 13:
		simulated_annealing(sol);
		break;
	case 14:
		VNS(sol);
		break;
	case 15:
		greedy_opt_search(sol);
		break;
	default:
		print_error("Unknown model type!");
		break;
	}
}

static int CPXPUBLIC subtour_elimination_lazycb(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	data_wrap* dw = (data_wrap*)cbhandle;
	int ncols = CPXgetnumcols(env, dw->lp);

	// get solution xstar
	double* xstar = (double*)malloc(ncols * sizeof(double));
	if (CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, ncols - 1))
		return 1;

	// get some random information at the node (as an example)
	solution local_solution;
	memset(&local_solution, 0x00, sizeof(local_solution));
	double objval = CPX_INFBOUND;
	int mythread = -1;
	double zbest;
	CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);

	//apply cut separator and possibly add violated cuts
	CPXgetx(env, dw->lp, xstar, 0, ncols - 1);
	identify_connected_components(xstar, &local_solution, dw->sol->inst->nnodes);
	free(xstar);

	if (local_solution.ncomp > 1)
	{
		add_sec_cb(&local_solution, dw->sol->inst->nnodes, env, cbdata, wherefrom);
		if (VERB_PRINT_COMP_CUT && !benchmark)
			printf(" callback: cut %i components\n", local_solution.ncomp);
		*useraction_p = CPX_CALLBACK_SET;
	}
	close_solution(&local_solution);
	return 0;
}

static inline void swap(int* a, int* b)
{
	if (*a == *b)
		return;
	int temp = *a;
	*a = *b;
	*b = temp;
}

void multi_start(solution * sol)
{
	solution* best_solution = (solution*)calloc(1, sizeof(solution));
	solution* current_sol = (solution*)calloc(1, sizeof(solution));
	current_sol->inst = sol->inst;
	printf("Total number of trials: %i\n", sol->trials);

	double min_cost = DBL_MAX;
	for (int i = 0; i < sol->trials; i++)
	{
		sol->heuristic_method(current_sol);
		if (current_sol->latest_cost < min_cost)
		{
			min_cost = current_sol->latest_cost;
			best_solution->inst = NULL;
			close_solution(best_solution);
			*best_solution = *current_sol;
			printf(" +++ Current best solution updated at trial %i: %.2f\n", i + 1, min_cost);
		}
		else
		{
			// Temporarily detaching the instance from 'current_solution' 
			// to safely close the solution without closing the instance.
			current_sol->inst = NULL;
			close_solution(current_sol);
			current_sol->inst = sol->inst;
		}
		if (!(i % 10))
			printf("Computed %i/%i trials\n", i + 1, sol->trials);
	}
	*sol = *best_solution;
	free(current_sol);
	free(best_solution);
	printf("Final best solution cost: %.2f\n", min_cost);
}

void multi_start_parallel(solution * sol)
{
	double min_cost = DBL_MAX;
	int trials = sol->trials;
	#pragma omp parallel
	{
		#pragma omp master
		{
			printf("Total number of trials: %i\n", sol->trials);
			printf("Using %i threads\n", omp_get_num_threads());
		}

		solution* best_solution = (solution*)calloc(1, sizeof(solution));
		solution* current_sol = (solution*)calloc(1, sizeof(solution));
		current_sol->inst = sol->inst;

		srand(sol->inst->random_seed + omp_get_thread_num());

		double local_min_cost = DBL_MAX;
		int i;
		#pragma omp for
		for (i = 0; i < trials; i++)
		{
			sol->heuristic_method(current_sol);
			if (current_sol->latest_cost < local_min_cost)
			{
				local_min_cost = current_sol->latest_cost;
				best_solution->inst = NULL;
				close_solution(best_solution);
				*best_solution = *current_sol;
				printf("T%02i: Found local candidate solution: %.2f\n", omp_get_thread_num(), local_min_cost);
			}
			else
			{
				// Temporarily detaching the instance from 'current_solution' 
				// to safely close the solution without closing the instance.
				current_sol->inst = NULL;
				close_solution(current_sol);
				current_sol->inst = sol->inst;
			}
		}

		// Each thread generated a best_solution so we have to select the best among them
		#pragma omp critical
		{
			if (local_min_cost < min_cost)
			{
				*sol = *best_solution;
				min_cost = local_min_cost;
			}
			else
			{
				best_solution->inst = NULL;
				close_solution(best_solution);
			}
		}
		free(current_sol);
		free(best_solution);
	}

	printf("Final best solution cost: %.2f\n", min_cost);
}

double obtain_SA_mean_iteration_per_sec(solution * sol, double initial_temperature, double final_temperature, double normalizing_factor)
{
	const int nnodes = sol->inst->nnodes;
	const int T = (int)1e6;
	double delta, sol_cost, temperature;
	const double cooling = pow(final_temperature / initial_temperature, 1.0 / (double)T);
	temperature = initial_temperature;
	int a, b;

	instance * inst = sol->inst;
	int* tour = (int*)malloc((nnodes + 1) * sizeof(int));
	memcpy(tour, sol->tour, (nnodes + 1) * sizeof(int));
	sol_cost = sol->latest_cost;

	double timer = -clock();
	// simulating simulated annealing
	for (long long t = 0; t < T; t++)
	{
		// Selecting randomly two indexes (a,b) both in [0; nnodes-1] such that a<b
		while ((a = (int)((double)rand() / (double)RAND_MAX * (nnodes - 1))) >=
			(b = (int)((double)rand() / (double)RAND_MAX * (nnodes - 1))) ||
			b - a >= nnodes - 2);

		// calculating the result of applying the 'inverse' operator on the tour cost
		delta =
			-dist(tour[(a - 1 + nnodes) % nnodes], tour[a], inst)
			- dist(tour[b], tour[(b + 1 + nnodes) % nnodes], inst)
			+ dist(tour[(a - 1 + nnodes) % nnodes], tour[b], inst)
			+ dist(tour[a], tour[(b + 1 + nnodes) % nnodes], inst);

		if (delta <= 0 ||
			(((double)rand() / (double)RAND_MAX) < exp(-delta * normalizing_factor / temperature)))
		{
			sol_cost += delta;

			// creating the tour i.e. applying the 'inverse' operator to the tour
			for (int i = 0; i < (b - a + 1) / 2; i++)
				swap(&tour[a + i], &tour[b - i]);
			tour[nnodes] = tour[0];
		}

		temperature *= cooling;
	}

	timer += clock();
	timer /= CLOCKS_PER_SEC;

	free(tour);
	// restoring random seed as all of this never happened ;)
	srand(sol->inst->random_seed);

	return (double)T / timer;
}

int TSPopt(solution * sol)
{
	// Setting up CPLEX environment
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	if (VERB_PRINT_CPLEX_INFOS && !benchmark)
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);

	// Setting user parameters for cplex environment
	CPXsetdblparam(env, CPXPARAM_TimeLimit, sol->inst->time_limit);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, sol->inst->random_seed);

	//Building model and solving problem
	build_model(sol->inst, env, lp);
	solve_lp(sol, env, lp);

	// Getting infos about final solution
	get_sol(sol, env, lp);

	if (VERB_PRINT_CPLEX_SOLUTION && !benchmark)
		print_cplex_solution(sol->xstar, sol);

	if (!benchmark)
		build_sol(sol->xstar, sol);

	// Closing CPLEX utilites 
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0;
}

// Variable Neighborhood Search
void VNS(solution* sol)
{
	const long long nnodes = sol->inst->nnodes;
	const double* xcoord = sol->inst->xcoord;
	const double* ycoord = sol->inst->ycoord;
	double residual_time = sol->inst->time_limit;
	
	int* best_succ = (int*)malloc(nnodes * sizeof(int));
	int* diversificated_succ = (int*)malloc(nnodes * sizeof(int));
	int min_kick = 31;
	int max_kick = 39;
	int current_kick = min_kick;

	double best_cost = DBL_MAX;
	double initial_cost, refined_cost = DBL_MAX;
	double latest_refined_cost, iteration_time;
	double diversificated_cost;
	int iterations = 0;

	double timer = -clock();
	
	do
	{
		iteration_time = -clock();
		if (iterations == 0)
		{
			printf("Getting initial solution...\n");
			first_nearest_neighbor(sol);
			greedy_refinement_2opt(sol);
			best_cost = initial_cost = sol->latest_cost;
			latest_refined_cost = -1;
			memcpy(best_succ, sol->succ, nnodes * sizeof(int));
		}
		else
		{
			// Diversification phase
			diversificate(sol, current_kick);
			memcpy(diversificated_succ, sol->succ, nnodes * sizeof(int));
			diversificated_cost = sol->latest_cost;
		}

		// Intensification phase
		refine(sol);

		// Printing info about incumbent solution
		refined_cost = sol->latest_cost;
		printf("Cost of refined solution: %.8lf (iter %d) \n", refined_cost, iterations + 1);

		if (best_cost - refined_cost > 1e-3)
		{
			best_cost = sol->latest_cost;
			memcpy(best_succ, sol->succ, nnodes * sizeof(int));
			printf(" +++ Current best solution updated at trial %i: %.8f\n", iterations + 1, best_cost);
			current_kick = min_kick;
		}
		else
		{
			if (current_kick == max_kick)
				current_kick = min_kick;
			else
				current_kick += 2;
		}

		iteration_time += clock();
		residual_time -= iteration_time / CLOCKS_PER_SEC;
		if (residual_time < 0)
			break;
		iterations++;
	} while (residual_time > 0);

	memcpy(sol->succ, best_succ, nnodes * sizeof(int));
	sol->latest_cost = best_cost;
	free(best_succ);
	free(diversificated_succ);

	printf("Computed %i iterations\n", iterations);
	printf("\nCost of initial solution: %.4lf \n", initial_cost);
	printf("Cost of best solution found: %.4lf \n", best_cost);
	printf("Real cost %lf\n", get_cost(sol));
	double rel_improv = (initial_cost - best_cost) / initial_cost;
	printf("Initial solution improved by %.2f %% \n\n", 100 * rel_improv);
	sol->latest_improvement = 100 * rel_improv;
}

// Intensification phase using 2-opt
void refine(solution * sol)
{
	// Getting pointers from solution
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;

	int* succ = sol->succ;

	int improvable = 1;
	int swapped = 0;
	const int print_step = nnodes / 100;

	while (improvable)
	{
		double min_delta = DBL_MAX;
		double current_delta;

		int best_node_a = -1;
		int best_node_b = -1;
		int best_node_c = -1;
		int best_node_d = -1;

		for (int i = 0; i < nnodes - 1; i++)
		{
			int node_a = i;
			int node_b = succ[i];
			for (int j = i + 1; j < nnodes; j++)
			{
				int node_c = j;
				int node_d = succ[j];

				// Not interested in adjacent arches
				if (node_a == node_d || node_c == node_b)
					continue;

				// Computing distances 
				double cost_ab = sqrt((xcoord[node_a] - xcoord[node_b]) * (xcoord[node_a] - xcoord[node_b]) +
					(ycoord[node_a] - ycoord[node_b]) * (ycoord[node_a] - ycoord[node_b]));
				double cost_cd = sqrt((xcoord[node_c] - xcoord[node_d]) * (xcoord[node_c] - xcoord[node_d]) +
					(ycoord[node_c] - ycoord[node_d]) * (ycoord[node_c] - ycoord[node_d]));

				double cost_ac = sqrt((xcoord[node_a] - xcoord[node_c]) * (xcoord[node_a] - xcoord[node_c]) +
					(ycoord[node_a] - ycoord[node_c]) * (ycoord[node_a] - ycoord[node_c]));
				double cost_bd = sqrt((xcoord[node_b] - xcoord[node_d]) * (xcoord[node_b] - xcoord[node_d]) +
					(ycoord[node_b] - ycoord[node_d]) * (ycoord[node_b] - ycoord[node_d]));

				// Variation of the obj func when operating inspected swap
				current_delta = cost_ac + cost_bd - (cost_ab + cost_cd);

				// Updating best swap
				if (current_delta < min_delta)
				{
					min_delta = current_delta;
					best_node_a = node_a;
					best_node_b = node_b;
					best_node_c = node_c;
					best_node_d = node_d;
				}
			}
		}

		// Improving swap found
		if (min_delta < 0)
		{
			improvable = 1;

			// Adding new arches and inverting the sense of a part of the tour
			// When getting info, succ has the meaning of prec!
			succ[best_node_a] = best_node_c; // new arch ac

			int to = best_node_b; // node b
			int from = succ[best_node_b]; // prec of node b
			int fromfrom = succ[from]; // prec of prec of node b

			succ[to] = best_node_d; // new arch bd
			while (to != best_node_c)
			{
				succ[from] = to;
				to = from;
				from = fromfrom;
				fromfrom = succ[fromfrom];
			}

			// Updating value of cost function
			sol->latest_cost += min_delta;
			//printf("\r2-opt refining: %4i tuples swapped", ++swapped);
		}
		// End of the local search
		else
			improvable = 0;
	}
	printf("\n");
}

// Diversification phase using KICK_SIZE number of arches
void diversificate(solution * sol, int kick_size)
{
	// Getting pointers from solution
	long long nnodes = sol->inst->nnodes;
	double* xcoord = sol->inst->xcoord;
	double* ycoord = sol->inst->ycoord;

	int* succ = sol->succ;

	typedef struct
	{
		int from;
		int to;
	} arc;

	arc* selection = (arc*)malloc(kick_size * sizeof(arc));

	// Selecting KICK_SIZE arches
	int from = 0;
	int steps = 0;
	for (int i = 0; i < kick_size; i++)
	{
		// Jump between lower and upper is num = (rand() % (upper – lower + 1)) + lower
		int jump = (rand() % ((int)(nnodes - steps) / (kick_size - i) - 4 + 1)) + 2;
		steps += (jump + 1);

		// Looking for next arc
		for (int j = 0; j < jump; j++)
			from = succ[from];

		selection[i].from = from;
		selection[i].to = succ[from];
	}

	double delta = 0;

	// Kicking 
	for (int i = 0; i < kick_size - 1; i++)
	{
		// Updating delta (removal of old arch)
		delta -= sqrt(INLINE_SQUARED_DISTANCE(xcoord[selection[i].from], xcoord[selection[i].to],
			ycoord[selection[i].from], ycoord[selection[i].to]));

		// Adding new arch and updating delta
		succ[selection[i].from] = selection[i + 1].to;
		delta += sqrt(INLINE_SQUARED_DISTANCE(xcoord[selection[i].from], xcoord[selection[i + 1].to],
			ycoord[selection[i].from], ycoord[selection[i + 1].to]));
	}
	// The last selected arch need to be reconnected with the first one! 
	// Removal
	delta -= sqrt(INLINE_SQUARED_DISTANCE(xcoord[selection[kick_size - 1].from], xcoord[selection[kick_size - 1].to],
		ycoord[selection[kick_size - 1].from], ycoord[selection[kick_size - 1].to]));
	// Adding
	succ[selection[kick_size - 1].from] = selection[0].to;
	delta += sqrt(INLINE_SQUARED_DISTANCE(xcoord[selection[kick_size - 1].from], xcoord[selection[0].to],
		ycoord[selection[kick_size - 1].from], ycoord[selection[0].to]));

	// Updating value of cost function
	sol->latest_cost += delta;
	free(selection);
}

void greedy_opt_search(solution * sol)
{
	printf("Total number of trials: %i\n", sol->trials);
	printf("Maximum checks per refinement function: ");
	if (sol->checks == CHECKS_NO_LIMITS)
		printf("UNLIMITED\n");
	else
		printf("%iM\n", sol->checks);

	double current_cost, cost_3opt, initial_cost;
	printf("Getting initial solution ...\n");
	first_nearest_neighbor(sol);
	int checks = sol->checks;
	sol->checks = CHECKS_NO_LIMITS;
	greedy_refinement_2opt(sol);
	sol->checks = checks;
	initial_cost = sol->latest_cost;
	printf("Initial cost: %lf\n", sol->latest_cost);

	int i = 0;
	do
	{
		printf("\n *** Trial %i\n", i + 1);
		current_cost = sol->latest_cost;
		greedy_refinement_2opt(sol);
		cost_3opt = sol->latest_cost;
		printf(" # Improvement: %.2lf%%\n", (current_cost - cost_3opt) / current_cost * 100.0);
		greedy_refinement_3opt(sol);
		printf(" # Improvement: %.2lf%%\n", (cost_3opt - sol->latest_cost) / cost_3opt * 100.0);
		i++;
	} while (sol->latest_cost != current_cost && i < sol->trials);

	printf("\n ### Total improvement: %.2lf%%\n\n", (initial_cost - sol->latest_cost) / initial_cost * 100.0);
}

// Symmetric model
long xpos(long i, long j, long n)
{
	if (i == j)
		print_error(" i == j in xpos");
	if (i > j)
		return j * n + i - ((j + 1L) * (j + 2L)) / 2L;
	return i * n + j - ((i + 1L) * (i + 2L)) / 2L;
}

// Asymmetric model (all positions in the matrix)
int xpos_asymm(int i, int j, int n)
{
	return (i * n + j);
}

// Position of the additional variables u in MTZ model
int upos_mtz(int i, int ustart)
{
	return (ustart + i);
}

// Position of the additional variables y in GG model
int ypos_gg(int i, int j, int n, int ystart)
{
	return (i * n + j + ystart);
}


