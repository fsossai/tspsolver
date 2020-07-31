#include "tsp.h"
#include <windows.h>
#include <time.h>
#include <string.h>

void parse_command_line(int argc, char* argv[], solution *sol);
void print_help();
void print_sol_info(solution* sol);
void read_input(instance* inst);

#ifndef TESTING
int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		print_help();
		return 1;
	}
	if (VERB_PRINT_CLINE_ARGUMENTS)
		for (int i = 1; i < argc; i++)
			printf("arg %i\t: %s\n", i, argv[i]);

	solution sol;
	setup_solution(&sol);

	parse_command_line(argc, argv, &sol);
	read_input(sol.inst);

	printf("Model: %s \n\n", get_model_name(sol.inst->model_type));
	fflush(stdout);

	double timer = -clock();
	TSPopt(&sol);
	timer += clock();

	sol.execution_time = timer / CLOCKS_PER_SEC;

	if(VERB_PRINT_SOL_INFO)
		print_sol_info(&sol);

	if (GNUPLOT_OUTPUT)
		gnuplot_plt(&sol);

	close_solution(&sol);
	return 0;
}
#endif

void print_sol_info(solution* sol)
{
	printf("Cost of best integer solution: %.2lf \n", sol->latest_cost);
	printf("Value of best bound: %.2lf \n", sol->best_lb);
	printf("TSP problem solved in %02i:%02i:%02i (%.3lf s)\n",
		(int)(sol->execution_time / 3600),
		(int)((int)sol->execution_time % 3600) / 60,
		(int)sol->execution_time % 60,
		sol->execution_time);
}

void parse_command_line(int argc, char* argv[], solution * sol)
{
	// Setting default values for env parameters   
	instance* inst = sol->inst;
	inst->model_type = 0;
	strcpy(inst->input_file, "NULL");
	inst->cplex_write_model = 0;
	inst->cplex_model_file = NULL;
	inst->integer_costs = 0;
	inst->time_limit = (double)INT_MAX;
	inst->random_seed = time(0);
	inst->neighbors = 1;
	inst->lambda = 1.5;
	sol->trials = 1;
	sol->best_lb = -1;
	sol->checks = CHECKS_NO_LIMITS;

	int help = 0; if (argc < 1) help = 1;
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			
		if (strcmp(argv[i], "-model_type") == 0) { inst->model_type = atoi(argv[++i]); continue; } 	
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 
		if (strcmp(argv[i], "-time_limit") == 0) { inst->time_limit = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-random_seed") == 0) { inst->random_seed = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-neighbors") == 0) { inst->neighbors = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-trials") == 0) { sol->trials = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-checks") == 0) { sol->checks = atoi(argv[++i]); continue; }
		if (strcmp(argv[i], "-lambda") == 0) { inst->lambda = atof(argv[++i]); continue; }
		if (strcmp(argv[i], "-write_model") == 0)	
		{
			inst->cplex_write_model = 1;
			inst->cplex_model_file = (char*)malloc(STRING_MAX_SIZE);
			strcpy(inst->cplex_model_file, argv[++i]);
			continue;
		}
		help = 1;
	}

	if (help)
		print_help();
	srand(inst->random_seed);
}

void read_input(instance * inst) // simplified TSP parser
{
	FILE* fin = fopen(inst->input_file, "r");
	if (fin == NULL)
		print_error(" input file not found!");

	inst->nnodes = -1;

	char line[STRING_MAX_SIZE + 1];
	char* par_name;
	char* token1;
	char* token2;

	int active_section = 0; // =1 NODE_COORD_SECTION

	while (fgets(line, sizeof(line), fin) != NULL)
	{
		if (strlen(line) <= 1)
			continue; // skip empty lines
		par_name = strtok(line, " :");

		if (strncmp(par_name, "NAME", 4) == 0)
		{
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "COMMENT", 7) == 0)
		{
			active_section = 0;
			token1 = strtok(NULL, ":");
			printf(" Problem description\t: %s", token1);
			printf(" Model type\t\t: %i\n", inst->model_type);
			continue;
		}

		if (strncmp(par_name, "TYPE", 4) == 0)
		{
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "TSP", 3) != 0)
				print_error(" TYPE format error\t: unknown TYPE\n");
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "DIMENSION", 9) == 0)
		{
			if (inst->nnodes >= 0)
				print_error(" Format error\t: multiple DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			printf(" Number of nodes\t: %d\n", inst->nnodes);
			inst->xcoord = (double*)calloc(inst->nnodes, sizeof(double));
			inst->ycoord = (double*)calloc(inst->nnodes, sizeof(double));
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0)
		{
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "EUC_2D", 6) && strncmp(token1, "ATT", 3))
				print_error(" Format error\t: unknown EDGE_WEIGHT_TYPE");
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0)
		{
			if (inst->nnodes <= 0)
				print_error(" Format error\t: DIMENSION section must appear before NODE_COORD_SECTION section");
			active_section = 1;
			continue;
		}

		if (strncmp(par_name, "EOF", 3) == 0)
		{
			active_section = 0;
			break;
		}

		if (active_section == 1)
		{
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= inst->nnodes)
				print_error(" Format error\t: unknown NODE_COORD_SECTION section");
			token1 = strtok(NULL, " ");
			token2 = strtok(NULL, " ");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if (VERB_PRINT_FILE_INPUT_NODES)
				printf(" Node\t#%4i\t= (%9.3f\t, %9.3f\t)\n", i + 1, inst->xcoord[i], inst->ycoord[i]);
			continue;
		}

		print_error(" Format error\t: unknown statements");
	}

	fclose(fin);
}

void print_help()
{
	printf("\n\nAvailable parameters (vers.july-2020)--------------------------------------------\n");
	printf("-input\n");
	printf("-model_type [ 0 ... 15 ]\n");
	for (int i = 0; i<=15; i++)
		printf("\t%2i\t: %s\n", i, get_model_name(i));
	printf("-trials\t\t(for Simulated Annealing are expressed in millions\n");
	printf("-random_seed\n");
	printf("-time_limit\n");
	printf("\nenter --help for help\n");
	printf("-----------------------------------------------------------------------------------------\n");
}