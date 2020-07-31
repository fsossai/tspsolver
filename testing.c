#include "tsp.h"

#ifdef TESTING

#include "testing.h"
#include <Windows.h>
#include <time.h>

extern char benchmark;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		print_help();
		return 1;
	}
	if (VERB_PRINT_CLINE_ARGUMENTS)
		for (int i = 1; i < argc; i++)
			printf("arg %i\t: %s\n", i, argv[i]);

	benchmark = 1;
	testset ts;
	memset(&ts, 0, sizeof(testset));
	parse_test_command_line(argc, argv, &ts);
	init_testset(&ts);
	if (ts.comparing_heuristics)
	{
		read_opt_costs(&ts);
		compare_heuristics(&ts);
	}
	else if (ts.comparing_times)
	{
		compare_exec_time(&ts);
	}
	else
		test_all(&ts);

	close_testset(&ts);
	return 0;
}

void close_testset(testset * ts)
{
	for (int i = 0; i < ts->size; i++)
		close_instance(&ts->insts[i]);
	free(ts->insts);
	if (ts->comparing_heuristics)
		free(ts->opt_costs);
}

void init_testset(testset * ts)
{
	WIN32_FIND_DATA fdFile;
	HANDLE hFind = NULL;

	char path[STRING_MAX_SIZE + 1];
	char** inst_names = (char**)calloc(TESTSET_MAX_SIZE, sizeof(char**));
	int tsize = 0;

	sprintf(path, "%s\\%s", ts->directory, ts->testset_filename);
	hFind = FindFirstFile(path, &fdFile);
	if (hFind != INVALID_HANDLE_VALUE)
	{
		FILE* f = fopen(path, "r");
		char line[STRING_MAX_SIZE + 1];
		int last;
		if (f != NULL)
		{
			while (fgets(line, sizeof(line), f) != NULL)
			{
				last = strlen(line) - 1;
				if (line[0] == '#') // comments starts with '#'
					continue;
				// removing CR, LF and blanks from the tail
				while (line[last] == '\r' ||
					line[last] == '\n' ||
					line[last] == ' ')
					line[last--] = '\0';
				if (last != -1) // found a new valid file name
				{
					inst_names[tsize] = (char*)calloc(STRING_MAX_SIZE + 1, sizeof(char));
					sprintf(inst_names[tsize], "%s\\%s", ts->directory, line);
					tsize++;
				}
			}
			fclose(f);
		}
		else
			print_error("can't read \"testset.txt\"");
	}

	FindClose(hFind);

	if (tsize == 0)
	{
		sprintf(path, "%s\\*.tsp", ts->directory);
		hFind = FindFirstFile(path, &fdFile);
		if (hFind == INVALID_HANDLE_VALUE)
			print_error("Invalid test directory");
		do
		{
			if (!(fdFile.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
			{
				inst_names[tsize] = (char*)calloc(STRING_MAX_SIZE + 1, sizeof(char));
				sprintf(inst_names[tsize], "%s\\%s", ts->directory, fdFile.cFileName);
				tsize++;
			}
		} while (FindNextFile(hFind, &fdFile));

		FindClose(hFind);
	}

	ts->size = tsize;
	ts->insts = (instance*)calloc(tsize, sizeof(instance));
	for (int i = 0; i < tsize; i++)
	{
		strcpy(ts->insts[i].input_file, inst_names[i]);
		free(inst_names[i]);
	}
	free(inst_names);
}

void read_opt_costs(testset* ts)
{
	WIN32_FIND_DATA fdFile;
	HANDLE hFind = NULL;

	char path[STRING_MAX_SIZE + 1];
	double* opt_costs = (char**)calloc(ts->size, sizeof(double*));
	int tsize = 0;

	sprintf(path, "%s\\%s", ts->directory, ts->opt_filename);
	hFind = FindFirstFile(path, &fdFile);
	if (hFind != INVALID_HANDLE_VALUE)
	{
		FILE* f = fopen(path, "r");
		char line[STRING_MAX_SIZE + 1];
		int last;
		if (f != NULL)
		{
			while (fgets(line, sizeof(line), f) != NULL && tsize < ts->size)
			{
				last = strlen(line) - 1;
				if (line[0] == '#') // comments starts with '#'
					continue;
				// removing CR, LF and blanks from the tail
				while (line[last] == '\r' ||
					line[last] == '\n' ||
					line[last] == ' ')
					line[last--] = '\0';
				if (last != -1) // found a new valid file name
				{
					// removing the name at the beginning of the line
					int i = 0;
					while (line[i++] != ',')
						if (i >= STRING_MAX_SIZE)
							print_error("invalid opt cost file");
					opt_costs[tsize] = atof(&line[i]);
					tsize++;
				}
			}
			fclose(f);
		}
		else
			print_error("can't read \"opt.txt\"");
	}

	FindClose(hFind);
	
	ts->opt_costs = opt_costs;
	if (tsize < ts->size)
		printf("WARNING: opt costs file doesn't contains all of the optimal costs!\n");
}

void parse_test_command_line(int argc, char* argv[], testset * ts)
{
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-test_dir") == 0)
		{
			strcpy(ts->directory, argv[++i]);
		}
		else if (strcmp(argv[i], "-test_models") == 0)
		{
			strcpy(ts->model_types, argv[++i]);
		}
		else if (strcmp(argv[i], "-heuristics") == 0)
		{
			ts->comparing_heuristics = 1;
		}
		else if (strcmp(argv[i], "-times") == 0)
		{
			ts->comparing_times = 1;
		}
		else if (strcmp(argv[i], "-ts_file") == 0)
		{
			strcpy(ts->testset_filename, argv[++i]);
		}
		else if (strcmp(argv[i], "-opt_file") == 0)
		{
			strcpy(ts->opt_filename, argv[++i]);
		}
		else if (strcmp(argv[i], "-seed") == 0)
		{
			ts->random_seed = atoi(argv[++i]);
		}
	}
}

void test_all(testset * ts)
{
	double timer, sec;
	int* models = (int*)calloc(strlen(ts->model_types), sizeof(int));
	int nmodels = 0;
	for (int m = 0; m < strlen(ts->model_types); m++) // select which model types to test with
	{
		if (ts->model_types[m] >= '0' && ts->model_types[m] <= '9')
			models[nmodels++] = ts->model_types[m] - '0';
	}
	FILE* output = fopen("benchmark.txt", "w");
	if (output == NULL)
		print_error("can't create benchmark file");

	// header: column number and names
	fprintf(output, "%i,", nmodels);
	for (int m = 0; m < nmodels; m++)
		fprintf(output, "%s%s", get_model_name(models[m]), (m + 1 < nmodels) ? "," : "\n");

	solution sol;
	for (int i = 0; i < ts->size; i++)
	{
		printf("\nTesting problem %03i\t: %s\n", i, ts->insts[i].input_file);
		memset(&sol, 0x00, sizeof(solution));
		read_input(&ts->insts[i]);
		ts->insts[i].random_seed = ts->random_seed;
		ts->insts[i].time_limit = 3600;
		fprintf(output, "%s, ", ts->insts[i].input_file);

		for (int m = 0; m < nmodels; m++) // for each model type
		{
			sol.inst = &ts->insts[i];
			sol.inst->model_type = models[m];
			printf("  using model %i %-20s ... ", models[m], get_model_name(models[m]));
			fflush(stdout);

			timer = -clock();
			if (TSPopt(&sol))
				print_error(" error within TSPopt()");
			timer += clock();

			sec = timer / CLOCKS_PER_SEC;
			printf("\tsolved in %02i:%02i:%02i (%.3lf s)\n",
				(int)(sec / 3600),
				(int)((int)sec % 3600) / 60,
				(int)sec % 60,
				sec);

			fprintf(output, "%10.3lf%s", sec, (m + 1 < nmodels) ? "," : "\n");
			fflush(output);
			if (sol.comp)
				free(sol.comp);
			if (sol.succ)
				free(sol.succ);
			if (sol.nodes_order)
				free(sol.nodes_order);
			memset(&sol, 0x00, sizeof(solution));
		}
	}

	fclose(output);
	free(models);
}

void compare_heuristics(testset *ts)
{
	solution sol;
	memset(&sol, 0x00, sizeof(solution));
	const int M_FIRST_NN = 9;
	const int M_GRASP = 11;
	const int M_INSERTION = 12;
	const int M_HARD_FIXING = 7;
	const int M_LOCAL_BRANCHING = 8;
	const int M_SIMULATED_ANNEALING = 13;
	const int M_VNS = 14;

	FILE* output = fopen("heuristics_comparison.txt", "w");
	if (output == NULL)
		print_error("can't create comparison file");

	double time;
	for (int i = 0; i<ts->size; i++)
	{
		solution sol;
		memset(&sol, 0x00, sizeof(solution));
		read_input(&ts->insts[i]);
		sol.inst = &ts->insts[i];
		printf("\nTesting problem %03i\t: %s\n", i, ts->insts[i].input_file);
		
		fprintf(output, "%s ", sol.inst->input_file);

		printf("\t%s\t\t...\n", get_model_name(M_VNS));
		sol.inst->time_limit = 3600;
		sol.inst->model_type = M_VNS;
		double time = -clock();

		VNS(&sol);

		time += clock();
		printf("apx: %.2lf%%\n", (sol.latest_cost / ts->opt_costs[i] - 1.0) * 100.0);
		printf(" TIME %7.3lf \n", time / (double)CLOCKS_PER_SEC);
		fprintf(output, "%7.3lf ", (sol.latest_cost / ts->opt_costs[i] - 1.0) * 100.0);
		fprintf(output, "%7.3lf ", time / (double)CLOCKS_PER_SEC);
		fprintf(output, "%7.3lf \n", sol.latest_improvement);
		fflush(output);
		sol.inst = NULL;
		close_solution(&sol);
		sol.inst = &ts->insts[i];
	}

	fclose(output);
}

void compare_exec_time(testset* ts)
{
	solution sol;
	memset(&sol, 0x00, sizeof(solution));
	const int M_FIRST_NN = 9;
	const int M_GRASP = 11;
	const int M_INSERTION = 12;
	

	FILE* output = fopen("time_comparison.txt", "w");
	if (output == NULL)
		print_error("can't create comparison file");

	double time;
	for (int i = 0; i < ts->size; i++)
	{
		solution sol;
		memset(&sol, 0x00, sizeof(solution));
		read_input(&ts->insts[i]);
		sol.inst = &ts->insts[i];
		int* aux = (int*)malloc(sol.inst->nnodes * sizeof(int));
		printf("\nTesting problem %03i\t: %s\n", i, ts->insts[i].input_file);

		fprintf(output, "%s ", sol.inst->input_file);

		// ---> REFINEMENT TEST
		printf("\t%s\t\t...\t", get_model_name(M_FIRST_NN));
		first_nearest_neighbor(&sol);
		double initial_cost = sol.latest_cost;
		memcpy(aux, sol.succ, sol.inst->nnodes * sizeof(int));
		double time;

		time = -clock();
		refine(&sol);
		time += clock();
		fprintf(output, "%7.3lf ", (initial_cost - sol.latest_cost) / initial_cost * 100.0);
		fprintf(output, "%.3lf ", time / (double)CLOCKS_PER_SEC);
		printf("%7.3lf%% ", (initial_cost - sol.latest_cost) / initial_cost * 100.0);
		printf("%.3lf sec\n", time / (double)CLOCKS_PER_SEC);

		memcpy(sol.succ, aux, sol.inst->nnodes * sizeof(int)); //restoring solution

		time = -clock();
		sol.latest_cost = initial_cost;
		greedy_refinement_2opt(&sol);
		time += clock();
		fprintf(output, "%7.3lf ", (initial_cost - sol.latest_cost) / initial_cost * 100.0);
		fprintf(output, "%.3lf\n", time / (double)CLOCKS_PER_SEC);
		printf("%7.3lf%% ", (initial_cost - sol.latest_cost) / initial_cost * 100.0);
		printf("%.3lf sec\n", time / (double)CLOCKS_PER_SEC);
		sol.inst = NULL;
		close_solution(&sol);

		printf("\n");
		free(aux);
	}

	fclose(output);
}

#endif