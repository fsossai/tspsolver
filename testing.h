#ifndef _TESTING_H
#define _TESTING_H

#include "tsp.h"

typedef struct
{
	int size;
	instance* insts;
	char directory[STRING_MAX_SIZE + 1];
	char model_types[STRING_MAX_SIZE + 1];
	char testset_filename[STRING_MAX_SIZE + 1];
	char opt_filename[STRING_MAX_SIZE + 1];
	char comparing_heuristics;
	char comparing_times;
	double* opt_costs;
	int random_seed;
} testset;

void close_testset(testset* ts);
void parse_test_command_line(int argc, char* argv[], testset* ts);
void init_testset(testset* ts);
void test_all(testset* ts);
void compare_heuristics(testset* ts);
void compare_exec_time(testset* ts);
void read_opt_costs(testset* ts);

#endif
