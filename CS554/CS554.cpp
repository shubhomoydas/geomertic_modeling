// CS554.cpp : Defines the entry point for the console application.
//

#include "glut.h"
#include "graphutils.h"
#include <map>
#include "stdafx.h"
#include "mathfuncs.h"
#include "unittests.h"

int learnply_main(CommandOptions *options);

void test_tensors(int argc, char **argv);

int test_texture(int argc, char** argv);

void test_rndperm();

int _tmain(int argc, _TCHAR* argv[])
{
	//std::cout << "Hello World" << std::endl;
	//test_image_load(); exit(0);
	//test_sparse_vertex_weights(); exit(0);
	//test_sparse_solver(); exit(0);
	//test_sparse_matrix(); exit(0);

	char ** argn = allocate_argn(argc, argv);

	//test_texture(argc, argn); exit(0);
	//test_tensors(argc, argn); exit(0);
	//test_rndperm(); exit(0);

	CommandOptions *options = parse_command_line_options(argc, argn);
	if (options) {
		options->print_options();

		learnply_main(options);
	} else {
		printf("Invalid Commandline options.");
	}
	release_argn(argc, argn);

	return 0;
}

