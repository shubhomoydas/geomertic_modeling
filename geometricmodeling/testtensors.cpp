#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <map>
#include "glut.h"
#include <string.h>
#include <fstream>
#include "ply.h"
#include "icVector.H"
#include "icMatrix.H"
#include "learnply.h"
#include "trackball.h"
#include "tmatrix.h"
#include "mathfuncs.h"
#include "graphutils.h"
#include "utils.h"
#include "CornerPoints.h"
#include "graphcut.h"

void test_solve();

void test_tensors(int argc, char **argv) {

	char *progname;
	int num = 1;
	FILE *this_file;
	float rotmat[4][4];

	Polyhedron *poly;
	CommandOptions *options = parse_command_line_options(argc, argv);
	if (options) {
		options->print_options();
	} else {
		printf("Invalid Commandline options.");
		exit(-1);
	}

	progname = options->argv[0];

	options->file = "C:/smd/classes/CS554/HW1/project1/tempmodels/bunny.ply";
	//options->file = "C:/smd/classes/CS554/HW1/project1/tempmodels/icosahedron.ply";
	this_file = fopen((options->file).c_str(), "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	vertex_tensor_list *vtl = compute_curvature_tensors(poly, WeightScheme::MEAN_CURV, 10, 0.1);
	tris_tensor_list *ttl = compute_tris_curvature_tensors(poly, vtl);

	segments *segs = segment_polyhedron(poly, ttl, options);

	delete segs;
	delete ttl;
	delete vtl;
	//exit(0);
	//poly->finalize();  // finalize everything
}

void test_solve() {
	icMatrix3x3 a(4, 1, -2, 2, -3, 3, -6, -2, 1);
	icVector3 b(0,9,0);
	double sol[3];
	solve_linear(a, b, sol); // correct sol: (3/4, -2, 1/2)
	printf("solution: x=%4.3f, y=%4.3f, z=%4.3f\n",sol[0],sol[1],sol[2]);
	exit(0);
}

