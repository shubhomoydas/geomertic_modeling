#include <map>
#include <queue>
#include "stdafx.h"

#include "mathfuncs.h"

#include "stb_image_aug.h"
#include "unittests.h"
#include "segment.h"

void test_rndperm() {
	rndperm r(5);
	std::priority_queue< valpair*, std::vector<valpair*>, compareValpair > pq;
	int i = 0;
	while ((i = r.next()) >= 0) {
		printf("%d\n",i);
		pq.push(new valpair(i,rand()));
	}
	while (!pq.empty()) {
		valpair *v = pq.top();
		printf("Retrieved in order: %d, %f\n", v->i, v->val);
		pq.pop();
	}
}

void test_image_load() {
	int x, y, channels, force_channel=0;
	unsigned char *img = stbi_load( \
		"C:\\smd\\classes\\CS554\\HW3\\sample_textures\\3.jpg", \
		&x, &y, &channels, force_channel);
	printf ("x: %d, y: %d, channels: %d\n", x, y, channels);
	stbi_image_free(img);
}

void test_sparse_vertex_weights() {

	Vertex* a = new Vertex(0,0,0); a->index = 0;
	vert_weights* wa = new vert_weights(a, 2);
	wa->vidxs[0] = 1; wa->vidxs[1] = 3;
	wa->wts[0] = 0.5; wa->wts[1] = 0.5;

	Vertex* b = new Vertex(0,0,0); b->index = 1;
	vert_weights* wb = new vert_weights(b, 2);
	wb->vidxs[0] = 0; wb->vidxs[1] = 2;
	wb->wts[0] = 0.5; wb->wts[1] = 0.5;

	Vertex* c = new Vertex(0,0,0); c->index = 2;
	vert_weights* wc = new vert_weights(c, 2);
	wc->vidxs[0] = 1; wc->vidxs[1] = 3;
	wc->wts[0] = 0.5; wc->wts[1] = 0.5;

	Vertex* d = new Vertex(0,0,0); d->index = 3;
	vert_weights* wd = new vert_weights(d, 2);
	wd->vidxs[0] = 0; wd->vidxs[1] = 2;
	wd->wts[0] = 0.5; wd->wts[1] = 0.5;

	vert_weights_list* vwts = new vert_weights_list(4);
	vwts->wts[0] = wa;
	vwts->wts[1] = wb;
	vwts->wts[2] = wc;
	vwts->wts[3] = wd;

	std::map<int, double>* predefs = new std::map<int,double>();
	predefs->insert(std::pair<int,double>(0,1.0));
	predefs->insert(std::pair<int,double>(1,0.0));
	linear_system* ls = create_sparse_laplace_system(vwts, predefs);
	ls->a->print();
	
	double **m = new double*[2];
	m[0] = new double[2];
	m[1] = new double[2];
	ls->a->dense(m);
	printf("Matrix A:\n");
	print_matrix(m,2);
	printf("Matrix b:\n");
	print_matrix(&ls->b,1, 2);
	ls->solve();
	printf("Matrix x:\n");
	print_matrix(&ls->x,1, 2);

	double *ret = new double[4];
	ls->get_soln_with_predef(predefs, ret);
	printf("Matrix x with predefined vals:\n");
	print_matrix(&ret,1, 4);

	delete ret;
	delete ls;
	delete predefs;
	delete vwts; // deletes the vert_weights as well
	delete a; delete b; delete c; delete d;

}

void test_sparse_solver_2x2a() {
	double _a[2][2] = {
		{ 1, 2 },
		{ 4, 5 }
	};
	double *a[2] = {_a[0], _a[1]};
	double _b[2]  = { 3, 6 }, *b = _b;
	double _x[2] = { 0, 0 }, *x = _x;
	sparse s(a, 2);
	s.linear(b,x);
	s.print();
	printf("A=\n");
	print_matrix(a, 2);
	printf("b=\n");
	print_matrix(&b, 1, 2);
	printf("x=\n");
	print_matrix(&x, 1, 2); // [-1, 2]
}

void test_sparse_solver_2x2b() {
	double _a[2][2] = {
		{ 1, -0.5 },
		{ -0.5, 1 }
	};
	double *a[2] = {_a[0], _a[1]};
	double _b[2]  = { 0, 0.5 }, *b = _b;
	double _x[2] = { 0, 0 }, *x = _x;
	sparse s(a, 2);
	s.linear(b,x);
	s.print();
	printf("A=\n");
	print_matrix(a, 2);
	printf("b=\n");
	print_matrix(&b, 1, 2);
	printf("x=\n");
	print_matrix(&x, 1, 2); // [-1, 2]
}

void test_sparse_solver_3x3() {
	double _a[3][3] = {
		{ 1, 1, 1 },
		{ 0, 2, 2 },
		{ 0, 0, 2 }
	};
	double *a[3] = {_a[0], _a[1], _a[2]};
	double _b[3]  = { 2, -2, 2 }, *b = _b;
	double _x[3] = { 0, 0, 0 }, *x = _x;
	sparse s(a, 3);
	s.linear(b,x);
	s.print();
	printf("A=\n");
	print_matrix(a, 3);
	printf("b=\n");
	print_matrix(&b, 1, 3);
	printf("x=\n");
	print_matrix(&x, 1, 3); // [3, -2, 1]
}

void test_sparse_solver() {
	printf("Testing 2x2a\n");
	test_sparse_solver_2x2a();
	printf("Testing 3x3\n");
	test_sparse_solver_3x3();
	printf("Testing 2x2b\n");
	test_sparse_solver_2x2b();
}

void test_sparse_matrix() {
	double _a[5][5] = {
		{ 3.0, 0.0, 1.0, 0.0, 0.0 },
		{ 0.0, 4.0, 0.0, 0.0, 0.0 },
		{ 0.0, 7.0, 5.0, 9.0, 0.0 },
		{ 0.0, 0.0, 0.0, 0.0, 2.0 },
		{ 0.0, 0.0, 0.0, 6.0, 5.0 }
	};
	double *a[5] = {_a[0], _a[1], _a[2], _a[3], _a[4]};

	double **b = new double*[5];
	b[0] = new double[5];
	b[1] = new double[5];
	b[2] = new double[5];
	b[3] = new double[5];
	b[4] = new double[5];

	print_matrix(a,5);
	sparse s(a, 5);
	s.print();
	
	s.dense(b);

	print_matrix(b, 5);
}

void print_matrix(double **a, int n) {
	print_matrix(a, n, n);
}

void print_matrix(double **a, int nrow, int ncol) {
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			printf(" %+4.3f",a[i][j]);
		}
		printf("\n");
	}
}
