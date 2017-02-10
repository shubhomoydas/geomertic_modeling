#include <stdlib.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <time.h>
#include "redsvd.hpp"
#include "mathfuncs.h"
#include "graphcut.h"
#include "nrutil.h"
#include "CornerPoints.h"
#include "utils.h"

typedef Eigen::Triplet<float> FTriplet;

inline float euc_sq_dist_tris(float *a1, float *a2, int l) {
	float dist = 0;
	for (int i = 0; i < 6; i++) {
		float a = (a2[i]-a1[i])/(i < 6 ? 1 : 10);
		dist += a*a;
	}
	return dist;
	//float d[5] = {a1[0]-a2[0], a1[1]-a2[1], a1[2]-a2[2], (a1[3]-a2[3])/10, (a1[4]-a2[4])/10};
	//return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2] + d[3]*d[3] + d[4]*d[4]);
}

void checkMatrix(REDSVD::SMatrixXf& mat, int mcols);

int bi_partition_with_graph_cut(
	Polyhedron *poly, tris_tensor_list *ttl, segments *segs, CommandOptions *options) {

	time_t start, end;

	int ntris = poly->ntris;

	int* partition = segs->segidxs;

	for (int i = 0; i < ntris; i++) {
		tris_tensor *tt = ttl->tt[i];
		Triangle *t = tt->t;
		icVector3 v1(t->verts[0]->x, t->verts[0]->y, t->verts[0]->z);
		icVector3 v2(t->verts[1]->x, t->verts[1]->y, t->verts[1]->z);
		icVector3 v3(t->verts[2]->x, t->verts[2]->y, t->verts[2]->z);
		tt->median = (v1+v2+v3)*(1/3);
	}

	std::vector<FTriplet> fv;
	fv.reserve(ntris*20);

	time(&start);
	std::set<int> bk;
	float *rowsums = new float[ntris];
	for (int i = 0; i < ntris; i++) rowsums[i] = 0;
	for (int i = 0; i < ntris; i++) {
		bk.clear();
		bk.insert(i); // triangle must not be self-neighbor
		Triangle *t = poly->tlist[i];
		tris_tensor *tt = ttl->tt[i];
		float d1[8] = {tt->median.x,tt->median.y,tt->median.z,t->normal.x,t->normal.y,t->normal.z,tt->kappas.x,tt->kappas.y};
		fv.push_back(FTriplet(i,i,0));
		for (int j = 0; j < t->nverts; j++) {
			Vertex *v = t->verts[j];
			for (int k = 0; k < v->ntris; k++) {
				int idx = v->tris[k]->index;
				// this will be a symmetric matrix
				if (idx > i && bk.count(idx) == 0) {
					bk.insert(idx);
					tris_tensor *ttx = ttl->tt[idx];
					Triangle *tn = ttx->t;
					float d2[8] = {ttx->median.x,ttx->median.y,ttx->median.z,tn->normal.x,tn->normal.y,tn->normal.z,ttx->kappas.x,ttx->kappas.y};
					//float dist = length(ttx->kappas - tt->kappas);
					float dist = euc_sq_dist_tris(d1,d2,8);
					dist = exp(-dist/10);
					rowsums[i  ] += dist;
					rowsums[idx] += dist;
					fv.push_back(FTriplet(i,idx,dist));
					fv.push_back(FTriplet(idx,i,dist));
				}
			}
		}
	}

	time(&end);
	printf("Time to fill Triplets: %f\n", ((float)(end-start)));

	time(&start);
	REDSVD::SMatrixXf adj(ntris, ntris);
	adj.setFromTriplets(fv.begin(), fv.end());
	printf("Time to fill sparse: %f\n", ((float)(end-start)));
	time(&end);

	time(&start);
	std::vector<FTriplet>::iterator it = fv.begin();
	for (int i=0; it != fv.end(); it++,i++) {
		//printf("%5.4f ",rowsums[i]);
		float val = it->value() / rowsums[it->row()];
		adj.coeffRef(it->row(), it->col()) = val;
		//printf("val(%4d, %4d)=%5.4f -> (%5.4f / %5.4f)\n", it->row(), it->col(), val, it->value(), rowsums[it->row()]);
	}
	printf("\n");
	time(&end);
	printf("Time to modify sparse: %f\n", ((float)(end-start)));

	printf("Triangles: %d; Adj: %d, %d\n",ntris, adj.rows(), adj.cols());
	checkMatrix(adj, ntris);
	
	time(&start);
	REDSVD::RedSymEigen eig(adj, IMIN(adj.rows(), 100));
	time(&end);
	printf("Time to compute eigen: %f\n", ((float)(end-start)));

	printf("Returned matrix: [%d, %d]\n", eig.eigenVectors().rows(), eig.eigenVectors().cols());
	
	int cn[2]; cn[0] = cn[1] = 0;
	for (int i = 0; i < ntris; i++) {
		float val = eig.eigenVectors()(i, 0);
		if (val < 0) {
			cn[0]++; 
			partition[i] = 0;
		} else {
			cn[1]++;
			partition[i] = 1;
		}
		if (options->debug)
			printf("%5.4f\n", val);
	}
	printf("#left: %d, #right: %d\n", cn[0], cn[1]);
	
	if (cn[0] > 0 && cn[1] > 0) segs->maxsegid++;

	//delete nn;
	delete rowsums;

	return 1;

}

void checkMatrix(REDSVD::SMatrixXf& mat, int mcols) {

	if (mat.rows() > 100) return;

	time_t start, end;

	Eigen::MatrixXf A(mat.rows(), mat.cols()); // = mat;
	printf("Mat(%d, %d); A(%d, %d)\n", mat.rows(), mat.cols(), A.rows(), A.cols());
	A.setIdentity();
	A = A * mat;
	int rows = A.rows(), cols = A.cols();
	//int rows = mat.rows(), cols = mat.cols();
	printf("Mat: %d,%d; print: %d,%d\n",rows,cols,rows,mcols);

	time(&start);
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(A);
	time(&end);
	printf("Time to compute eigen with dense: %f\n", ((float)(end-start)));
	if (eigensolver.info() != Eigen::Success) {
		printf("Could not find the values...\n");
	} else {
		printf("1st Eigen Value: %f\n",eigensolver.eigenvalues()(0,0));
	}

	if (A.rows() > 100) return;

	// DONOT call the below for big matrices
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < mcols; j++) {
			printf("%5.4f ",A(i,j));
		}
		printf("\n");
	}
}