#include <stdlib.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include "mathfuncs.h"
#include "graphcut.h"
#include "nrutil.h"
#include "CornerPoints.h"
#include "utils.h"

void get_vertex_edges(CornerPoint *c, icVector3 *pts);
icMatrix3x3 get_rotation_to_global(Vertex *v, icVector3 *edge);
double compute_normal_curvature(icVector3 &uv, icVector3 &n);

void prepare_linear_system_kappa(Vertex *v, icVector3 *edges, double *kappa, double **x, icMatrix3x3 &M, int valence);
void prepare_linear_system_xTx(double *kappa, double **x, icMatrix3x3 &xTx, icVector3 &b, int valence);

void smooth_tensor(vertex_tensor **vt, int wtscheme, int iters, double lambda);
void smooth_tensor(vertex_tensor_list *vtl, vert_weights_list *vwts, int wtscheme, int iters, double lambda);
icVector3* copy_tensor_list(vertex_tensor_list *vtl);
void update_tensors(vertex_tensor_list *vtl, vert_weights_list *vwts, icVector3 *lmns, double lambda);
void compute_weighted_neighbor_tensors(vert_weights *vwts, icVector3 *lmns, double *dest);
void print_tensor_lmn(vertex_tensor_list *vtl);
double compute_gaussian_curvature(CornerPoint &o);
void set_plot_params(vertex_tensor_list *vtl);

void print_vector_lengths(icVector3 *vecs, int n);

double compute_gaussian_curvature(CornerPoint &c) {
	CornerPoint* o = &c;
	double curv = 0;
	do {
		curv = curv + o->angle;
		o = o->p->o->p;
	} while (c.index != o->index);
	return (2*PI - curv);
}

double compute_mean_curvature(CornerPoint &c) {
	CornerPoint* o = &c;
	double curv = 0;
	do {
		curv = curv + (COTANGENT(o->p->angle) + COTANGENT(o->p->o->angle))/2;
		o = o->p->o->p;
	} while (c.index != o->index);
	return curv;
}

tris_tensor_list *compute_tris_curvature_tensors(Polyhedron *poly, vertex_tensor_list *vtl) {
	tris_tensor_list *ttl = new tris_tensor_list(poly->ntris);
	for (int i = 0; i < ttl->length; i++) {
		Triangle *t = poly->tlist[i];
		tris_tensor *tt = new tris_tensor();
		tt->t = t;
		ttl->tt[t->index] = tt;

		vertex_tensor *v1 = vtl->vt[t->verts[0]->index];
		vertex_tensor *v2 = vtl->vt[t->verts[1]->index];
		vertex_tensor *v3 = vtl->vt[t->verts[2]->index];

		// centroid of the triangle
		tt->p.x = (v1->v->x + v2->v->x + v3->v->x)/3;
		tt->p.y = (v1->v->y + v2->v->y + v3->v->y)/3;
		tt->p.z = (v1->v->z + v2->v->z + v3->v->z)/3;

		// compute the reference axis at the face center.
		tt->ref_x.x = v2->v->x - v1->v->x; 
		tt->ref_x.y = v2->v->y - v1->v->y; 
		tt->ref_x.z = v2->v->z - v1->v->z;
		normalize(tt->ref_x);

		tt->ref_z.x = t->normal.x; tt->ref_z.y = t->normal.y; tt->ref_z.z = t->normal.z;

		icVector3 _y = cross(tt->ref_z, tt->ref_x);
		tt->ref_y.x = _y.x, tt->ref_y.y = _y.y; tt->ref_y.z = _y.z;

		// setup translation matrices
		tt->iM.entry[0][0] = tt->ref_x.x; tt->iM.entry[0][1] = tt->ref_y.x; tt->iM.entry[0][2] = tt->ref_z.x;
		tt->iM.entry[1][0] = tt->ref_x.y; tt->iM.entry[1][1] = tt->ref_y.y; tt->iM.entry[1][2] = tt->ref_z.y;
		tt->iM.entry[2][0] = tt->ref_x.z; tt->iM.entry[2][1] = tt->ref_y.z; tt->iM.entry[2][2] = tt->ref_z.z;

		//printf("Inverse Transform Matrix: %d\n",t->index);
		//print_3x3_matrix(tt->iM);

		// Tensor at point face center will be average of the vertex tensors.
		icMatrix2x2 &T = tt->T;
		T.entry[0][0] = (v1->T.entry[0][0] + v2->T.entry[0][0] + v3->T.entry[0][0])/3;
		T.entry[0][1] = (v1->T.entry[0][1] + v2->T.entry[0][1] + v3->T.entry[0][1])/3;
		T.entry[1][0] = T.entry[0][1];
		T.entry[1][1] = (v1->T.entry[1][1] + v2->T.entry[1][1] + v3->T.entry[1][1])/3;

		tt->computeLocalTensor();
	}
	return ttl;
}

double get_total_angle_deficit(vertex_tensor_list *vtl) {
	double deficit = 0;
	for (int i = 0; i < vtl->length; i++) {
		deficit += vtl->vt[i]->gauss_curv;
	}
	return deficit;
}

vertex_tensor_list *compute_curvature_tensors(Polyhedron *poly, int wtscheme, int iters, double lambda) {

	CornerPoints*cp = new CornerPoints(poly);
	CornerPoint *corners = cp->corners;

	vertex_tensor_list *vtl = new vertex_tensor_list(poly->nverts);
	vertex_tensor **vt = vtl->vt;

	cp->setValences();
	//cp->printValences();

	int npts = 0;
	icVector3 *edges = new icVector3[cp->maxvalence];

	double min_gauss_curv = FLT_MAX, max_gauss_curv = -FLT_MAX;
	double min_mean_curv = FLT_MAX, max_mean_curv = -FLT_MAX;

	std::vector<double> mean_curves;

	for (int i = 0; i < cp->length; i++) {

		if (!corners[i].mincornerforvertex) continue; // process one vertex only once

		get_vertex_edges(&corners[i], edges);

		//printf("Vertex: %d normal:\n",corners[i].v->index);
		//print_vector(corners[i].v->normal);

		vt[corners[i].v->index] = new vertex_tensor(corners[i].v, corners[i].valence);

		double gauss_curv = compute_gaussian_curvature(corners[i]);
		vt[corners[i].v->index]->gauss_curv = gauss_curv;
		min_gauss_curv = DMIN(min_gauss_curv, gauss_curv);
		max_gauss_curv = DMAX(max_gauss_curv, gauss_curv);

		vtl->min_gauss_curv = min_gauss_curv;
		vtl->max_gauss_curv = max_gauss_curv;

		double mean_curv = compute_mean_curvature(corners[i]);
		vt[corners[i].v->index]->mean_curv = mean_curv;
		min_mean_curv = DMIN(min_mean_curv, mean_curv);
		max_mean_curv = DMAX(max_mean_curv, mean_curv);
		mean_curves.push_back(mean_curv);

		//printf("Curvature of vertex %d: Gauss=%4.3f, Mean=%4.3f\n", vt[corners[i].v->index]->v->index, gauss_curv, mean_curv);

		vt[corners[i].v->index]->computeGlobalTensor(edges);

		//printf("Rotation Matrix:\n",corners[i].v->index);
		//print_3x3_matrix(corners[i].mv);
		//printf("\n");

	}
	
	//printf("Gauss: max: %4.3f min: %4.3f\n",min_gauss_curv,max_gauss_curv);
	//printf("Mean:  max: %4.3f min: %4.3f\n",min_mean_curv,max_mean_curv);
	std::sort(mean_curves.begin(), mean_curves.end());
	//printf("Sorted...\n");
	int max_pos = mean_curves.size();
	if (max_pos > 20) 
		max_pos = (int)((float)mean_curves.size()*0.95);
	max_mean_curv = mean_curves.at(max_pos-1);
	//printf("Mean:  max: %4.3f min: %4.3f\n",min_mean_curv,max_mean_curv);

	// Normalize the mean and gaussian curvature
	double gauss_range = max_gauss_curv - min_gauss_curv;
	double mean_range = max_mean_curv - min_mean_curv;
	for (int i = 0; i < vtl->length; i++) {
		vertex_tensor *v = vtl->vt[i];
		v->norm_gauss_curv = (v->gauss_curv - min_gauss_curv) / gauss_range;
		v->norm_mean_curv = (v->mean_curv - min_mean_curv) / mean_range;
		if (v->norm_mean_curv > 1.0) v->norm_mean_curv = 1.0;
		if (v->norm_mean_curv < 0.0) v->norm_mean_curv = 0.0;
	}

	// smooth
	vert_weights_list* vwts = get_vertex_weights(poly, cp, wtscheme);

	//printf("Before Smoothing: \n");
	//print_tensor_lmn(vtl);

	smooth_tensor(vtl, vwts, wtscheme, iters, lambda);

	//printf("After Smoothing: \n");
	//print_tensor_lmn(vtl);

	for (int i = 0; i < vtl->length; i++) {
		vertex_tensor *lt = vt[corners[i].v->index];
		lt->computeLocalTensor();
		//printf("Principal curvature axes:\n");
		//print_2x2_matrix(lt->curvAxis);
		//printf("principal curvatures: %4.3f, %4.3f\n", lt->kappas.entry[0], lt->kappas.entry[1]);
	}

	set_plot_params(vtl);

	delete edges;
	if (cp) delete cp;

	return vtl;

}

void set_plot_params(vertex_tensor_list *vtl) {
	std::vector<double> kappas;
	double edge_length = 0.0;
	vtl->edge_length = 0.0;
	for (int i = 0; i < vtl->length; i++) {
		edge_length += vtl->vt[i]->edge_length;
		kappas.push_back(abs(vtl->vt[i]->kappas.x));
		kappas.push_back(abs(vtl->vt[i]->kappas.y));
	}
	vtl->edge_length = edge_length / vtl->length;
	std::sort(kappas.begin(), kappas.end());
	if (kappas.size() > 20) {
		// keep the 90th  percentile
		int pos = (int)(((float)kappas.size())*0.90);
		vtl->max_kappa = kappas.at(pos-1);
	} else {
		vtl->max_kappa = kappas.at(kappas.size()-1);
	}
}

void _vertex_tensor::computeGlobalTensor(icVector3 *edges) {

	// compute the three orthogonal coordinates ox, oy, oz
	// ox -> with first edge as reference -> tangent to surface
	// oy -> perpendicular to normal and ox/tangent
	// oz -> vertex normal
	icVector3 tangent;
	for (int i = 0; i < valence; i++) {
		//printf("Edge:\n");
		//print_vector(edges[0]);
		tangent = edges[i] - dot(v->normal, edges[i])*(v->normal);
		//printf("Length tangent: %9.8f, valence: %d\n",length(tangent), valence);
		edge_length = length(tangent);
		if (edge_length > DBL_EPSILON) break;
		//printf("Trying next edge %d\n",i+1);
		//print_vector(edges[i]);
	}

	normalize(tangent);
	icVector3 third = cross(v->normal, tangent);
	//icVector3 third = cross(tangent, v->normal);
	normalize(third);

	//double c1,c2,c3;
	//c1 = dot(tangent, third);
	//c2 = dot(tangent, v->normal);
	//c3 = dot(v->normal, third);
	//printf("Checks: %4.3f, %4.3f, %4.3f\n", c1, c2, c3);

	iM = icMatrix3x3(
		tangent.x, third.x, v->normal.x, 
		tangent.y, third.y, v->normal.y, 
		tangent.z, third.z, v->normal.z
	);
	M = inverse(iM);

	ref_x.x = tangent.x;   ref_x.y = tangent.y;   ref_x.z = tangent.z;
	ref_y.x = third.x;     ref_y.y = third.y;     ref_y.z = third.z;
	ref_z.x = v->normal.x; ref_z.y = v->normal.y; ref_z.z = v->normal.z;

	// no need to normalize the ref_* axis since they are already normalized

	//printf("iM:\n");
	//print_3x3_matrix(iM);

	//icMatrix3x3 tM = transpose(M);
	//icMatrix3x3 check = multiply(tM, M);
	//printf("M^t * M:\n");
	//print_3x3_matrix(check);

	double *kappa = new double[valence];
	double **x = new double*[valence];
	for (int i = 0; i < valence; i++) {
		x[i] = new double[3];
	}
	//print_vector_lengths(edges, valence);
	prepare_linear_system_kappa(v, edges, kappa, x, M, valence);
	//printf("kappa for vertex %d\n",v->index);
	//printArray<double>(kappa, valence);

	// compute 3x3 matrix (x^T)x
	icMatrix3x3 xTx;
	icVector3 b;
	prepare_linear_system_xTx(kappa, x, xTx, b, valence);

	double sol[3];
	solve_linear(xTx, b, sol);
	T.entry[0][0] = sol[0];                 // l
	T.entry[0][1] = T.entry[1][0] = sol[1]; // m
	T.entry[1][1] = sol[2];                 // n

	//printf("Global Coordinate Tensor:\n");
	//print_2x2_matrix(T);

	//printf("inverted xTx:\n");
	//icMatrix3x3 ixTx = inverse(xTx);
	//print_3x3_matrix(ixTx);
	//printf("b:\n");
	//print_vector(b);

	for (int i = 0; i < valence; i++) delete x[i];
	delete x;

	delete kappa;

}

void print_vector_lengths(icVector3 *vecs, int n) {
	printf("len: ");
	for (int i = 0; i < n; i++) {
		printf(" %5.4f",length(vecs[i]));
	}
	printf("\n");
}

void _vertex_tensor::computeLocalTensor() {
	icMatrix3x3 _lt(0,0,0,0,0,0,0,0,0);
	_lt.entry[0][0] = T.entry[0][0];
	_lt.entry[1][1] = T.entry[1][1];
	_lt.entry[0][1] = _lt.entry[1][0] = T.entry[0][1];
	icMatrix3x3 iMt = transpose(iM);
	//printf("Transpose of iM:\n");
	//print_3x3_matrix(iMt);
	icMatrix3x3 m1 = multiply(iMt, _lt);
	//printf("Intermediate Coordinate Tensor:\n");
	//print_3x3_matrix(m1);

	icMatrix3x3 m2 = multiply(m1,iM);
	locT.entry[0][0] = m2.entry[0][0];
	locT.entry[0][1] = m2.entry[0][1];
	locT.entry[1][0] = m2.entry[1][0];
	locT.entry[1][1] = m2.entry[1][1];

	//printf("Local Coordinate Tensor:\n");
	//printf("A = [%4.3f,%4.3f;%4.3f,%4.3f]\n",locT.entry[0][0],locT.entry[0][1],locT.entry[1][0],locT.entry[1][1]);
	//print_2x2_matrix(locT);
	//print_3x3_matrix(m2);

	compute_eigen_symmetric_2x2(locT, curvAxis, kappas);

	//printf("Local Tensor Eigen:\n");
	//print_2x2_matrix(curvAxis);
	//printf("kappas: %4.3f, %4.3f\n",kappas.x,kappas.y);

	// first column of the principal curvature direction matrix
	minor = ref_x * curvAxis.entry[0][0] + ref_y * curvAxis.entry[1][0];
	// second column of the principal curvature direction matrix
	major = ref_x * curvAxis.entry[0][1] + ref_y * curvAxis.entry[1][1];

}

void print_tensor_lmn(vertex_tensor_list *vtl) {
	for (int i = 0; i < vtl->length; i++) {
		vertex_tensor *vt = vtl->vt[i];
		printf("%4.3f %4.3f %4.3f\n", vt->T.entry[0][0], vt->T.entry[0][1], vt->T.entry[1][1]);
	}
}

icVector3* copy_tensor_list(vertex_tensor_list *vtl) {
	icVector3 *copy = new icVector3[vtl->length];
	for (int i = 0; i < vtl->length; i++) {
		vertex_tensor *vt = vtl->vt[i];
		icMatrix2x2 &T = vt->T;
		copy[i].x = T.entry[0][0];
		copy[i].y = T.entry[0][1];
		copy[i].z = T.entry[1][1];
	}
	return copy;
}

void compute_weighted_neighbor_tensors(vert_weights *vwts, icVector3 *lmns, double *dest) {
	dest[0] = dest[1] = dest[2] = 0.0;
	icVector3 &cv = lmns[vwts->v->index];
	for (int i = 0; i < vwts->valence; i++) {
		double wt = vwts->wts[i];
		icVector3 &v = lmns[vwts->vidxs[i]];
		dest[0] += wt*(v.x - cv.x);
		dest[1] += wt*(v.y - cv.y);
		dest[2] += wt*(v.z - cv.z);
	}
}

void update_tensors(vertex_tensor_list *vtl, vert_weights_list *vwts, icVector3 *lmns, double lambda) {
	int nverts = vtl->length;
	double vdelta[3];
	for (int i = 0; i < nverts; i++) {
		compute_weighted_neighbor_tensors(vwts->wts[i], lmns, vdelta);
		vertex_tensor* vt = vtl->vt[i];
		vt->T.entry[0][0] += lambda * vdelta[0];
		vt->T.entry[0][1] += lambda * vdelta[1];
		vt->T.entry[1][1] += lambda * vdelta[2];
		vt->T.entry[1][0] = vt->T.entry[0][1];
	}
}

void smooth_tensor(vertex_tensor_list *vtl, vert_weights_list *vwts, int wtscheme, int iters, double lambda) {
	for (int iter = 0; iter < iters; iter++) {
		icVector3 *lmns = copy_tensor_list(vtl);
		update_tensors(vtl, vwts, lmns, lambda);
		delete lmns;
	}
}

void prepare_linear_system_xTx(double *kappa, double **x, icMatrix3x3 &xTx, icVector3 &b, int valence) {
	for (int i = 0; i < 3; i++) {
		b.entry[i] = 0;
		for (int j = 0; j < 3; j++) {
			xTx.entry[i][j] = 0;
			for (int k = 0; k < valence; k++) {
				xTx.entry[i][j] += x[k][i]*x[k][j];
			}
		}
	}
	for (int i = 0; i < 3; i++) {
		b.entry[i] = 0;
		for (int l = 0; l < valence; l++) {
			b.entry[i] += x[l][i]*kappa[l];
		}
	}
}

void prepare_linear_system_kappa(Vertex *v, icVector3 *edges, double *kappa, double **x, icMatrix3x3 &M, int valence) {
	for (int i = 0; i < valence; i++) {
		double a,b,c;
		kappa[i] = compute_normal_curvature(edges[i], v->normal);
		icVector3 ab_ = edges[i] - dot(v->normal, edges[i])*(v->normal);
		if (length(ab_) == 0) printf("WARN::Zero length edge found.\n");
		normalize(ab_);
		// rotate the edge so as to align with the global coordinates
		a = M.entry[0][0]*ab_.x + M.entry[0][1]*ab_.y + M.entry[0][2]*ab_.z;
		b = M.entry[1][0]*ab_.x + M.entry[1][1]*ab_.y + M.entry[1][2]*ab_.z;
		c = M.entry[2][0]*ab_.x + M.entry[2][1]*ab_.y + M.entry[2][2]*ab_.z; // this should be zero
		if (abs(c) > FLT_EPSILON) {
			printf("WARN::Translated edge: (%5.4f,%5.4f,%14.13f)\n",a,b,c);
		}
		x[i][0] = a*a;
		x[i][1] = 2*a*b;
		x[i][2] = b*b;
	}
}

double compute_normal_curvature(icVector3 &uv, icVector3 &n) {
	double l = length(uv);
	return -2*dot(uv,n)/(l*l); // negated since we want xi-xj and _not_ xj-xi
}

//_vertex_tensor::~_vertex_tensor() {
//	//if (edges) delete edges;
//}

// Matrix of the form:
// | a  b |
// | b  c |
void compute_eigen_symmetric_2x2(icMatrix2x2 &m, icMatrix2x2 &evec, icVector2 &eval) {
	if (abs(m.entry[0][1]-m.entry[1][0]) > FLT_EPSILON) {
		printf("WARN::2x2 matrix not symmetric in evaluating eigen\n");
		print_2x2_matrix(m);
	}
	double a = m.entry[0][0], b = m.entry[0][1], c = m.entry[1][1];
	//double det = determinant(m);
	double eval1 = ((a+c) + sqrt((a+c)*(a+c) - 4*(a*c - b*b)))/2;
	double eval2 = ((a+c) - sqrt((a+c)*(a+c) - 4*(a*c - b*b)))/2;
	double e1x = b, e1y = -(a - eval1); // corresponding to eigen val1
	double e2x = b, e2y = -(a - eval2); // corresponding to eigen val2
	double e1l = sqrt(e1x*e1x + e1y*e1y), e2l = sqrt(e2x*e2x + e2y*e2y);
	if (eval2 > eval1) {
		eval.entry[0] = eval1; eval.entry[1] = eval2;
		evec.entry[0][0] = e1l > FLT_EPSILON ? e1x/e1l : 0;
		evec.entry[1][0] = e1l > FLT_EPSILON ? e1y/e1l : 0;
		evec.entry[0][1] = e2l > FLT_EPSILON ? e2x/e2l : 0;
		evec.entry[1][1] = e2l > FLT_EPSILON ? e2y/e2l : 0;
	} else {
		eval.entry[0] = eval2; eval.entry[1] = eval1;
		evec.entry[0][0] = e2l > FLT_EPSILON ? e2x/e2l : 0;
		evec.entry[1][0] = e2l > FLT_EPSILON ? e2y/e2l : 0;
		evec.entry[0][1] = e1l > FLT_EPSILON ? e1x/e1l : 0;
		evec.entry[1][1] = e1l > FLT_EPSILON ? e1y/e1l : 0;
	}
}

void solve_linear(icMatrix3x3 &A, icVector3 &b, double *x) {
	icMatrix3x3 ia = inverse(A);
	x[0] = ia.entry[0][0]*b.x + ia.entry[0][1]*b.y + ia.entry[0][2]*b.z;
	x[1] = ia.entry[1][0]*b.x + ia.entry[1][1]*b.y + ia.entry[1][2]*b.z;
	x[2] = ia.entry[2][0]*b.x + ia.entry[2][1]*b.y + ia.entry[2][2]*b.z;
}

void print_vector(icVector3 &m) {
	printf("%4.3f %4.3f %4.3f\n",m.x,m.y,m.z);
}
void print_3x3_matrix(icMatrix3x3 &m) {
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			printf(" %4.3f",m.entry[i][j]);
		}
		printf("\n");
	}
}
void print_2x2_matrix(icMatrix2x2 &m) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			printf(" %4.3f",(float)m.entry[i][j]);
		}
		printf("\n");
	}
}

void get_vertex_edges(CornerPoint *c, icVector3 *pts) {
	CornerPoint* o = c;
	Vertex *v = c->v;
	int i = 0;
	do {
		Vertex *n = o->n->v;
		pts[i].x = n->x - v->x;
		pts[i].y = n->y - v->y;
		pts[i].z = n->z - v->z;
		//pts[i].x = v->x - n->x;
		//pts[i].y = v->y - n->y;
		//pts[i].z = v->z - n->z;
		i++;
		o = o->p->o->p;
	} while (c->index != o->index);
}
