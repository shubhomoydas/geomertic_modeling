#include <stdlib.h>
#include <map>
#include "mathfuncs.h"
#include "nrutil.h"
#include "CornerPoints.h"
#include "utils.h"

void sort_vals(float* vals, int* dest, int num);
float** compute_covariance_from_aggregates(float** moments, int nverts);
void compute_eig_min_max(Vertex** vlist, int nverts, float** eigen_vec, float* min, float* max);
void compute_eigen(vertex_list_summary* summary);
void inverse_eigen(float** bounding_box, int n, float** eigen_vec);
void print_array(float** arr, int rows, int cols);

bool is_vertex_corner_min(CornerPoint& c);
int compute_valence(CornerPoint& c);
void populate_vertex_weights(CornerPoint& c, vert_weights* vwts, int scheme);
void smooth_polynomial(Polyhedron* poly, CornerPoints*cp, vert_weights_list* vwts, int wtscheme, int iters, double lambda);
Vertex** copy_vertex_list(Polyhedron* poly);
void compute_weighted_neighbor_vertex_coords(vert_weights* vwts, Vertex** verts, double* dest);
void update_poly_coords_for_smoothing(Polyhedron* poly, vert_weights_list* vwts, Vertex** vcopy, double lambda);
void update_poly_after_smoothing(Polyhedron* poly);

void compute_bounding_box(Vertex** vlist, vertex_list_summary* summary);
void setFacesAndCornersOfBoundingBox(vertex_list_summary* summary, float* min, float* max);
int compute_saddles(CornerPoint& c, double *diffusion, int *minverts, int *maxverts);
void compute_saddles(Polyhedron *poly, CornerPoints *cp, double *diffusion_colors, int *saddleverts, \
					 int *saddles, int nsaddles, int *minverts, int *maxverts);

#define INIT_3D(var, x, y, z) var = new float[3]; var[0]=x; var[1]=y; var[2]=z;
#define FACE_3D(var, a, b, c, d) var = new int[4]; var[0]=a; var[1]=b; var[2]=c; var[3]=d;

void compute_saddles(Polyhedron *poly, CornerPoints *cp, double *diffusion_colors, int *saddleverts, \
					 int *saddles, int nsaddles, int *minverts, int *maxverts) {
	CornerPoint* corners = cp->corners;
	*minverts = 0;
	*maxverts = 0;
	for (int i = 0; i < nsaddles; i++) saddles[i] = 0;
	for (int i = 0; i < cp->length; i++) {
		// check if this is the lowest indexed corner for the vertex
		// to ensure that each vertex is processed only once
		if (!is_vertex_corner_min(corners[i])) continue;
		Vertex* v = corners[i].v;
		CornerPoint& c = corners[i];
		int vidx = v->index;
		int nsaddle = compute_saddles(c, diffusion_colors, minverts, maxverts);
		if (nsaddle > 1) saddleverts[vidx] = 1;
		//printf("Saddle for %d: %d\n",vidx, nsaddle);
		if (nsaddle >= nsaddles) {
			printf("WARN:: Larger number of saddles found (%d) that requested (%d).\n",nsaddle,nsaddles);
		} else {
			saddles[nsaddle]++;
		}
	}
}

int compute_saddles(CornerPoint& c, double *diffusion, int *minverts, int *maxverts) {
	CornerPoint* o = &c;
	int direction = 0, olddirection = 0, nsaddle = 0;
	double currdiffusion = diffusion[o->v->index];
	// set the min and max to the first neighbor's values.
	double mi, mx; mi = mx = diffusion[o->n->v->index];
	do {
		double mn = diffusion[o->n->v->index];
		mi = FMIN(mi, mn);
		mx = FMAX(mx, mn);
		if (mn > currdiffusion)
			direction = 1;
		else
			direction = -1;
		if (olddirection == 1 && direction == -1)
			nsaddle++;
		olddirection = direction;
		o = o->p->o->p;
	} while (c.index != o->index);

	// process the last transition...
	if (diffusion[o->n->v->index] > currdiffusion)
		direction = 1;
	else
		direction = -1;
	if (olddirection == 1 && direction == -1)
		nsaddle++;

	if (currdiffusion < mi) (*minverts)++;
	if (currdiffusion > mx) (*maxverts)++;
	return nsaddle;
}

void set_diffusion_colors(Polyhedron *poly, double *diffusion_colors, int *saddleverts, std::map<int,double> *predefs, int wtscheme) {
	bool default_val = (predefs->size() < 2);
	for (int i = 0; i < poly->nverts; i++)
		saddleverts[i] = 0;
	if (default_val) {
		for (int i = 0; i < poly->nverts; i++)
			diffusion_colors[i] = 0.0;
		return;
	}
	CornerPoints*cp = new CornerPoints(poly);
	vert_weights_list* vwts = get_vertex_weights(poly, cp, wtscheme);
	linear_system* ls = create_sparse_laplace_system(vwts, predefs);
	//printf("in create sa[1]=%5.4f\n",ls->a->sa[1]);
	//ls->a->print();
	ls->solve();
	ls->get_soln_with_predef(predefs, diffusion_colors);
	int saddles[11], minverts, maxverts;
	compute_saddles(poly, cp, diffusion_colors, saddleverts, saddles, 11, &minverts, &maxverts);
	int m = minverts + maxverts;
	for (int i = 2; i < 11; i++) {
		printf("#%d-Saddles: %d\n",i-1,saddles[i]);
		m = m - saddles[i]*(i-1);
	}
	printf("#Min: %d, #Max: %d, M: %d\n",minverts, maxverts, m);
	//printf("Diffusion Colors: \n");
	//for (int i = 0; i < poly->nverts; i++) {
	//	printf("%5.4f\n",diffusion_colors[i]);
	//}
	//exit(0);
	delete ls;
	delete vwts;
	delete cp;
}

linear_system* create_sparse_laplace_system(vert_weights_list* vwts, std::map<int,double>* predefs) {
	int npre = predefs->size();
	int n = vwts->length - npre;
	int nmax = vwts->total_valence() + 1 + n;

	sparse* s = new sparse(n, nmax);
	linear_system* ls = new linear_system(s);
	double* x = ls->x, *b = ls->b;

	unsigned long* ija = s->ija; // IMPORTANT: sa, ija user indexing 1...nmax
	double* sa = s->sa;
	ija[1] = n + 2;
	sa[n+1] = 0; // arbitrary value
	int k = n + 1;
	vert_weights** wts = vwts->wts;

	// form a translation array to help ignore pre-defined values
	int* pre = new int[vwts->length];
	for (int i = 0, _i = 0; i < vwts->length; i++) {
		if (predefs->count(i) > 0) {
			pre[i] = -1;
		} else {
			pre[i] = _i++;
		}
		//printf("%d -> %d\n",i,pre[i]);
	}

	// index for sa incremented by 1 since sa indexes relative 1, NOT 0
	for (int i = 0; i < n; i++) sa[i+1] = 1.0;
	for (int i = 0; i < vwts->length; i++) {
		int idx = pre[i]; // translated index ignoring predefined variables
		if (idx < 0) continue; // found pre-defined variable
		b[idx] = 0;
		for (int j = 0; j < wts[i]->valence; j++) {
			int vidx = wts[i]->vidxs[j];
			if (pre[vidx] < 0) {
				b[idx] += predefs->find(vidx)->second * wts[i]->wts[j];
			} else {
				k++;
				sa[k] = -wts[i]->wts[j];
				// index incremented by 1 since ija indexes relative 1, NOT 0
				ija[k] = pre[vidx]+1;
			}
		}
		ija[idx+1+1] = k+1;
		//if (sa[1] < 1.0) printf("Error::sa[1]=%5.4f at i=%d, k=%d\n",sa[1], i, k);
	}
	delete pre;
	return ls;
}

void update_poly_after_smoothing(Polyhedron* poly) {
	for (int i=0; i < poly->nedges; i++) {
		delete(poly->elist[i]->tris);
		delete(poly->elist[i]);
	}
	delete poly->elist;
	poly->elist = NULL;

	// re-create all edges and other details
	poly->initialize();
	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();
}

void smooth_polynomial(Polyhedron* poly, int wtscheme, int iterscheme, int iters, double lambda) {
	if (iterscheme == IterScheme::INITIAL) {
		smooth_polynomial(poly, wtscheme, iters, lambda);
	} else if (iterscheme == IterScheme::RECOMPUTE) {
		for (int i = 0; i < iters; i++) {
			smooth_polynomial(poly, wtscheme, 1, lambda);
		}
	} else {
		printf("ERR:: Invalid iteration scheme %d\n",iterscheme);
	}
}

void smooth_polynomial(Polyhedron* poly, int wtscheme, int iters, double lambda) {
	if (iters <= 0) return;
	CornerPoints*cp = new CornerPoints(poly);
	//poly->print_vertices();
	//cp->print_all_corners();
	vert_weights_list* vwts = get_vertex_weights(poly, cp, wtscheme);
	//printf("#Vertex Weights: %d, Total Valence: %d\n",vwts->length,vwts->total_valence());
	//vwts->print();
	smooth_polynomial(poly, cp, vwts, wtscheme, iters, lambda);
	printf("Completed Smoothing...\n");
	//poly->print_vertices();
	delete vwts;
	delete cp;
	//exit(0);
}

void smooth_polynomial(Polyhedron* poly, CornerPoints*cp, vert_weights_list* vwts, int wtscheme, int iters, double lambda) {
	int nverts = poly->nverts;
	for (int iter = 0; iter < iters; iter++) {
		Vertex** vcopy = copy_vertex_list(poly);
		update_poly_coords_for_smoothing(poly, vwts, vcopy, lambda);
		release_objects<Vertex>(vcopy, poly->nverts);
	}
	update_poly_after_smoothing(poly);
}

void update_poly_coords_for_smoothing(Polyhedron* poly, vert_weights_list* vwts, Vertex** vcopy, double lambda) {
	int nverts = poly->nverts;
	double vdelta[3];
	for (int i = 0; i < nverts; i++) {
		compute_weighted_neighbor_vertex_coords(vwts->wts[i], vcopy, vdelta);
		Vertex* v = poly->vlist[i];
		v->x += lambda * vdelta[0];
		v->y += lambda * vdelta[1];
		v->z += lambda * vdelta[2];
	}
}

void compute_weighted_neighbor_vertex_coords(vert_weights* vwts, Vertex** verts, double* dest) {
	dest[0] = dest[1] = dest[2] = 0.0;
	Vertex* cv = vwts->v;
	for (int i = 0; i < vwts->valence; i++) {
		double wt = vwts->wts[i];
		Vertex *v = verts[vwts->vidxs[i]];
		dest[0] += wt*(v->x - cv->x);
		dest[1] += wt*(v->y - cv->y);
		dest[2] += wt*(v->z - cv->z);
	}
}

vert_weights_list* get_vertex_weights(Polyhedron* poly, CornerPoints* cp, int scheme) {
	vert_weights_list* vwtslist = new vert_weights_list(poly->nverts);
	vert_weights** vwts = vwtslist->wts;
	CornerPoint* corners = cp->corners;
	for (int i = 0; i < cp->length; i++) {
		// check if this is the lowest indexed corner for the vertex
		// to ensure that each vertex is processed only once
		if (!is_vertex_corner_min(corners[i])) continue;
		Vertex* v = corners[i].v;
		CornerPoint& c = corners[i];
		int vidx = v->index;
		int valence = compute_valence(c);
		//printf("Valence for %d: %d\n",vidx, valence);
		vwts[vidx] = new vert_weights(v, valence);
		populate_vertex_weights(c, vwts[vidx], scheme);
		vwts[vidx]->normalize();
		//vwts[vidx]->print(); printf("\n");
	}
	return vwtslist;
}

int order_indexes(vert_weights* vwts, int index) {
	int i = vwts->valence - 2;
	while (i >= 0 && index < vwts->vidxs[i]) {
		if (vwts->vidxs[i] < INT_MAX) {
			vwts->vidxs[i+1] = vwts->vidxs[i];
			vwts->wts[i+1] = vwts->wts[i];
		}
		i--;
	}
	return i+1;
}

void populate_vertex_weights(CornerPoint& c, vert_weights* vwts, int scheme) {
	CornerPoint* o = &c;
	int i = 0;
	do {
		i = order_indexes(vwts, o->n->v->index);
		vwts->vidxs[i] = o->n->v->index;
		switch(scheme) {
		case WeightScheme::UNIF:
			vwts->wts[i] = 1;
			break;
		case WeightScheme::CORD:
			// prev corner is opposite to edge connecting c to next corner
			vwts->wts[i] = 1. / o->p->e->length;
			break;
		case WeightScheme::MEAN_CURV:
			vwts->wts[i] = (COTANGENT(o->p->angle) + COTANGENT(o->p->o->angle))/2;
			break;
		case WeightScheme::MEAN_VAL:
			vwts->wts[i] = (tan(o->angle / 2) + tan(o->p->o->p->angle / 2))/2;
			break;
		default:
			printf("WARN:: Unknown Weight Scheme: %d\n", scheme);
			break;
		}
		o = o->p->o->p;
	} while (c.index != o->index);
}

Vertex** copy_vertex_list(Polyhedron* poly) {
	int nverts = poly->nverts;
	Vertex** copy = new Vertex*[nverts];
	for (int i = 0; i < nverts; i++) {
		copy[i] = new Vertex(poly->vlist[i]->x, poly->vlist[i]->y, poly->vlist[i]->z);
	}
	return copy;
}

bool is_vertex_corner_min(CornerPoint& c) {
	CornerPoint* o = c.p->o->p;
	while (c.index < o->index) {
		o = o->p->o->p;
	}
	return (c.index == o->index);
}

int compute_valence(CornerPoint& c) {
	CornerPoint* o = &c;
	int valence = 0;
	do {
		valence++;
		o = o->p->o->p;
	} while (c.index != o->index);
	return valence;
}

void sort_vals(float* vals, int* dest, int num) {
	for (int i = 0; i < num; i++) dest[i] = i;
	// bubble sort
	for (int i = num-1; i > 0; i--) {
		for (int j = 0; j < i; j++) {
			if (vals[dest[j]] > vals[dest[j+1]]) {
				int t = dest[j];
				dest[j] = dest[j+1];
				dest[j+1] = t;
			}
		}
	}
}

void compute_bounding_box(Vertex** vlist, vertex_list_summary* summary) {
	float min[3], max[3];
	compute_eig_min_max(vlist, summary->nverts, summary->eigen_vec, min, max);
	setFacesAndCornersOfBoundingBox(summary, min, max);
}

void setFacesAndCornersOfBoundingBox(vertex_list_summary* summary, float* min, float* max) {
	float** bounding_box = new float*[8];
	INIT_3D(bounding_box[0], min[0], min[1], min[2]); // 0
	INIT_3D(bounding_box[1], min[0], min[1], max[2]); // 1
	INIT_3D(bounding_box[2], min[0], max[1], max[2]); // 2
	INIT_3D(bounding_box[3], min[0], max[1], min[2]); // 3
	INIT_3D(bounding_box[4], max[0], min[1], min[2]); // 4
	INIT_3D(bounding_box[5], max[0], min[1], max[2]); // 5
	INIT_3D(bounding_box[6], max[0], max[1], min[2]); // 6
	INIT_3D(bounding_box[7], max[0], max[1], max[2]); // 7
	int** faces = new int*[6];
	FACE_3D(faces[0], 1, 2, 3, 0); // green
	FACE_3D(faces[1], 2, 7, 6, 3); // blue
	FACE_3D(faces[2], 4, 6, 7, 5); // green
	FACE_3D(faces[3], 0, 4, 5, 1); // blue
	FACE_3D(faces[4], 0, 3, 6, 4); // red
	FACE_3D(faces[5], 2, 1, 5, 7); // red
	inverse_eigen(bounding_box, 8, summary->eigen_vec);
	summary->bounding_box = bounding_box;
	summary->box_faces = faces;
}

void inverse_eigen(float** bounding_box, int n, float** eigen_vec) {
	float t[3];
	for (int i = 0; i < n; i++) {
		t[0]=bounding_box[i][0]; t[1]=bounding_box[i][1]; t[2]=bounding_box[i][2];
		bounding_box[i][0] = eigen_vec[0][0]*t[0]+eigen_vec[0][1]*t[1]+eigen_vec[0][2]*t[2];
		bounding_box[i][1] = eigen_vec[1][0]*t[0]+eigen_vec[1][1]*t[1]+eigen_vec[1][2]*t[2];
		bounding_box[i][2] = eigen_vec[2][0]*t[0]+eigen_vec[2][1]*t[1]+eigen_vec[2][2]*t[2];
	}
}

void compute_eig_min_max(Vertex** vlist, int nverts, float** eigen_vec, float* min, float* max) {
	float t[3];
	for (int v = 0; v < nverts; v++) {
		Vertex* vert = vlist[v];
		t[0] = (float)(vert->x*eigen_vec[0][0]+vert->y*eigen_vec[1][0]+vert->z*eigen_vec[2][0]);
		t[1] = (float)(vert->x*eigen_vec[0][1]+vert->y*eigen_vec[1][1]+vert->z*eigen_vec[2][1]);
		t[2] = (float)(vert->x*eigen_vec[0][2]+vert->y*eigen_vec[1][2]+vert->z*eigen_vec[2][2]);
		if (v == 0) {
			min[0] = t[0]; min[1] = t[1]; min[2] = t[2];
			max[0] = t[0]; max[1] = t[1]; max[2] = t[2];
		} else {
			min[0] = FMIN(min[0],t[0]); min[1] = FMIN(min[1],t[1]); min[2] = FMIN(min[2],t[2]);
			max[0] = FMAX(max[0],t[0]); max[1] = FMAX(max[1],t[1]); max[2] = FMAX(max[2],t[2]);
		}
	}
}

void print_array(float** arr, int rows, int cols) {
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			printf("%8.7f ", arr[i][j]);
		}
		printf("\n");
	}
}

/*
Returns the aggregates (sum, squared sum, pairwise product sum, min, max) 
and covariances between the three dimensions.
*/
vertex_list_summary* compute_vertex_list_summary(Vertex** vlist, int nverts, bool withnorm) {
	float** aggregates = 0;
	if (withnorm)
		aggregates = compute_vertex_norm_aggregates(vlist, nverts);
	else
		aggregates = compute_vertex_aggregates(vlist, nverts);
	float** cov = compute_covariance_from_aggregates(aggregates, nverts);
	vertex_list_summary* summary = new vertex_list_summary(aggregates, cov, nverts);
	compute_eigen(summary);
	return summary;
}

/*
Returns components for calculating moments and variances for 3d vertices.
The output is a 3x3 matrix having the following values:
sum(x)   sum(y)   sum(z)
sum(x^2) sum(y^2) sum(z^2)
sum(x*y) sum(x*z) sum(y*z)
max(x)   max(y)   max(z)
min(x)   min(y)   min(z)
*/
float** compute_vertex_aggregates(Vertex** vlist, int nverts) {
	if (nverts == 0) {return 0;}
	float** aggregates = new float*[5];
	for (int i = 0; i < 5; i++) {aggregates[i] = new float[3];}
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 3; j++) {
			aggregates[i][j] = 0.0;
		}
	}
	float minx = (float)vlist[0]->x, miny = (float)vlist[0]->y, minz = (float)vlist[0]->z;
	float maxx = (float)vlist[0]->x, maxy = (float)vlist[0]->y, maxz = (float)vlist[0]->z;
	for (int i = 0; i < nverts; i++) {
		aggregates[0][0] += (float)vlist[i]->x;
		aggregates[0][1] += (float)vlist[i]->y;
		aggregates[0][2] += (float)vlist[i]->z;

		aggregates[1][0] += (float)(vlist[i]->x * vlist[i]->x);
		aggregates[1][1] += (float)(vlist[i]->y * vlist[i]->y);
		aggregates[1][2] += (float)(vlist[i]->z * vlist[i]->z);

		aggregates[2][0] += (float)(vlist[i]->x * vlist[i]->y);
		aggregates[2][1] += (float)(vlist[i]->x * vlist[i]->z);
		aggregates[2][2] += (float)(vlist[i]->y * vlist[i]->z);
		
		minx = FMIN(minx, (float)vlist[i]->x);
		miny = FMIN(miny, (float)vlist[i]->y);
		minz = FMIN(minz, (float)vlist[i]->z);

		maxx = FMAX(minx, (float)vlist[i]->x);
		maxy = FMAX(miny, (float)vlist[i]->y);
		maxz = FMAX(minz, (float)vlist[i]->z);
	}
	aggregates[3][0] = minx;
	aggregates[3][1] = miny;
	aggregates[3][2] = minz;
	aggregates[4][0] = maxx;
	aggregates[4][1] = maxy;
	aggregates[4][2] = maxz;
	return aggregates;
}

/*
Returns components for calculating moments and variances for 3d vertex norms.
The output is a 3x3 matrix having the following values:
sum(x)   sum(y)   sum(z)
sum(x^2) sum(y^2) sum(z^2)
sum(x*y) sum(x*z) sum(y*z)
max(x)   max(y)   max(z)
min(x)   min(y)   min(z)
*/
float** compute_vertex_norm_aggregates(Vertex** vlist, int nverts) {
	if (nverts == 0) {return 0;}
	float** aggregates = new float*[5];
	for (int i = 0; i < 5; i++) {aggregates[i] = new float[3];}
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 3; j++) {
			aggregates[i][j] = 0.0;
		}
	}
	float minx = (float)vlist[0]->normal.x, miny = (float)vlist[0]->normal.y, minz = (float)vlist[0]->normal.z;
	float maxx = (float)vlist[0]->normal.x, maxy = (float)vlist[0]->normal.y, maxz = (float)vlist[0]->normal.z;
	for (int i = 0; i < nverts; i++) {
		aggregates[0][0] += (float)vlist[i]->normal.x;
		aggregates[0][1] += (float)vlist[i]->normal.y;
		aggregates[0][2] += (float)vlist[i]->normal.z;

		aggregates[1][0] += (float)(vlist[i]->normal.x * vlist[i]->normal.x);
		aggregates[1][1] += (float)(vlist[i]->normal.y * vlist[i]->normal.y);
		aggregates[1][2] += (float)(vlist[i]->normal.z * vlist[i]->normal.z);

		aggregates[2][0] += (float)(vlist[i]->normal.x * vlist[i]->normal.y);
		aggregates[2][1] += (float)(vlist[i]->normal.x * vlist[i]->normal.z);
		aggregates[2][2] += (float)(vlist[i]->normal.y * vlist[i]->normal.z);
		
		minx = FMIN(minx, (float)vlist[i]->normal.x);
		miny = FMIN(miny, (float)vlist[i]->normal.y);
		minz = FMIN(minz, (float)vlist[i]->normal.z);

		maxx = FMAX(minx, (float)vlist[i]->normal.x);
		maxy = FMAX(miny, (float)vlist[i]->normal.y);
		maxz = FMAX(minz, (float)vlist[i]->normal.z);
	}
	aggregates[3][0] = minx;
	aggregates[3][1] = miny;
	aggregates[3][2] = minz;
	aggregates[4][0] = maxx;
	aggregates[4][1] = maxy;
	aggregates[4][2] = maxz;
	return aggregates;
}

void compute_eigen(vertex_list_summary* summary) {
	if (!summary->cov) return;
	int nrot = 10;
	float** cov = new float*[3]; // temporary array that gets destroyed in eigen computation
	float** eigen_vec = new float*[3];
	float* eigen_val = new float[3];
	for (int i = 0; i < 3; i++) {
		cov[i] = new float[3];
		eigen_vec[i] = new float[3];
		eigen_val[i] = 0.0;
	}
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cov[i][j] = summary->cov[i][j];
			eigen_vec[i][j] = 0.0;
		}
	}
	for (int i = 0; i < 3; i++) {cov[i]--; eigen_vec[i]--;} // adjust for index-1 addressing
	jacobi(cov-1, 3, eigen_val-1, eigen_vec-1, &nrot);
	for (int i = 0; i < 3; i++) {cov[i]++; eigen_vec[i]++;} // readjust for index-0 addressing

	// sort according to eigen values
	int p[3];
	sort_vals(eigen_val, p, 3);
	float** sorted_eigen_vec = new float*[3];
	float* sorted_eigen_val = new float[3];
	for (int i = 0; i < 3; i++) {
		sorted_eigen_val[i] = eigen_val[p[i]];
		sorted_eigen_vec[i] = new float[3];
		for (int j = 0; j < 3; j++) {
			sorted_eigen_vec[i][j] = eigen_vec[i][p[j]];
		}
	}

	summary->eigen_vec = sorted_eigen_vec;
	summary->eigen_val = sorted_eigen_val;

	delete eigen_val;
	delete_pointer_array((void**)eigen_vec, 3);
	delete_pointer_array((void**)cov, 3);
}

float** compute_covariance_from_aggregates(float** moments, int nverts) {
	float** cov = new float*[3];
	cov[0] = new float[3];
	cov[1] = new float[3];
	cov[2] = new float[3];

	// Compute diagonal elements
	// var X = sum(X^2) - (sum(X)^2)/n
	cov[0][0] = moments[1][0] - (moments[0][0]*moments[0][0])/(float)nverts;
	cov[1][1] = moments[1][1] - (moments[0][1]*moments[0][1])/(float)nverts;
	cov[2][2] = moments[1][2] - (moments[0][2]*moments[0][2])/(float)nverts;
	
	// Compute off-diagonal elements
	// cov XY = sum(XY) - sum(X)*sum(Y)/n
	cov[0][1] = cov[1][0] = moments[2][0] - (moments[0][0]*moments[0][1])/(float)nverts;
	// cov XZ = sum(XZ) - sum(X)*sum(Z)/n
	cov[0][2] = cov[2][0] = moments[2][1] - (moments[0][0]*moments[0][2])/(float)nverts;
	// cov YZ = sum(YZ) - sum(Y)*sum(Z)/n
	cov[1][2] = cov[2][1] = moments[2][2] - (moments[0][1]*moments[0][2])/(float)nverts;

	return cov;
}
float** compute_vertex_covariance(Vertex** vlist, int nverts) {
	float** aggregates = compute_vertex_aggregates(vlist, nverts);
	float** cov = compute_covariance_from_aggregates(aggregates, nverts);
	delete_pointer_array((void**)aggregates, 3);
	return cov;
}

float** vertex_to_array(Vertex** vlist, int nverts) {
	float** arr = new float*[nverts];
	for (int i = 0; i < nverts; i++) {
		arr[i] = new float[3];
		arr[i][0] = (float)vlist[i]->x;
		arr[i][1] = (float)vlist[i]->y;
		arr[i][2] = (float)vlist[i]->z;
	}
	return arr;
}

/*
Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. On
output, elements of a above the diagonal are destroyed. d[1..n] returns the eigenvalues of a.
v[1..n][1..n] is a matrix whose columns contain, on output, the normalized eigenvectors of
a. nrot returns the number of Jacobi rotations that were required.
*/
void jacobi(float **a, int n, float d[], float **v, int *nrot)
{
	int j,iq,ip,i;
	float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
	b=vector(1,n);
	z=vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) { 
			for (iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=0.2*sm/(n*n);
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=100.0*fabs(a[ip][iq]);
				if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((float)(fabs(h)+g) == (float)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}

int* permute(int n) {
	int* shuffle = new int[n];
	for (int i = 0; i < n; i++) shuffle[i] = i;
	for (int i = 1, j = n; i < n; i++, j--) {
		int pos = rand() % j;
		int t = shuffle[pos];
		shuffle[pos] = shuffle[j-1];
		shuffle[j-1] = t;
	}
	return shuffle;
}

int compute_3d_color_interval_length(int num) {
	return (int)(ceil(pow((double)num, (double)1/3)));
}

/* Creates a randomized unique color map */
float** create_rand_3d_color_map(int num) {
	int p = compute_3d_color_interval_length(num);
	int maxcol = p*p*p;
	float** color_map = new float* [maxcol];
	int* shuffle = permute(maxcol);
	int l = 0;
	for (int i = 0; i < p; i++) {
		for (int j = 0; j < p; j++) {
			for (int k = 0; k < p; k++) {
				int s = shuffle[l];
				color_map[s] = new float[3];
				color_map[s][0] = (double)i/(double)p;
				color_map[s][1] = (double)j/(double)p;
				color_map[s][2] = (double)k/(double)p;
				l++;
			}
		}
	}
	delete shuffle;
	return color_map;
}

void delete_color_map(float** color_map, int num) {
	int p = compute_3d_color_interval_length(num);
	int maxcol = p*p*p;
	delete_pointer_array((void**)color_map, num);
}

void delete_pointer_array(void** ptrs, int num) {
	for (int i = 0; i < num; i++) {
		delete ptrs[i];
	}
	delete ptrs;
}

float vertex_list_summary::getSmoothingThreshold(float pct) {
	if (!aggregates) return -1;
	float thres = FMAX( abs(aggregates[SUMMARY_MAX_IDX][0]-aggregates[SUMMARY_MIN_IDX][0]), \
		FMAX(abs(aggregates[SUMMARY_MAX_IDX][1]-aggregates[SUMMARY_MIN_IDX][1]), \
		abs(aggregates[SUMMARY_MAX_IDX][2]-aggregates[SUMMARY_MIN_IDX][2])));
	return thres*pct;
}

void vertex_list_summary::print_vertex_sumary() {
	printf("nverts: %d\n", nverts);
	if (aggregates) {
		printf("Aggregates: \n");
		print_array(aggregates, 5, 3);
	}
	if (cov) {
		printf("Covariances: \n");
		print_array(cov, 3, 3);
	}
	if (eigen_vec) {
		printf("Eigen Vectors: \n");
		print_array(eigen_vec, 3, 3);
	}
	if (eigen_val) {
		printf("Eigen Values: \n");
		print_array(&(eigen_val), 1, 3);
	}
	if (bounding_box) {
		printf("Bounding Box: \n");
		print_array(bounding_box, 8, 3);
	}
}

