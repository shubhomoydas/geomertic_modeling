
#include <stdio.h>
#include <math.h>
#include <map>
#include "icVector.H"
#include "icMatrix.H"
#include "tmatrix.h"
#include "learnply.h"
#include "CornerPoints.h"
#include "nrutil.h"
#include "stb_image_aug.h"

#ifndef __MATHFUNCS_H__

#define __MATHFUNCS_H__

namespace WeightScheme {
	enum WEIGHT_SCHEME {UNIF=0, CORD=1, MEAN_CURV=2, MEAN_VAL=3};
}

namespace IterScheme {
	enum ITER_SCHEME {INITIAL=0, RECOMPUTE=1};
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau); a[k][l]=h+s*(g-h*tau);

#define COTANGENT(x) tan((PI/2) - (x))

#define EPS 1.0e-14

enum{SUMMARY_MIN_IDX=3,SUMMARY_MAX_IDX=4};

void delete_pointer_array(void** ptrs, int num);

typedef struct _vert_weights {
	Vertex* v;
	int valence;
	int *vidxs;
	double *wts;
	_vert_weights(Vertex* _v, int _valence): v(_v), valence(_valence) {
		vidxs = new int[_valence];
		wts = new double[_valence];
		for (int i = 0; i < _valence; i++) vidxs[i] = INT_MAX;
	}
	_vert_weights(): v(0), valence(0), vidxs(0), wts(0) {}
	void normalize() {
		double tot = 0.0;
		for (int i = 0; i < valence; i++) {
			tot += wts[i];
		}
		for (int i = 0; i < valence; i++) {
			wts[i] = wts[i]/tot;
		}
	}
	void print() {
		printf("v(%3d):",v->index);
		for (int i = 0; i < valence; i++) {
			printf(" %3d:%5.4f",vidxs[i],wts[i]);
		}
	}
	~_vert_weights() {
		if (vidxs) delete vidxs;
		if (wts) delete wts;
	}
} vert_weights;

typedef struct _vert_weights_list {
	int length;
	vert_weights** wts;
	_vert_weights_list(): length(0), wts(0){}
	_vert_weights_list(int _length): length(_length) {
		wts = new vert_weights*[_length];
	}
	int total_valence() {
		int tot = 0;
		for (int i=0; i < length; i++) {
			tot += wts[i]->valence;
		}
		return tot;
	}
	void print() {
		for (int i = 0; i < length; i++) {
			wts[i]->print(); printf("\n");
		}
	}
	~_vert_weights_list() {
		if (!wts) return;
		for (int i = 0; i < length; i++) {
			delete wts[i];
		}
	}
} vert_weights_list;

typedef struct _vertex_list_summary {
	float** aggregates;
	float** cov;
	float** eigen_vec;
	float* eigen_val;
	float** bounding_box;
	int** box_faces;
	int nverts;
	_vertex_list_summary(): \
		aggregates(0), cov(0), eigen_vec(0), eigen_val(0), bounding_box(0), box_faces(0), nverts(0) {}
	_vertex_list_summary(float** _aggregates, float** _cov, int _nverts): \
		aggregates(_aggregates), cov(_cov), \
		eigen_vec(0), eigen_val(0), bounding_box(0), box_faces(0), nverts(_nverts) {}
	~_vertex_list_summary() {
		if (aggregates) {
			delete_pointer_array((void**)aggregates, 5);
			aggregates = 0;
		}
		if (cov) {
			delete_pointer_array((void**)cov,3);
			cov = 0;
		}
		if (eigen_vec) {
			delete_pointer_array((void**)eigen_vec,3);
			eigen_vec = 0;
		}
		if (eigen_val) {
			delete eigen_val;
			eigen_val = 0;
		}
		if (bounding_box) {
			delete_pointer_array((void**)bounding_box,8);
			bounding_box = 0;
		}
		if (box_faces) {
			delete_pointer_array((void**)box_faces,8);
			box_faces = 0;
		}
	}
	void replace_eigen(float** _eigen_vec, float* _eigen_val) {
		if (eigen_vec) {
			delete_pointer_array((void**)eigen_vec,3);
			eigen_vec = 0;
		}
		if (eigen_val) {
			delete eigen_val;
			eigen_val = 0;
		}
		eigen_vec = _eigen_vec;
		eigen_val = _eigen_val;
	}
	float getSmoothingThreshold(float pct);
	void print_vertex_sumary();
} vertex_list_summary;

/* Code mostly borrowed from Numerical Recipies in C by William H. Press.
 * 
 * The internal represenation of this data structure uses array indexing
 * relative 1, NOT 0.
 */
typedef struct _sparse {
	// IMPORTANT: sa, ija user indexing 1...nmax
	double* sa;
	unsigned long* ija;
	unsigned long n;
	unsigned long nmax;
	_sparse(): sa(0), ija(0), n(0), nmax(0) {}
	_sparse(unsigned long _n, unsigned long _nmax): n(_n), nmax(_nmax) {
		//ija = new unsigned long[_nmax];
		//sa = new double[_nmax];
		ija = lvector(1, _nmax);
		sa = dvector(1, _nmax);
	}

	/*
	 * Converts a square matrix a[0..n-1][0..n-1] into row-indexed sparse storage mode. Only elements
	 * of a with magnitude =thresh are retained. Output is in two linear arrays with dimension
	 * nmax (an input parameter): sa[1..] contains array values, indexed by ija[1..]. The
	 * number of elements filled of sa and ija on output are both ija[ija[1]-1]-1 (see text).
	 */
	_sparse(double **a, unsigned long _n, double thresh) {
		init(a, _n, thresh);
	}

	_sparse(double **a, unsigned long _n) {
		init(a, _n, EPS);
	}

	void init(double **a, unsigned long _n, double thresh) {
		n = _n;
		nmax = n*n + 1;
		//printf("Init sparse: n=%d, nmax=%d\n",n,nmax);
		ija = lvector(1, nmax);
		sa = dvector(1, nmax);
		sprsin(a, thresh);
	}

	~_sparse() {
		//if (sa) delete sa;
		//if (ija) delete ija;
		if (sa) free_dvector(sa, 1, nmax);
		if (ija) free_lvector(ija, 1, nmax);
	}

	/*
	 * b[0..n-1] is the column matrix in A.x = b.
	 * x[0..n-1] is the matrix that will contain the solution as output. This should be
	 *   initialized to some random values or zeros.
	 */
	void linear(double b[], double x[]) {
		int iter;
		int itol = 3;
		double err;
		// convert to a ralative-1 indexing.
		linbcg(b-1, x-1, itol, EPS, 1000, &iter, &err);
		printf("Solved linear system in %d iters, %8.7f err.\n",iter, err);
	}

	// dense() converts the sparse representation to a dense representation
	// a[0..n-1][0..n-1] is output matrix where the dense matrix should be stored
	void dense(double** a) {
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a[i][j] = 0.0;
		for (int i = 0; i < n; i++) a[i][i] = sa[i+1];
		for (int i = 0; i < n; i++) {
			for (int k = ija[i+1]; k < ija[i+2]; k++) {
				a[i][ija[k]-1] = sa[k];
			}
		}
	}

	void print() {
		printf("Sparse n = %d, tot = %d, nmax = %d\n", n, ija[n+1], nmax);
		for (int i = 1; i < ija[n+1]; i++) {
			printf("%5d %5d %4.3f\n",i,ija[i],sa[i]);
		}
		//for (int i = 1; i < ija[n+1]; i++) {
		//	printf("%5d ",i);
		//}
		//printf("\nija:\n");
		//for (int i = 1; i < ija[n+1]; i++) {
		//	printf("%5d ",ija[i]);
		//}
		//printf("\nsa:\n");
		//for (int i = 1; i < ija[n+1]; i++) {
		//	printf("%4.3f ",sa[i]);
		//}
		//printf("\n");
	}

private:
	/*
	 * Solves A.x = b for x[1..n], given b[1..n], by the iterative biconjugate gradient method.
	 * On input x[1..n] should be set to an initial guess of the solution (or all zeros); itol is 1,2,3,
	 * or 4, specifying which convergence test is applied (see text); itmax is the maximum number
	 * of allowed iterations; and tol is the desired convergence tolerance. On output, x[1..n] is
	 * reset to the improved solution, iter is the number of iterations actually taken, and err is the
	 * estimated error. The matrix A is referenced only through the user-supplied routines atimes,
	 * which computes the product of either A or its transpose on a vector; and asolve, which solves
	 * A~.x = b or A^T.x = b for some preconditioner matrix A~ (possibly the trivial diagonal part of A).
	 */
	void linbcg(double b[], double x[], int itol, double tol,
		int itmax, int *iter, double *err) {
		unsigned long j;
		double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
		double *p,*pp,*r,*rr,*z,*zz; // Double precision is a good idea in this routine.

		p=dvector(1,n);
		pp=dvector(1,n);
		r=dvector(1,n);
		rr=dvector(1,n);
		z=dvector(1,n);
		zz=dvector(1,n);

		// Calculate initial residual.
		*iter=0;
		atimes(n,x,r,0); // Input to atimes is x[1..n], output is r[1..n];
		// the final 0 indicates that the matrix (not its transpose) is to be used.
		for (j=1;j<=n;j++) {
			r[j]=b[j]-r[j];
			rr[j]=r[j];
		}
		/* atimes(n,r,rr,0); */ // Uncomment this line to get the “minimum residual”
		if (itol == 1) { // variant of the algorithm.
			bnrm=snrm(n,b,itol);
			asolve(n,r,z,0); // Input to asolve is r[1..n], output is z[1..n];
			// the final 0 indicates that the matrix A (not its transpose) is to be used.
		}
		else if (itol == 2) {
			asolve(n,b,z,0);
			bnrm=snrm(n,z,itol);
			asolve(n,r,z,0);
		}
		else if (itol == 3 || itol == 4) {
			asolve(n,b,z,0);
			bnrm=snrm(n,z,itol);
			asolve(n,r,z,0);
			znrm=snrm(n,z,itol);
		} else nrerror("illegal itol in linbcg");
		while (*iter <= itmax) { // Main loop.
			++(*iter);
			asolve(n,rr,zz,1); // Final 1 indicates use of transpose matrix AT.
			for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
			// Calculate coefficient bk and direction vectors p and pp.
			if (*iter == 1) {
				for (j=1;j<=n;j++) {
					p[j]=z[j];
					pp[j]=zz[j];
				}
			}
			else {
				bk=bknum/bkden;
				for (j=1;j<=n;j++) {
					p[j]=bk*p[j]+z[j];
					pp[j]=bk*pp[j]+zz[j];
				}
			}
			bkden=bknum; // Calculate coefficient ak, newi terate x, and new
			atimes(n,p,z,0); // residuals r and rr.

			for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
			ak=bknum/akden;
			atimes(n,pp,zz,1);
			for (j=1;j<=n;j++) {
				x[j] += ak*p[j];
				r[j] -= ak*z[j];
				rr[j] -= ak*zz[j];
			}
			asolve(n,r,z,0); // Solve A . z = r and check stopping criterion.
			if (itol == 1)
				*err=snrm(n,r,itol)/bnrm;
			else if (itol == 2)
				*err=snrm(n,z,itol)/bnrm;

			else if (itol == 3 || itol == 4) {
				zm1nrm=znrm;
				znrm=snrm(n,z,itol);
				if (fabs(zm1nrm-znrm) > EPS*znrm) {
					dxnrm=fabs(ak)*snrm(n,p,itol);
					*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
				} else {
					*err=znrm/bnrm; // Error may not be accurate, so loop again.
					continue;
				}
				xnrm=snrm(n,x,itol);
				if (*err <= 0.5*xnrm) *err /= xnrm;
				else {
					*err=znrm/bnrm; // Error may not be accurate, so loop again.
					continue;
				}
			}
			//printf("iter=%4d err=%12.6f\n",*iter,*err);
			if (*err <= tol) break;
		}
		free_dvector(p,1,n);
		free_dvector(pp,1,n);
		free_dvector(r,1,n);
		free_dvector(rr,1,n);
		free_dvector(z,1,n);
		free_dvector(zz,1,n);
	}

	/*
	 * Converts a square matrix a[0..n-1][0..n-1] into row-indexed sparse storage mode. Only elements
	 * of a with magnitude =thresh are retained. Output is in two linear arrays with dimension
	 * nmax (an input parameter): sa[1..] contains array values, indexed by ija[1..]. The
	 * number of elements filled of sa and ija on output are both ija[ija[1]-1]-1 (see text).
	 */
	void sprsin(double **a, double thresh) {
		int i,j;
		unsigned long k;
		for (j=1;j<=n;j++) sa[j]=a[j-1][j-1]; // Store diagonal elements.
		ija[1]=n+2; // Index to 1st rowoff- diagonal element, if any.
		k=n+1;
		for (i=1;i<=n;i++) { // Loop over rows.
			for (j=1;j<=n;j++) { // Loop over columns.
				if (fabs(a[i-1][j-1]) >= thresh && i != j) {
					if (++k > nmax) nrerror("sprsin: nmax too small");
					sa[k]=a[i-1][j-1]; // Store off-diagonal elements and their columns.
					ija[k]=j;
				}
			}
			ija[i+1]=k+1; // As each row is completed, store index to next.
		}
		sa[n+1] = 0.0; // an arbitrary value
	}

	/* 
	 * Multiply a matrix in row-index sparse storage arrays sa and ija by a vector x[1..n], 
	 * giving a vector b[1..n].
	 */
	void sprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n) {
		unsigned long i,k;
		if (ija[1] != n+2) nrerror("sprsax: mismatched vector and matrix");
		for (i=1;i<=n;i++) {
			b[i]=sa[i]*x[i];                 // Start with diagonal term.
			for (k=ija[i];k<=ija[i+1]-1;k++) // Loop over off-diagonal terms.
				b[i] += sa[k]*x[ija[k]];
		}
	}

	/* 
	 * Multiply the transpose of a matrix in row-index sparse storage arrays sa and ija by a vector
	 * x[1..n], giving a vector b[1..n].
	 */
	void sprstx(double sa[], unsigned long ija[], double x[], double b[],
		unsigned long n)
	{
		unsigned long i,j,k;
		if (ija[1] != n+2) nrerror("mismatched vector and matrix in sprstx");
		for (i=1;i<=n;i++) b[i]=sa[i]*x[i]; // Start with diagonal terms.
		for (i=1;i<=n;i++) {                // Loop over off-diagonal terms.
			for (k=ija[i];k<=ija[i+1]-1;k++) {
				j=ija[k];
				b[j] += sa[k]*x[i];
			}
		}
	}

	/*
	 * Compute one of two norms for a vector sx[1..n], as signaled by itol. Used by linbcg.
	 */
	double snrm(unsigned long n, double sx[], int itol)
	{
		unsigned long i,isamax;
		double ans;
		if (itol <= 3) {
			ans = 0.0;
			for (i=1;i<=n;i++) ans += sx[i]*sx[i]; // Vector magnitude norm.
			return sqrt(ans);
		} else {
			isamax=1;
			for (i=1;i<=n;i++) {                   // Largest component norm.
				if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
			}
			return fabs(sx[isamax]);
		}
	}

	void atimes(unsigned long n, double x[], double r[], int itrnsp)
	{
		// These are double versions of sprsax and sprstx.
		if (itrnsp) sprstx(sa,ija,x,r,n);
		else sprsax(sa,ija,x,r,n);
	}

	void asolve(unsigned long n, double b[], double x[], int itrnsp)
	{
		unsigned long i;
		for(i=1;i<=n;i++) x[i]=(sa[i] != 0.0 ? b[i]/sa[i] : b[i]);
		// The matrix A is the diagonal part of A, stored in the first n elements of sa. Since the
		// transpose matrix has the same diagonal, the flag itrnsp is not used.
	}

} sparse;

typedef struct _linear_system {
	sparse *a;
	double *x;
	double *b;
	_linear_system(): a(0), x(0), b(0) {}
	_linear_system(sparse* _a): a(_a) {
		if (!_a) { x = 0; b = 0; return; }
		x = new double[_a->n];
		b = new double[_a->n];
		for (int i = 0; i < _a->n; i++) x[i] = 0.0;
	}
	void solve() {
		if (!a) return;
		a->linear(b,x);
	}
	void get_soln_with_predef(std::map<int,double>* predefs, double *ret) {
		int len = a->n+predefs->size();
		for (int i = 0, j = 0; i < len && j < a->n; i++) {
			std::map<int,double>::iterator it = predefs->find(i);
			if (it == predefs->end()) {
				ret[i] = x[j];
				j++;
			} else {
				ret[i] = it->second;
			}
		}
	}
	~_linear_system() {
		if (a) delete a;
		if (x) delete x;
		if (b) delete b;
	}
} linear_system;

typedef struct _vertex_tensor {

	Vertex *v;
	int valence;
	double gauss_curv;
	double norm_gauss_curv;
	double mean_curv;
	double norm_mean_curv;
	double edge_length; // length of one reference edge
	icMatrix3x3 M; // converts from local to global
	icMatrix3x3 iM; // inverse M
	icMatrix2x2 T; // Global coordinate Tensor
	icMatrix2x2 locT; // Local coordinate Tensor

	icMatrix2x2 curvAxis;
	icVector2 kappas;

	icVector3 ref_x, ref_y, ref_z;

	icVector3 minor, major;

	_vertex_tensor(): v(0), valence(0), \
		gauss_curv(0), norm_gauss_curv(0), \
		mean_curv(0), norm_mean_curv(0), edge_length(0) {}
	_vertex_tensor(Vertex *_v, int _valence): v(_v), valence(_valence), \
		gauss_curv(0), norm_gauss_curv(0), mean_curv(0), norm_mean_curv(0), edge_length(0) {}

	void computeGlobalTensor(icVector3 *edges);
	void _vertex_tensor::computeLocalTensor();

	//~_vertex_tensor();

} vertex_tensor;

typedef struct _vertex_tensor_list {
	int length;
	vertex_tensor **vt;
	double edge_length;
	double max_kappa;
	double min_gauss_curv, max_gauss_curv;
	_vertex_tensor_list(): length(0), vt(0), edge_length(0), max_kappa(0.0), \
		min_gauss_curv(0.0), max_gauss_curv(0.0) {}
	_vertex_tensor_list(int _length): length(_length), edge_length(0), max_kappa(0.0), \
		min_gauss_curv(0.0), max_gauss_curv(0.0) {
		vt = new vertex_tensor*[length];
	}
	~_vertex_tensor_list() {
		if (vt) {
			for (int i = 0; i < length; i++) {
				delete vt[i];
			}
			delete vt;
		}
	}
} vertex_tensor_list;

typedef struct _tris_tensor : vertex_tensor {
	Triangle *t;
	icVector3 p;
	icVector3 median;
	_tris_tensor(): t(0) {v = 0;}
} tris_tensor;

typedef struct _tris_tensor_list {
	int length;
	tris_tensor **tt;
	_tris_tensor_list(): length(0), tt(0) {}
	_tris_tensor_list(int _length) {
		length = _length;
		tt = new tris_tensor*[length];
	}
	~_tris_tensor_list() {
		if (tt) {
			for (int i = 0; i < length; i++) {
				delete tt[i];
			}
			delete tt;
		}
	}
} tris_tensor_list;

void jacobi(float **a, int n, float d[], float **v, int *nrot);

int* permute(int n);
float** compute_vertex_aggregates(Vertex** vlist, int nverts);
float** compute_vertex_norm_aggregates(Vertex** vlist, int nverts);
float** compute_vertex_covariance(Vertex** vlist, int nverts);
vertex_list_summary* compute_vertex_list_summary(Vertex** vlist, int nverts, bool withnorm);
void compute_bounding_box(Vertex** vlist, vertex_list_summary* summary);

linear_system* create_sparse_laplace_system(vert_weights_list* vwts, std::map<int,double>* predefs);
void set_diffusion_colors(Polyhedron *poly, double *diffusion_colors, int *saddleverts, std::map<int,double> *predefs, int wtscheme);

vert_weights_list* get_vertex_weights(Polyhedron* poly, CornerPoints* cp, int scheme);
void smooth_polynomial(Polyhedron* poly, int wtscheme, int iterscheme, int iters, double lambda);
void smooth_polynomial(Polyhedron* poly, int wtscheme, int iters, double lambda);

void sort_vals(float* vals, int* dest, int num);

// solves for x in Ax = b and stores result in x
void solve_linear(icMatrix3x3 &A, icVector3 &b, double *x);

float** create_rand_3d_color_map(int num);
void delete_color_map(float** color_map, int num);

void print_vector(icVector3 &m);
void print_3x3_matrix(icMatrix3x3 &m);
void print_2x2_matrix(icMatrix2x2 &m);

vertex_tensor_list *compute_curvature_tensors(Polyhedron *poly, int wtscheme, int iters, double lambda);
tris_tensor_list *compute_tris_curvature_tensors(Polyhedron *poly, vertex_tensor_list *vtl);
double get_total_angle_deficit(vertex_tensor_list *vtl);

void compute_eigen_symmetric_2x2(icMatrix2x2 &m, icMatrix2x2 &evec, icVector2 &eval);

typedef struct _texture_map {
	int n;
	std::map<int,double> *f_predefs;
	std::map<int,double> *g_predefs;
	double *f, *g;
	int *f_saddles, *g_saddles;
	unsigned char *texture;
	int x, y, channels;
	_texture_map(): n(0), f_predefs(0), g_predefs(0), f(0), g(0), f_saddles(0), g_saddles(0), \
		texture(0), x(0), y(0), channels(0) {}
	_texture_map(int _n, std::string &texturefile): n(_n) {
		f = new double[_n];
		g = new double[_n];
		f_saddles = new int[_n];
		g_saddles = new int[_n];
		f_predefs = new std::map<int,double>();
		g_predefs = new std::map<int,double>();
		//texture = stbi_load(texturefile.c_str(), &x, &y, &channels, 1); // grey scale only for now...
		if (!texture) {
			printf("Error: Could not load texture file %s\n", texturefile.c_str());
		} else {
			printf ("Texture Loaded: x: %d, y: %d, channels: %d\n", x, y, channels);
		}
	}
	void set_diffusion_colos(Polyhedron *poly, int wtscheme) {
		set_diffusion_colors(poly, f, f_saddles, f_predefs, wtscheme);
		set_diffusion_colors(poly, g, g_saddles, g_predefs, wtscheme);
	}
	void get_color(int vidx, double *color) {
		int row = FMIN(y-1, floor(f[vidx]*y));
		int col = FMIN(x-1, floor(g[vidx]*x));
		int coloridx = row*x + col;
		if (!texture) return;
		//if (vidx < 20) {
		//	printf("vidx: %d, row = %d, col = %d, coloridx = %d\n", vidx, row, col, coloridx);
		//}
		color[0] = color[1] = color[2] = ((double)texture[coloridx]) / 256.0;
	}
	void print() {
		for (int i = 0; i < n; i++) {
			printf("%d\t%5.4f\t%5.4f\n",i,f[i],g[i]);
		}
	}
	~_texture_map() {
		if (f) delete f;
		if (g) delete f;
		if (f_saddles) delete f_saddles;
		if (g_saddles) delete g_saddles;
		if (f_predefs) delete f_predefs;
		if (g_predefs) delete g_predefs;
		if (texture) stbi_image_free(texture);
	}
} texture_map;

#endif
