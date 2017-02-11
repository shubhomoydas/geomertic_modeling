#include <stdlib.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include "mathfuncs.h"
#include "CornerPoints.h"
#include "utils.h"

#ifndef _CS554_SEGMENT_H_
#define _CS554_SEGMENT_H_

struct rndperm {
	int n;
	int p;
	int *list;
	rndperm(int _n): n(_n) {
		p = n;
		list = new int[n];
		for (int i = 0; i < n; i++) {
			list[i] = i;
		}
	}
	int next() {
		if (p == 0) return -1;
		int l = rand() % p;
		int t = list[l];
		list[l] = list[p-1];
		list[p-1] = t;
		p--;
		return t;
	}
	~rndperm() {
		if (list) delete list;
	}
};

struct valpair {
	int i;
	double val;
	valpair(): i(0), val(0) {}
	valpair(int _i, double _val): i(_i), val(_val) {}
};

class compareValpair {
public:
	// dist1 has highest prio than dist2 if dist1.dist is smaller than dist2.dist
	bool operator()(const valpair* dist1, const valpair* dist2) {
	   return (dist1->val > dist2->val);
	}
};

struct transform_data {
	Triangle *face;
	icVector3 face_center;
	//icVector3 new_center;
	icMatrix3x3 M;
	icMatrix3x3 iM;
	transform_data(tris_tensor *t) {
		face = t->t;
		face_center = t->p;
		//iM.entry[0][0] = t->minor.x; iM.entry[0][1] = t->major.x; iM.entry[0][2] = t->ref_z.x;
		//iM.entry[1][0] = t->minor.y; iM.entry[1][1] = t->major.y; iM.entry[1][2] = t->ref_z.y;
		//iM.entry[2][0] = t->minor.z; iM.entry[2][1] = t->major.z; iM.entry[2][2] = t->ref_z.z;
		iM.entry[0][0] = t->ref_x.x; iM.entry[0][1] = t->ref_y.x; iM.entry[0][2] = t->ref_z.x;
		iM.entry[1][0] = t->ref_x.y; iM.entry[1][1] = t->ref_y.y; iM.entry[1][2] = t->ref_z.y;
		iM.entry[2][0] = t->ref_x.z; iM.entry[2][1] = t->ref_y.z; iM.entry[2][2] = t->ref_z.z;
		M = inverse(iM);
		//new_center = transform(face_center);
		double err = dot(t->minor, t->major)+dot(t->major,t->ref_z)+dot(t->minor,t->ref_z);
		if (abs(err) > FLT_EPSILON) {
			printf("Incorrect orthogonalization! err = %f\n",err);
			exit(0);
		}

		// check transforms
		icVector3 vp(t->t->verts[0]->x,t->t->verts[0]->y,t->t->verts[0]->z);
		//icVector3 c_x = transform(v);
		icVector3 p = vp - t->p;
		icVector3 c_x = transform(p);
		if (abs(c_x.z) > DBL_EPSILON) {
			printf("Incorrect computation! (%5.4f,%5.4f,%5.4f)\n",c_x.x,c_x.y,c_x.z);
			exit(0);
		}

	}
	icVector3 transform(Vertex *v) {
		icVector3 vp(v->x,v->y,v->z);
		icVector3 r = transform(vp - face_center);
		return r;
		//return r - new_center;
	}
	icVector3 transform(icVector3 &a) {
		icVector3 r;
		r.x = M.entry[0][0]*a.x + M.entry[0][1]*a.y + M.entry[0][2]*a.z;
		r.y = M.entry[1][0]*a.x + M.entry[1][1]*a.y + M.entry[1][2]*a.z;
		r.z = M.entry[2][0]*a.x + M.entry[2][1]*a.y + M.entry[2][2]*a.z;
		return r;
	}
};

struct segments {
	int length;
	int *segidxs;
	std::vector< std::set<int>* > *segs;
	std::vector< int >* centers;
	std::vector< transform_data* > *tds;
	int maxsegid;
	float **colormap;
	icVector2 **params;
	segments(): length(0), segidxs(0), segs(0), centers(0), tds(0), maxsegid(-1), colormap(0), params(0) {}
	segments(int _length): length(_length) { //, segidxs(0), segids(0), segs(0), maxsegid(0){}
		segidxs = new int[length];
		maxsegid = -1;
		centers = new std::vector<int>();
		segs = new std::vector< std::set<int>* >();
		tds = new std::vector< transform_data* >();
		for (int i = 0; i < length; i++) {
			segidxs[i] = -1;
		}
	}
	void createColorMap() {
		colormap = create_rand_3d_color_map(maxsegid+1);
	}
	void deleteColorMap() {
		if (colormap) {
			delete_color_map(colormap, maxsegid+1);
		}
	}
	~segments() {
		if (centers) delete centers;
		if (segs) release_vector_elements< std::set<int> >(segs);
		if (tds) release_vector_elements< transform_data >(tds);
		delete segs;
		delete tds;
		if (segidxs) delete segidxs;
		deleteColorMap();
		if (params) {
			for (int i = 0; i < length; i++) {
				delete params[i];
			}
			delete params;
		}
	}
};

segments *segment_polyhedron(Polyhedron *poly, tris_tensor_list *ttl, CommandOptions *options);

void save_segments(segments *segs, std::string &file);
segments* load_segments(std::string &file, tris_tensor_list *ttl);
void prepare_parameterizations(Polyhedron *poly, segments *segs);

#endif
