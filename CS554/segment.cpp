#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <algorithm>
#include <time.h>
#include "mathfuncs.h"
#include "nrutil.h"
#include "CornerPoints.h"
#include "utils.h"
#include "segment.h"

#define MAX_LEN 1.0
#define MAX_TEX_X 0.5
#define MAX_TEX_Y 0.5

typedef std::priority_queue< valpair*, std::vector<valpair*>, compareValpair > facequeue;

void clear_queue(facequeue &pq);

struct tr_queue {
	icVector3 median;
	facequeue pq;
	tr_queue(tris_tensor *t) {
		median = t->median;
	}
	void add(tris_tensor *t) {
		double dist = length(median - t->median);
		int sz = pq.size()+1;
		median = ((sz*median) + t->median)*(1/(sz+1));
		valpair *v = new valpair(t->t->index, dist);
		pq.push(v);
	}
	void reorder(tris_tensor_list *ttl) {
		facequeue tpq;
		tris_tensor **tt = ttl->tt;
		valpair *v = 0;
		for (; !pq.empty(); pq.pop()) {
			v = pq.top();
			v->val = length(median - tt[v->i]->median);
			tpq.push(v);
		}
		pq = tpq;
	}
	bool empty() {
		return pq.empty();
	}
	int top() {
		if (pq.empty()) return -1;
		valpair *v = pq.top(); pq.pop();
		int tidx = v->i;
		delete v;
		return tidx;
	}
	void clear() {
		clear_queue(pq);
	}
	~tr_queue() {clear();}
};

void compute_triangle_medians(tris_tensor_list *ttl);
int get_unsegmented_face(rndperm &r, segments *segs);
int segment_greedy(Polyhedron *poly, tris_tensor_list *ttl, segments *segs, CommandOptions *options);
int grow_patch(int c, CornerPoints *cp, tris_tensor_list *ttl, segments *segs, CommandOptions *options);
void add_unsegmented_neighbors_to_queue(int ref, CornerPoint &p, segments *segs, tris_tensor_list *ttl, tr_queue &pq, std::set<int> &queued);
bool check_validity(CornerPoint &p, segments *segs, transform_data &td, int segid, float maxlen);

void prepare_parameterizations(Polyhedron *poly, segments *segs) {
	icVector2 **params = new icVector2*[segs->length];
	for (int i = 0; i < segs->segs->size(); i++) {
		int center = segs->centers->at(i);
		std::set<int> *segids = segs->segs->at(i);
		std::set<int>::iterator segitr = segids->begin();
		double min_x=DBL_MAX, max_x=-DBL_MAX, min_y=DBL_MAX, max_y=-DBL_MAX;
		transform_data *td = segs->tds->at(i);
		for(; segitr != segids->end(); segitr++) {
			params[*segitr] = new icVector2[3];
			Triangle *t = poly->tlist[*segitr];
			//if (i < 1) printf("%4d:",i);
			for(int n = 0; n < 3; n++) {
				icVector3 tmp = td->transform(t->verts[n]);
				params[*segitr][n].x = FMAX(FMIN(tmp.x + MAX_TEX_X, 1),0);
				params[*segitr][n].y = FMAX(FMIN(tmp.y + MAX_TEX_Y, 1),0);
				max_x = DMAX(max_x, params[*segitr][n].x);
				min_x = DMIN(min_x, params[*segitr][n].x);
				max_y = DMAX(max_y, params[*segitr][n].y);
				min_y = DMIN(min_y, params[*segitr][n].y);
				//if (i < 1) printf("[%5.4f,%5.4f] ",params[*segitr][n].x,params[*segitr][n].y);
				//if (i < 1) printf("[%4.3f,%4.3f,%4.3f]",tmp.x,tmp.y,tmp.z);
			}
			//if (i < 1) printf("\n");
		}
		double diff_x = max_x - min_x;
		double diff_y = max_y - min_y;
		bool prt = i < 10;
		if (prt) printf("%4d:",i);
		segitr = segids->begin();
		for(; segitr != segids->end(); segitr++) {
			for(int n = 0; n < 3; n++) {
				if (abs(diff_x) > FLT_EPSILON) params[*segitr][n].x = (params[*segitr][n].x - min_x)/diff_x;
				if (abs(diff_y) > FLT_EPSILON) params[*segitr][n].y = (params[*segitr][n].y - min_y)/diff_y;
				if (prt) printf("[%5.4f,%5.4f] ",params[*segitr][n].x,params[*segitr][n].y);
			}
			if (prt) printf("\n");
		}
	}
	segs->params = params;
}

segments *segment_polyhedron(Polyhedron *poly, tris_tensor_list *ttl, CommandOptions *options) {
	segments *segs = new segments(poly->ntris);
	int ret = segment_greedy(poly, ttl, segs, options);
	printf("Segments Added: %d\n",segs->maxsegid+1);
	return segs;
}

// This is going to be a very greedy method and the logic will be:
// Start at random un segmented corner point and create a segment for it
// Add the bordering faces not yet processed to a priority queue.
// Priority is based on distance from the starting face.
// Also maintain a black-list for faces that should not be added to current segment.
int segment_greedy (
	Polyhedron *poly, tris_tensor_list *ttl, segments *segs, CommandOptions *options) {

	srand(10987); // so that results can be duplicated.

	CornerPoints *cp = new CornerPoints(poly);
	rndperm r(poly->ntris);

	compute_triangle_medians(ttl);
	int segsadded = 0;
	for (int c = get_unsegmented_face(r, segs); c > -1; \
		c = get_unsegmented_face(r, segs)) {
		segs->maxsegid++;
		segs->segidxs[c] = segs->maxsegid;
		segs->centers->push_back(c);
		std::set<int> *segtris = new std::set<int>();
		segtris->insert(c);
		segs->segs->push_back(segtris);
		segsadded++;
		int patchsize = grow_patch(c, cp, ttl, segs, options);
		printf("Growing patch at %d of size %d\n", c, patchsize);
	}

	return segsadded;

}

int grow_patch(int c, CornerPoints *cp, tris_tensor_list *ttl, segments *segs, CommandOptions *options) {
	std::set<int>* faceset = segs->segs->at(segs->maxsegid);
	tr_queue pq(ttl->tt[c]);
	transform_data *td = new transform_data(ttl->tt[c]);
	segs->tds->push_back(td);
	std::set<int> queued;
	int nfaces = 1; // includes the starting face
	CornerPoint &p = cp->corners[c*3]; // there are 3 corner points for every triangle
	add_unsegmented_neighbors_to_queue(p.t->index, p, segs, ttl, pq, queued);
	while (!pq.empty() && nfaces < options->segmentmaxfaces) {
		if (nfaces %10) pq.reorder(ttl);
		int tidx = pq.top();
		if (queued.count(tidx) > 0) queued.erase(queued.find(tidx));
		CornerPoint &np = cp->corners[tidx*3];
		bool valid = check_validity(np, segs, *td, segs->maxsegid, options->segmentmaxlen);
		if (valid) {
			segs->segidxs[tidx] = segs->maxsegid;
			faceset->insert(tidx);
			nfaces++;
			add_unsegmented_neighbors_to_queue(p.t->index, np, segs, ttl, pq, queued);
		}
	}
	return nfaces;
}

inline int get_segid(int tidx, int refidx, int* segidxs, int ref) {
	return tidx == refidx ? ref : segidxs[tidx];
}

inline bool valid_face_add(CornerPoint *p, int* segidxs, int ref) {
	int changes_from_ref_to_diff = 0;
	CornerPoint *c = p;
	int refidx = p->t->index;
	int prevsegid = -1;
	do {
		int currsegid = get_segid(c->t->index, refidx, segidxs, ref);
		if (prevsegid == ref && currsegid != ref) changes_from_ref_to_diff++;
		if (changes_from_ref_to_diff > 1) return false;
		prevsegid = currsegid;
		c = c->p->o->p;
	} while(c->index != p->index);
	return true;
}

bool check_boundary(transform_data &td, Vertex *v, float maxlen) {
	icVector3 nv = td.transform(v);
	//printf("nv: %5.4f,%5.4f,%5.4f\n",nv.x,nv.y,nv.z);
	if (nv.x > maxlen || nv.y > maxlen) {
		//printf("FALSE+ nv: %5.4f,%5.4f\n",nv.x,nv.y);
		return false;
	}
	if (nv.x < -maxlen || nv.y < -maxlen) {
		//printf("FALSE- nv: %5.4f,%5.4f\n",nv.x,nv.y);
		return false;
	}
	return true;
}

// check for self-intersections
bool check_validity(CornerPoint &p, segments *segs, transform_data &td, int segid, float maxlen) {
	int curridx = p.t->index;
	CornerPoint *c = &p;
	do {
		// TODO: check if within texture boundary
		// if the opposite face has the same seg id, then this vertex is projecting out
		// of the closed region. Therefore, check if this is going beyond the texture area
		if (dot(td.face->normal, c->t->normal) < 0) return false;
		if (segs->segidxs[c->o->t->index] == segid) {
			if (!check_boundary(td, c->v, maxlen)) return false;
		}
		if (!valid_face_add(c, segs->segidxs, segid)) return false;
		c = c->n;
	} while(c->index != p.index);
	return true;
}

void add_unsegmented_neighbors_to_queue(int ref, CornerPoint &p, segments *segs, tris_tensor_list *ttl, tr_queue &pq, std::set<int> &queued) {
	CornerPoint *c = &p;
	icVector3 &refT = ttl->tt[ref]->median;
	do {
		int tidx = c->o->t->index;
		if (segs->segidxs[tidx] == -1 && queued.count(tidx) == 0) {
			pq.add(ttl->tt[tidx]);
			queued.insert(tidx);
		}
		c = c->n;
	} while(c->index != p.index);
}

void compute_triangle_medians(tris_tensor_list *ttl) {
	int ntris = ttl->length;
	for (int i = 0; i < ntris; i++) {
		tris_tensor *tt = ttl->tt[i];
		Triangle *t = tt->t;
		icVector3 v1(t->verts[0]->x, t->verts[0]->y, t->verts[0]->z);
		icVector3 v2(t->verts[1]->x, t->verts[1]->y, t->verts[1]->z);
		icVector3 v3(t->verts[2]->x, t->verts[2]->y, t->verts[2]->z);
		tt->median = (v1+v2+v3)*(1/3);
	}
}

void clear_queue(facequeue &pq) {
	for (; !pq.empty(); pq.pop()) {
		valpair *v = pq.top();
		delete v;
	}
}

int get_unsegmented_face(rndperm &r, segments *segs) {
	int c;
	do {
		c = r.next();
	} while (c > -1 && segs->segidxs[c] != -1);
	return c;
}

void save_segments(segments *segs, std::string &file) {
	std::ofstream File(file.c_str());
	if (File) {
		File << segs->length << std::endl;
		File << segs->maxsegid << std::endl;
		for (int i = 0; i < segs->length; i++) {
			File << segs->segidxs[i] << std::endl;
		}
		std::vector<int>::iterator it;
		for (it = segs->centers->begin(); it != segs->centers->end(); it++) {
			File << *it << std::endl;
		}
		File.close();
	}
}

segments* load_segments(std::string &file, tris_tensor_list *ttl) {
	std::string line;
	segments *segs = 0;
	std::ifstream File(file.c_str());
	if (File) {
		int nline = 0;
		int idx = 0;
		while(getline(File,line)) {
			nline++;
			if (nline == 1) {
				int length = atoi(line.c_str());
				segs = new segments(length);
			} else if (nline == 2) {
				int maxsegid = atoi(line.c_str());
				segs->maxsegid = maxsegid;
				for (int i = 0; i < maxsegid+1; i++) {
					segs->segs->push_back(new std::set<int>());
				}
			} else if (nline <= segs->length+2) {
				int seg = atoi(line.c_str());
				segs->segidxs[idx] = seg;
				segs->segs->at(seg)->insert(idx);
				idx++;
			} else {
				int center = atoi(line.c_str());
				segs->centers->push_back(center);
				segs->tds->push_back(new transform_data(ttl->tt[center]));
			}
		}
		File.close();
	}
	return segs;
}
