#include <vector>
#include <algorithm>
#include <stdlib.h>
#include "stdafx.h"
#include "CornerPoints.h"

void CornerPoint::setAngle() {
	// The triangle vertices are ABC where this vertex is A.
	// B -> next vertex, C -> prev vertex
	// opposite edge lengths are a,b,c resp.
	// angle A = arccos((b^2 + c^2 - a^2) / (a*b*c))
	double a = e->length;
	double b = p->e->length;
	double c = n->e->length;
	if (a < DBL_EPSILON || b < DBL_EPSILON || c < DBL_EPSILON) {
		std::printf("WARN:: Corner %4d has length zero!\n", index);
		angle = 0;
		return;
	}
	angle = acos((b*b + c*c - a*a)/(2*b*c));
	//printf("%4d: acos((%5.4f + %5.4f - %5.4f)/%5.4f) = %5.4f\n",index,b*b,c*c,a*a,2*b*c,angle);
}

int CornerPoint::setValence() {
	CornerPoint* o = this;
	int minindex = index;
	mincornerforvertex = false;
	valence = 0;
	do {
		minindex = IMIN(minindex,o->index);
		valence++;
		o = o->p->o->p;
	} while (index != o->index);
	if (index == minindex) mincornerforvertex = true;
	return valence;
}

bool corner_comp(CornerPoint* c1, CornerPoint* c2) {
	int l1 = IMIN(c1->p->v->index, c1->n->v->index);
	int r1 = IMAX(c1->p->v->index, c1->n->v->index);
	int l2 = IMIN(c2->p->v->index, c2->n->v->index);
	int r2 = IMAX(c2->p->v->index, c2->n->v->index);
	return l1 < l2 || (l1 == l2 && r1 < r2);
}

Edge* getEdge(Triangle* t, Vertex* v) {
	Edge* e = t->edges[0];
	if (t->edges[0]->verts[0]->index == v->index || t->edges[0]->verts[1]->index == v->index) {
		if (t->edges[1]->verts[0]->index == v->index || t->edges[1]->verts[1]->index == v->index) {
			e = t->edges[2];
		} else {
			e = t->edges[1];
		}
	}
	return e;
}

CornerPoints::CornerPoints(Polyhedron *poly)
{
	length = 3 * poly->ntris;
	corners = new CornerPoint[length];
	int ntris = poly->ntris;
	std::vector<CornerPoint*> cv;
	for (int i = 0; i < ntris; i++) {
		corners[i*3 + 0].v = poly->tlist[i]->verts[0];
		corners[i*3 + 1].v = poly->tlist[i]->verts[1];
		corners[i*3 + 2].v = poly->tlist[i]->verts[2];
		corners[i*3 + 0].t = poly->tlist[i];
		corners[i*3 + 1].t = poly->tlist[i];
		corners[i*3 + 2].t = poly->tlist[i];
		corners[i*3 + 0].p = &corners[i*3 + 2];
		corners[i*3 + 1].p = &corners[i*3 + 0];
		corners[i*3 + 2].p = &corners[i*3 + 1];
		corners[i*3 + 0].n = &corners[i*3 + 1];
		corners[i*3 + 1].n = &corners[i*3 + 2];
		corners[i*3 + 2].n = &corners[i*3 + 0];
		corners[i*3 + 0].index = i*3 + 0;
		corners[i*3 + 1].index = i*3 + 1;
		corners[i*3 + 2].index = i*3 + 2;
		cv.push_back(&corners[i*3 + 0]);
		cv.push_back(&corners[i*3 + 1]);
		cv.push_back(&corners[i*3 + 2]);
		// c.e = edge that does not contain the vertex c.v
		corners[i*3 + 0].e = getEdge(poly->tlist[i], corners[i*3 + 0].v);
		corners[i*3 + 1].e = getEdge(poly->tlist[i], corners[i*3 + 1].v);
		corners[i*3 + 2].e = getEdge(poly->tlist[i], corners[i*3 + 2].v);
	}
	std::sort(cv.begin(), cv.end(), corner_comp);
	int l, r, pl = -1, pr = -1;
	CornerPoint* prev = 0;
	std::vector<CornerPoint*>::iterator it;
	for (it = cv.begin(); it != cv.end(); it++) {
		l = IMIN((*it)->p->v->index, (*it)->n->v->index);
		r = IMAX((*it)->p->v->index, (*it)->n->v->index);
		//printf("Corner Index %4d, %4d\n", l, r);
		if (pl == l && pr == r) {
			prev->o = *it;
			(*it)->o = prev;
			//printf(">> Corner Index %4d, %4d\n", l, r);
		}
		pl = l;
		pr = r;
		prev = *it;
	}
	// set all angles
	for (int i = 0; i < length; i++) {
		corners[i].setAngle();
	}
}

void CornerPoints::setValences() {
	maxvalence = -1;
	for (int i = 0; i < length; i++) {
		int val = corners[i].setValence();
		maxvalence = IMAX(val, maxvalence);
	}
}

void CornerPoints::printValences() {
	for (int i = 0; i < length; i++) {
		printf("Corner %d: %d %s\n", corners[i].index,corners[i].valence, \
			corners[i].mincornerforvertex?"->min":"");
	}
	printf("Max Valence: %d\n",maxvalence);
}

void CornerPoints::validate() {
	for (int i = 0; i < length; i++) {
		if (!(corners[i].p->n->v == corners[i].n->p->v)) {
			printf("Error!! Failed assertion p.n == n.p\n");
			break;
		}
		if (!(corners[i].n->n->n->v == corners[i].v)) {
			printf("Error!! Failed assertion n.n.n == this\n");
			break;
		}
		if (!(corners[i].p->p->p->v == corners[i].v)) {
			printf("Error!! Failed assertion p.p.p == this\n");
			break;
		}
		if (corners[i].o) {
			if (!(corners[i].o->n->v == corners[i].p->v)) {
				printf("Error!! Failed assertion o.n == p\n");
				//if (!corners[i].o) printf("o is NULL\n");
				//printf("o.n: %4d, p: %4d\n", corners[i].o->v->index, corners[i].p->v->index);
				break;
			}
			if (!(corners[i].o->p->v == corners[i].n->v)) {
				printf("Error!! Failed assertion o.p == n\n");
				break;
			}
			if (corners[i].t->index == corners[i].o->t->index) {
				printf("Error!! Failed assertion opposite faces must be different\n");
				break;
			}
		}
		if (corners[i].p->o) {
			if (!(corners[i].p->o->p->v == corners[i].v)) {
				printf("Error!! Failed assertion p.o.p == this\n");
				break;
			}
			if (corners[i].p->o->p->p->o) {
				if (!(corners[i].p->o->p->p->o->p->v == corners[i].v)) {
					printf("Error!! Failed assertion p.o.p == this\n");
					break;
				}
			}
		}
		if (!(corners[i].n->v == corners[i].e->verts[0] || corners[i].n->v == corners[i].e->verts[1])) {
			printf("Error!! Failed assertion n.v{%4d} == e.v[0]{%4d}\n", corners[i].n->v->index, corners[i].e->verts[1]->index);
			break;
		}
	}
}


void CornerPoints::print_all_corners() {
	for (int i = 0; i < length; i++) {
		printf("Corner #%3d Face: %3d Vertex: %3d Opp: %3d OppFace: %3d, angle: %5.4f, e.len: %5.4f\n", \
			corners[i].index, corners[i].t->index, corners[i].v->index, \
			corners[i].o->index, corners[i].o->t->index, corners[i].angle, corners[i].e->length);
	}
}

CornerPoints::~CornerPoints(void)
{
	//for (int i = 0; i < length; i++) {
	//	if (corners[i].mv) delete corners[i].mv;
	//}
	delete corners;
}
