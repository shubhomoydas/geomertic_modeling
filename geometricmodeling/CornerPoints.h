#pragma once

#include <stdio.h>
#include "learnply.h"

typedef struct _CornerPoint {
	_CornerPoint* n;
	_CornerPoint* p;
	_CornerPoint* o;
	Vertex* v;
	Triangle* t;
	Edge* e;
	int index;
	double angle; // in radians
	int valence; // number of connected vertices
	bool mincornerforvertex;
	//icMatrix3x3 mv; // translation matrix that converts from local to global co-ordinates
	_CornerPoint(): n(0), p(0), o(0), v(0), t(0), e(0), index(-1), \
		angle(0), valence(0), mincornerforvertex(false) {}
	void setAngle();
	int setValence();
} CornerPoint;

class CornerPoints
{
public:
	CornerPoint *corners;
	int length;
	int maxvalence;
	CornerPoints(Polyhedron *poly);
	void validate();
	void print_all_corners();
	void setValences();
	void printValences();
	~CornerPoints(void);
};

