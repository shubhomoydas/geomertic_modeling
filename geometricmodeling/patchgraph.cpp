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

void display_patches (
	GLenum mode, 
	Polyhedron *this_poly, 
	vertex_tensor_list *vtl, 
	tris_tensor_list *ttl,
	segments *segs,
	int display_mode, 
	GLuint texture,
	CommandOptions* options
) {

	unsigned int i, j;
	GLfloat mat_diffuse[4];
	double cycle_color[3];

	double a = 1;

	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	std::set<int> centers;
	for (int i = 0; i < segs->centers->size(); i++) {
		centers.insert(segs->centers->at(i));
	}

	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
			glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];
		int seg = segs->segidxs[i];
		icVector2 *param = 0; 
		if(options->segmenttexture) param = segs->params[temp_t->index];

		bool color_vertex = true; //(centers.count(temp_t->index) > 0);

		switch (display_mode) {
		case 0:
			glBegin(GL_POLYGON);
			if(color_vertex) {
				mat_diffuse[0] = segs->colormap[seg][0];
				mat_diffuse[1] = segs->colormap[seg][1];
				mat_diffuse[2] = segs->colormap[seg][2];
			} else {
				mat_diffuse[0] = 1;
				mat_diffuse[1] = 1;
				mat_diffuse[2] = 1;
			}
			mat_diffuse[3] = 1.0;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
			glEdgeFlag(TRUE);
			for (j=0; j<3; j++) {
				Vertex *v = temp_t->verts[j];
				glNormal3d(v->normal.entry[0], v->normal.entry[1], v->normal.entry[2]);
				glVertex3d(v->x, v->y, v->z);
				if(options->segmenttexture) glTexCoord2d(param[j].x,param[j].y);
			}
			glEnd();
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				Vertex *v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(v->x, v->y, v->z);
			}
			glEnd();
			break;

		case 10:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				mat_diffuse[0] = 1.0;
				mat_diffuse[1] = 0.0;
				mat_diffuse[2] = 0.0;
				mat_diffuse[3] = 1.0;

				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);

				glColor3f(1.0, 0.0, 0.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;
		}
	}
	if(options->segmenttexture) glDisable(GL_TEXTURE_2D);
}

