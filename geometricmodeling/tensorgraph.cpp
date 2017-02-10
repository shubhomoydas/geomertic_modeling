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

void display_curv (
	GLenum mode, 
	Polyhedron *this_poly, 
	vertex_tensor_list *vtl, 
	tris_tensor_list *ttl,
	int display_mode, 
	CommandOptions* options
) {

	unsigned int i, j;
	GLfloat mat_diffuse[4];
	double cycle_color[3];

	double a = 1;
	//if (options->curvtensor) a = 0.5;
	//if (options->curvtype == 2) a = 0.0;

	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	for (i=0; i<this_poly->ntris && options->curvtype != 3; i++) {
		if (mode == GL_SELECT)
			glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];
		bool setpolycolor = true;

		switch (display_mode) {
		case 0:
			glBegin(GL_POLYGON);
			glEdgeFlag(TRUE);
			for (j=0; j<3; j++) {
				Vertex *v = temp_t->verts[j];
				glNormal3d(v->normal.entry[0], v->normal.entry[1], v->normal.entry[2]);
				bool setvertcolor = true;
				cycle_color[0] = 0;
				cycle_color[1] = 0;
				cycle_color[2] = 0;
				switch (options->curvtype) {
				case 0:
					color_mapping(1 - vtl->vt[v->index]->norm_gauss_curv, cycle_color);
					break;
				case 1:
					color_mapping(vtl->vt[v->index]->norm_mean_curv, cycle_color);
					break;
				case 2:
					cycle_color[0] = 1;
					cycle_color[1] = 1;
					cycle_color[2] = 1;
					break;
				default:
					color_mapping(1 - vtl->vt[v->index]->norm_mean_curv, cycle_color);
					break;
				}
				mat_diffuse[0] = cycle_color[0];
				mat_diffuse[1] = cycle_color[1];
				mat_diffuse[2] = cycle_color[2];
				mat_diffuse[3] = a;
				if (setvertcolor) glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
				glVertex3d(v->x, v->y, v->z);
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
	if (options->curvtensor) 
		draw_tensor_edges(vtl, options);
	else if (options->curvstream) 
		draw_streamlines(vtl, ttl,options);
}

void draw_streamlines(vertex_tensor_list *vtl, tris_tensor_list *ttl, CommandOptions* options) {
	GLfloat mat_diffuse[4];
	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 0.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;

	icVector3 lightpos(5.5,0.0,0.0);
	normalize(lightpos);
	bool penink = (options->curvaxis == 3);
	double *lighting = new double[ttl->length];
	for (int i = 0; i < ttl->length; i++) {
		icVector3 lnorm = ttl->tt[i]->p - lightpos;
		lighting[i] = -dot(lnorm, ttl->tt[i]->t->normal);
	}

	glDisable(GL_BLEND);
	//printf("Draw stream lines...\n");
	for (int i = 0; i < ttl->length; i++) {
		tris_tensor *tt = ttl->tt[i];
		Triangle *t = tt->t;

		double k1 = tt->kappas.x;
		double k2 = tt->kappas.y;

		double l1 = DMIN(vtl->edge_length*0.5,0.1), l2 = DMIN(vtl->edge_length*0.5,0.1);
		icVector3 minor = tt->minor*l1;
		icVector3 major = tt->major*l2;

		bool drawmajor = (options->curvaxis == 0 || options->curvaxis == 1);
		bool drawminor = (options->curvaxis == 0 || options->curvaxis == 2);

		if (options->curvaxis == 3) {
			if (lighting[i] > 0.95) {
				drawmajor = false;
				drawminor = false;
			} else if (lighting[i] > 0.4) {
				drawmajor = true;
				drawminor = false;
			} else {
				drawmajor = true;
				drawminor = true;
			}
		}

		glLineWidth(1.0);
		glEnable(GL_LINE_SMOOTH);
		glBegin(GL_LINES);
		
		icVector3 &p = tt->p;
		if (drawmajor) {
			mat_diffuse[0] = penink ? 0 : 1;
			mat_diffuse[1] = penink ? 0 : 0;
			mat_diffuse[2] = penink ? 0 : 0;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
			//glColor3d(mat_diffuse[0], mat_diffuse[1], mat_diffuse[2]);
			glVertex3f(p.x - major.x, p.y - major.y, p.z - major.z);
			glVertex3f(p.x + major.x, p.y + major.y, p.z + major.z);
		}

		if (drawminor) {
			mat_diffuse[0] = penink ? 0 : 0;
			mat_diffuse[1] = penink ? 0 : 1;
			mat_diffuse[2] = penink ? 0 : 0;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
			//glColor3d(mat_diffuse[0], mat_diffuse[1], mat_diffuse[2]);
			glVertex3f(p.x - minor.x, p.y - minor.y, p.z - minor.z);
			glVertex3f(p.x + minor.x, p.y + minor.y, p.z + minor.z);
		}

		glEnd();

	}

	delete lighting;

}

void draw_tensor_edges(vertex_tensor_list *vtl, CommandOptions* options) {
	GLfloat mat_diffuse[4];
	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 0.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;
	glDisable(GL_BLEND);
	for (int i = 0; i < vtl->length; i++) {

		vertex_tensor *vt = vtl->vt[i];
		Vertex *v = vt->v;

		double k1 = vt->kappas.x;
		double k2 = vt->kappas.y;

		//double k1 = DMIN(vtl->max_kappa, abs(vt->kappas.x))/vtl->max_kappa;
		//double k2 = DMIN(vtl->max_kappa, abs(vt->kappas.y))/vtl->max_kappa;

		double k = k1+k2;

		//if (k < FLT_EPSILON) continue;
		//if (DMAX(k1,k2) < 0.0) continue;

		//printf("Kappas: %4.3f, %4.3f\n",vt->kappas.x,vt->kappas.y);

		//k1 = vt->edge_length*1*k1/k; k2 = vt->edge_length*1*k2/k;
		//k1 = DMIN(vt->edge_length*1*k1,0.1); k2 = DMIN(vt->edge_length*1*k2,0.1);
		double l1 = DMIN(vtl->edge_length*0.5,0.1), l2 = DMIN(vtl->edge_length*0.5,0.1);
		//k1 = DMIN(vtl->edge_length*2*k1,0.1); k2 = DMIN(vtl->edge_length*2*k2,0.1);

		//icVector3 minor = vt->minor*k2; // k2 by intent - the more curvature, shorter length because it curves fast
		//icVector3 major = vt->major*k1; // k1 by intent - the less curvature, longer length because it is flat

		icVector3 minor = vt->minor*l1;
		icVector3 major = vt->major*l2;

		bool drawmajor = (options->curvaxis == 0 || options->curvaxis == 1);
		bool drawminor = (options->curvaxis == 0 || options->curvaxis == 2);

		if (options->curvaxis == 3) {
			if (abs(k2) < abs(k1)) drawminor = true;
			else drawmajor = true;
		} else if (options->curvaxis == 4) {
			if (abs(k2) > abs(k1)) drawminor = true;
			else drawmajor = true;
		}

		//if (options->curvaxis > 2) {
		//	if (abs(k1)-abs(k2) < FLT_EPSILON) {
		//		drawminor = true;
		//		drawmajor = true;
		//	}
		//}

		//if (abs(dot(major,minor)) > FLT_EPSILON) printf("WARN::Axis not perpendicular!!");

		glLineWidth(1.0);
		glEnable(GL_LINE_SMOOTH);
		glBegin(GL_LINES);
		
		if (drawmajor) {
			mat_diffuse[0] = 1; //k1 < k2 ? 1 : 0;
			mat_diffuse[1] = 0; //k1 < k2 ? 0 : 1;
			mat_diffuse[2] = 0;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
			//glColor3d(mat_diffuse[0], mat_diffuse[1], mat_diffuse[2]);
			glVertex3f(v->x - major.x, v->y - major.y, v->z - major.z);
			glVertex3f(v->x + major.x, v->y + major.y, v->z + major.z);
		}

		if (drawminor) {
			mat_diffuse[0] = 0; //k1 < k2 ? 0 : 1;
			mat_diffuse[1] = 1; //k1 < k2 ? 1 : 0;
			mat_diffuse[2] = 0;
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
			//glColor3d(mat_diffuse[0], mat_diffuse[1], mat_diffuse[2]);
			glVertex3f(v->x - minor.x, v->y - minor.y, v->z - minor.z);
			glVertex3f(v->x + minor.x, v->y + minor.y, v->z + minor.z);
		}

		glEnd();

	}
}
