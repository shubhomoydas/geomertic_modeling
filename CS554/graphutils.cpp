#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include "glut.h"
#include "mathfuncs.h"
#include "utils.h"
#include "CornerPoints.h"
#include "stb_image_aug.h"
#include "SOIL.h"

void line_from_and_direction(float x1, float y1, float z1, float a, float b, float c, float scale, GLfloat* rgb);
void draw_bounding_box(float** bounding_box, int** faces);
void drawSilhouetteEdges(CornerPoints *cp, icVector3& r);
void drawEdges(std::vector<Edge*>* sil);
void drawSilhouetteEdgesLevelset(CornerPoints *cp, icVector3& r);
void drawLevelsetEdges(std::vector<Edge*>* sil);

Vertex* new_middle_vertex(Vertex* a, Vertex* b, Vertex* c, Vertex* d);
Vertex* new_middle_vertex(CornerPoint* q, std::vector<Vertex*>& new_verts, std::map<int, Vertex*>& corner_vert);
CornerPoint* corner_opp_largest_edge(CornerPoint* c);
Triangle* new_triangle(Vertex* v1, Vertex* v2, Vertex* v3);
void update_polyhedron(Polyhedron* poly, std::vector<Vertex*>& new_verts, std::vector<Triangle*>& new_tris, \
		std::vector<int>& del_edge, std::set<int>& del_tris);

// load a RGB .JPG file as a texture
GLuint loadTextureJPG( const char * filename, int wrap ) {

	/* load an image file directly as a new OpenGL texture */
	GLuint tex_2d = SOIL_load_OGL_texture
		(
			filename,
			SOIL_LOAD_AUTO,
				SOIL_CREATE_NEW_ID,
				SOIL_FLAG_POWER_OF_TWO
				| SOIL_FLAG_MIPMAPS
				| SOIL_FLAG_MULTIPLY_ALPHA
				| SOIL_FLAG_COMPRESS_TO_DXT
				| SOIL_FLAG_DDS_LOAD_DIRECT
				| SOIL_FLAG_INVERT_Y
		);
	return tex_2d;

}

GLuint loadTextureJPG_( const char * filename, int wrap ) {

    GLuint texture;

	int x, y, channels, force_channel=1;
	unsigned char *img = stbi_load(filename, &x, &y, &channels, force_channel);
	printf ("Loaded Texture Image %s\nx: %d, y: %d, channels: %d\n", filename, x, y, channels);
	
    // allocate a texture name
    glGenTextures( 1, &texture );
	printf("Texture name: %d\n",texture);
	if (texture == GL_INVALID_VALUE) {
		printf("Invalid Texture name\n",texture);
		exit(0);
	}

    // select our current texture
    glBindTexture( GL_TEXTURE_2D, texture );

    // select modulate to mix texture with color for shading
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );

    // when texture area is small, bilinear filter the closest mipmap
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                     GL_LINEAR_MIPMAP_NEAREST );
    // when texture area is large, bilinear filter the first mipmap
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    // if wrap is true, the texture wraps over at the edges (repeat)
    //       ... false, the texture ends at the edges (clamp)
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                     wrap ? GL_REPEAT : GL_CLAMP );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                     wrap ? GL_REPEAT : GL_CLAMP );

    // build our texture mipmaps
	gluBuild2DMipmaps( GL_TEXTURE_2D, force_channel ? force_channel : channels, x, y,
                       GL_RGB, GL_UNSIGNED_BYTE, img );
	
    // free buffer
	stbi_image_free(img);
	
    return texture;
}

double findLargestEdgeLength(Polyhedron* poly) {
	double len = 0.0;
	for (int i = 0; i < poly->nedges; i++) {
		if (poly->elist[i]->length > len) {
			len = poly->elist[i]->length;
		}
	}
	return len;
}

void divide_one_edge(CornerPoint* q, Polyhedron* poly, std::vector<Vertex*>& new_verts, std::vector<Triangle*>& new_tris, \
		std::vector<int>& del_edge, std::set<int>& del_tris, std::map<int, Vertex*>& corner_vert) {
	del_tris.insert(q->t->index);
	del_edge.push_back(q->e->index);
	Vertex* vert = new_middle_vertex(q, new_verts, corner_vert);
	new_tris.push_back(new_triangle(q->v, q->n->v, vert));
	new_tris.push_back(new_triangle(q->v, vert, q->p->v));
}

void divide_two_edges(CornerPoint* q1, CornerPoint* q2, \
		Polyhedron* poly, std::vector<Vertex*>& new_verts, std::vector<Triangle*>& new_tris, \
		std::vector<int>& del_edge, std::set<int>& del_tris, std::map<int, Vertex*>& corner_vert) {
	del_tris.insert(q1->t->index);
	del_edge.push_back(q1->e->index);
	del_edge.push_back(q2->e->index);
	Vertex* v1 = new_middle_vertex(q1, new_verts, corner_vert);
	Vertex* v2 = new_middle_vertex(q2, new_verts, corner_vert);
	new_tris.push_back(new_triangle(q1->v, q2->v, v1));
	new_tris.push_back(new_triangle(q1->v, v1, v2));
	new_tris.push_back(new_triangle(v1, q1->p->v, v2));
}

void divide_three_edges(CornerPoint* q1, CornerPoint* q2, CornerPoint* q3, \
		Polyhedron* poly, std::vector<Vertex*>& new_verts, std::vector<Triangle*>& new_tris, \
		std::vector<int>& del_edge, std::set<int>& del_tris, std::map<int, Vertex*>& corner_vert) {
	del_tris.insert(q1->t->index);
	del_edge.push_back(q1->e->index);
	del_edge.push_back(q2->e->index);
	del_edge.push_back(q3->e->index);
	Vertex* v1 = new_middle_vertex(q1, new_verts, corner_vert);
	Vertex* v2 = new_middle_vertex(q2, new_verts, corner_vert);
	Vertex* v3 = new_middle_vertex(q3, new_verts, corner_vert);
	new_tris.push_back(new_triangle(q1->v, v3, v2));
	new_tris.push_back(new_triangle(q2->v, v1, v3));
	new_tris.push_back(new_triangle(q3->v, v2, v1));
	new_tris.push_back(new_triangle(v1, v2, v3));
}

int find_edges_to_divide(CornerPoint* c, CornerPoint** r, float thres) {
	int l = 0;
	CornerPoint *q = c;
	r[0]=r[1]=r[2]=0;
	if (q->e->length > thres) r[l++] = q;
	q = q->n;
	if (q->e->length > thres) r[l++] = q;
	q = q->n;
	if (q->e->length > thres) r[l++] = q;
	// order the corners
	if (l == 2 && r[1]->n->index == r[0]->index) {
		CornerPoint* t = r[0];
		r[0] = r[1]; r[1] = t;
	}
	return l;
}

void subdivide_irregular(Polyhedron* poly, CornerPoints* cp, float thres) {

	CornerPoint* c = cp->corners;
	std::set<int> del_tris;
	std::vector<int> del_edge;
	std::vector<Vertex*> new_verts;
	std::vector<Triangle*> new_tris;
	std::map<int, Vertex*> corner_vert;

	for (int i = 0; i < cp->length; i++) {

		if (c[i].index > c[i].n->index || c[i].index > c[i].p->index) continue; // only one of the vertices of a triangle

		if (!c[i].o) { printf("WARN:: Unclosed edge!!!\n"); continue; }

		CornerPoint *r[3];
		int l = find_edges_to_divide(&c[i], r, thres);
		//printf("Num divides for %3d: %3d\n", c[i].index, l);
		if (l == 1) 
			divide_one_edge(r[0], poly, new_verts, new_tris, del_edge, del_tris, corner_vert);
		else if (l == 2) 
			divide_two_edges(r[0], r[1], poly, new_verts, new_tris, del_edge, del_tris, corner_vert);
		else if (l == 3) 
			divide_three_edges(r[0], r[1], r[2], poly, new_verts, new_tris, del_edge, del_tris, corner_vert);

	}

	//printf("Before update # Vertices: %d\n", poly->nverts);
	//printf("Before update # Faces: %d\n", poly->ntris);
	//printf("# Triangle sub divides: %d / %d\n", del_tris.size(), poly->ntris);
	//printf("# Edge sub divides: %d / %d\n", del_edge.size(), poly->nedges);
	//printf("# New Faces: %d + %d - %d\n", new_tris.size(), poly->ntris, del_tris.size());

	// important:: The input CornerPoints structure will be obsolete with update_polyhedron()
	update_polyhedron(poly, new_verts, new_tris, del_edge, del_tris);

}

Vertex* new_middle_vertex(Vertex* a, Vertex* b, Vertex* c, Vertex* d) {
	Vertex* vert = new Vertex((a->x+b->x)*3/8+(c->x+d->x)*1/8, \
		(a->y+b->y)*3/8+(c->y+d->y)*1/8, \
		(a->z+b->z)*3/8+(c->z+d->z)*1/8);
	return vert;
}

Vertex* new_middle_vertex(CornerPoint* q, std::vector<Vertex*>& new_verts, std::map<int, Vertex*>& corner_vert) {
	Vertex* vert = 0;
	std::map<int, Vertex*>::iterator it = corner_vert.find(q->index);
	if (it != corner_vert.end()) 
		vert = it->second;
	else {
		Vertex *v1=q->e->verts[0], *v2=q->e->verts[1];
		vert = new Vertex((v1->x+v2->x)/2, (v1->y+v2->y)/2, (v1->z+v2->z)/2);
		new_verts.push_back(vert);
		corner_vert.insert(std::pair<int, Vertex*>(q->index, vert));
		if (q->o) corner_vert.insert(std::pair<int, Vertex*>(q->o->index, vert));
	}
	return vert;
}

CornerPoint* corner_opp_largest_edge(CornerPoint* c) {
	CornerPoint *q = c->n, *r = c;
	if (q->e->length > r->e->length) r = q;
	q = q->n;
	if (q->e->length > r->e->length) r = q;
	return r;
}

Triangle* new_triangle(Vertex* v1, Vertex* v2, Vertex* v3) {
	Triangle* tris = new Triangle;
	tris->nverts = 3;
	tris->verts[0] = v1;
	tris->verts[1] = v2;
	tris->verts[2] = v3;
	return tris;
}

void update_polyhedron(Polyhedron* poly, std::vector<Vertex*>& new_verts, std::vector<Triangle*>& new_tris, \
		std::vector<int>& del_edge, std::set<int>& del_tris) {
	if (new_verts.size() > 0) {
		int max_verts = poly->nverts + new_verts.size();
		Vertex** new_vlist = new Vertex*[max_verts];
		for (int i = 0; i < poly->nverts; i++) {
			new_vlist[i] = poly->vlist[i];
			delete poly->vlist[i]->tris;
			poly->vlist[i]->tris = NULL;
		}
		std::vector<Vertex*>::iterator it = new_verts.begin();
		for (int i = poly->nverts; i < max_verts; i++, it++) {
			new_vlist[i] = *it;
		}
		delete poly->vlist;
		poly->max_verts = poly->nverts = max_verts;
		poly->vlist = new_vlist;
	}
	if (new_tris.size() > 0 || del_tris.size() > 0) {
		int max_tris = poly->ntris + new_tris.size() - del_tris.size();
		Triangle** new_tlist = new Triangle*[max_tris];
		int j = 0;
		for (int i = 0; i < poly->ntris; i++) {
			if (del_tris.find(poly->tlist[i]->index) == del_tris.end())
				new_tlist[j++] = poly->tlist[i];
			else {
				//free(poly->tlist[i]->other_props);
				delete poly->tlist[i];
			}
		}
		std::vector<Triangle*>::iterator it = new_tris.begin();
		for (; j < max_tris; j++, it++) {
			new_tlist[j] = *it;
		}
		delete poly->tlist;
		poly->max_tris = poly->ntris = max_tris;
		poly->tlist = new_tlist;
	}
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

void subdivide(Polyhedron* poly, CornerPoints* cp, float thres) {
	CornerPoint* c = cp->corners;
	std::set<int> del_tris;
	std::vector<int> del_edge;
	std::vector<Vertex*> new_verts;
	std::vector<Triangle*> new_tris;
	for (int i = 0; i < cp->length; i++) {
		if (!c[i].o) { printf("WARN:: Unclosed edge!!!\n"); continue; }
		if (c[i].o && c[i].index > c[i].o->index) continue; // only one of the opposite points
		if (del_tris.find(c[i].t->index) != del_tris.end()) continue; // face already processed
		if (del_tris.find(c[i].o->t->index) != del_tris.end()) continue; // face already processed
		CornerPoint* q = corner_opp_largest_edge(&c[i]);
		if (del_tris.find(q->o->t->index) != del_tris.end()) q = &c[i]; // use the current corner
		bool divide = false;
		if (q->e->length > thres) {
			del_tris.insert(q->t->index);
			del_tris.insert(q->o->t->index);
			//printf("Deleting corner %3d faces: %3d, %3d\n", q->index, q->t->index, q->o->t->index);
			del_edge.push_back(q->e->index);
			Vertex* vert = new_middle_vertex(q->e->verts[0], q->e->verts[1], q->v, q->o->v);
			new_verts.push_back(vert);
			new_tris.push_back(new_triangle(q->v, q->n->v, vert));
			new_tris.push_back(new_triangle(q->v, vert, q->p->v));
			new_tris.push_back(new_triangle(q->o->v, q->o->n->v, vert));
			new_tris.push_back(new_triangle(q->o->v, vert, q->o->p->v));
			divide = true;
		}
	}

	//printf("Before update # Vertices: %d\n", poly->nverts);
	//printf("Before update # Faces: %d\n", poly->ntris);
	//printf("# Triangle sub divides: %d / %d\n", del_tris.size(), poly->ntris);
	//printf("# Edge sub divides: %d / %d\n", del_edge.size(), poly->nedges);
	//printf("# New Faces: %d + %d - %d\n", new_tris.size(), poly->ntris, del_tris.size());

	// important:: The input CornerPoints structure will be obsolete with update_polyhedron()
	update_polyhedron(poly, new_verts, new_tris, del_edge, del_tris);

}

void drawSilhouetteEdgesLevelset(CornerPoints *cp, icVector3& r) {
	CornerPoint* c = cp->corners;
	GLfloat mat_diffuse[4];
	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 0.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;
	glDisable(GL_BLEND);
	//printf("Silhouette Edges\n");
	for (int i = 0; i < cp->length; i++) {
		// process only one triangle... alternatively, i = i + 3.
		if (c[i].index > c[i].n->index || c[i].index > c[i].p->index) continue;
		CornerPoint* q = &c[i];
		int k = 0;
		Vertex v1(0,0,0), v2(0,0,0);
		Vertex* pv[2]; pv[0] = &v1; pv[1] = &v2;
		for (int j = 0; j < 3; j++, q = q->n) {
			double d1 = dot(q->e->verts[0]->normal, r);
			double d2 = dot(q->e->verts[1]->normal, r);
			double max1=std::max<double>(d1, d2), min1=std::min<double>(d1,d2);
			if (max1 < 0 && min1 < 0) continue;
			if (max1 > 0 && min1 > 0) continue;
			if (max1 >= 0.0 && min1 <= 0.0) {
				double f = max1/(max1-min1);
				Vertex *v = q->e->verts[0], *w = q->e->verts[1];
				if (d1 > d2) { v = q->e->verts[1]; w = q->e->verts[0]; };
				pv[k]->x = f*v->x + (1-f)*(w->x);
				pv[k]->y = f*v->y + (1-f)*(w->y);
				pv[k]->z = f*v->z + (1-f)*(w->z);
				k++;
				//printf(">Edge: %d, %4.3f\n",c[i].e->index, f);
			}
			if (k == 2) {
				glLineWidth(3.0);
				glEnable(GL_LINE_SMOOTH);
				glBegin(GL_LINES);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
				glVertex3f(v1.x, v1.y, v1.z);
				glVertex3f(v2.x, v2.y, v2.z);
				glEnd();
				break;
			}
		}
	}
}

void drawSilhouetteEdges(CornerPoints *cp, icVector3& r) {
	std::vector<Edge*>* sil = new std::vector<Edge*>();
	CornerPoint* c = cp->corners;
	//printf("Silhouette Edges\n");
	for (int i = 0; i < cp->length; i++) {
		if (!c[i].o) continue;
		if (c[i].index > c[i].o->index) continue;
		double d1 = dot(c[i].t->normal, r);
		double d2 = dot(c[i].o->t->normal, r);
		if (std::max<double>(d1, d2) > 0.0 && std::min<double>(d1,d2) <= 0.0) {
			sil->push_back(c[i].e);
			//printf("Edge: %d\n",c[i].e->index);
		}
	}
	drawEdges(sil);
	delete sil;
}

void drawEdges(std::vector<Edge*>* sil) {
	GLfloat mat_diffuse[4];
	mat_diffuse[0] = 1.0;
	mat_diffuse[1] = 0.0;
	mat_diffuse[2] = 0.0;
	mat_diffuse[3] = 1.0;
	std::vector<Edge*>::iterator it;
	glDisable(GL_BLEND);
	for (it = sil->begin(); it != sil->end(); it++) {
		Edge* e = *it;
		glLineWidth(3.0);
		glEnable(GL_LINE_SMOOTH);
		glBegin(GL_LINES);
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glVertex3f(e->verts[0]->x, e->verts[0]->y, e->verts[0]->z);
		glVertex3f(e->verts[1]->x, e->verts[1]->y, e->verts[1]->z);
		glEnd();
	}
}

void display_silhouette(GLenum mode, Polyhedron *this_poly, CornerPoints *cp, float** color_map, int display_mode, CommandOptions* options)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];

	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	GLdouble mvmatrix[16];
	glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);

	//GLdouble projmatrix[16];
	//glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);

	icVector3 r(mvmatrix[2], mvmatrix[6], mvmatrix[10]);

	bool setpolycolor = (options->coloring_type_poly > 0);
	bool setvertcolor = (options->coloring_type_vert > 0);

	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
			glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];

		if (setpolycolor) {
			double g = dot(r, temp_t->normal);
			//printf("Face: %d, g=%f", temp_t->index, g);
			switch (options->coloring_type_poly) {
			case 1:
				if (g >= 0) {
					mat_diffuse[0] = 1.0;
					mat_diffuse[1] = 1.0;
					mat_diffuse[2] = 1.0;
					mat_diffuse[3] = 1.0;
				} else {
					mat_diffuse[0] = 0.0;
					mat_diffuse[1] = 0.0;
					mat_diffuse[2] = 0.0;
					mat_diffuse[3] = 1.0;
				}
				if (false) {
					mat_diffuse[0] = fabs(temp_t->normal.entry[0]);
					mat_diffuse[1] = fabs(temp_t->normal.entry[1]);
					mat_diffuse[2] = fabs(temp_t->normal.entry[2]);
					mat_diffuse[3] = 1.0;
				}
				break;
			case 2:
				mat_diffuse[0] = color_map[i][0];
				mat_diffuse[1] = color_map[i][1];
				mat_diffuse[2] = color_map[i][2];
				mat_diffuse[3] = 1.0;
				break;
			}
			//glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		}
		glBegin(GL_POLYGON);
		for (j=0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
			if (setvertcolor) {
				switch (options->coloring_type_vert) {
				case 1:
					mat_diffuse[0] = fabs(temp_t->normal.entry[0]);
					mat_diffuse[1] = fabs(temp_t->normal.entry[1]);
					mat_diffuse[2] = fabs(temp_t->normal.entry[2]);
					mat_diffuse[3] = 1.0;
					break;
				case 2:
					mat_diffuse[0] = 1.0 - ((int)(floor(temp_v->x/options->l)) % 2);
					mat_diffuse[1] = 1.0 - ((int)(floor(temp_v->y/options->l)) % 2);
					mat_diffuse[2] = 1.0 - ((int)(floor(temp_v->z/options->l)) % 2);
					mat_diffuse[3] = 1.0;
					break;
				}
				glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			}
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}

	if (options->silhouette == 1)
		drawSilhouetteEdges(cp, r);
	else if (options->silhouette == 2)
		drawSilhouetteEdgesLevelset(cp, r);

}

void draw_axes_bbox(vertex_list_summary* summary, bool axes, bool bbox) {

	if (!(axes || bbox)) return;

	float center[3];
	center[0] = summary->aggregates[0][0]/(float)summary->nverts;
	center[1] = summary->aggregates[0][1]/(float)summary->nverts;
	center[2] = summary->aggregates[0][2]/(float)summary->nverts;

	float* min = summary->aggregates[SUMMARY_MIN_IDX];
	float* max = summary->aggregates[SUMMARY_MAX_IDX];
	float** eig = summary->eigen_vec;

	// all axis
	GLfloat axiscolors[] = {0.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0};
	if (axes) {
		for (int i=0; i < 3; i++) {
			line_from_and_direction(center[0], center[1], center[2], eig[0][i], eig[1][i], eig[2][i],  2.0, axiscolors+i*3);
			line_from_and_direction(center[0], center[1], center[2], eig[0][i], eig[1][i], eig[2][i], -2.0, axiscolors+i*3);
		}
	}
	
	if (bbox) {
		if (!summary->bounding_box) return;
		draw_bounding_box(summary->bounding_box, summary->box_faces);
	}
	
}

void draw_bounding_box(float** bounding_box, int** faces) {
	GLfloat mat_diffuse[4];
	glEnable(GL_DEPTH_TEST);
	for (int i = 0; i < 6; i++) { // faces
		switch (i) {
		case 4:
		case 5:
			mat_diffuse[0] = 1.0;
			mat_diffuse[1] = 0.0;
			mat_diffuse[2] = 0.0;
			break;
		case 1:
		case 3:
			mat_diffuse[0] = 0.0;
			mat_diffuse[1] = 0.0;
			mat_diffuse[2] = 1.0;
			break;
		case 0:
		case 2:
			mat_diffuse[0] = 0.0;
			mat_diffuse[1] = 1.0;
			mat_diffuse[2] = 0.0;
			break;
		}
		mat_diffuse[3] = 0.3;
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
		glBegin(GL_POLYGON);
		for (int j = 0; j < 4; j++) {
			float* vert = bounding_box[faces[i][j]];
			glVertex3f(vert[0], vert[1], vert[2]);
		}
		glEnd();
	}
}

void line_from_and_direction(float x1, float y1, float z1, float a, float b, float c, float scale, GLfloat* rgb) {
	glColor3f(rgb[0], rgb[1], rgb[2]);
	//glColor3f(1.0,0.0,0.0);
	glLineWidth(3.0);
	glEnable(GL_LINE_SMOOTH);
	glBegin(GL_LINES);
	glVertex3f(x1, y1, z1);
	glVertex3f(x1+a*scale, y1+b*scale, z1+c*scale);
	glEnd();
}
