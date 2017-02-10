/*

Functions for learnply

Eugene Zhang, 2005
*/

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
#include "segment.h"

unsigned char orientation;  // 0=ccw, 1=cw

FILE *this_file;
const int win_width=1024;
const int win_height=1024;

double radius_factor = 0.9;

int display_mode = 0; 
double error_threshold = 1.0e-13;
char reg_model_name[128];
FILE *f;	
int ACSIZE = 1; // for antialiasing
int view_mode=0;  // 0 = othogonal, 1=perspective
float s_old, t_old;
float rotmat[4][4];
static Quaternion rvec;

int mouse_mode = -2;  // -2=no action, -1 = down, 0 = zoom, 1 = rotate x, 2 = rotate y, 3 = tranlate x, 4 = translate y, 5 = cull near 6 = cull far
int mouse_button = -1; // -1=no button, 0=left, 1=middle, 2=right
int last_x, last_y;

struct jitter_struct{
	double x;
	double y;
} jitter_para;

jitter_struct ji1[1] = {{0.0, 0.0}};
jitter_struct ji16[16] = {{0.125, 0.125}, {0.375, 0.125}, {0.625, 0.125}, {0.875, 0.125}, 
{0.125, 0.375}, {0.375, 0.375}, {0.625, 0.375}, {0.875, 0.375}, 
{0.125, 0.625}, {0.375, 0.625}, {0.625, 0.625}, {0.875, 0.625}, 
{0.125, 0.875}, {0.375, 0.875}, {0.625, 0.875}, {0.875, 0.875}, };

Polyhedron *poly;

void init(void);
void set_view(GLenum mode, Polyhedron *poly);
void set_scene(GLenum mode, Polyhedron *poly);
void keyboard(unsigned char key, int x, int y);
void motion(int x, int y);
void display(void);
void mouse(int button, int state, int x, int y);

void display_original(GLenum mode, Polyhedron *this_poly);
void display_shape(GLenum mode, Polyhedron *poly);

/* MOY:: Added for HW1 Problem 2 */
CommandOptions *options;
float** color_map;
vertex_list_summary* poly_summary = 0;
CornerPoints *cp = 0;

std::map<int,double> predefs;
double *diffusion_colors;
int *saddleverts;

texture_map *tm = 0;
GLuint texture = 0;

vertex_tensor_list *vtl = 0;
tris_tensor_list *ttl = 0;

segments *segs = 0;

/******************************************************************************
Main program.
******************************************************************************/

int learnply_main(CommandOptions *_options)
{
	char *progname;
	int num = 1;
	FILE *this_file;

	options = _options;
	progname = options->argv[0];

	this_file = fopen((options->file).c_str(), "r");
	poly = new Polyhedron (this_file);
	fclose(this_file);
	mat_ident( rotmat );

	poly->initialize(); // initialize everything

	poly->calc_bounding_sphere();
	poly->calc_face_normals_and_area();
	poly->average_normals();

	if (options->euler) {
		int euler = poly->nverts - poly->nedges + poly->ntris;
		printf("Euler Characteristic: %4d - %4d + %4d = %d\n",poly->nverts, poly->nedges, poly->ntris, euler);
		printf("Handles: %d\n",(1-euler/2));
		if (!(options->axes || options->bbox || options->subdivide || options->silhouette)) {
			//poly->finalize();
			return 0;
		}
	}

	orientation = poly->orientation;

	bool withnorm = options->norm;
	
	//printf("Vertices: %d\n",poly->nverts);
	if (options->subdivide || options->axes || options->bbox) {
		poly_summary = compute_vertex_list_summary(poly->vlist, poly->nverts, false);
	}

	if (options->subdivide > 0) {

		double thres = poly_summary->getSmoothingThreshold(options->t);
		double largest = findLargestEdgeLength(poly);
		int ndivs = ceil(log(largest/thres)/log(2.0));
		ndivs = min(options->r, ndivs);
		if (options->debug) {
			printf("Threshold = %8.7f, largest = %8.7f, ndivs = %3d, r = %3d\n", thres, largest, ndivs, options->r);
		}

		if (options->debug) {
			printf("Before subdivide # Vertices: %d\n", poly->nverts);
			printf("Before subdivide # Edge: %d\n", poly->nedges);
			printf("Before subdivide # Faces: %d\n", poly->ntris);
		}
	
		for (int divs = 0; divs < ndivs; divs++) {
			if (options->debug) printf("Subdivision iteration: %3d/%3d\n", divs+1, ndivs);

			cp = new CornerPoints(poly);

			if (options->subdivide == 1)
				subdivide(poly, cp, thres);
			else if (options->subdivide == 2)
				subdivide_irregular(poly, cp, thres);

			delete cp; // cp is obsolete...
			cp = 0;
		}

		if (options->debug) {
			printf("After subdivide # Vertices: %d\n", poly->nverts);
			printf("After subdivide # Edge: %d\n", poly->nedges);
			printf("After subdivide # Faces: %d\n", poly->ntris);
		}

		if (options->euler) {
			int euler = poly->nverts - poly->nedges + poly->ntris;
			printf("Euler Characteristic: %4d - %4d + %4d = %d\n", \
				poly->nverts, poly->nedges, poly->ntris, euler);
			printf("Handles: %d\n",(1-euler/2));
		}

	}

	if (options->smooth) {
		smooth_polynomial(poly, options->wtscheme, options->iterscheme, options->nsmooth, options->lambda);
		//if (true) return 0;
	}

	if (options->heatdiff) {
		diffusion_colors = new double[poly->nverts];
		saddleverts = new int[poly->nverts];

		//dragon sample points for DEBUG
		//predefs.insert(std::pair<int,double>(7605,0.0));
		//predefs.insert(std::pair<int,double>(5313,1.0));

		set_diffusion_colors(poly, diffusion_colors, saddleverts, &predefs, options->wtscheme);
	}

	if (options->curv) {
		vtl = compute_curvature_tensors(poly, options->wtscheme, options->curvnsmooth, options->lambda);
		printf("Computed vertex tensors...\n");
		double deficit = get_total_angle_deficit(vtl);
		printf("Deficit: %5.4f\n", (deficit/(2*PI)));
		//for (int i = 0; i < vtl->length; i++) {
		//	vertex_tensor *vt = vtl->vt[i];
		//	printf("Vertex: %d: Gauss: %4.3f, Mean: %4.3f\n",vt->v->index, vt->norm_gauss_curv, vt->norm_mean_curv);
		//}
		ttl = compute_tris_curvature_tensors(poly, vtl);
		printf("Computed face tensors...\n");
	}

	if (options->segment) {
		vtl = compute_curvature_tensors(poly, options->wtscheme, options->curvnsmooth, options->lambda);
		printf("Computed vertex tensors...\n");
		double deficit = get_total_angle_deficit(vtl);
		printf("Deficit: %5.4f\n", (deficit/(2*PI)));
		ttl = compute_tris_curvature_tensors(poly, vtl);
		printf("Computed face tensors...\n");
		if (!options->segmentloadfile.empty()) {
			segs = load_segments(options->segmentloadfile, ttl);
		} else {
			segs = segment_polyhedron(poly, ttl, options);
		}
		if(options->segmenttexture) prepare_parameterizations(poly, segs);
		//for (int i = 0; i < poly->ntris; i++) {
		//	for (int n = 0; n < 3; n++) {
		//		printf("%5d: [%5.4f,%5.4f] ",i,segs->params[i][n].x,segs->params[i][n].y);
		//	}
		//	printf("\n");
		//}
		if (!options->segmentsavefile.empty()) {
			save_segments(segs, options->segmentsavefile);
		}
		if (options->debug) exit(0);
		segs->createColorMap();
	}

	if (options->texture) {

		//if (options->debug) printf("Loading texture file %s\n",options->texturefile.c_str());
		//texture = loadTextureJPG(options->texturefile.c_str(), 0);
		//if (options->debug) printf("Texture file Loaded.\n");
		//exit(0);

		tm = new texture_map(poly->nverts, options->texturefile);

		// test data for bunny for DEBUG only
		tm->f_predefs->insert(std::pair<int,double>(4556, 0.0));
		tm->f_predefs->insert(std::pair<int,double>(1577, 1.0));
		tm->f_predefs->insert(std::pair<int,double>(2083, 1.0));
		tm->f_predefs->insert(std::pair<int,double>(165, 0.0));

		tm->g_predefs->insert(std::pair<int,double>(4556, 0.0));
		tm->g_predefs->insert(std::pair<int,double>(1656, 1.0));
		tm->g_predefs->insert(std::pair<int,double>(3133, 1.0));
		tm->g_predefs->insert(std::pair<int,double>(165, 0.0));

		tm->set_diffusion_colos(poly, options->wtscheme);

		if (options->debug) {
			tm->print();
			exit(0);
		}
	}

	if (options->silhouette) {
		cp = new CornerPoints(poly);
	}

	if (options->axes || options->bbox) {
		if (withnorm) {
			vertex_list_summary* norm_summary = \
				compute_vertex_list_summary(poly->vlist, poly->nverts, true);
			poly_summary->replace_eigen(norm_summary->eigen_vec, norm_summary->eigen_val);
			norm_summary->eigen_vec = 0;
			norm_summary->eigen_val = 0;
			delete norm_summary;
		}
		compute_bounding_box(poly->vlist, poly_summary);
		if (options->debug) poly_summary->print_vertex_sumary();
	}

	color_map = create_rand_3d_color_map(poly->ntris);

	glutInit(&options->argc, options->argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowPosition (20, 20);
	glutInitWindowSize (win_width, win_height); 
	glutCreateWindow ("Geometric Modeling");
	init ();
	glutKeyboardFunc (keyboard);
	glutDisplayFunc(display); 
	glutMotionFunc (motion);
	glutMouseFunc (mouse);
	glutMainLoop(); 
	poly->finalize();  // finalize everything

	delete_color_map(color_map, poly->ntris);
	if (poly_summary) delete poly_summary;

	if (cp) delete cp;

	if (texture) glDeleteTextures( 1, &texture );

	if (vtl) delete vtl;
	if (ttl) delete ttl;
	if (segs) delete segs;

	return 0;    /* ANSI C requires main to return int. */
}

void display_shape(GLenum mode, Polyhedron *this_poly)
{
	if (options->silhouette > 0)
		display_silhouette(mode, this_poly, cp, color_map, display_mode, options);
	else if (options->curv)
		display_curv(mode, this_poly, vtl, ttl, display_mode, options);
	else if (options->segment) {
		if(options->segmenttexture) {
			if (!texture) texture = loadTextureJPG(options->texturefile.c_str(), 0);
			glEnable( GL_TEXTURE_2D );
			//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
			glBindTexture( GL_TEXTURE_2D, texture );
		}
		display_patches(mode, this_poly, vtl, ttl, segs, display_mode, texture, options);
	} else
		display_original(mode, this_poly);
}

void display_original(GLenum mode, Polyhedron *this_poly)
{
	unsigned int i, j;
	GLfloat mat_diffuse[4];
	double cycle_color[3];

	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1., 1.);

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);

	if (options->texture) {
		if (!texture) texture = loadTextureJPG(options->texturefile.c_str(), 0);

		glEnable( GL_TEXTURE_2D );
		//glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
		glBindTexture( GL_TEXTURE_2D, texture );
		GLboolean istex = glIsTexture(texture);
		//printf("Texture name: %d, proper: %s\n",texture,(istex==GL_TRUE ? "true" : "false"));
	}

	for (i=0; i<this_poly->ntris; i++) {
		if (mode == GL_SELECT)
			glLoadName(i+1);

		Triangle *temp_t=this_poly->tlist[i];
		bool setpolycolor = true;

		switch (display_mode) {
		case 0:
			switch (options->coloring_type_poly) {
			case 1:
				mat_diffuse[0] = fabs(temp_t->normal.entry[0]);
				mat_diffuse[1] = fabs(temp_t->normal.entry[1]);
				mat_diffuse[2] = fabs(temp_t->normal.entry[2]);
				mat_diffuse[3] = 1.0;
				break;
			case 2:
				mat_diffuse[0] = color_map[i][0];
				mat_diffuse[1] = color_map[i][1];
				mat_diffuse[2] = color_map[i][2];
				mat_diffuse[3] = 1.0;
				break;
			default:
				setpolycolor = false;
			}
			if (setpolycolor) glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
			glBegin(GL_POLYGON);
			glEdgeFlag(TRUE);
			for (j=0; j<3; j++) {

				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_v->normal.entry[0], temp_v->normal.entry[1], temp_v->normal.entry[2]);
				bool setvertcolor = true;
				switch (options->coloring_type_vert) {
				case 1:
					mat_diffuse[0] = fabs(temp_t->normal.entry[0]);
					mat_diffuse[1] = fabs(temp_t->normal.entry[1]);
					mat_diffuse[2] = fabs(temp_t->normal.entry[2]);
					mat_diffuse[3] = 1.0;
					//glColor3d(fabs(temp_v->normal.entry[0]), \
					//	fabs(temp_v->normal.entry[1]), fabs(temp_v->normal.entry[2]));
					break;
				case 2:
					mat_diffuse[0] = 1.0 - ((int)(floor(temp_v->x/options->l)) % 2);
					mat_diffuse[1] = 1.0 - ((int)(floor(temp_v->y/options->l)) % 2);
					mat_diffuse[2] = 1.0 - ((int)(floor(temp_v->z/options->l)) % 2);
					mat_diffuse[3] = 1.0;
					//glColor3d(1 - ((int)(floor(temp_v->x/options->l)) % 2), \
					//	1 - ((int)(floor(temp_v->y/options->l)) % 2), \
					//	1 - ((int)(floor(temp_v->z/options->l)) % 2));
					break;
				case 3: // heat diffusion
					color_mapping(diffusion_colors[temp_v->index], cycle_color);
					mat_diffuse[0] = cycle_color[0];
					mat_diffuse[1] = cycle_color[1];
					mat_diffuse[2] = cycle_color[2];
					mat_diffuse[3] = 1.0;
					break;
				case 4: // mark only saddlepoints
					//color_mapping(diffusion_colors[temp_v->index], cycle_color);
					mat_diffuse[0] = saddleverts[temp_v->index] ? 1 : 1;
					mat_diffuse[1] = saddleverts[temp_v->index] ? 0 : 1;
					mat_diffuse[2] = saddleverts[temp_v->index] ? 0 : 1;
					mat_diffuse[3] = 1.0;
					break;
				case 15: // texture
					tm->get_color(temp_v->index, cycle_color);
					mat_diffuse[0] = cycle_color[0];
					mat_diffuse[1] = cycle_color[1];
					mat_diffuse[2] = cycle_color[2];
					mat_diffuse[3] = 1.0;
					setvertcolor = false;
					break;
				default:
					setvertcolor = false;
					break;
				}
				if (setvertcolor) glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
				if (options->texture) {
					glTexCoord2d(tm->f[temp_v->index],tm->g[temp_v->index]);
				}
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
			}
			glEnd();
			break;

		case 6:
			glBegin(GL_POLYGON);
			for (j=0; j<3; j++) {
				Vertex *temp_v = temp_t->verts[j];
				glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
				glColor3f(1.0, 1.0, 1.0);
				glVertex3d(temp_v->x, temp_v->y, temp_v->z);
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
	if (options->texture) {
		glDisable(GL_TEXTURE_2D);
	}
}

void display(void)
{
	GLint viewport[4];
	int jitter;

	glClearColor (1.0, 1.0, 1.0, 1.0);  // background for rendering color coding and lighting
	glGetIntegerv (GL_VIEWPORT, viewport);

	glClear(GL_ACCUM_BUFFER_BIT);
	for (jitter = 0; jitter < ACSIZE; jitter++) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		set_view(GL_RENDER, poly);
		glPushMatrix ();
		switch(ACSIZE){
		case 1:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;

		case 16:
			glTranslatef (ji16[jitter].x*2.0/viewport[2], ji16[jitter].y*2.0/viewport[3], 0.0);
			break;

		default:
			glTranslatef (ji1[jitter].x*2.0/viewport[2], ji1[jitter].y*2.0/viewport[3], 0.0);
			break;
		}
		set_scene(GL_RENDER, poly);

		display_shape(GL_RENDER, poly);

		// Below line is for CS 554 hw1
		draw_axes_bbox(poly_summary, options->axes, options->bbox);

		glPopMatrix ();
		glAccum(GL_ACCUM, 1.0/ACSIZE);
	}
	glAccum (GL_RETURN, 1.0);
	glFlush();
	glutSwapBuffers();
	glFinish();
}

void init(void) {
	/* select clearing color */ 

	glClearColor (0.0, 0.0, 0.0, 0.0);  // background
	glShadeModel (GL_FLAT);
	glPolygonMode(GL_FRONT, GL_FILL);

	glDisable(GL_DITHER);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	// may need it
	glPixelStorei(GL_PACK_ALIGNMENT,1);
	glEnable(GL_NORMALIZE);
	if (orientation == 0) 
		glFrontFace(GL_CW);
	else 
		glFrontFace(GL_CCW);
}


/******************************************************************************
Process a keyboard action.  In particular, exit the program when an
"escape" is pressed in the window.
******************************************************************************/

void keyboard(unsigned char key, int x, int y) {
	int i;

	/* set escape key to exit */
	switch (key) {
	case 27:
		poly->finalize();  // finalize_everything
		exit(0);
		break;

	case '0':
		display_mode = 0;
		display();
		break;

	case '1':
		display_mode = 0;
		display();
		break;

	case '2':
		display_mode = 0;
		display();
		break;

	case '3':
		display_mode = 3;
		display();
		break;

	case '4':
		display_mode = 4;
		display();
		break;

	case '5':
		display_mode = 5;
		display();
		break;

	case '6':
		display_mode = 6;
		display();
		break;

	case '7':
		display_mode = 7;
		display();
		break;

	case '8':
		display_mode = 8;
		display();
		break;

	case '9':
		display_mode = 9;
		display();
		break;

	case 'x':
		switch(ACSIZE){
		case 1:
			ACSIZE = 16;
			break;

		case 16:
			ACSIZE = 1;
			break;

		default:
			ACSIZE = 1;
			break;
		}
		fprintf(stderr, "ACSIZE=%d\n", ACSIZE);
		display();
		break;

	case '|':
		this_file = fopen("rotmat.txt", "w");
		for (i=0; i<4; i++) 
			fprintf(this_file, "%f %f %f %f\n", rotmat[i][0], rotmat[i][1], rotmat[i][2], rotmat[i][3]);
		fclose(this_file);
		break;

	case '^':
		this_file = fopen("rotmat.txt", "r");
		for (i=0; i<4; i++) 
			fscanf(this_file, "%f %f %f %f ", (&rotmat[i][0]), (&rotmat[i][1]), (&rotmat[i][2]), (&rotmat[i][3]));
		fclose(this_file);
		display();
		break;

	}
}

void multmatrix(const Matrix m)
{ 
	int i,j, index = 0;

	GLfloat mat[16];

	for ( i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			mat[index++] = m[i][j];

	glMultMatrixf (mat);
}

void set_view(GLenum mode, Polyhedron *poly)
{
	icVector3 up, ray, view;
	GLfloat light_ambient0[] = { 0.3, 0.3, 0.3, 1.0 };
	GLfloat light_diffuse0[] = { 0.7, 0.7, 0.7, 1.0 };
	GLfloat light_specular0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse1[] = { 0.5, 0.5, 0.5, 1.0 };
	GLfloat light_specular1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_ambient2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_diffuse2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position[] = { 0.0, 0.0, 0.0, 1.0 };

	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient0);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse0);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular0);
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular1);


	glMatrixMode(GL_PROJECTION);
	if (mode == GL_RENDER)
		glLoadIdentity();

	if (view_mode == 0)
		glOrtho(-radius_factor, radius_factor, -radius_factor, radius_factor, 0.0, 40.0);
	else
		gluPerspective(45.0, 1.0, 0.1, 40.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	light_position[0] = 5.5;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	light_position[0] = -0.1;
	light_position[1] = 0.0;
	light_position[2] = 0.0;
	glLightfv(GL_LIGHT2, GL_POSITION, light_position);
}

void set_scene(GLenum mode, Polyhedron *poly)
{
	glTranslatef(0.0, 0.0, -3.0);
	multmatrix( rotmat );

	glScalef(1.0/poly->radius, 1.0/poly->radius, 1.0/poly->radius);
	glTranslatef(-poly->center.entry[0], -poly->center.entry[1], -poly->center.entry[2]);
}

void motion(int x, int y) {
	float r[4];
	float xsize, ysize, s, t;

	switch(mouse_mode){
	case -1:

		xsize = (float) win_width;
		ysize = (float) win_height;

		s = (2.0 * x - win_width) / win_width;
		t = (2.0 * (win_height - y) - win_height) / win_height;

		if ((s == s_old) && (t == t_old))
			return;

		mat_to_quat( rotmat, rvec );
		trackball( r, s_old, t_old, s, t );
		add_quats( r, rvec, rvec );
		quat_to_mat( rvec, rotmat );

		s_old = s;
		t_old = t;

		display();
		break;
	}
}

int processHits(GLint hits, GLuint buffer[])
{
	unsigned int i, j;
	GLuint names, *ptr;
	double smallest_depth=1.0e+20, current_depth;
	int seed_id=-1; 
	unsigned char need_to_update;

	//printf("hits = %d\n", hits);
	ptr = (GLuint *) buffer;
	for (i = 0; i < hits; i++) {  /* for each hit  */
		need_to_update = 0;
		names = *ptr;
		ptr++;

		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		current_depth = (double) *ptr/0x7fffffff;
		if (current_depth < smallest_depth) {
			smallest_depth = current_depth;
			need_to_update = 1;
		}
		ptr++;
		for (j = 0; j < names; j++) {  /* for each name */
			if (need_to_update == 1)
				seed_id = *ptr - 1;
			ptr++;
		}
	}
	//printf("triangle id = %d\n", seed_id);
	return seed_id;
}

void mouse(int button, int state, int x, int y) {
	//if (button == GLUT_LEFT_BUTTON || button == GLUT_RIGHT_BUTTON) {
	if (button == GLUT_LEFT_BUTTON) {
		switch(mouse_mode) {
		case -2:  // no action
			if (state == GLUT_DOWN) {
				float xsize = (float) win_width;
				float ysize = (float) win_height;

				float s = (2.0 * x - win_width) / win_width;
				float t = (2.0 * (win_height - y) - win_height) / win_height;

				s_old = s;
				t_old = t;

				mouse_mode = -1;  // down
				mouse_button = button;
				last_x = x;
				last_y = y;
			}
			break;

		default:
			if (state == GLUT_UP) {
				button = -1;
				mouse_mode = -2;
			}
			break;
		}
	//} else if (button == GLUT_MIDDLE_BUTTON) {
	} else if (button == GLUT_RIGHT_BUTTON) {
		if (state == GLUT_DOWN) {  // build up the selection feedback mode

			GLuint selectBuf[win_width];
			GLint hits;
			GLint viewport[4];

			glGetIntegerv(GL_VIEWPORT, viewport);

			glSelectBuffer(win_width, selectBuf);
			(void) glRenderMode(GL_SELECT);

			glInitNames();
			glPushName(0);

			glMatrixMode(GL_PROJECTION);
			glPushMatrix();
			glLoadIdentity();
			/*  create 5x5 pixel picking region near cursor location */
			gluPickMatrix((GLdouble) x, (GLdouble) (viewport[3] - y),
				1.0, 1.0, viewport);

			set_view(GL_SELECT, poly);
			glPushMatrix ();
			set_scene(GL_SELECT, poly);
			display_shape(GL_SELECT, poly);
			glPopMatrix();
			glFlush();

			hits = glRenderMode(GL_RENDER);
			//poly->seed = processHits(hits, selectBuf);
			int tid = processHits(hits, selectBuf);
			int vid = poly->tlist[tid]->verts[0]->index;
			double val = 0.0;
			if (predefs.size() % 2 > 0) val = 1.0; // even places are 0.0, odd places 1.0
			printf("Triangle: %d, Vertex: %d\n", tid, vid);
			predefs.insert(std::pair<int,double>(vid,val));
			set_diffusion_colors(poly, diffusion_colors, saddleverts, &predefs, options->wtscheme);
			display();
		}
	}
}

void display_object()
{
	unsigned int i, j;
	Polyhedron *the_patch = poly;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	for (i=0; i<poly->ntris; i++) {
		Triangle *temp_t=poly->tlist[i];
		glBegin(GL_POLYGON);
		GLfloat mat_diffuse[] = {1.0, 1.0, 1.0, 1.0};

		glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);

		glColor3f(1.0, 1.0, 1.0);
		glNormal3d(temp_t->normal.entry[0], temp_t->normal.entry[1], temp_t->normal.entry[2]);
		for (j=0; j<3; j++) {
			Vertex *temp_v = temp_t->verts[j];
			glVertex3d(temp_v->x, temp_v->y, temp_v->z);
		}
		glEnd();
	}
}

