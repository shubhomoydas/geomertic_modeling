#include "glut.h"
#include "mathfuncs.h"
#include "utils.h"
#include "CornerPoints.h"
#include "graphcut.h"

#ifndef __GRAPHUTILS_H__

#define __GRAPHUTILS_H__

GLuint loadTextureJPG( const char * filename, int wrap );
void draw_axes_bbox(vertex_list_summary* summary, bool axes, bool bbox);
void display_silhouette(GLenum mode, Polyhedron *this_poly, CornerPoints *cp, float** color_map, int display_mode, CommandOptions* options);
void subdivide(Polyhedron* poly, CornerPoints* cp, float thres);
void subdivide_irregular(Polyhedron* poly, CornerPoints* cp, float thres);
double findLargestEdgeLength(Polyhedron* poly);

void display_curv (
	GLenum mode, 
	Polyhedron *this_poly, 
	vertex_tensor_list *vtl,
	tris_tensor_list *ttl,
	int display_mode, 
	CommandOptions* options
);

void display_patches (
	GLenum mode, 
	Polyhedron *this_poly, 
	vertex_tensor_list *vtl, 
	tris_tensor_list *ttl,
	segments *segs,
	int display_mode, 
	GLuint texture,
	CommandOptions* options
);

void draw_tensor_edges(vertex_tensor_list *vtl, CommandOptions* options);
void draw_streamlines(vertex_tensor_list *vtl, tris_tensor_list *ttl, CommandOptions* options);

#endif
