#include <string>
#include <vector>
#include <stdio.h>

#ifndef _CS554_UTILS_H_
#define _CS554_UTILS_H_

typedef struct _CommandOptions {
	std::string file;
	int coloring_type_poly;
	int coloring_type_vert;
	float l; // checkerboard divider
	float t; // threshold for sub-division in terms of percent of max bounding box length
	bool norm;
	bool axes;
	bool bbox;
	bool euler;
	bool debug;
	int subdivide;
	int silhouette;
	int r;
	bool smooth;
	int wtscheme;
	int iterscheme;
	double lambda;
	int nsmooth;
	bool heatdiff;
	bool saddle;
	bool texture;
	std::string texturefile;
	bool curv;
	int curvtype;
	bool curvtensor;
	int curvaxis;
	int curvnsmooth;
	bool curvstream;
	bool segment;
	int segmentdepth;
	int segmentmaxfaces;
	std::string segmentloadfile;
	std::string segmentsavefile;
	bool segmenttexture;
	float segmentmaxlen;
	int argc;
	char **argv;
	_CommandOptions(): \
		file(std::string()), coloring_type_poly(0), coloring_type_vert(0), l(1.0), t(0.1), \
		norm(false), axes(false), bbox(false), euler(false), debug (false), subdivide(false), \
		silhouette(0), r(1), smooth(false), wtscheme(0), iterscheme(0), lambda(0.1), nsmooth(10), \
		heatdiff(false), saddle(false), texture(false), \
		texturefile(std::string("sample_textures\\checkerboard.jpg")), \
		curv(false), curvtype(0), curvtensor(false), curvaxis(0), curvnsmooth(0), curvstream(false), \
		segment(false), segmentdepth(1), segmentmaxfaces(40), \
		segmentloadfile(std::string()), segmentsavefile(std::string()), segmenttexture(false), segmentmaxlen(0.5), \
		argc(0), argv(0) {}
	void print_options() {
		printf("file=%s\n",file.c_str());
		printf("coloring_type_poly=%d\n",coloring_type_poly);
		printf("coloring_type_vert=%d\n",coloring_type_vert);
		printf("l=%f\n",l);
		printf("t=%f\n",t);
		printf("norm=%s\n",(norm ? "true" : "false"));
		printf("axes=%s\n",(axes ? "true" : "false"));
		printf("bbox=%s\n",(bbox ? "true" : "false"));
		printf("euler=%s\n",(euler ? "true" : "false"));
		printf("debug=%s\n",(debug ? "true" : "false"));
		printf("subdivide=%d\n",subdivide);
		printf("silhouette=%d\n",silhouette);
		printf("r=%d\n",r);
		printf("smooth=%s\n",(smooth ? "true" : "false"));
		printf("wtscheme=%d\n",wtscheme);
		printf("iterscheme=%d\n",iterscheme);
		printf("lambda=%f\n",((float)lambda));
		printf("nsmooth=%d\n",nsmooth);
		printf("heatdiff=%s\n",(heatdiff ? "true" : "false"));
		printf("saddle=%s\n",(saddle ? "true" : "false"));
		printf("texture=%s\n",(texture ? "true" : "false"));
		printf("texturefile=%s\n",texturefile.c_str());
		printf("curv=%s\n",(curv ? "true" : "false"));
		printf("curv.type=%d\n",curvtype);
		printf("curv.tensor=%s\n",(curvtensor ? "true" : "false"));
		printf("curv.axis=%d\n",curvaxis);
		printf("curv.nsmooth=%d\n",curvnsmooth);
		printf("curv.stream=%s\n",(curvstream ? "true" : "false"));
		printf("segment=%s\n",(segment ? "true" : "false"));
		printf("segment.depth=%d\n",segmentdepth);
		printf("segment.maxfaces=%d\n",segmentmaxfaces);
		printf("segment.loadfile=%s\n",segmentloadfile.c_str());
		printf("segment.savefile=%s\n",segmentsavefile.c_str());
		printf("segment.texture=%s\n",(segmenttexture ? "true" : "false"));
		printf("segment.maxlen=%5.4f\n",segmentmaxlen);
	}
} CommandOptions;

CommandOptions *parse_command_line_options(int argc, char **argv);
void color_mapping(double percentage, double col[3]);

template<class T>
void release_vector_elements(std::vector< T* > *v) {
	if (!v) return;
	for (std::vector< T* >::iterator it = v->begin(); it != v->end(); it++) {
		delete *it;
	}
}

template<class T>
void release_objects(T** arr, int len) {
	if (!arr) return;
	for (int i = 0; i < len; i++) {
		delete arr[i];
	}
}

template<class T>
void printArray(T* array, int len) {
	for (int i = 0; i < len; i++) {
		printf("%4.3f ",(float)array[i]);
	}
	printf("\n");
}

#endif
