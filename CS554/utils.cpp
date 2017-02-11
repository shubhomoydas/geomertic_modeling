#include "stdafx.h"

CommandOptions *parse_command_line_options(int argc, char **argv) {
	CommandOptions *options = new CommandOptions();
	options->file = "ply/tetrahedron.ply";
	//options->file = "ply/icosahedron.ply";
	options->argc = argc;
	options->argv = argv;
	bool returnOptions = true;
    for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];
		//printf("ARG:::%s\n",arg.c_str());
		if ((arg == "-h") || (arg == "-help")) {
			returnOptions = false;
			break;
		} else if ((arg == "-color_poly")) {
			if (i + 1 < argc) {
				options->coloring_type_poly = atoi(argv[++i]);
			} else {
				printf("-color_poly option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-color_vert")) {
			if (i + 1 < argc) {
				options->coloring_type_vert = atoi(argv[++i]);
				options->coloring_type_vert = IMAX(0, options->coloring_type_vert);
				options->coloring_type_vert = IMIN(2, options->coloring_type_vert);
			} else {
				printf("-color_vert option requires one argument between 0 and 2.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-l")) {
			if (i + 1 < argc) {
				options->l = (float)atof(argv[++i]);
			} else {
				printf("-l option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-t")) {
			if (i + 1 < argc) {
				options->t = (float)atof(argv[++i]);
			} else {
				printf("-t option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-f") || arg == "-file") {
			if (i + 1 < argc) {
				options->file = argv[++i];
			} else {
				printf("-file option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-norm")) {
			options->norm = true;
		} else if ((arg == "-axes")) {
			options->axes = true;
		} else if ((arg == "-bbox")) {
			options->bbox = true;
		} else if ((arg == "-euler")) {
			options->euler = true;
		} else if ((arg == "-debug")) {
			options->debug = true;
		} else if ((arg == "-smooth")) {
			options->smooth = true;
		} else if ((arg == "-heatdiff")) {
			options->heatdiff = true;
		} else if ((arg == "-saddle")) {
			options->saddle = true;
		} else if ((arg == "-texture")) {
			options->texture = true;
		} else if ((arg == "-curv")) {
			options->curv = true;
		} else if ((arg == "-curv.tensor")) {
			options->curvtensor = true;
		} else if ((arg == "-curv.stream")) {
			options->curvstream = true;
		} else if ((arg == "-segment")) {
			options->segment = true;
		} else if ((arg == "-segment.texture")) {
			options->segmenttexture = true;
		} else if ((arg == "-subdivide")) {
			if (i + 1 < argc) {
				options->subdivide = atoi(argv[++i]);
			} else {
				printf("-subdivide option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-silhouette")) {
			if (i + 1 < argc) {
				options->silhouette = atoi(argv[++i]);
			} else {
				printf("-silhouette option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-r")) {
			if (i + 1 < argc) {
				options->r = atoi(argv[++i]);
				options->r = IMAX(1,options->r);
			} else {
				printf("-r option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-wtscheme")) {
			if (i + 1 < argc) {
				options->wtscheme = atoi(argv[++i]);
				options->wtscheme = IMAX(0,options->wtscheme);
			} else {
				printf("-wtscheme option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-iterscheme")) {
			if (i + 1 < argc) {
				options->iterscheme = atoi(argv[++i]);
				options->iterscheme = IMAX(0,options->iterscheme);
			} else {
				printf("-iterscheme option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-lambda")) {
			if (i + 1 < argc) {
				options->lambda = atof(argv[++i]);
			} else {
				printf("-lambda option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-nsmooth")) {
			if (i + 1 < argc) {
				options->nsmooth = atoi(argv[++i]);
				options->nsmooth = IMAX(0,options->nsmooth);
			} else {
				printf("-nsmooth option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if (arg == "-texturefile") {
			if (i + 1 < argc) {
				options->texturefile = argv[++i];
			} else {
				printf("-texturefile option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-curv.type")) {
			if (i + 1 < argc) {
				options->curvtype = atoi(argv[++i]);
				options->curvtype = IMAX(0,options->curvtype);
			} else {
				printf("-curv.type option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-curv.axis")) {
			if (i + 1 < argc) {
				options->curvaxis = atoi(argv[++i]);
				options->curvaxis = IMAX(0,options->curvaxis);
			} else {
				printf("-curv.axis option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-curv.nsmooth")) {
			if (i + 1 < argc) {
				options->curvnsmooth = atoi(argv[++i]);
				options->curvnsmooth = IMAX(0,options->curvnsmooth);
			} else {
				printf("-curv.nsmooth option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-segment.depth")) {
			if (i + 1 < argc) {
				options->segmentdepth = atoi(argv[++i]);
				options->segmentdepth = IMAX(0,options->segmentdepth);
			} else {
				printf("-segment.depth option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-segment.maxfaces")) {
			if (i + 1 < argc) {
				options->segmentmaxfaces = atoi(argv[++i]);
				options->segmentmaxfaces = IMAX(0,options->segmentmaxfaces);
			} else {
				printf("-segment.maxfaces option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if (arg == "-segment.loadfile") {
			if (i + 1 < argc) {
				options->segmentloadfile = argv[++i];
			} else {
				printf("-segment.loadfile option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if (arg == "-segment.savefile") {
			if (i + 1 < argc) {
				options->segmentsavefile = argv[++i];
			} else {
				printf("-segment.savefile option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else if ((arg == "-segment.maxlen")) {
			if (i + 1 < argc) {
				options->segmentmaxlen = (float)atof(argv[++i]);
				options->segmentmaxlen = abs(options->segmentmaxlen);
			} else {
				printf("-segment.maxlen option requires one argument.\n");
				returnOptions = false;
				break;
			}
		} else {
			printf("WARN: unrecognized argument '%s' will be ignored.\n",arg);
		}
	}
	if (options->curvstream) {
		options->curv = true;
		options->curvtensor = false;
	}
	if (options->curvtensor) {
		options->curv = true;
	}
	if (options->curv) {
		options->silhouette = false;
		options->axes = false;
		options->bbox = false;
		options->heatdiff = false;
	}
	if (options->texture) {
		options->axes = false;
		options->bbox = false;
		options->heatdiff = false;
		options->coloring_type_vert = 5;
		options->coloring_type_poly = 0;
	}
	if (options->heatdiff) {
		options->axes = false;
		options->bbox = false;
		options->coloring_type_vert = 3;
		options->coloring_type_poly = 0;
	}
	if (options->heatdiff && options->saddle) {
		options->coloring_type_vert = 4;
		options->coloring_type_poly = 0;
	}
	if (options->iterscheme < 0 || options->iterscheme > 1) {
		printf("-iterscheme supports only values 0 or 1.\n");
		returnOptions = false;
	}
	if (options->silhouette > 0) {
		//options->subdivide = 0;
		options->axes = false;
		options->bbox = false;
		if (options->silhouette > 2) {
			printf("-silhouette supports only options 1 or 2.\n");
			returnOptions = false;
		}
	}
	if (options->subdivide < 0 || options->subdivide > 2) {
		printf("-subdivide supports only options 1 or 2.\n");
		returnOptions = false;
	}
	if (!returnOptions) {
		delete options;
		return 0;
	}
	if (options->coloring_type_poly == 0 && options->coloring_type_vert == 0)
		options->coloring_type_poly = 1;
	return options;
}

void color_mapping(double percentage, double col[3])
{
	if (percentage == 0.0){
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 1.0;
	}
	else if (percentage <= 1.0/3){
		col[0] = 1.0-percentage*3.0;
		col[1] = 1.0;
		col[2] = 1.0-percentage*3.0;
	}
	else if (percentage <= 2.0/3){
		col[0] = percentage*3.0-1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
	else if (percentage <= 3.0/3){
		col[0] = 1.0;
		col[1] = 3.0-percentage*3.0;
		col[2] = 0.0;
	}
	else {
		col[0] = 1.0;
		col[1] = 1.0;
		col[2] = 0.0;
	}
}

/* Redundant functions copied over from learnply.cpp */
void sort(unsigned int *A, unsigned int *B, unsigned int *C, unsigned int sid, unsigned int eid){
	unsigned int i;
	unsigned int *tempA, *tempB, *tempC;
	unsigned int current1, current2, current0;

	if (sid>=eid)
		return;
	sort(A, B, C, sid, (sid+eid)/2);
	sort(A, B, C, (sid+eid)/2+1, eid);
	tempA = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempB = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	tempC = (unsigned int *)malloc(sizeof(unsigned int)*(eid-sid+1));
	for (i=0; i<eid-sid+1; i++){
		tempA[i] = A[i+sid];
		tempB[i] = B[i+sid];
		tempC[i] = C[i+sid];
	}
	current1 = sid;
	current2 = (sid+eid)/2+1;
	current0 = sid;
	while ((current1<=(sid+eid)/2) && (current2<=eid)){
		if (tempA[current1-sid] < tempA[current2-sid]) {
			A[current0] = tempA[current1-sid];
			B[current0] = tempB[current1-sid];
			C[current0] = tempC[current1-sid];
			current1++;		
		}
		else if (tempA[current1-sid] > tempA[current2-sid]){
			A[current0] = tempA[current2-sid];
			B[current0] = tempB[current2-sid];
			C[current0] = tempC[current2-sid];
			current2++;		
		}
		else {
			if (tempB[current1-sid] < tempB[current2-sid]) {
				A[current0] = tempA[current1-sid];
				B[current0] = tempB[current1-sid];
				C[current0] = tempC[current1-sid];
				current1++;		
			} else {
				A[current0] = tempA[current2-sid];
				B[current0] = tempB[current2-sid];
				C[current0] = tempC[current2-sid];
				current2++;		
			}
		}
		current0++;
	}
	if (current1<=(sid+eid)/2){
		for (i=current1; i<=(sid+eid)/2; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}
	if (current2<=eid){
		for (i=current2; i<=eid; i++){
			A[current0] = tempA[i-sid];
			B[current0] = tempB[i-sid];
			C[current0] = tempC[i-sid];
			current0++;
		}
	}

	free(tempA);
	free(tempB);
	free(tempC);
}

