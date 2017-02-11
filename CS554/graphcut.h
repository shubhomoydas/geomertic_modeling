#include <stdlib.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include "mathfuncs.h"
#include "CornerPoints.h"
#include "utils.h"
#include "segment.h"

#ifndef _CS554_GRAPHCUT_H_
#define _CS554_GRAPHCUT_H_

int bi_partition_with_graph_cut(
	Polyhedron *poly, tris_tensor_list *ttl, segments *segs, CommandOptions *options);

#endif

