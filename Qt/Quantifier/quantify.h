#include <bitset>
#include <stdio.h>

#include "network.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct segment_str {
	POINT end1, end2;
};
typedef segment_str SEGMENT_TYPE;

#define STR_LEN 128
#define BIG 1.0e6

FILE *fperr=NULL, *fpout=NULL;

int nxc, nyc, nzc, nx8, ny8, nz8, nbytes;
float voxelsize[3];
unsigned char *closedata = NULL;
int nsegments;
SEGMENT_TYPE *segment = NULL;
NETWORK *NP0 = NULL;

bool is_setup = false;
