#ifndef QUANTIFY_H
#define QUANTIFY_H

#include <bitset>
#include <stdio.h>

#include "network.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define PI 3.14159
#define ALPHAMAX PI/12

struct segment_str {
	POINT end1, end2;
    double diam;
    double len;
};
typedef segment_str SEGMENT_TYPE;

#define STR_LEN 128
#define BIG 1.0e6

FILE *fperr=NULL;

int nxc, nyc, nzc, nx8, ny8, nz8, nbytes;
double voxelsize[3];
unsigned char *closedata = NULL;
int nsegments;
SEGMENT_TYPE *segment = NULL;
NETWORK *NP0 = NULL;

bool is_setup = false;

#endif
