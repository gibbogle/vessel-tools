#ifndef NETWORK_H
#define NETWORK_H

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct point_str
{
	float x,y,z;
	float d;
	bool used;
};
typedef point_str POINT;	
// A POINT has position (x,y,z), diameter (d), and usage flag (used)

struct vertex_str
{
	POINT point;
	int nlinks;
	int *edge;	// list of connected edges
	int pt[10];
	bool used;
	int netID;	// needs to be initialised to 0
	int point_index;
};
typedef vertex_str VERTEX;	
// A VERTEX has a POINT (point), number of linked segments (nlinks)

struct edge_str
{
	int npts, npts_used;
	int *pt;
	int *pt_used;
	int vert[2];
	float segavediam;
	float length_vox;	// voxel-voxel length
	bool used;
	int netID;	// connected network ID, initially 0
	float length_um;
};
typedef edge_str EDGE;	
// An EDGE has two VERTEX indices (vert[2]) to the array vertex[], 
// and a list of npts POINT indices (pt[]) to the array point[]

struct network_str
{
	int np;
	int nv;
	int ne;
	POINT *point;
	VERTEX *vertex;
	EDGE *edgeList;
};
typedef network_str NETWORK;

#define PI 3.14159
#define STR_LEN 128
#define NEMAX 1000
#define NBOX 500

float dist(NETWORK *net, int k1, int k2);
int WriteCmguiData(char *basename, NETWORK *net, float origin_shift[]);
int EdgeDimensions(EDGE *edges, POINT *points, int ne);
int CreateDistributions(NETWORK *net);
bool EqualPoints(POINT p1, POINT p2);
int WriteAmiraFile(char *amFileOut, char *amFileIn, NETWORK *net, float origin_shift[]);
int ReadAmiraFile(char *amFile, NETWORK *net);

#endif // NETWORK_H