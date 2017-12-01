#ifndef NETWORK2_H
#define NETWORK2_H

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MAXNBRS 30
#define MAXEDGES 1000000
#define STR_LEN 128

struct voxel_str {
	int pos[3];
	int initial_pos[3];
	int nbrs;
	int nbr[MAXNBRS];
	double diameter;
	int vertex_num; 
};
typedef voxel_str VOXEL;

struct edge_str
{
	int vert[2];
	int npts;
//	int npts_used;
	int *pt;
//	int *pt_used;
	double segavediam;
	double volume;
//	double segmindiam;
	float length_um;	// voxel-voxel length
	bool used;
//	bool ok;
};
typedef edge_str EDGE;

struct network_str
{
	int np;
	int nv;
	int ne;
//	POINT *point;
	VOXEL *point;
//	VERTEX *vertex;
	EDGE *edgeList;
};
typedef network_str NETWORK;

struct vertex_str
{
	int voxel_index;
	float diameter;
};
typedef vertex_str VERTEX;	


/*
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

float dist(NETWORK *net, int k1, int k2);
int EdgeDimensions(EDGE *edges, POINT *points, int ne);
int CreateDistributions(NETWORK *net);
bool EqualPoints(POINT p1, POINT p2);
*/
int WriteCmguiData(char *basename, NETWORK *net, float origin_shift[]);
//int WriteAmiraFile(char *amFileOut, char *amFileIn, NETWORK *net, float origin_shift[]);
//int ReadAmiraFile(char *amFile, NETWORK *net);

#endif // NETWORK2_H