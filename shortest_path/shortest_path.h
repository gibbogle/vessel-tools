#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct point_str
{
	float x,y,z;
	float d;
	float distance[2];
	bool used;
	int nb[2];
};
typedef point_str POINT;

struct vertex_str
{
	POINT point;
	int pt[10];
	bool used;
	float distance[2];	// [0] for artery, [1] for vein
	// Info for junctions, derived from edge[].vert[] and .segavediam
	int nlinks;
	int nb[2];
	int edge_num[10];
	float edge_diam[10];
};
typedef vertex_str VERTEX;

struct edge_str
{
	int npts, npts_used;
	int *pt;
	int *pt_used;
	int vert[2];
	float segavediam;
	float length_vox;	// voxel-voxel length
//	double length_jun;	// junction-junction length
	bool used;
	int dist10;
	int nb[2];
};
typedef edge_str EDGE;

struct shnode_str
{
	int id;
	float val;
	float diam;
	float dist;
};
typedef shnode_str SHNODE;
