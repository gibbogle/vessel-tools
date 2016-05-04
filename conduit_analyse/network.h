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
	int pt[10];
	bool used;
	int nf;			// added for fibredir
	int fib[10];	// added for fibredir
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
	float length_um;	// voxel-voxel length
//	double length_jun;	// junction-junction length
	bool used;
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

