#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct point_str
{
	float x,y,z;
	float d;
	bool used;
	int iedge;		// ADDED edge that the point is on
	int nearestpt;	// nearest point on another edge (actually pt index in net->point[])
	float d2near;	// distance^2 to nearest point
};
typedef point_str POINT;	
// A POINT has position (x,y,z), diameter (d), and usage flag (used)

struct vertex_str
{
	POINT point;
	int nlinks;		// number of linked edges
	int link[10];	// list of linked edges
	int pt;			// point index = vertex index because vertexes are the first entries in the point list
	bool used;
	int kvnearest;
	float dnearest;
//	int nf;			// added for fibredir etc.
//	int fib[10];	// added for fibredir
};
typedef vertex_str VERTEX;	
// A VERTEX kv has a POINT (point) with index pt = kv in net->point[], number of linked segments (nlinks)
// i.e. net->vertex[kv].point = net->point[kv]
// NOTE: therefore can remove struct member pt, in fact member point can also be removed.

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
	int jumpable_pt;	// pt index (0 - npts-1) of pt on the edge nearest to another edge
	float dnearest;		// distance to the nearest edge
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

