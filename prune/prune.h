#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct point_str
{
	float x,y,z;
	float d;
	bool used;
};
typedef point_str POINT;

struct vertex_str
{
	POINT point;
	int nlinks;
	int pt[10];
	bool used;
};
typedef vertex_str VERTEX;

struct edge_str
{
	int npts, npts_used;
	int *pt;
	int *pt_used;
	int vert[2];
	float segavediam;
	float length_um;	// length in um
	bool used;
};
typedef edge_str EDGE;


