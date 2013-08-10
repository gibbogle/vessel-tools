// shunt.cpp
// Process a Spatialgraph file to detect early shunts

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

#include "shortest_path.h"

/*
struct pair_str
{
	int nd_from, nd_to;
};
typedef pair_str PAIR;
*/
extern "C" int dikbd_run(char *, int *, int *);

/*
struct loosend_str
{
	int iedge;			// edge
	int iend;			// free end (0, 1)
	double dir[3];
	double len;
	double dave;
	int joined;
};
typedef loosend_str LOOSEND;

struct conn_str
{
	int n;
	int v[5];
};
typedef conn_str CONN;
*/


int WriteCmguiData(char *basename);
float dist(int k1, int k2);

int nv, ne, np;
int nv_used, ne_used, np_used;
int vertex_case;	// 1 = artery, 2 = vein, 3 = artery + vein
bool use_nb;
int nbmax[2];
float artery_col[3], vein_col[3];
float distmin[2], distmax[2];
EDGE *edgeList;
VERTEX *vertex;
POINT *point;

FILE *fperr, *fpout;
FILE *exelem, *exnode;
char output_basename[128];
float ratio_limit;

#define NBOX 400
#define STR_LEN 128
#define PI 3.14159

/*
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float dist(int k1, int k2)
{
	float dx = point[k2].x - point[k1].x;
	float dy = point[k2].y - point[k1].y;
	float dz = point[k2].z - point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}
*/

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void showEdges(void)
{
	int i, k, j;
	EDGE edge;

	fprintf(fperr,"Edge points:\n");
	for (i=0; i<100; i++) {
		edge = edgeList[i];
		fprintf(fperr,"edge: %6d %6d\n",i,edge.npts);
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			fprintf(fperr,"%6d %6d  %6.1f %6.1f %6.1f\n",k,j,point[j].x,point[j].y,point[j].z);
		}
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int showEdge(int ie)
{
	int i, kv0, kv1;
	EDGE edge;

	printf("Edge: %6d Used: %d Points: ",ie,edgeList[ie].used);
	for (i=0; i<edgeList[ie].npts; i++) 
		printf("%6d ",edgeList[ie].pt[i]);
	printf("\n");
	kv0 = edgeList[ie].vert[0];
	kv1 = edgeList[ie].vert[1];
	if (kv0 != edgeList[ie].pt[0]) {
		printf("kv0 not pt[0]: %d %d\n",kv0,edgeList[ie].pt[0]);
		exit(1);
	}
	if (kv1 != edgeList[ie].pt[edgeList[ie].npts-1]) {
		printf("kv1 not pt[n-1]: %d %d\n",kv1,edgeList[ie].pt[edgeList[ie].npts-1]);
		exit(1);
	}
	for (i=0;i<ne;i++) {
		if (i == ie) continue;
		edge = edgeList[i];
		if (!edge.used) continue;
		if (edge.vert[0] == kv0 || edge.vert[1] == kv0)
			printf("  vert[0]: %d connected to edge: %d\n",kv0,i);
		if (edge.vert[0] == kv1 || edge.vert[1] == kv1)
			printf("  vert[1]: %d connected to edge: %d\n",kv1,i);
	}
	return 0;
}


//-----------------------------------------------------------------------------------------------------
// Write Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile(char *amFileOut, char *amFileIn)
{
	int i, k, j, npts, npts_used;
	EDGE edge;

	fprintf(fpout,"\nWriteAmiraFile: %s\n",amFileOut);
	npts = 0;
	npts_used = 0;
	for (i=0;i<ne;i++) {
		npts += edgeList[i].npts;
		npts_used += edgeList[i].npts_used;
	}

	FILE *fpam = fopen(amFileOut,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
	fprintf(fpam,"# Created by prune.exe from: %s\n",amFileIn);
	fprintf(fpam,"\n");
	fprintf(fpam,"define VERTEX %d\n",nv);
	fprintf(fpam,"define EDGE %d\n",ne);
	fprintf(fpam,"define POINT %d\n",npts);
	fprintf(fpam,"\n");
	fprintf(fpam,"Parameters {\n");
	fprintf(fpam,"    ContentType \"HxSpatialGraph\"\n");
	fprintf(fpam,"}\n");
	fprintf(fpam,"\n");
	fprintf(fpam,"VERTEX { float[3] VertexCoordinates } @1\n");
	fprintf(fpam,"EDGE { int[2] EdgeConnectivity } @2\n");
	fprintf(fpam,"EDGE { int NumEdgePoints } @3\n");
	fprintf(fpam,"POINT { float[3] EdgePointCoordinates } @4\n");
	fprintf(fpam,"POINT { float thickness } @5\n");
	fprintf(fpam,"\n");
	fprintf(fpam,"\n@1\n");

	for (i=0;i<nv;i++) {
		fprintf(fpam,"%6.1f %6.1f %6.1f\n",vertex[i].point.x,vertex[i].point.y,vertex[i].point.z);
	}
	fprintf(fpam,"\n@2\n");
	for (i=0;i<ne;i++) {
		edge = edgeList[i];
		fprintf(fpam,"%d %d\n",edge.vert[0],edge.vert[1]);
	}
	fprintf(fpam,"\n@3\n");
	for (i=0;i<ne;i++) {
		edge = edgeList[i];
		fprintf(fpam,"%d\n",edge.npts);
	}
	fprintf(fpam,"\n@4\n");
	for (i=0;i<ne;i++) {
		edge = edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
			fprintf(fpam,"%6.1f %6.1f %6.1f\n",point[j].x,point[j].y,point[j].z);
		}
	}
	fprintf(fpam,"\n@5\n");
	for (i=0;i<ne;i++) {
		edge = edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
			fprintf(fpam,"%6.2f\n",point[j].d);
		}
	}
	fclose(fpam);
	printf("Completed WriteAmiraFile\n");
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
// Note that Amira uses zero-based numbering
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile(char *amFile)
{
	int i, j, k, kp, npts, nvv, nee, npp;
	EDGE edge;
	char line[STR_LEN];

//	fprintf(fpout,"c Amira file: %s\n",amFile);
	FILE *fpam = fopen(amFile,"r");

	npts = 0;
	kp = 0;
	k = 0;
	while (k < 3) {
		fgets(line, STR_LEN, fpam);		// reads until newline character
		printf("%s\n",line);
		if (strncmp(line,"define VERTEX",13) == 0) {
			sscanf(line+13,"%d",&nv);
			nvv = nv+1;		// +1 to allow for 1-based distance[]
			k++;
		}
		if (strncmp(line,"define EDGE",11) == 0) {
			sscanf(line+11,"%d",&ne);
			nee = ne;
			k++;
		}
		if (strncmp(line,"define POINT",12) == 0) {
			sscanf(line+12,"%d",&np);
			npp = np;
			k++;
		}
	}

	vertex = (VERTEX *)malloc(nvv*sizeof(VERTEX));
	edgeList = (EDGE *)malloc(nee*sizeof(EDGE));
	point = (POINT *)malloc(npp*sizeof(POINT));	

	// Initialize
	for (i=0; i<nee; i++) {
		edgeList[i].used = false;
	}
	for (i=0; i<npp; i++) {
		point[i].used = false;
	}

	while (1) {
		if (fgets(line, STR_LEN, fpam) == NULL) {
			printf("Finished reading SpatialGraph file\n\n");
			fclose(fpam);
			break;
		}
		if (line[0] == '@') {
			sscanf(line+1,"%d",&k);
			if (k == 1) {
				for (i=0;i<nv;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @1\n");
						return 1;
					}
					sscanf(line,"%f %f %f\n",&vertex[i].point.x,&vertex[i].point.y,&vertex[i].point.z);
					kp = i;
					vertex[i].point.d = 0;
//					vertex[i].point.used = true;
					point[kp] = vertex[i].point;
				}
				kp++;
			} else if (k == 2) {
				for (i=0;i<ne;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @2\n");
						return 1;
					}
					sscanf(line,"%d %d",&edgeList[i].vert[0],&edgeList[i].vert[1]);
					edgeList[i].used = true;
				}
				printf("Got edge vertex indices\n");
			} else if (k == 3) {
				for (i=0;i<ne;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @3\n");
						return 1;
					}
					sscanf(line,"%d",&edgeList[i].npts);
					if (edgeList[i].npts < 1) {
						printf("ReadAmiraFile: i: %d npts: %d\n",i,edgeList[i].npts);
						return 1;
					}
					edgeList[i].npts_used = edgeList[i].npts;
					edgeList[i].pt = (int *)malloc(edgeList[i].npts*sizeof(int));
					edgeList[i].pt_used = (int *)malloc(edgeList[i].npts*sizeof(int));
					npts += edgeList[i].npts;
					edgeList[i].pt[0] = edgeList[i].vert[0];
					edgeList[i].pt[edgeList[i].npts-1] = edgeList[i].vert[1];
				}
				printf("Got edge npts, total: %d\n",npts);
			} else if (k == 4) {
				for (i=0;i<ne;i++) {
					edge = edgeList[i];
					float len = 0;
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
						if (k > 0 && k<edge.npts-1) {
							sscanf(line,"%f %f %f",&point[kp].x,&point[kp].y,&point[kp].z);
							edgeList[i].pt[k] = kp;
							edgeList[i].pt_used[k] = kp;
							kp++;
						}
						if (k > 0) {
							len = len + dist(edgeList[i].pt[k-1],edgeList[i].pt[k]);
						}
					}
					edgeList[i].length_vox = len;
				}
			} else if (k == 5) {
				for (i=0;i<ne;i++) {
					edge = edgeList[i];
					float dave = 0;
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @5\n");
							return 1;
						}
						j = edge.pt[k];
						sscanf(line,"%f",&point[j].d);
						if (point[j].d == 0) {
							printf("Error: ReadAmiraFile: zero diameter: i: %d npts: %d k: %d j: %d\n",i,edge.npts,k,j);
							return 1;
						}
						if (j < nv) {		// because the first nv points are vertices
							vertex[j].point.d = point[j].d;
						}
						dave += point[j].d;
						edgeList[i].segavediam = dave/edge.npts;
					}
				}
				printf("Got point thicknesses\n");
			}
		}
	}
	// Flag used points
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			point[j].used = true;
		}
	}
	fclose(fpam);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int CreateVertexLinks()
{
	int kv, ie, kv0, kv1, nlinks, nlmax;
	EDGE edge;

	nlmax = 0;
	for (kv=0; kv<nv; kv++) {
		nlinks = 0;
		for (ie=0; ie<ne; ie++) {
			kv0 = edgeList[ie].vert[0];
			kv1 = edgeList[ie].vert[1];
			if (kv == kv0 || kv == kv1) {
				vertex[kv].edge_num[nlinks] = ie;
				vertex[kv].edge_diam[nlinks] = edgeList[ie].segavediam;
				nlinks++;
			}
		}
		vertex[kv].nlinks = nlinks;
		nlmax = MAX(nlmax,nlinks);
//		if (nlinks == 0) {
//			printf("Error: CreateVertexLinks: nlinks: %d %d\n",kv,nlinks);
//			fprintf(fpout,"Error: CreateVertexLinks: nlinks: %d %d\n",kv,nlinks);
//		}
	}
	printf("Did CreateVertexLinks: max nlinks: %d\n",nlmax);
	fprintf(fpout,"Did CreateVertexLinks: max nlinks: %d\n",nlmax);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double dotproduct(double v1[], double v2[])
{
	double sum = 0;
	for (int i=0; i<3; i++) {
		sum += v1[i]*v2[i];
	}
	return sum;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double angle(double *v1, double *v2)
{
	double a1, a2, cosa;

	a1 = dotproduct(v1,v1);
	a2 = dotproduct(v2,v2);
	cosa = dotproduct(v1,v2)/sqrt(a1*a2);
	return acos(cosa);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double fang(double ang)
{
	return (1 + cos(ang))/2;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double fdist(double d)
{
	return 1/d;
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int inlist(int *list, int n, int k)
{
	int i;
	if (n == 0) return 0;
	for (i=0; i<n; i++) {
		if (list[i] == k) return i;
	}
	return 0;
}





//-----------------------------------------------------------------------------------------------------
// Reorder element and node identifiers
//-----------------------------------------------------------------------------------------------------
int squeezer(void)
{
	int i, j, k;
	int ne_x, nv_x, np_x, knew, i_x;
	EDGE edge;
	EDGE *edgeList_x;
	VERTEX *vertex_x;
	POINT *point_x;
	int *oldpt;

	printf("squeezer\n");
	vertex_x = (VERTEX *)malloc(nv*sizeof(VERTEX));
	edgeList_x = (EDGE *)malloc(2*ne*sizeof(EDGE));		// 2* for added joining edges
	point_x = (POINT *)malloc(2*np*sizeof(POINT));
	oldpt = (int *)malloc(2*np*sizeof(int));

	ne_x = 0;
	nv_x = 0;
	np_x = 0;
	// First process the vertices
	for (i=0; i<ne; i++) {
		if (!edgeList[i].used) continue;
		edge = edgeList_x[ne_x];
		edge.used = true;
		for (j=0; j<2; j++) {
			knew = inlist(oldpt,np_x,edgeList[i].vert[j]);
			if (knew == 0) {
				oldpt[np_x] = edgeList[i].vert[j];
				vertex_x[np_x].point = point[oldpt[np_x]];
				point_x[np_x] = point[oldpt[np_x]];
				knew = np_x;
				np_x++;
			}
			edge.vert[j] = knew;
		}
		edgeList_x[ne_x] = edge;
		ne_x++;
	}
	nv_x = np_x;

	// Now add the edge interior points to the list.  These occur once only.
	i_x = 0;
	for (i=0; i<ne; i++) {
//		printf("edge: %d used: %d npts:%d\n",i,edgeList[i].used,edgeList[i].npts);
		if (!edgeList[i].used) continue;
		int npts = edgeList[i].npts;
		if (npts < 1) {
			printf("squeezer: i: %d npts: %d\n",i,npts);
			return 1;
		}
		edgeList_x[i_x].pt = (int *)malloc(npts*sizeof(int));
		edge = edgeList_x[i_x];
		edge.npts = npts;
		edge.pt[0] = edge.vert[0];
		for (k=1; k<npts-1; k++) {
			j = edgeList[i].pt[k];
			oldpt[np_x] = j;
			point_x[np_x] = point[j];
			edge.pt[k] = np_x;
			np_x++;
		}
		edge.pt[npts-1] = edge.vert[1];
		edgeList_x[i_x] = edge;
		i_x++;
	}
	printf("Added interior edge points\n");
	printf("ne, ne_x: %d %d  nv, nv_x: %d %d  np, np_x: %d %d\n",ne,ne_x,nv,nv_x,np,np_x);

	// Now copy the revised data back into the original arrays
	nv = nv_x;
	for (i=0; i<nv; i++) {
		vertex[i] = vertex_x[i];
	}
	ne = ne_x;
	for (i=0; i<ne; i++) {
		edgeList[i] = edgeList_x[i];
	}
	np = np_x;
	for (i=0; i<np; i++) {
		point[i] = point_x[i];
	}
	free(vertex_x);
	free(edgeList_x); 
	free(point_x);
	free(oldpt);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// For edge i, computes the distance from one vertex to the other.  This distance is muliplied by 10
// and stored as an integer value edgeList[i].dist10
// The maximum node number is also determined
//-----------------------------------------------------------------------------------------------------
int getdist(int i, int *ndmax)
{
	EDGE edge;
	POINT p;
	int npts, k, j, kd, nd1, nd2;
	float d,x,y,z,xp,yp,zp,dx,dy,dz;

	edge = edgeList[i];
	if (!edge.used) {
		edgeList[i].dist10 = -1;
		 return 1;
	}
	npts = edge.npts;
	d = 0;
	for (k=0; k<npts; k++) {
		j = edge.pt[k];
		p = point[j];
		x = p.x;
		y = p.y;
		z = p.z;
		if (k > 0) {
			dx = x - xp;
			dy = y - yp;
			dz = z - zp;
			d += sqrt(dx*dx + dy*dy + dz*dz);
		}
		xp = x;
		yp = y;
		zp = z;
	}
	kd = 10*d;
	nd1 = edge.pt[0] + 1;		// Note "+ 1" for SP compatibility
	nd2 = edge.pt[npts-1] + 1;
	if (nd1 != nd2) {
//		fprintf(fpout,"a %6d %6d %4d\n",nd1,nd2,kd);
//		fprintf(fpout,"a %6d %6d %4d\n",nd2,nd1,kd);
		edgeList[i].dist10 = kd;
		if (nd1 > *ndmax) *ndmax = nd1;
		if (nd2 > *ndmax) *ndmax = nd2;
		return 0;
	} else {
		edgeList[i].dist10 = -1;
		return 1;
	}
}

//-----------------------------------------------------------------------------------------------------
// Creates the SP input file containing for each edge with non-negative .dist10 a pair of entries:
// nd1, nd1, kd
// nd2, nd1, kd
// where nd1 qnd nd2 are the vertex node numbers and kd = dist10 (= 10*edge length)
//-----------------------------------------------------------------------------------------------------
int	WriteSPfile(char *spfile, int iv0, int ndmax, int narcs)
{
	EDGE edge;
	int npts, i, kd, nd1, nd2;
	FILE *fpsp;

	fpsp = fopen(spfile,"w");
	fprintf(fpsp,"p sp %8d %8d\n",ndmax,narcs);
	fprintf(fpsp,"c\n");
	fprintf(fpsp,"n %8d\n",iv0);
	fprintf(fpsp,"c\n");
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		if (!edge.used) continue;
		kd = edge.dist10;
		if (kd < 0) continue;
		npts = edge.npts;
		nd1 = edge.pt[0] + 1;		// The SP code uses 1-based node numbering, like cmgui
		nd2 = edge.pt[npts-1] + 1;
		fprintf(fpsp,"a %6d %6d %4d\n",nd1,nd2,kd);
		fprintf(fpsp,"a %6d %6d %4d\n",nd2,nd1,kd);
	}
	fclose(fpsp);
	printf("Created SP file: %s\n",spfile);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Descending sort
//-----------------------------------------------------------------------------------------------------
int shcmp(const void *v1, const void *v2)
 {
	 SHNODE *vv1, *vv2;

	 vv1 = (SHNODE *)v1;
	 vv2 = (SHNODE *)v2;
	 double temp = vv1->val - vv2->val;
     if (temp > 0)
		 return -1;
	 else if (temp < 0)
		 return 1;
	 else
		 return 0;
 }

//-----------------------------------------------------------------------------------------------------
// Check out the kth entry in shlist[]
// First need to determine the number of edges at this vertex.
//-----------------------------------------------------------------------------------------------------
int explore(int i0)
{
	int i, nev, ev[8], npts;
	EDGE edge;

	nev = 0;
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		npts = edge.npts;
		if (edge.vert[0] == i0) {
			ev[nev] = i;
			nev++;
		} else if (edge.vert[1] == i0) {
			ev[nev] = -i;
			nev++;
		}
	}
	printf("vertex: %d nev: %d\n",i0,nev);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int removeEdge(int nd1, int nd2)
{
	int i;
	EDGE edge;

	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		if ((edge.vert[0] == nd1 && edge.vert[1] == nd2) 
		 || (edge.vert[1] == nd1 && edge.vert[0] == nd2)) {
			edgeList[i].used = false;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Returns the index of the edge with vertices (kv0,kv1)
//-----------------------------------------------------------------------------------------------------
int getedge(int kv0, int kv1, bool *ok)
{
	int nlinks, k, ie;
	float d, dmin, dmax;
	EDGE edge;

	dmin = 1.0e10;
	dmax = 0;
	nlinks = vertex[kv0].nlinks;
	for (k=0; k<nlinks; k++) {
		d = vertex[kv0].edge_diam[k];
		dmin = MIN(dmin,d);
		dmax = MAX(dmax,d);
	}
	if (ratio_limit == 0 || dmax/dmin < ratio_limit)
		*ok = true;
	else
		*ok = false;

	for (k=0; k<nlinks; k++) {
		ie = vertex[kv0].edge_num[k];
		edge = edgeList[ie];
		if (edge.vert[0] == kv0 && edge.vert[1] == kv1) return ie;
		if (edge.vert[0] == kv1 && edge.vert[1] == kv0) return ie;
	}
	printf("kv0 edges: %d\n",nlinks);
	nlinks = vertex[kv0].nlinks;
	for (k=0; k<nlinks; k++) {
		ie = vertex[kv0].edge_num[k];
		edge = edgeList[ie];
		printf("%d %d\n",edge.vert[0],edge.vert[1]);
	}
	return -1;
}

#define MAX_BRANCH 100
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void fixEdges(int ivtype)
{
	int ie, k, kv, nb, ip, kp;
	float d;
	EDGE edge;

	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (edge.nb[ivtype] != 0) continue;
		nb = 1000;
		d = 1.0e10;
		for (k=0; k<2; k++) {
			kv = edge.vert[k];
			nb = MIN(nb,vertex[kv].nb[ivtype]);
			d = MIN(d,vertex[kv].distance[ivtype]);
		}
		if (nb == 0) {
			printf("fixEdges: nb=0: ie: %d vert: %d %d\n",ie,edge.vert[0],edge.vert[1]);
			exit(1);
		}
		edgeList[ie].nb[ivtype] = MIN(nb,MAX_BRANCH-1);
//		edgeList[ie].distance[ivtype] = d);
		for (ip=1; ip<edgeList[ie].npts; ip++) {
			kp = edgeList[ie].pt[ip];
			point[kp].nb[ivtype] = nb;
			point[kp].distance[ivtype] = d;
		}
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int quantify(int *path, int ivtype, int vertex_no)
{
	int i, kv, kp, ie, ie0, ip, n, nmax, nb, totcnt;
	int *vlist, *elist;
	int branch_cnt[MAX_BRANCH];
	bool ok;
	float d;

	for (i=0; i<MAX_BRANCH; i++) {
		branch_cnt[i] = 0;
	}
	distmin[ivtype] = 1.0e10;
	distmax[ivtype] = 0;
	nmax = 0;
	nbmax[ivtype] = 0;
	totcnt = 0;
	fprintf(fpout,"ivtype, vertex_no, d: %d %d %f\n",ivtype,vertex_no,vertex[vertex_no].distance[ivtype]);
	for (kv=0; kv<nv; kv++) {
		if (kv == vertex_no) continue;
		if (vertex[kv].nlinks == 0) {
//			printf("vertex with nlinks = 0: kv: %d\n",kv);
			continue;
		}
//		printf("kv: %d\n",kv);
		d = vertex[kv].distance[ivtype];
		distmin[ivtype] = MIN(distmin[ivtype],d);
		distmax[ivtype] = MAX(distmax[ivtype],d);
		kp = kv;
		n = 1;
		for (;;) {
			kp = path[kp];
			n++;
			if (kp == vertex_no) break;
		}
		nmax = MAX(nmax,n);
		vlist = (int *)malloc(n*sizeof(int));
		elist = (int *)malloc(n*sizeof(int));
		kp = kv;
		vlist[0] = kp;
		n = 1;
		nb = 0;
		for (;;) {
			kp = path[kp];
			vlist[n] = kp;
			ie = getedge(vlist[n-1],vlist[n],&ok);
			if (ie < 0) {
				printf("Error: quantify: edge not found: %d %d\n",vlist[n-1],vlist[n]);
				fprintf(fpout,"Error: quantify: edge not found: %d %d\n",vlist[n-1],vlist[n]);
				return 1;
			}
//			if (!ok) {
//				printf("Failed to find edge: %d: %d %d\n",n,vlist[n-1],vlist[n]);
//				return 1;
//			}
			if (nb == 0) ie0 = ie;	// first edge
			if (ok) nb++;
			elist[n-1] = ie;
			n++;
			if (kp == vertex_no) break;
		}

		free(vlist);
		free(elist);
		nb = MIN(nb,MAX_BRANCH-1);
		nbmax[ivtype] = MAX(nbmax[ivtype],nb);
		branch_cnt[nb] = branch_cnt[nb] + 1;
		totcnt++;
		if (nb == 0) {
			printf("quantify: edge with nb=0: %d\n",ie0);
		}
		edgeList[ie0].nb[ivtype] = nb;
		point[edgeList[ie0].vert[0]].nb[ivtype] = nb;
		point[edgeList[ie0].vert[1]].nb[ivtype] = nb;
		point[edgeList[ie0].vert[0]].distance[ivtype] = d;
		point[edgeList[ie0].vert[1]].distance[ivtype] = d;
		vertex[kv].nb[ivtype] = nb;
		vertex[kv].distance[ivtype] = d;
		for (ip=1; ip<edgeList[ie0].npts; ip++) {
			kp = edgeList[ie0].pt[ip];
			point[kp].nb[ivtype] = nb;
			point[kp].distance[ivtype] = d;
		}
	}
	fixEdges(ivtype);
	printf("Min, max distance: %f %f\n",distmin[ivtype],distmax[ivtype]);
	fprintf(fpout,"Min, max distance: %f %f\n",distmin[ivtype],distmax[ivtype]);
	printf("quantify: ratio_limit: %f max branches: %d\n",ratio_limit,nbmax[ivtype]);
	fprintf(fpout,"quantify: ratio_limit: %f max branches: %d\n",ratio_limit,nbmax[ivtype]);
	fprintf(fpout,"Distribution of number of branches\n");
	for (i=0; i<MAX_BRANCH; i++) {
		fprintf(fpout,"%4d %6d %8.4f\n",i,branch_cnt[i],float(branch_cnt[i])/totcnt);
		if (branch_cnt[i] == 0 && i > 10) break;
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err, i, k, i0;
	char *input_amfile;
	char spfile[] = "sp_arc.dat";
	char drive[32], dir[128],filename[1024], ext[32], outfilename[1024];
	char errfilename[1024], output_amfile[1024], result_file[128];
	int artery_vertex_no, vein_vertex_no, use_branch;
	int ivtype;		// 0 = artery, 1 = vein
//	int artery_vertex, vein_vertex;
	int *distance10 = NULL;
//	bool artery = true;
//	bool vein = false;
	int *path;

	if (argc != 13) {
		printf("Usage: short input_amfile output_file artery_vertex aR aG aB vein_vertex vR vG vB use_branch ratio_limit\n");
		printf("       artery_vertex is the main artery vertex number (1-based)\n");
		printf("       aR aG aB are the RGB colour components for artery distance (0-1)\n");
		printf("       vein_vertex is the main vein vertex number (1-based)\n");
		printf("       vR vG vB are the RGB colour components for vein distance (0-1)\n");
		printf("       use_branch = 1 if the colouring is to be based on # of branches, = 0 if based on distance\n");
		printf("       ratio_limit is the threshold value of max_diam/min_diam at a branch\n");
		fperr = fopen("shunt_error.log","w");
		fprintf(fperr,"Usage: short input_amfile output_file artery_vertex aR aG aB vein_vertex vR vG vB use_branch ratio_limit\n");
		fprintf(fperr,"       artery_vertex is the main artery vertex number (1-based)\n");
		fprintf(fperr,"       aR aG aB are the RGB colour components for artery distance (0-1)\n");
		fprintf(fperr,"       vein_vertex is the main vein vertex number (1-based)\n");
		fprintf(fperr,"       vR vG vB are the RGB colour components for vein distance (0-1)\n");
		fprintf(fperr,"       use_branch = 1 if the colouring is to be based on # of branches, = 0 if based on distance\n");
		fprintf(fperr,"       ratio_limit is the threshold value of max_diam/min_diam at a branch\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	sscanf(argv[3],"%d",&artery_vertex_no);
	sscanf(argv[4],"%f",&artery_col[0]);
	sscanf(argv[5],"%f",&artery_col[1]);
	sscanf(argv[6],"%f",&artery_col[2]);
	sscanf(argv[7],"%d",&vein_vertex_no);
	sscanf(argv[8],"%f",&vein_col[0]);
	sscanf(argv[9],"%f",&vein_col[1]);
	sscanf(argv[10],"%f",&vein_col[2]);
	sscanf(argv[11],"%d",&use_branch);
	sscanf(argv[12],"%f",&ratio_limit);

	strcpy(outfilename,argv[2]);
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);

	sprintf(errfilename,"%s_short.log",output_basename);
	fperr = fopen(errfilename,"w");
	fpout = fopen(outfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	use_nb = (use_branch == 1);
	err = ReadAmiraFile(input_amfile);
	if (err != 0) return 2;
	err = CreateVertexLinks();
	if (err != 0) return 3;

	int ndmax = 0;
	int narcs = 0;
	for (int ie=0; ie<ne; ie++) {
		if (getdist(ie, &ndmax) == 0)
			narcs += 2;
	}
	printf("\nndmax: %d nv: %d\n",ndmax,nv);
	fprintf(fpout,"\nndmax: %d nv: %d\n",ndmax,nv);
	path = (int *)malloc(ndmax*sizeof(int));
	distance10 = (int *)malloc((ndmax+1)*sizeof(int));

	nbmax[0] = 0;
	nbmax[1] = 0;
	vertex_case = 0;
	if (artery_vertex_no != 0) {
		ivtype = 0;
		vertex_case += 1;
		WriteSPfile(spfile,artery_vertex_no,ndmax,narcs);	// shortest-path code is 1-based
		dikbd_run(spfile, distance10, path);
		for (i=1; i<=ndmax; i++) {
			vertex[i-1].distance[ivtype] = distance10[i]/10.0;
		}
		err = quantify(path, ivtype, artery_vertex_no-1);
		if (err != 0) return 4;
	}
	if (vein_vertex_no != 0) {
		ivtype = 1;
		vertex_case += 2;
		WriteSPfile(spfile,vein_vertex_no,ndmax,narcs);	// shortest-path code is 1-based
		dikbd_run(spfile, distance10, path);
		for (i=1; i<=ndmax; i++) {
			vertex[i-1].distance[ivtype] = distance10[i]/10.0;
		}
		err = quantify(path, ivtype, vein_vertex_no-1);
		if (err != 0) return 4;
	}
	fprintf(fpout,"nbmax: %3d %3d\n",nbmax[0],nbmax[1]);
	err = WriteCmguiData(output_basename);
	if (err != 0) {
		fprintf(fpout,"Error writing CMGUI files: %d\n",err);
		return 5;
	}

	fclose(fpout);

	return 0;

	/*
	if (vein) {
		// Vein distances
		WriteSPfile(spfile,vein_vertex,ndmax,narcs);
		dikbd_run(spfile, distance, path);
		for (i=1; i<=ndmax; i++) {
			vertex[i-1].distance[1] = distance[i]/10.;
		}
	}
	*/
	SHNODE *shlist = (SHNODE *)malloc((ndmax+1)*sizeof(SHNODE));

/*
	float val;
	for (i=1; i<=ndmax; i++) {
		k = i-1;
		fprintf (fpout,"%8d  %8.1f\n", k, vertex[k].distance);
		shlist[k].id = i;
		shlist[k].val = vertex[k].distance;
		shlist[k].diam = point[i].d;
		shlist[k].dist = vertex[k].distance;
//		shlist[k].dist[1] = vertex[k].distance[1];
	}
*/
	// Sort the array into ascending order.

    qsort(shlist, ndmax, sizeof(shlist[0]), shcmp);

	// The following code was an unsuccessful attempt to invent a way to find shunts
    // Display the sorted array.
	// Need to explore the node with highest val, to determine if it is indeed a shunt, then
	// seek a way to deshunt automatically.
	// A node with 4 links is probably the intersection of 2 vessels, should be able to determine the
	// links that are part of a vessel from the link directions.

	 /*
     puts("\nPress a key to continue.");
     getch();

     // Enter a search key.

     printf("Enter a value to search for: ");
     scanf("%d", &key);

     // Perform the search.

     ptr = (int *)bsearch(&key, arr, MAX, sizeof(arr[0]),intcmp);

     if ( ptr != NULL )
         printf("%d found at arr[%d].", key, (ptr - arr));
     else
         printf("%d not found.", key);
	 */


	// A potential spurious shunt is a vertex with a large diameter, especially if
	// the sum of the distances is small.

	/*
	err = WriteAmiraFile(output_amfile,input_amfile);
	if (err != 0) return 7;
	err = CreateDistributions();
	if (err != 0) return 8;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename);
		if (err != 0) return 9;
	}
	*/
	return 0;
}