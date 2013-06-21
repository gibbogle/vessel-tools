// Create an Amira spatialgraph file by joining pieces (blocks).
// Cutting planes are YZ, XZ, XY
// Each block is made up of all the vessel segments with a vertex within the block ranges.

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

#include "network.h"

#define STR_LEN 128
#define BIG 1.0e6

// These are guessed values for the maximum array sizes - may need to be increased
#define MAX_VERTEX 100000
#define MAX_EDGE 100000
#define MAX_POINT 1000000

#define EPSILON 0.051
#define EQU(a,b) (fabs((a)-(b)) < EPSILON)
#define EQUIV(a,b) (fabs((a)[0]-(b)[0]) < EPSILON && fabs((a)[1]-(b)[1]) < EPSILON && fabs((a)[2]-(b)[2]) < EPSILON)

FILE *fperr, *fpout;
EDGE *elist;
int *ivlist;

bool dbug = false;

int ivlistAdd(int iv, int *ivlist, int *nv);
int WriteCmguiData(char *basename, NETWORK *net);

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float dist(NETWORK *net, int k1, int k2)
{
	float dx = net->point[k2].x - net->point[k1].x;
	float dy = net->point[k2].y - net->point[k1].y;
	float dz = net->point[k2].z - net->point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}


//-----------------------------------------------------------------------------------------------------
// We hope that MAX_VERTEX, MAX_EDGE and MAX_POINT are all big enough.
//-----------------------------------------------------------------------------------------------------
int InitNetwork(NETWORK *net)
{
	net->vertex = (VERTEX *)malloc(MAX_VERTEX*sizeof(VERTEX));
	net->edgeList = (EDGE *)malloc(MAX_EDGE*sizeof(EDGE));
	net->point = (POINT *)malloc(MAX_POINT*sizeof(POINT));
	printf("Allocated arrays: %d %d %d\n",MAX_VERTEX,MAX_EDGE,MAX_POINT);
	net->nv = 0;
	net->ne = 0;
	net->np = 0;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Vessels of network net1 are added to network net0
// A check is made for duplicate edges by comparing the vertex positions.  If a vessel is already in the 
// edgeList nothing is added, otherwise new vertex number(s) is(are) generated for the added vertex(es), 
// and the points are all added.
//-----------------------------------------------------------------------------------------------------
int AddNetwork(NETWORK *net0, NETWORK *net1)
{
	int ie1, ie0, new_kv0, new_kv1, add_kv0, add_kv1, vmap[2], edir, ne0, nv1, new_ne;
	bool ehit, vhit00, vhit01, vhit10, vhit11;
	float new_v0[3], new_v1[3];
	float old_v0[3], old_v1[3];
	EDGE edge0, edge1;

	ne0 = net0->ne;
	new_ne = 0;
	printf("net0: ne, nv, np: %d %d %d\n",ne0,net0->nv,net0->np);
	fprintf(fperr,"net0: ne, nv, np: %d %d %d\n",ne0,net0->nv,net0->np);
	ivlist = (int *)malloc(2*MAX_VERTEX*sizeof(int));
	nv1 = 0;
	for (ie1=0; ie1<net1->ne; ie1++) {
		// Look for vertex matches
		edge1 = net1->edgeList[ie1];
		new_kv0 = edge1.vert[0];
		new_kv1 = edge1.vert[1];
		if (dbug) {
			printf("\nnew edge: %d vertex: %d %d\n",ie1,new_kv0,new_kv1);
			fprintf(fperr,"\nnew edge: %d vertex: %d %d\n",ie1,new_kv0,new_kv1);
		}
		ehit = false;	// true if new_v0 = old_v0, new_v1 = old_v1 (edir = 1), or  new_v0 = old_v1, new_v1 = old_v0 (edir = -1)
		vmap[0] = vmap[1] = -1;
		new_v0[0] = net1->vertex[edge1.vert[0]].point.x;
		new_v0[1] = net1->vertex[edge1.vert[0]].point.y;
		new_v0[2] = net1->vertex[edge1.vert[0]].point.z;
		new_v1[0] = net1->vertex[edge1.vert[1]].point.x;
		new_v1[1] = net1->vertex[edge1.vert[1]].point.y;
		new_v1[2] = net1->vertex[edge1.vert[1]].point.z;
		for (ie0=0; ie0<ne0; ie0++) {
			edge0 = net0->edgeList[ie0];
			old_v0[0] = net0->vertex[edge0.vert[0]].point.x;
			old_v0[1] = net0->vertex[edge0.vert[0]].point.y;
			old_v0[2] = net0->vertex[edge0.vert[0]].point.z;
			old_v1[0] = net0->vertex[edge0.vert[1]].point.x;
			old_v1[1] = net0->vertex[edge0.vert[1]].point.y;
			old_v1[2] = net0->vertex[edge0.vert[1]].point.z;
			vhit00 = vhit01 = vhit10 = vhit11 = false;
//			if (new_v0[0]==old_v0[0]&&new_v0[1]==old_v0[1]&&new_v0[2]==old_v0[2]) {
//			if (EQU(new_v0[0],old_v0[0])&&EQU(new_v0[1],old_v0[1])&&EQU(new_v0[2],old_v0[2])) {
			if (EQUIV(new_v0,old_v0)) {
				vhit00 = true;
				vmap[0] = edge0.vert[0];
//				printf("vhit00: %d\n",vmap[0]);
//			} else if (new_v0[0]==old_v1[0]&&new_v0[1]==old_v1[1]&&new_v0[2]==old_v1[2]) {
//			} else if (EQU(new_v0[0],old_v1[0])&&EQU(new_v0[1],old_v1[1])&&EQU(new_v0[2],old_v1[2])) {
			} else if (EQUIV(new_v0,old_v1)) {
				vhit01 = true;
				vmap[0] = edge0.vert[1];
//				printf("vhit01: %d\n",vmap[0]);
			}
//			if (new_v1[0]==old_v0[0]&&new_v1[1]==old_v0[1]&&new_v1[2]==old_v0[2]) {
//			if (EQU(new_v1[0],old_v0[0])&&EQU(new_v1[1],old_v0[1])&&EQU(new_v1[2],old_v0[2])) {
			if (EQUIV(new_v1,old_v0)) {
				vhit10 = true;
				vmap[1] = edge0.vert[0];
//				printf("vhit10: %d\n",vmap[1]);
//			} else if (new_v1[0]==old_v1[0]&&new_v1[1]==old_v1[1]&&new_v1[2]==old_v1[2]) {
//			} else if (EQU(new_v1[0],old_v1[0])&&EQU(new_v1[1],old_v1[1])&&EQU(new_v1[2],old_v1[2])) {
			} else if (EQUIV(new_v1,old_v1)) {
				vhit11 = true;
				vmap[1] = edge0.vert[1];
//				printf("vhit11: %d\n",vmap[1]);
			}
			if (vhit00 && vhit11) {
				ehit = true;
				edir = 1;
			}
			if (vhit01 && vhit10) {
				ehit = true;
				edir = -1;
			}
		}
		if (!ehit) {		// add the vessel edge1
			if (dbug) {
				printf("No hit: vmap: %d %d\n",vmap[0],vmap[1]);
				fprintf(fperr,"No hit: vmap: %d %d\n",vmap[0],vmap[1]);
			}
			if (vmap[0] >= 0 && vmap[1] >= 0) {
//				printf("AddNetwork: !ehit vmap: %d %d\n",vmap[0],vmap[1]);
//				return 1;
			}
			ie0 = net0->ne;	// new edge #
			net0->ne++;
			if (vmap[0] < 0) {
				add_kv0 = net0->nv + ivlistAdd(new_kv0,ivlist,&nv1);
			} else {
				add_kv0 = vmap[0];	// vertex in net0
			}
			if (vmap[1] < 0) {
				add_kv1 = net0->nv + ivlistAdd(new_kv1,ivlist,&nv1);
			} else {
				add_kv1 = vmap[1];	// vertex in net0
			}
			if (dbug) {
				printf("add_kv0, add_kv1: %d %d\n",add_kv0,add_kv1);
				fprintf(fperr,"add_kv0, add_kv1: %d %d\n",add_kv0,add_kv1);
			}
			net0->edgeList[ie0].vert[0] = add_kv0;
			net0->vertex[add_kv0].point = net1->vertex[new_kv0].point;
			net0->edgeList[ie0].vert[1] = add_kv1;
			net0->vertex[add_kv1].point = net1->vertex[new_kv1].point;
			net0->edgeList[ie0].npts = edge1.npts;
			net0->edgeList[ie0].used = true;
			net0->edgeList[ie0].pt = (int *)malloc(net0->edgeList[ie0].npts*sizeof(int));
			if (dbug) {
				printf("Add pts: %d\n",edge1.npts);
				fprintf(fperr,"Add pts: %d %d\n",net0->np,edge1.npts);
			}
			int kp = net0->np;
			for (int i=0; i<edge1.npts; i++) {
				net0->edgeList[ie0].pt[i] = edge1.pt[i];
				net0->point[kp] = net1->point[edge1.pt[i]];
				net0->point[kp].used = true;
				kp++;
			}
			net0->np = kp;
		} else {
//			printf("Hit\n");
//			fprintf(fperr,"Hit\n");
		}
	}
	net0->nv += nv1;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Write Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile(char *amFileOut, char *amFileIn, NETWORK *net)
{
	int i, k, j, npts;
	EDGE edge;

	fprintf(fpout,"\nWriteAmiraFile: %s\n",amFileOut);
	npts = 0;
	for (i=0;i<net->ne;i++) {
		npts += net->edgeList[i].npts;
	}

	FILE *fpam = fopen(amFileOut,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
	fprintf(fpam,"# Created by zoom.exe from: %s\n",amFileIn);
	fprintf(fpam,"\n");
	fprintf(fpam,"define VERTEX %d\n",net->nv);
	fprintf(fpam,"define EDGE %d\n",net->ne);
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

	for (i=0;i<net->nv;i++) {
//		fprintf(fpam,"%6.1f %6.1f %6.1f\n",net->vertex[i].point.x,net->vertex[i].point.y,net->vertex[i].point.z);
		fprintf(fpam,"%6.1f %6.1f %6.1f\n",
			net->vertex[i].point.x,
			net->vertex[i].point.y,
			net->vertex[i].point.z);
	}
	fprintf(fpam,"\n@2\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		fprintf(fpam,"%d %d\n",edge.vert[0],edge.vert[1]);
	}
	fprintf(fpam,"\n@3\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		fprintf(fpam,"%d\n",edge.npts);
	}
	fprintf(fpam,"\n@4\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
//			fprintf(fpam,"%6.1f %6.1f %6.1f\n",net->point[j].x,net->point[j].y,net->point[j].z);
			fprintf(fpam,"%6.1f %6.1f %6.1f\n",
				net->point[j].x,
				net->point[j].y,
				net->point[j].z);
		}
	}
	fprintf(fpam,"\n@5\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
			fprintf(fpam,"%6.2f\n",net->point[j].d);
		}
	}
	fclose(fpam);
	printf("Completed WriteAmiraFile: %s\n",amFileOut);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile(char *amFile, NETWORK *net)
{
	int i, j, k, kp, npts;
	int np_used, ne_used;
	EDGE edge;
	char line[STR_LEN];

	printf("ReadAmiraFile: %s\n",amFile);
	fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
	FILE *fpam = fopen(amFile,"r");

	npts = 0;
	kp = 0;
	k = 0;
	while (k < 3) {
		fgets(line, STR_LEN, fpam);		// reads until newline character
//		printf("%s\n",line);
		if (strncmp(line,"define VERTEX",13) == 0) {
			sscanf(line+13,"%d",&net->nv);
			k++;
		}
		if (strncmp(line,"define EDGE",11) == 0) {
			sscanf(line+11,"%d",&net->ne);
			k++;
		}
		if (strncmp(line,"define POINT",12) == 0) {
			sscanf(line+12,"%d",&net->np);
			k++;
		}
	}

	net->vertex = (VERTEX *)malloc(net->nv*sizeof(VERTEX));
	net->edgeList = (EDGE *)malloc(net->ne*sizeof(EDGE));
	net->point = (POINT *)malloc(net->np*sizeof(POINT));
	printf("Allocated arrays: ne,nv,np: %d %d %d\n",net->ne,net->nv,net->np);
	fprintf(fperr,"Allocated arrays: ne,nv,np: %d %d %d\n",net->ne,net->nv,net->np);

	// Initialize
	for (i=0; i<net->ne; i++) {
		net->edgeList[i].used = false;
	}
	for (i=0; i<net->np; i++) {
		net->point[i].used = false;
	}
//	printf("Initialised\n");

	while (1) {
		if (fgets(line, STR_LEN, fpam) == NULL) {
			printf("Finished reading SpatialGraph file\n\n");
			fclose(fpam);
			break;
		}
		if (line[0] == '@') {
			sscanf(line+1,"%d",&k);
			if (k == 1) {
				for (i=0;i<net->nv;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @1\n");
						return 1;
					}
					sscanf(line,"%f %f %f\n",&(net->vertex[i].point.x),&(net->vertex[i].point.y),&(net->vertex[i].point.z));
					kp = i;
					net->vertex[i].point.d = 0;
					net->point[kp] = net->vertex[i].point;
				}
				kp++;
//				printf("Got vertices\n");
			} else if (k == 2) {
				for (i=0;i<net->ne;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @2\n");
						return 1;
					}
					sscanf(line,"%d %d",&net->edgeList[i].vert[0],&net->edgeList[i].vert[1]);
					net->edgeList[i].used = true;
				}
//				printf("Got edge vertex indices\n");
			} else if (k == 3) {
				for (i=0;i<net->ne;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @3\n");
						return 1;
					}
					sscanf(line,"%d",&net->edgeList[i].npts);
					if (net->edgeList[i].npts < 1) {
						printf("ReadAmiraFile: i: %d npts: %d\n",i,net->edgeList[i].npts);
						return 1;
					}
					net->edgeList[i].npts_used = net->edgeList[i].npts;
					net->edgeList[i].pt = (int *)malloc(net->edgeList[i].npts*sizeof(int));
					net->edgeList[i].pt_used = (int *)malloc(net->edgeList[i].npts*sizeof(int));
					npts += net->edgeList[i].npts;
					net->edgeList[i].pt[0] = net->edgeList[i].vert[0];
					net->edgeList[i].pt[net->edgeList[i].npts-1] = net->edgeList[i].vert[1];
				}
//				printf("Got edge npts, total: %d\n",npts);
			} else if (k == 4) {
				for (i=0;i<net->ne;i++) {
					edge = net->edgeList[i];
					float len = 0;
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
						if (k > 0 && k<edge.npts-1) {
							sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[k] = kp;
							net->edgeList[i].pt_used[k] = kp;
							kp++;
						}
						if (k > 0) {
							len = len + dist(net,net->edgeList[i].pt[k-1],net->edgeList[i].pt[k]);
						}
					}
					net->edgeList[i].length_vox = len;
				}
			} else if (k == 5) {
				for (i=0;i<net->ne;i++) {
					edge = net->edgeList[i];
					float dave = 0;
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @5\n");
							return 1;
						}
						j = edge.pt[k];
						sscanf(line,"%f",&net->point[j].d);
						if (net->point[j].d == 0) {
							printf("Error: ReadAmiraFile: zero diameter: i: %d npts: %d k: %d j: %d\n",i,edge.npts,k,j);
							return 1;
						}
						if (j < net->nv) {		// because the first nv points are vertices
							net->vertex[j].point.d = net->point[j].d;
						}
						dave += net->point[j].d;
						net->edgeList[i].segavediam = dave/edge.npts;
					}
				}
//				printf("Got point thicknesses\n");
			}
		}
	}
	// Flag used points
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			net->point[j].used = true;
		}
	}
	fclose(fpam);
	np_used = 0;
	for (j=0; j<net->np; j++) {
		if (net->point[j].used) np_used++;
	}
//	printf("Points: np: %d np_used: %d\n",net->np,np_used);
	ne_used = 0;
	for (j=0; j<net->ne; j++) {
		if (net->edgeList[j].used) ne_used++;
	}
//	printf("Edges: ne: %d ne_used: %d\n",net->ne,ne_used);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void free_network(NETWORK *net)
{
//	free(elist);
	free(ivlist);
	free(net->vertex);
	free(net->edgeList);
	free(net->point);
}

//-----------------------------------------------------------------------------------------------------
// Check for vertex index iv in ivlist.  If it exists, return the index. 
// Otherwise add it to the list, increment nv, return the index.
//-----------------------------------------------------------------------------------------------------
int ivlistAdd(int iv, int *ivlist, int *nv)
{
	if (*nv == 0) {
		ivlist[0] = iv;
		*nv = 1;
		return 0;
	}
	for (int i=0; i<*nv; i++) {
		if (iv == ivlist[i]) {
			return i;
		}
	}
	ivlist[*nv] = iv;
	(*nv)++;
	return *nv-1;
}

//-----------------------------------------------------------------------------------------------------
// Add the vessels of net0 to net1, avoiding duplicates
//-----------------------------------------------------------------------------------------------------
int PrevAddNetwork(NETWORK *net0, NETWORK *net1)
{
	int i, j, ie, iv, kp, kv, ne, nv;
	float xv, yv, zv;
	float x1=0, x2=0, y1=0, y2=0, z1=0, z2=0;
	bool in;
	VERTEX v;

	ne = 0;
	for (i=0; i<net0->ne; i++) {
		in = false;
		for (j=0; j<2; j++) {
			iv = net0->edgeList[i].vert[j];
			v = net0->vertex[iv];
			xv = v.point.x;
			yv = v.point.y;
			zv = v.point.z;
			if (x1 <= xv && xv <= x2 && y1 <= yv && yv <= y2 && z1 <= zv && zv <= z2) {
				in = true;
			}
		}
		if (in) ne++;
	}
	printf("ne: %d\n",ne);
	// Create a list of edges, and a list of vertices
	// The index number of a vertex in the list replaces vert[] in the elist entry.
	elist = (EDGE *)malloc(ne*sizeof(EDGE));
	ivlist = (int *)malloc(2*ne*sizeof(int));
	nv = 0;
	ne = 0;
	for (i=0; i<net0->ne; i++) {
		in = false;
		for (j=0; j<2; j++) {
			iv = net0->edgeList[i].vert[j];
			v = net0->vertex[iv];
			xv = v.point.x;
			yv = v.point.y;
			zv = v.point.z;
			if (x1 <= xv && xv <= x2 && y1 <= yv && yv <= y2 && z1 <= zv && zv <= z2) {
				in = true;
			}
		}
		if (in) {
			elist[ne] = net0->edgeList[i];
			for (j=0; j<2; j++) {
				iv = net0->edgeList[i].vert[j];
				kv = ivlistAdd(iv,ivlist,&nv);
				elist[ne].vert[j] = kv;
			}
			ne++;
		}
	}
	printf("ne: %d  nv: %d\n",ne,nv);
	// Now from this list of edges we need to create the network net1.
	// Note that the vertex indices in ivlist[] are for the net0 vertex list.
	net1->vertex = (VERTEX *)malloc(nv*sizeof(VERTEX));
	net1->edgeList = (EDGE *)malloc(ne*sizeof(EDGE));
	net1->point = (POINT *)malloc(net0->np*sizeof(POINT));
	// Set up net1 vertices
	for (iv=0; iv<nv; iv++) {
		net1->vertex[iv] = net0->vertex[ivlist[iv]];
	}
	// Set up net1 edges
	for (i=0; i<ne; i++) {
		net1->edgeList[i] = elist[i];
	}
	// Copy net0 points to net1, initially set all points unused
	for (i=0; i<net0->np; i++) {
		net1->point[i] = net0->point[i];
		net1->point[i].used = false;
	}
	// Now flag those points that are used.
	for (ie=0; ie<ne; ie++) {
		for (kp=0; kp<elist[ie].npts; kp++) {
			i = elist[ie].pt[kp];
			net1->point[i].used = true;
		}
	}

	net1->ne = ne;
	net1->np = net0->np;
	net1->nv = nv;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *output_amfile;
	char drive[32], dir[128],filename[1024], ext[32];
	char errfilename[1024], input_amfile[1024], result_file[1024];
	char basename[64], input_basename[1024], output_basename[1024];
	char xstr[4], ystr[4], zstr[4];
	int ixcut, iycut, izcut, cmgui_flag;
	int ix, iy, iz;
	NETWORK *NP0, *NP1;

	if (argc != 7) {
		printf("Usage: am_join output_amfile basename ixcut iycut izcut cmgui_flag\n");
		fperr = fopen("am_join_error.log","w");
		printf("Usage: am_join output_amfile basename ixcut iycut izcut cmgui_flag\n");
		fprintf(fperr,"Usage: am_join output_amfile basename\n");
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	output_amfile = argv[1];
	strcpy(basename,argv[2]);
	sscanf(argv[3],"%d",&ixcut);
	sscanf(argv[4],"%d",&iycut);
	sscanf(argv[5],"%d",&izcut);
	sscanf(argv[6],"%d",&cmgui_flag);
	_splitpath(output_amfile,drive,dir,filename,ext);
	strcpy(input_basename,drive);
	strcat(input_basename,dir);
	strcat(input_basename,basename);
	sprintf(errfilename,"%s_am_join.log",input_basename);
	sprintf(result_file,"%s_am_join.out",input_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	fpout = fopen(result_file,"w");	
	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	InitNetwork(NP0);
//	err = ReadAmiraFile(input_amfile,NP0);
//	if (err != 0) return 2;
	NP1 = (NETWORK *)malloc(sizeof(NETWORK));
	for (ix=1; ix<=2; ix++) {
		if (ixcut == 0) {
			if (ix == 1) {
				strcpy(xstr,"_x1");
			} else {
				break;
			}
		} else {
			if (ix == 1) {
				strcpy(xstr,"_x1");
			} else {
				strcpy(xstr,"_x2");
			}
		}
		for (iy=1; iy<=2; iy++) {
			if (iycut == 0) {
				if (iy == 1) {
					strcpy(ystr,"_y1");
				} else {
					break;
				}
			} else {
				if (iy == 1) {
					strcpy(ystr,"_y1");
				} else {
					strcpy(ystr,"_y2");
				}
			}
			for (iz=1; iz<=2; iz++) {
				if (izcut == 0) {
					if (iz == 1) {
						strcpy(zstr,"_z1");
					} else {
						break;
					}
				} else {
					if (iz == 1) {
						strcpy(zstr,"_z1");
					} else {
						strcpy(zstr,"_z2");
					}
				}
				// make input_amfile
				strcpy(input_amfile,input_basename);
				strcat(input_amfile,xstr);
				strcat(input_amfile,ystr);
				strcat(input_amfile,zstr);
				strcat(input_amfile,".am");
				err = ReadAmiraFile(input_amfile,NP1);
				if (err != 0) return 2;
//				printf("Ranges: x: %6.1f %6.f  y: %6.1f %6.1f  z: %6.1f %6.1f\n",x1,x2,y1,y2,z1,z2);
//				err = CreateZoomNet(NP0,NP1,x1,x2,y1,y2,z1,z2);
				err = AddNetwork(NP0,NP1);
				if (err != 0) return 3;
				free_network(NP1);
			}
		}
	}

	// make output_amfile
	err = WriteAmiraFile(output_amfile,basename,NP0);
	if (err != 0) return 4;
	if (cmgui_flag == 1) {
		strcpy(output_basename,drive);
		strcat(output_basename,dir);
		strcat(output_basename,filename);
		err = WriteCmguiData(output_basename,NP0);
		if (err != 0) return 5;
	}
	return 0;
}