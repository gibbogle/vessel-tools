// Specify a subregion, either a block (range of x, y, z) or a sphere (centre, radius)
// and generate a network comprising all vessels that have an end point within the region.
// Optionally, the largest connected network is chosen.
// Working with the Amira spatialgraph file format.

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

int WriteCmguiData(char *basename, NETWORK *net, float origin_shift[]);

#define PI 3.14159
#define STR_LEN 128
#define NEMAX 1000
#define NBOX 500

FILE *fperr, *fpout;
int nconnected;
int ne_net[NEMAX];
float xmin, xmax, ymin, ymax, zmin,zmax, centre[3], radius, origin_shift[3];
bool use_sphere;

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float dist(NETWORK *net, int k1, int k2)
{
	float dx = net->point[k2].x - net->point[k1].x;
	float dy = net->point[k2].y - net->point[k1].y;
	float dz = net->point[k2].z - net->point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------
// mode = 'N' for numbers
//        'D' for diameters
//        'P' for positions
//-----------------------------------------------------------------------------------
void showedge(NETWORK *net, int ie, char mode)
{
	int k, kp;

	printf("edge: %d npts: %d\n",ie,net->edgeList[ie].npts);
	for (k=0; k<net->edgeList[ie].npts; k++) {
		kp = net->edgeList[ie].pt[k];
		if (mode == 'N')
			printf("%7d ",kp);
		else if (mode == 'D')
			printf("%7.1f ",net->point[kp].d);
		else 
			printf("%7.1f %7.1f %7.1f   ",net->point[kp].x,net->point[kp].y,net->point[kp].z);
		printf("\n");
	}
}

//-----------------------------------------------------------------------------------------------------
// This code is faulty - because points can appear multiple times, there are multiple subtractions.
// NOT USED
//-----------------------------------------------------------------------------------------------------
int ShiftOrigin(NETWORK *net, float origin_shift[])
{
	int i, k, j;
	EDGE edge;

	for (i=0;i<net->nv;i++) {
		net->vertex[i].point.x -= origin_shift[0];
		net->vertex[i].point.y -= origin_shift[1];
		net->vertex[i].point.z -= origin_shift[2];
	}
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
			net->point[j].x -= origin_shift[0];
			net->point[j].y -= origin_shift[1];
			net->point[j].z -= origin_shift[2];
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Write Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile(char *amFileOut, char *amFileIn, NETWORK *net, float origin_shift[])
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
			net->vertex[i].point.x - origin_shift[0],
			net->vertex[i].point.y - origin_shift[1],
			net->vertex[i].point.z - origin_shift[2]);
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
				net->point[j].x - origin_shift[0],
				net->point[j].y - origin_shift[1],
				net->point[j].z - origin_shift[2]);
		}
	}
	fprintf(fpam,"\n@5\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
			fprintf(fpam,"%6.2f\n",net->point[j].d);
			//if (net->point[j].d > diam_max) {
			//	printf("d > diam_max: edge: %d ne: %d k: %d d: %6.2f\n",i,net->ne,k,net->point[j].d);
			//	exit(1);
			//}
		}
	}
	fclose(fpam);
	printf("Completed WriteAmiraFile\n");
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

	fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
	FILE *fpam = fopen(amFile,"r");

	npts = 0;
	kp = 0;
	k = 0;
	while (k < 3) {
		fgets(line, STR_LEN, fpam);		// reads until newline character
		printf("%s\n",line);
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
	printf("Allocated arrays: %d %d %d\n",net->np,net->nv,net->ne);

	// Initialize
	for (i=0; i<net->ne; i++) {
		net->edgeList[i].used = false;
	}
	for (i=0; i<net->np; i++) {
		net->point[i].used = false;
	}
	printf("Initialised\n");

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
//					kp = i;
					net->vertex[i].point.d = 0;
//					net->point[kp] = net->vertex[i].point;
				}
//				kp++;
				printf("Got vertices\n");
			} else if (k == 2) {
				for (i=0;i<net->ne;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @2\n");
						return 1;
					}
					sscanf(line,"%d %d",&net->edgeList[i].vert[0],&net->edgeList[i].vert[1]);
					net->edgeList[i].used = true;
				}
				printf("Got edge vertex indices\n");
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
//					net->edgeList[i].pt[0] = net->edgeList[i].vert[0];
//					net->edgeList[i].pt[net->edgeList[i].npts-1] = net->edgeList[i].vert[1];
				}
				printf("Got edge npts, total: %d\n",npts);
			} else if (k == 4) {
				for (i=0;i<net->ne;i++) {
					edge = net->edgeList[i];
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
//						if (k > 0 && k<edge.npts-1) {
							sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[k] = kp;
							net->edgeList[i].pt_used[k] = kp;
							kp++;
//						}
					}
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
//						if (j < net->nv) {		// because the first nv points are vertices
//							net->vertex[j].point.d = net->point[j].d;
//						}
						dave += net->point[j].d;
						net->edgeList[i].segavediam = dave/edge.npts;
					}
				}
				printf("Got point thicknesses\n");
			}
		}
	}
	// Flag used points
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		net->edgeList[i].netID = 0;
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			net->point[j].used = true;
		}
	}
	for (i=0;i<net->nv;i++) {
		net->vertex[i].netID = 0;
	}
	fclose(fpam);
	np_used = 0;
	for (j=0; j<net->np; j++) {
		if (net->point[j].used) np_used++;
	}
	printf("Points: np: %d np_used: %d\n",net->np,np_used);
	ne_used = 0;
	for (j=0; j<net->ne; j++) {
		if (net->edgeList[j].used) ne_used++;
	}
	printf("Edges: ne: %d ne_used: %d\n",net->ne,ne_used);
	return 0;
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
//-----------------------------------------------------------------------------------------------------
int adjacentVertex(NETWORK *net, int kv, int ie)
{
	if (kv == net->edgeList[ie].vert[0])
		return net->edgeList[ie].vert[1];
	else
		return net->edgeList[ie].vert[0];
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int getStartingEdge(NETWORK *net)
{
	int ie;
	for (ie=0; ie<net->ne; ie++) {
		if (net->edgeList[ie].netID == 0) {
			return ie;
		}
	}
	return -1;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//procedure DFS(G,v):
int DFS(NETWORK *net, int kv, int netID)
{
	int i, ie, kw;

    //label v as explored
	net->vertex[kv].netID = netID;
    //for all edges e in G.adjacentEdges(v) do
	for (i=0; i<net->vertex[kv].nlinks;i++) {
		ie = net->vertex[kv].edge[i];
    //    if edge e is unexplored then
		if (net->edgeList[ie].netID == 0) {
    //        w <- G.adjacentVertex(v,e)
			kw = adjacentVertex(net,kv,ie);
    //        if vertex w is unexplored then
			if (net->vertex[kw].netID == 0) {
    //            label e as a discovery edge
				net->edgeList[ie].netID = netID;	// no distinction between discovery and back edges
				nconnected++;
    //            recursively call DFS(G,w)
				DFS(net,kw,netID);
    //        else
			} else {
    //            label e as a back edge
				net->edgeList[ie].netID = netID;	// no distinction between discovery and back edges
				nconnected++;
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// All edges connected to edge ie0 are tagged with netID.
// Returns the number of tagged edges as nconnected.
// This is the kernel of the program.
// Here is the pseudocode from Wikipedia Depth-first search
//procedure DFS(G,v):
//    label v as explored
//    for all edges e in G.adjacentEdges(v) do
//        if edge e is unexplored then
//            w <- G.adjacentVertex(v,e)
//            if vertex w is unexplored then
//                label e as a discovery edge
//                recursively call DFS(G,w)
//            else
//                label e as a back edge	
//-----------------------------------------------------------------------------------------------------
int ConnectEdges(NETWORK *net, int kv0, int netID, int *nconn)
{
	nconnected = 0;
	printf("kv0, netID: %d %d\n",kv0,netID);
	DFS(net,kv0,netID);
	*nconn = nconnected;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Network net1 is created from the edges of net0 that are tagged with netID
// Note the inefficiency in keeping all points, used or not.
//-----------------------------------------------------------------------------------------------------
int ExtractConnectedNet(NETWORK *net0, NETWORK *net1, int netID)
{
	int i, j, ie, iv, kp, kv, ne, nv;
	EDGE *elist;
	int *ivlist;

	ne = 0;
	for (i=0; i<net0->ne; i++) {
		if (net0->edgeList[i].netID == netID) ne++;
	}
	printf("Connected: ne: %d\n",ne);
	// Create a list of edges, and a list of vertices
	// The index number of a vertex in the list replaces vert[] in the elist entry.
	elist = (EDGE *)malloc(ne*sizeof(EDGE));
	ivlist = (int *)malloc(2*ne*sizeof(int));
	nv = 0;
	ne = 0;
	for (i=0; i<net0->ne; i++) {

		if (net0->edgeList[i].netID == netID) {
			elist[ne] = net0->edgeList[i];
			for (j=0; j<2; j++) {
				iv = net0->edgeList[i].vert[j];
				kv = ivlistAdd(iv,ivlist,&nv);
				elist[ne].vert[j] = kv;
			}
			ne++;
		}
	}
	printf("Connected: ne: %d  nv: %d\n",ne,nv);
	// Now from this list of edges we need to create the network net1.
	// Note that the vertex indices in ivlist[] are for the net0 vertex list.
	net1->vertex = (VERTEX *)malloc(nv*sizeof(VERTEX));
	net1->edgeList = (EDGE *)malloc(ne*sizeof(EDGE));
	net1->point = (POINT *)malloc(net0->np*sizeof(POINT));
	// Set up net1 vertices
	for (iv=0; iv<nv; iv++) {
		net1->vertex[iv] = net0->vertex[ivlist[iv]];
		net1->vertex[iv].netID = 0;
	}
	// Set up net1 edges
	for (i=0; i<ne; i++) {
		net1->edgeList[i] = elist[i];
		net1->edgeList[i].netID = 0;	// initially edges are not connected to a con_net
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
			//if (net1->point[i].d > 75) {
			//	printf("Error ExtractConnectedNet: edge: %d npts: %d pt: %d %d d: %6.1f\n",ie,elist[ie].npts,kp,i,net1->point[i].d);
			//	return 1;
			//}
		}
	}

	net1->ne = ne;
	net1->np = net0->np;
	net1->nv = nv;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int CreateLargestConnectedNet(NETWORK *net1, NETWORK *net2)
{
	int ie0, kv0, netID, nconn, err;
	int ne_max, netID_max;

	// First determine a net_ID for every edge
	netID = 0;
	for (;;) {
		// Start from any unconnected edge
		ie0 = getStartingEdge(net1);
		if (ie0 < 0) break;	// done
		netID++;
		if (netID > NEMAX) {
			printf("NEMAX exceeded\n");
			return 1;
		}
		kv0 = net1->edgeList[ie0].vert[0];
		printf("netID, ie0, kv0: %d %d %d\n",netID,ie0,kv0);
		ConnectEdges(net1,kv0,netID,&nconn);
		ne_net[netID] = nconn;
		printf("netID, nconn: %d %d\n",netID,nconn);
	}
	ne_max = 0;
	for (int i=1; i<=netID; i++) {
		if (ne_net[i] > ne_max) {
			ne_max = ne_net[i];
			netID_max = i;
		}
	}
	// Now create net2 from the largest subnetwork
	err = ExtractConnectedNet(net1,net2,netID_max);
	if (err != 0) return err;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Returns true if the point p is inside the specified region, false otherwise
//-----------------------------------------------------------------------------------------------------
bool inside(POINT p)
{
	float dx, dy, dz;

	if (use_sphere) {
		dx = p.x - centre[0];
		dy = p.y - centre[1];
		dz = p.z - centre[2];
		if (dx*dx+dy*dy+dz*dz < radius*radius)
			return true;
		else
			return false;
	} else {
		if (p.x < xmin || p.x > xmax || p.y < ymin || p.y > ymax || p.z < zmin || p.z > zmax)
			return false;
		else
			return true;
	}
}

//-----------------------------------------------------------------------------------------------------
// Compute ranges of the network
//-----------------------------------------------------------------------------------------------------
void getRanges(NETWORK *net)
{
	int ie, kp0, kp1;
	float rxmin, rxmax, rymin, rymax, rzmin, rzmax;
	POINT p;

	rxmin = rymin = rzmin = 1.0e10;
	rxmax = rymax = rzmax = -1.0e10;
	for (ie=0; ie<net->ne; ie++) {
		kp0 = net->edgeList[ie].pt[0];
		p = net->point[kp0];
		rxmin = MIN(rxmin,p.x);
		rxmax = MAX(rxmax,p.x);
		rymin = MIN(rymin,p.y);
		rymax = MAX(rymax,p.y);
		rzmin = MIN(rzmin,p.z);
		rzmax = MAX(rzmax,p.z);
		kp1 = net->edgeList[ie].pt[net->edgeList[ie].npts-1];
		p = net->point[kp1];
		rxmin = MIN(rxmin,p.x);
		rxmax = MAX(rxmax,p.x);
		rymin = MIN(rymin,p.y);
		rymax = MAX(rymax,p.y);
		rzmin = MIN(rzmin,p.z);
		rzmax = MAX(rzmax,p.z);
	}
	printf("Range of network: x: %6.1f - %6.1f y: %6.1f - %6.1f z: %6.1f - %6.1f\n",rxmin,rxmax,rymin,rymax,rzmin,rzmax);
	fprintf(fpout,"Range of network: x: %6.1f - %6.1f y: %6.1f - %6.1f z: %6.1f - %6.1f\n",rxmin,rxmax,rymin,rymax,rzmin,rzmax);
}

//-----------------------------------------------------------------------------------------------------
// Locate all edges with an end point within a specified region.
// Note the inefficiency in keeping all points, used or not.
//-----------------------------------------------------------------------------------------------------
int CreateSelectNet(NETWORK *net0, NETWORK *net1)
{
	int i, j, ie, iv, kp, kv, ne, nv, kp0, kp1;
	float d, dave;
	bool in;
	EDGE *elist;
	int *ivlist;

	ne = 0;
	for (ie=0; ie<net0->ne; ie++) {
		kp0 = net0->edgeList[ie].pt[0];
		kp1 = net0->edgeList[ie].pt[net0->edgeList[ie].npts-1];
		if (inside(net0->point[kp0]) || inside(net0->point[kp1])) ne++;
	}
	printf("ne: %d\n",ne);
	// Create a list of edges, and a list of vertices
	// The index number of a vertex in the list replaces vert[] in the elist entry.
	elist = (EDGE *)malloc(ne*sizeof(EDGE));
	ivlist = (int *)malloc(2*ne*sizeof(int));
	nv = 0;
	ne = 0;
	for (ie=0; ie<net0->ne; ie++) {
		kp0 = net0->edgeList[ie].pt[0];
		kp1 = net0->edgeList[ie].pt[net0->edgeList[ie].npts-1];
		if (inside(net0->point[kp0]) || inside(net0->point[kp1])) {
			elist[ne] = net0->edgeList[ie];
			for (j=0; j<2; j++) {
				iv = net0->edgeList[ie].vert[j];
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
	net1->point = (POINT *)malloc(net0->np*sizeof(POINT));	// keeping full original list of points
	// Set up net1 vertices
	for (iv=0; iv<nv; iv++) {
		net1->vertex[iv] = net0->vertex[ivlist[iv]];
		net1->vertex[iv].nlinks = 0;
	}
	// Set up net1 edges and count vertex edges
	for (i=0; i<ne; i++) {
		net1->edgeList[i] = elist[i];
		net1->edgeList[i].netID = 0;	// initially edges are not connected to a connected subnetwork
		for (int k=0; k<2;k++) {
			iv = net1->edgeList[i].vert[k];
			net1->vertex[iv].nlinks++;
		}
	}
	// Set up vertex edge lists
	for (iv=0; iv<nv; iv++) {
		net1->vertex[iv].edge = (int *)malloc(net1->vertex[iv].nlinks*sizeof(int));
		net1->vertex[iv].nlinks = 0;
	}
	for (i=0; i<ne; i++) {
		for (int k=0; k<2;k++) {
			iv = net1->edgeList[i].vert[k];
			net1->vertex[iv].edge[net1->vertex[iv].nlinks] = i;
			net1->vertex[iv].nlinks++;
		}
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

int myFactorial( int integer)
{
	if ( integer == 1)
		return 1;
	else {
		return (integer * (myFactorial(integer-1)));
	}
} 

//-----------------------------------------------------------------------------------------------------
// The average diameter of a vessel (edge) is now estimated by dividing the volume by the length.
// All distances in the .am file are in um.
//-----------------------------------------------------------------------------------------------------
int CreateDistributions(NETWORK *net)
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double lsegadbox[NBOX];
	double ad, len, ddiam, dlen, ltot, lsum, dsum, dvol, r2, r2prev, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam, ave_seg_diam;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot;
	double lenlimit = 3.0;
	EDGE edge;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lsegadbox[k] = 0;
		lvbox[k] = 0;
	}
	printf("Compute diameter distributions\n");
	fprintf(fperr,"Compute diameter distributions\n");
	// Diameters
	ddiam = 0.5;
	ndtot = 0;
	nsegdtot = 0;
	lsegdtot = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	volume = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
//		printf("ie: %d npts: %d\n",ie,edge.npts);
//		fprintf(fperr,"ie: %d npts: %d\n",ie,edge.npts);
//		fflush(fperr);
		bool dbug = false;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		dvol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = net->point[kp].d;
//			ad = avediameter[kp];
			ave_pt_diam += ad;
			if (dbug) {
				printf("%d  %d  %f  %f\n",ip,kp,ad,ddiam);
				fprintf(fperr,"%d  %d  %f  %f\n",ip,kp,ad,ddiam);
			}
			fflush(fperr);
//			dsum += ad;
			if (ad < 0.001) {
				printf("Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fperr,"Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				return 1;
			}
			ka = int(ad/ddiam + 0.5);
			if (ka >= NBOX) {
				printf("Vessel too wide (point): d: %f k: %d\n",ad,ka);
				fprintf(fperr,"Vessel too wide (point): d: %f k: %d\n",ad,ka);
				continue;
			}
			adbox[ka]++;
			ndtot++;
			if (ip > 0) {
				dlen = dist(net,kp,kpprev);
				r2 = ad*ad/4;
				dvol += PI*dlen*(r2 + r2prev)/2;
				lsum += dlen;
			}
			kpprev = kp;
			r2prev = r2;
		}
		net->edgeList[ie].length_um = float(lsum);
		volume += dvol;
		if (dbug) {
			printf("lsum: %f\n",lsum);
			fprintf(fperr,"lsum: %f\n",lsum);
			fflush(fperr);
			if (lsum == 0) return 1;
		}
//		ad = dsum/edge.npts;
		ad = 2*sqrt(dvol/(PI*lsum));	// segment diameter
		ave_seg_diam += ad;
		if (ad < 0.001) {
			printf("Zero segment diameter: edge: %d ad: %f\n",ie,ad);
			fprintf(fperr,"Zero segment diameter: edge: %d ad: %f\n",ie,ad);
			return 1;
		}
		ka = int(ad/ddiam + 0.5);
		if (ka >= NBOX) {
			printf("Vessel too wide (segment ave): d: %f k: %d\n",ad,ka);
			fprintf(fperr,"Vessel too wide (segment ave): d: %f k: %d\n",ad,ka);
			continue;
		}
		segadbox[ka]++;
		nsegdtot++;
		lsegadbox[ka] += lsum;
		lsegdtot += lsum;
	}
	// Determine d95, the diameter that >95% of points exceed.
	dsum = 0;
	for (k=0; k<NBOX; k++) {
		dsum += adbox[k]/float(ndtot);
		if (dsum > 0.05) {
			d95 = (k-1)*ddiam;
			break;
		}
	}
	printf("Compute length distributions: lower limit = %6.1f um\n",lenlimit);
	fprintf(fperr,"Compute length distributions: lower limit = %6.1f um\n",lenlimit);
	// Lengths
	dlen = 1;
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/dlen + 0.5);
		if (k*dlen <= lenlimit) continue;
		if (k >= NBOX) {
			printf("Edge too long: len: %d  %f  k: %d\n",ie,len,k);
			fprintf(fperr,"Edge too long: len: %d  %f  k: %d\n",ie,len,k);
			continue;
		}
		lvbox[k]++;
		ave_len += len;
		ltot++;
	}
	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	fprintf(fpout,"Total vertices: %d  points: %d\n",net->nv,net->np);
	fprintf(fpout,"Vessels: %d\n",net->ne);
	printf("Average pt diameter: %6.2f vessel diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	fprintf(fpout,"Average pt diameter: %6.2f vessel diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	printf("Average vessel length: %6.1f\n",ave_len/ltot);
	fprintf(fpout,"Average vessel length: %6.1f\n",ave_len/ltot);
	fprintf(fpout,"Volume: %10.0f\n\n",volume);

	//ndpts = 0;
	//for (k=0; k<NBOX; k++) {
	//	if (adbox[k]/float(ndtot) >= 0.0005) {
	//		ndpts = k+1;
	//	}
	//}
	for (k=NBOX-1; k>=0; k--) {
		if (segadbox[k] > 0) break;
	}
	ndpts = k+2;
	fprintf(fpout,"Vessel diameter distribution\n");
	fprintf(fpout,"   um    number  fraction    length  fraction\n");
	for (k=0; k<ndpts; k++) {
		fprintf(fpout,"%6.2f %8d %9.5f  %8.0f %9.5f\n",k*ddiam,segadbox[k],segadbox[k]/float(nsegdtot),
			lsegadbox[k],lsegadbox[k]/lsegdtot);
	}

	//nlpts = 0;
	//for (k=0; k<NBOX; k++) {
	//	if (lvbox[k]/ltot >= 0.0005) {
	//		nlpts = k+1;
	//	}
	//}
	for (k=NBOX-1; k>=0; k--) {
		if (lvbox[k] > 0) break;
	}
	nlpts = k+2;
	fprintf(fpout,"Vessel length distribution\n");
	fprintf(fpout,"   um    number  fraction\n");
	for (k=0; k<nlpts; k++) {
		fprintf(fpout,"%6.2f %8d %9.5f\n",k*dlen,lvbox[k],lvbox[k]/ltot);
	}
	return 0;
}


//-----------------------------------------------------------------------------------------------------
// If connect_flag==1, the largest connected network is found, other edges are dropped.
// Note: sphere region not implemented
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *input_amfile;
	char drive[32], dir[2048],filename[256], ext[32];
	char errfilename[2048], output_amfile[2048], result_file[2048];
	char output_basename[2048];
	int cmgui_flag, connect_flag;
	bool connect;
	NETWORK *NP0, *NP1, *NP2;

	if (argc == 11) {
		use_sphere = false;
		input_amfile = argv[1];
		strcpy(output_amfile,argv[2]);
		sscanf(argv[3],"%f",&xmin);
		sscanf(argv[4],"%f",&xmax);
		sscanf(argv[5],"%f",&ymin);
		sscanf(argv[6],"%f",&ymax);
		sscanf(argv[7],"%f",&zmin);
		sscanf(argv[8],"%f",&zmax);
		sscanf(argv[9],"%d",&connect_flag);
		sscanf(argv[10],"%d",&cmgui_flag);
	//} else if (argc == 9) {
	//	use_sphere = true;
	//	input_amfile = argv[1];
	//	strcpy(output_amfile,argv[2]);
	//	sscanf(argv[3],"%f",&centre[0]);
	//	sscanf(argv[4],"%f",&centre[1]);
	//	sscanf(argv[5],"%f",&centre[2]);
	//	sscanf(argv[6],"%f",&radius);
	//	sscanf(argv[7],"%d",&connect_flag);
	//	sscanf(argv[8],"%d",&cmgui_flag);
	} else {
		printf("Usage: am_block input_amfile output_amfile xmin xmax ymin ymax zmin zmax connect_flag cmgui_flag\n");
//		printf(" or\n");
//		printf("       am_block input_amfile output_amfile centre_x centre_y centre_z radius connect_flag cmgui_flag\n");
		fperr = fopen("select_error.log","w");
		fprintf(fperr,"Usage: am_block input_amfile output_amfile xmin xmax ymin ymax zmin zmax connect_flag cmgui_flag\n");
//		fprintf(fperr," or\n");
//		fprintf(fperr,"       am_block input_amfile output_amfile centre_x centre_y centre_z radius connect_flag cmgui_flag\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	connect = (connect_flag != 0);

	_splitpath(output_amfile,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_select.log",output_basename);
	sprintf(result_file,"%s_select.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	fpout = fopen(result_file,"w");	
	printf("Input .am file: %s\n",input_amfile);
	printf("connect_flag: %d\n",connect_flag);
	printf("cmgui_flag: %d\n",cmgui_flag);
	fprintf(fpout,"Input .am file: %s\n",input_amfile);
	fprintf(fpout,"connect_flag: %d\n",connect_flag);
	fprintf(fpout,"cmgui_flag: %d\n",cmgui_flag);
	if (use_sphere) {
		printf("sphere centre: %6.1f %6.1f %6.1f  radius: %6.1f\n",centre[0],centre[1],centre[2],radius);
		fprintf(fpout,"sphere centre: %6.1f %6.1f %6.1f  radius: %6.1f\n",centre[0],centre[1],centre[2],radius);
	} else {
		printf("block ranges: x: %6.1f - %6.1f y: %6.1f - %6.1f z: %6.1f - %6.1f\n",xmin,xmax,ymin,ymax,zmin,zmax);
		fprintf(fpout,"block ranges: x: %6.1f - %6.1f y: %6.1f - %6.1f z: %6.1f - %6.1f\n",xmin,xmax,ymin,ymax,zmin,zmax);
	}

	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;
	getRanges(NP0);
	NP1 = (NETWORK *)malloc(sizeof(NETWORK));
	err = CreateSelectNet(NP0,NP1);
	if (err != 0) return 3;

	printf("NP1: ne, nv, np: %d %d %d\n",NP1->ne,NP1->nv,NP1->np);

	if (connect) {
		NP2 = (NETWORK *)malloc(sizeof(NETWORK));
		err = CreateLargestConnectedNet(NP1,NP2);
		if (err != 0) return 4;
	}

	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;

	if (connect) {
		err = WriteAmiraFile(output_amfile,input_amfile,NP2,origin_shift);
		if (err != 0) return 5;
		if (cmgui_flag == 1) {
			err = WriteCmguiData(output_basename,NP2,origin_shift);
			if (err != 0) return 6;
		}
		err = CreateDistributions(NP2);
		if (err != 0) return 7;
	} else {
		err = WriteAmiraFile(output_amfile,input_amfile,NP1,origin_shift);
		if (err != 0) return 5;
		if (cmgui_flag == 1) {
			err = WriteCmguiData(output_basename,NP1,origin_shift);
			if (err != 0) return 6;
		}
		err = CreateDistributions(NP1);
		if (err != 0) return 7;
	}
	return 0;
}