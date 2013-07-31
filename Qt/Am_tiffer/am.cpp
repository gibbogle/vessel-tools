#include "mainwindow.h"

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>


#define STR_LEN 128
#define BIG 1.0e6

EDGE *elist;
int *ivlist;

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
// Write Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int MainWindow::writeAmiraFile(const char *amFileOut, const char *amFileIn, NETWORK *net)
{
	int i, k, j, npts;
	EDGE edge;

    printf("\nWriteAmiraFile: %s\n",amFileOut);
    fprintf(fpout,"\nWriteAmiraFile: %s\n",amFileOut);
    npts = 0;
	for (i=0;i<net->ne;i++) {
		npts += net->edgeList[i].npts;
	}

	FILE *fpam = fopen(amFileOut,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
    fprintf(fpam,"# Created by dead.exe from: %s\n",amFileIn);
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
int MainWindow::readAmiraFile(const char *amFile, NETWORK *net)
{
	int i, j, k, kp, npts;
	int np_used, ne_used;
	EDGE edge;
	char line[STR_LEN];

//    if (fpout == NULL) {
//        fpout = fopen("dead.out","w");
//    }
    printf("ReadAmiraFile: %s\n",amFile);
    fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
    FILE *fpam = fopen(amFile,"r");

	npts = 0;
	kp = 0;
	k = 0;
	while (k < 3) {
		fgets(line, STR_LEN, fpam);		// reads until newline character
		printf("%s\n",line);
        fprintf(fpout,"%s\n",line);
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
    net->point = (APOINT *)malloc(net->np*sizeof(APOINT));
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
					kp = i;
					net->vertex[i].point.d = 0;
					net->point[kp] = net->vertex[i].point;
				}
				kp++;
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
					net->edgeList[i].pt[0] = net->edgeList[i].vert[0];
					net->edgeList[i].pt[net->edgeList[i].npts-1] = net->edgeList[i].vert[1];
				}
				printf("Got edge npts, total: %d\n",npts);
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
				printf("Got point thicknesses\n");
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
	ne_used = 0;
	for (j=0; j<net->ne; j++) {
		if (net->edgeList[j].used) ne_used++;
	}
	printf("Edges: ne: %d ne_used: %d\n",net->ne,ne_used);
    printf("Vertices: nv: %d\n",net->nv);
    printf("Points: np: %d np_used: %d\n",net->np,np_used);
    fprintf(fpout,"Edges: ne: %d ne_used: %d\n",net->ne,ne_used);
    fprintf(fpout,"Vertices: nv: %d\n",net->nv);
    fprintf(fpout,"Points: np: %d np_used: %d\n",net->np,np_used);
    return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void MainWindow::freeNetwork(NETWORK *net)
{
    int ie;
//	free(elist);
//	free(ivlist);
    if (net->nv > 0) free(net->vertex);
    if (net->np > 0) free(net->point);
    if (net->ne > 0) {
        for (ie=0; ie<net->ne; ie++) {
            free(net->edgeList[ie].pt);
            free(net->edgeList[ie].pt_used);
        }
        free(net->edgeList);
    }
    net->ne = 0;
    net->nv = 0;
    net->np = 0;
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
// Locate all edges with at least one vertex within the block
//-----------------------------------------------------------------------------------------------------
int CreateBlockNet(NETWORK *net0, NETWORK *net1, float x1, float x2, float y1, float y2, float z1, float z2)
{
	int i, j, ie, iv, kp, kv, ne, nv;
	float xv, yv, zv;
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
    net1->point = (APOINT *)malloc(net0->np*sizeof(APOINT));
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
// Create a network made up of the dead-end vessels of the original network
//-----------------------------------------------------------------------------------------------------
int MainWindow::createNetwork(NETWORK *net, NETWORK *deadnet, DEADEND *deadlist, int ndead)
{
    int k, j, iv, ie, ip, iev;
    int *vert;
    EDGE edge;
    APOINT p;

    deadnet->ne = ndead;
    deadnet->edgeList = (EDGE *)malloc(deadnet->ne*sizeof(EDGE));
    deadnet->nv = 2*ndead;
    deadnet->vertex = (VERTEX *)malloc(deadnet->nv*sizeof(VERTEX));
    iv = 0;
    ie = 0;
    ip = 0;
    for (k=0; k<ndead; k++) {
        if (deadlist[k].intensity == 0) continue;
        iev = deadlist[k].ie;
        edge = net->edgeList[iev];
        vert = edge.vert;
        deadnet->vertex[iv].point = net->vertex[vert[0]].point;
        deadnet->vertex[iv].used = true;
        deadnet->edgeList[ie].vert[0] = iv;
        iv++;
        deadnet->vertex[iv].point = net->vertex[vert[1]].point;
        deadnet->vertex[iv].used = true;
        deadnet->edgeList[ie].vert[1] = iv;
        iv++;
        deadnet->edgeList[ie].npts = edge.npts;
        deadnet->edgeList[ie].npts_used = edge.npts;
        deadnet->edgeList[ie].pt = (int *)malloc(edge.npts*sizeof(int));
        ie++;
        ip += edge.npts;
    }
    deadnet->nv = iv;
    deadnet->ne = ie;
    deadnet->np = ip;
    deadnet->point = (APOINT *)malloc(deadnet->np*sizeof(APOINT));
    ip = 0;
    ie = 0;
    for (k=0; k<ndead; k++) {
        if (deadlist[k].intensity == 0) continue;
        iev = deadlist[k].ie;
        edge = net->edgeList[iev];
        for (j=0; j<edge.npts; j++) {
            p = net->point[edge.pt[j]];
            deadnet->edgeList[ie].pt[j] = ip;
            deadnet->point[ip] = p;
            ip++;
        }
        ie++;
    }
    printf("nv,ne,np: %d %d %d\n",deadnet->nv,deadnet->ne,deadnet->np);
    return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int MainWindow::readNetwork(NETWORK *net, const char *amfile)
{
	int err;

    if (net->ne != 0) freeNetwork(net);
    err = readAmiraFile(amfile,net);
    if (err != 0) return 1;
    am_read = true;
//    err = WriteAmiraFile(output_amfile,input_amfile,NP1);
//    if (err != 0) return 4;
//    free_network(NP1);

	return 0;
}

bool MainWindow::inSphere(APOINT p)
{
    float dx, dy, dz;

    dx = p.x - sphereCentre[0];
    dy = p.y - sphereCentre[1];
    dz = p.z - sphereCentre[2];
    if (dx*dx+dy*dy+dz*dz < sphereRadius*sphereRadius)
        return true;
    else
        return false;
}


