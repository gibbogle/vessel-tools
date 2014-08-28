#include <cstdio>
#include <string>

#include "network.h"

extern FILE *fperr, *fpout;

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
bool EqualPoints(POINT p1, POINT p2)
{
	if (p1.x != p2.x) return false;
	if (p1.y != p2.y) return false;
	if (p1.z != p2.z) return false;
	return true;
}

//-----------------------------------------------------------------------------------------------------
// Check for a point position repeated on an edge
//-----------------------------------------------------------------------------------------------------
int CheckNetwork(NETWORK *net, char *str)
{
	int ie, ip, ne, npts,kp0,kp1,kp, ie1, iv0, iv1;
	float dave, dave1;
	EDGE edge, edge1;
	POINT p0, p1, p;
	int repeated;

	printf("CheckNetwork: %s\n",str);
	fprintf(fpout,"CheckNetwork: %s\n",str);
	repeated = 0;
	ne = net->ne;
	for (ie=0; ie<ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		npts = edge.npts;
		kp0 = edge.pt[0];
		p0 = net->point[kp0];
		kp1 = edge.pt[npts-1];
		p1 = net->point[kp1];
		for (ip=1; ip<npts; ip++) {
			kp = edge.pt[ip];
			p = net->point[kp];
			if (EqualPoints(p0,p)) {
				printf("CheckNetwork: repeated (0): ie,npts,v0,ip,kp0,kp: %d %d %d %d %d %d\n",ie,npts,0,ip,kp0,kp);
				printf("removed\n");
				fprintf(fpout,"CheckNetwork: repeated (0): ie,npts,v0,ip,kp0,kp: %d %d %d %d %d %d\n",ie,npts,0,ip,kp0,kp);
				fprintf(fpout,"removed\n");
				repeated++;
				net->edgeList[ie].used = false;
			}
		}
		for (ip=0; ip<npts-1; ip++) {
			kp = edge.pt[ip];
			p = net->point[kp];
			if (EqualPoints(p,p1)) {
				printf("CheckNetwork: repeated (1): ie,npts,ip,v1,kp,kp1: %d %d %d %d %d %d\n",ie,npts,ip,npts-1,kp,kp1);
				printf("removed\n");
				fprintf(fpout,"CheckNetwork: repeated (1): ie,npts,ip,v1,kp,kp1: %d %d %d %d %d %d\n",ie,npts,ip,npts-1,kp,kp1);
				fprintf(fpout,"removed\n");
				repeated++;
				net->edgeList[ie].used = false;
			}
		}
	}
	// Now check for two edges connecting the same pair of vertices.
	for (ie=0; ie<ne; ie++) {
		edge =net->edgeList[ie];
		iv0 = edge.vert[0];
		iv1 = edge.vert[1];
		for (ie1=0; ie1<ne; ie1++) {
			if (ie1 == ie) continue;
			edge1 = net->edgeList[ie1];
			if ((edge1.vert[0]==iv0 && edge1.vert[1]==iv1) || (edge1.vert[0]==iv1 && edge1.vert[1]==iv0)) {
				// double connection between iv0 and iv1 - remove the thinnest
				dave = 0;
				for (ip=0; ip<edge.npts; ip++) {
					dave += net->point[edge.pt[ip]].d;
				}
				dave1 = 0;
				for (ip=0; ip<edge1.npts; ip++) {
					dave1 += net->point[edge1.pt[ip]].d;
				}
				dave /= edge.npts;
				dave1 /= edge1.npts;
				if (dave < dave1) {
					net->edgeList[ie].used = false;
					printf("CheckNetwork: double connecting edge removed: ie,iv0,iv1: %d %d %d diams: %f %f\n",ie,iv0,iv1,dave,dave1);
					fprintf(fpout,"CheckNetwork: double connecting edge removed: ie,iv0,iv1: %d %d %d diams: %f %f\n",ie,iv0,iv1,dave,dave1);
				} else {
					net->edgeList[ie1].used = false;
					printf("CheckNetwork: double connecting edge removed: ie,iv0,iv1: %d %d %d diams: %f %f\n",ie1,iv0,iv1,dave,dave1);
					fprintf(fpout,"CheckNetwork: double connecting edge removed: ie,iv0,iv1: %d %d %d diams: %f %f\n",ie1,iv0,iv1,dave,dave1);
				}
			}
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
					net->vertex[i].point.d = 0;
					net->vertex[i].used = true;
				}
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
//					fprintf(fpout,"i: %6d npts: %4d\n",i,net->edgeList[i].npts);
					if (net->edgeList[i].npts < 1) {
						printf("ReadAmiraFile: i: %d npts: %d\n",i,net->edgeList[i].npts);
						return 1;
					}
					net->edgeList[i].npts_used = net->edgeList[i].npts;
					net->edgeList[i].pt = (int *)malloc(net->edgeList[i].npts*sizeof(int));
					net->edgeList[i].pt_used = (int *)malloc(net->edgeList[i].npts*sizeof(int));
					npts += net->edgeList[i].npts;
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
						sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
						net->edgeList[i].pt[k] = kp;
						net->edgeList[i].pt_used[k] = kp;
						kp++;
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
		net->edgeList[j].netID = j;
		if (net->edgeList[j].used) ne_used++;
	}
	printf("Edges: ne: %d ne_used: %d\n",net->ne,ne_used);

	CheckNetwork(net, "ReadAmiraFile");

	return 0;
}

