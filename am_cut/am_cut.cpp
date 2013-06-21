// Cut an Amira spatialgraph file into pieces (blocks).
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

FILE *fperr, *fpout;
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
						if (j == 13) {
							printf("j=13: d: %f\n",net->point[j].d);
						}
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
	printf("Points: np: %d np_used: %d\n",net->np,np_used);
	ne_used = 0;
	for (j=0; j<net->ne; j++) {
		if (net->edgeList[j].used) ne_used++;
	}
	printf("Edges: ne: %d ne_used: %d\n",net->ne,ne_used);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void free_network(NETWORK *net)
{
	free(elist);
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
// This code assumes that the cutting planes are specified in um.
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *input_amfile;
	char drive[32], dir[128],filename[1024], ext[32];
	char errfilename[1024], output_amfile[1024], result_file[1024];
	char basename[64], output_basename[1024];
	char xstr[4], ystr[4], zstr[4];
	float xcut, ycut, zcut;
	int ix, iy, iz;
	float x1, x2, y1, y2, z1, z2;
	NETWORK *NP0, *NP1;

	if (argc != 6) {
		printf("Usage: am_cut input_amfile basename xcut ycut zcut\n");
		printf("       Note: if a cut value is 0 there will be no cut\n");
		fperr = fopen("am_cut_error.log","w");
		fprintf(fperr,"Usage: am_cut input_amfile basename xcut ycut zcut\n");
		fprintf(fperr,"       Note: if a cut value is 0 there will be no cut\n");
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(basename,argv[2]);
	sscanf(argv[3],"%f",&xcut);
	sscanf(argv[4],"%f",&ycut);
	sscanf(argv[5],"%f",&zcut);
	_splitpath(input_amfile,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,basename);
	sprintf(errfilename,"%s_am_cut.log",output_basename);
	sprintf(result_file,"%s_am_cut.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	fpout = fopen(result_file,"w");	
	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;
	NP1 = (NETWORK *)malloc(sizeof(NETWORK));
	for (ix=1; ix<=2; ix++) {
		if (xcut == 0) {
			if (ix == 1) {
				x1 = -BIG;
				x2 = BIG;
				strcpy(xstr,"_x1");
			} else {
				break;
			}
		} else {
			if (ix == 1) {
				x1 = -BIG;
				x2 = xcut;
				strcpy(xstr,"_x1");
			} else {
				x1 = xcut;
				x2 = BIG;
				strcpy(xstr,"_x2");
			}
		}
		for (iy=1; iy<=2; iy++) {
			if (ycut == 0) {
				if (iy == 1) {
					y1 = -BIG;
					y2 = BIG;
					strcpy(ystr,"_y1");
				} else {
					break;
				}
			} else {
				if (iy == 1) {
					y1 = -BIG;
					y2 = ycut;
					strcpy(ystr,"_y1");
				} else {
					y1 = ycut;
					y2 = BIG;
					strcpy(ystr,"_y2");
				}
			}
			for (iz=1; iz<=2; iz++) {
				if (zcut == 0) {
					if (iz == 1) {
						z1 = -BIG;
						z2 = BIG;
						strcpy(zstr,"_z1");
					} else {
						break;
					}
				} else {
					if (iz == 1) {
						z1 = -BIG;
						z2 = zcut;
						strcpy(zstr,"_z1");
					} else {
						z1 = zcut;
						z2 = BIG;
						strcpy(zstr,"_z2");
					}
				}
				printf("Ranges: x: %6.1f %6.f  y: %6.1f %6.1f  z: %6.1f %6.1f\n",x1,x2,y1,y2,z1,z2);
				err = CreateBlockNet(NP0,NP1,x1,x2,y1,y2,z1,z2);
				if (err != 0) return 3;
				// make output_amfile
				strcpy(output_amfile,output_basename);
				strcat(output_amfile,xstr);
				strcat(output_amfile,ystr);
				strcat(output_amfile,zstr);
				strcat(output_amfile,".am");
				err = WriteAmiraFile(output_amfile,input_amfile,NP1);
				if (err != 0) return 4;
				free_network(NP1);
			}
		}
	}

	return 0;
}