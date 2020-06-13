// translate.cpp
// To read an Amira file, compute distributions, and write CMGUI files

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

FILE *fperr, *fpout;
FILE *exelem, *exnode;

float ddiam, dlen;
bool use_len_limit, use_len_diam_limit;
float len_limit, len_diam_limit;

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
//-----------------------------------------------------------------------------------------------------
int inlist(int *list, int n, int k)
{
	int i;
	if (n == 0) return -1;
	for (i=0; i<n; i++) {
		if (list[i] == k) return i;
	}
	return -1;
}

//-----------------------------------------------------------------------------------------------------
// Look for (x,y,z) in the vertex position list.
// If found, if vertex[].point_index < 0 set it to kp, then set kpv = vertex[].point_index
//-----------------------------------------------------------------------------------------------------
void get_vertex_ptindex(NETWORK *net, float x, float y, float z, int kp, int *kpv)
{
	int kv, kvhit;
	float tol = 0.01;
	bool hit;

	hit = false;
	for (kv=0; kv<net->nv; kv++) {
		if (fabs(x - net->vertex[kv].point.x) < tol) {
			if (fabs(y - net->vertex[kv].point.y) < tol) {
				if (fabs(z - net->vertex[kv].point.z) < tol) {
					hit = true;
					kvhit = kv;
					break;
				}
			}
		}
	}
	if (hit) {
		if (net->vertex[kvhit].point_index < 0) {
			net->vertex[kvhit].point_index = kp;
		}
		*kpv = net->vertex[kvhit].point_index;
	} else {
		printf("Vertex (x,y,z) not found in vertex list, adding it: %f %f %f\n",x,y,z);
		net->vertex[net->nv].point.x = x;
		net->vertex[net->nv].point.y = y;
		net->vertex[net->nv].point.z = z;
		net->vertex[net->nv].point_index = kp;
		net->nv++;
	}
}

//-----------------------------------------------------------------------------------------------------
// Write Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile1(char *amFileOut, char *amFileIn, NETWORK *net, float origin_shift[])
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
		}
	}
	fclose(fpam);
	printf("Completed WriteAmiraFile\n");
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile1(char *amFile, NETWORK *net)
{
	int i, j, k, kp, kpv, npts, nee, npp;
	float x, y, z;
	EDGE edge;
	char line[STR_LEN];

	fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
	FILE *fpam = fopen(amFile,"r");

	npts = 0;
	kp = 0;
	k = 0;
	while (k < 3) {
		fgets(line, STR_LEN, fpam);		// reads until newline character
		if (strncmp(line,"define VERTEX",13) == 0) {
			sscanf(line+13,"%d",&net->nv);
			k++;
		}
		if (strncmp(line,"define EDGE",11) == 0) {
			sscanf(line+11,"%d",&net->ne);
			nee = 2*net->ne;
			k++;
		}
		if (strncmp(line,"define POINT",12) == 0) {
			sscanf(line+12,"%d",&net->np);
			npp = 2*net->np;
			k++;
		}
	}

	net->vertex = (VERTEX *)malloc(2*net->nv*sizeof(VERTEX));
	net->edgeList = (EDGE *)malloc(nee*sizeof(EDGE));	// 2* for added joining edges
	net->point = (POINT *)malloc(npp*sizeof(POINT));	// 2* for added joining edges

	// Initialize
	for (i=0; i<nee; i++) {
		net->edgeList[i].used = false;
	}
	for (i=0; i<npp; i++) {
		net->point[i].used = false;
	}

	while (1) {
		if (fgets(line, STR_LEN, fpam) == NULL) {
			printf("Finished reading SpatialGraph file\n\n");
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
					sscanf(line,"%f %f %f\n",&net->vertex[i].point.x,&net->vertex[i].point.y,&net->vertex[i].point.z);
					net->vertex[i].point.d = 0;
					net->vertex[i].point_index = -1;
				}
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
				}
				printf("Got edge npts, total: %d\n",npts);
			} else if (k == 4) {
				kp = 0;
				for (i=0;i<net->ne;i++) {
					edge = net->edgeList[i];
					float len = 0;
					printf("edge: %d  npts: %d\n",i,edge.npts);
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
						sscanf(line,"%f %f %f",&x,&y,&z);
						if (k == 0 || k == edge.npts-1) {
							get_vertex_ptindex(net,x,y,z,kp,&kpv);
							net->point[kpv].x = x;
							net->point[kpv].y = y;
							net->point[kpv].z = z;
							net->edgeList[i].pt[k] = kpv;
							net->edgeList[i].pt_used[k] = kpv;
						} else {
							net->point[kp].x = x;
							net->point[kp].y = y;
							net->point[kp].z = z;
							net->edgeList[i].pt[k] = kp;
							net->edgeList[i].pt_used[k] = kp;
						}
//						printf("point: %d  %f %f %f\n",kp,x,y,z);
						kp++;
						if (k > 0) {
							len = len + dist(net,net->edgeList[i].pt[k-1],net->edgeList[i].pt[k]);
						}
					}
					net->edgeList[i].length_um = len;
				}
				printf("got point positions\n");
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
	int err = 0;
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			net->point[j].used = true;
		}
		// Check for null edge
		if (edge.npts == 2) {
			if (edge.pt[0] == edge.pt[1]) {
				printf("Error: repeated end points on edge: %d  %d\n",i,edge.pt[0]);
				fprintf(fperr,"Error: repeated end points on edge: %d  %d\n",i,edge.pt[0]);
//				err = 1;
			}
			// Check for same end points with different indices
			int j1 = edge.pt[0];
			int j2 = edge.pt[1];
			if (net->point[j1].x == net->point[j2].x 
				&& net->point[j1].y == net->point[j2].y 
				&& net->point[j1].z == net->point[j2].z) {
				printf("edge: %d vertices: %d %d same pos: %f %f %f\n",i,j1,j2,net->point[j1].x,net->point[j1].y,net->point[j1].z);
				fprintf(fperr,"edge: %d vertices: %d %d same pos: %f %f %f\n",i,j1,j2,net->point[j1].x,net->point[j1].y,net->point[j1].z);
//				err = 1;
			}
		}
	}
	fclose(fpam);
	return err;
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int checkUnconnected(NETWORK *net)
{
	int i, kv0, kv1, n0, n1, ii;
	EDGE edge;

	// Step through edges, looking for unconnected edges.
	for (i=0; i<net->ne; i++) {
		if (!net->edgeList[i].used) continue;
		edge = net->edgeList[i];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		n0 = 0;
		n1 = 0;
		for (ii=0; ii<net->ne; ii++) {
			if (!net->edgeList[ii].used) continue;
			if (i == ii) continue;
			if (net->edgeList[ii].vert[0] == kv0 || net->edgeList[ii].vert[1] == kv0) n0=1;
			if (net->edgeList[ii].vert[0] == kv1 || net->edgeList[ii].vert[1] == kv1) n1=1;
		}
		if (n0+n1 == 0) {
			printf("Error: pruner: edge: %d is unconnected: npts: %d vert: %d %d\n",i,edge.npts,kv0,kv1);
			fprintf(fperr,"Error: pruner: edge: %d is unconnected: npts: %d vert: %d %d\n",i,edge.npts,kv0,kv1);
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *input_amfile;
	char drive[32], dir[128],filename[256], ext[32];
	char errfilename[256], output_amfile[256], outfilename[256],output_basename[256];
	int limit_mode, cmgui_flag;
	float limit_value, origin_shift[3];
	NETWORK *net;

	if (argc != 8) {
		printf("Usage: translate input_amfile output_file limit_mode limit_value ddiam dlen cmgui_flag\n");
		fperr = fopen("translate_error.log","w");
		fprintf(fperr,"Usage: translate input_amfile output_file limit_mode limit_value ddiam dlen cmgui_flag\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(outfilename,argv[2]);
	sscanf(argv[3],"%d",&limit_mode);
	sscanf(argv[4],"%f",&limit_value);
	sscanf(argv[5],"%f",&ddiam);
	sscanf(argv[6],"%f",&dlen);
	sscanf(argv[7],"%d",&cmgui_flag);
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_translate.log",output_basename);
	sprintf(output_amfile,"%s.am",output_basename);
//	sprintf(result_file,"%s_prune.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

//	fpout = fopen(result_file,"w");	
	fpout = fopen(outfilename,"w");	
	net = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,net);
	if (err != 0) return 2;

// Don't squeeze - best to keep the same numbering, to make debugging easier
//	err = squeezer();	// must squeeze, or SpatialGraph and CMGUI files are not consistent
//	if (err != 0) return 6;

	if (limit_mode == 0) {
		use_len_diam_limit = false;
		use_len_limit = false;
	} else 	if (limit_mode == 1) {
		use_len_diam_limit = false;
		use_len_limit = true;
		len_limit = limit_value;
	} else 	if (limit_mode == 2) {
		use_len_diam_limit = true;
		use_len_limit = false;
		len_diam_limit = limit_value;
	}

	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;

	err = WriteAmiraFile(output_amfile,input_amfile,net,origin_shift);
	if (err != 0) return 3;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,net,origin_shift);
		if (err != 0) return 4;
	}
	err = EdgeDimensions(net->edgeList,net->point,net->ne);
	if (err != 0) return 5;
	err = CreateDistributions(net);
	if (err != 0) return 6;
	return 0;
}