// Check a .am file

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

#define STR_LEN 128

FILE *fperr, *fpout;

#define EPSILON 0.001
#define PTEQU(a,b) (((a).x==(b).x)  && ((a).y==(b).y) && ((a).z==(b).z))
#define PTEQUIV(a,b) (fabs((a).x-(b).x) < EPSILON && fabs((a).y-(b).y) < EPSILON && fabs((a).z-(b).z) < EPSILON)

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
// This code is faulty - because points can appear multiple times, there are multiple subtractions.
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

	printf("\nWriteAmiraFile: %s\n",amFileOut);
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
	printf("did vertices\n");
	fprintf(fpam,"\n@2\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		fprintf(fpam,"%d %d\n",edge.vert[0],edge.vert[1]);
	}
	printf("did edge vert\n");
	fprintf(fpam,"\n@3\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		fprintf(fpam,"%d\n",edge.npts);
	}
	printf("did edge npts\n");
	fprintf(fpam,"\n@4\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
//		printf("edge: %d  npts: %d\n",i,edge.npts);
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
//			fprintf(fpam,"%6.1f %6.1f %6.1f\n",net->point[j].x,net->point[j].y,net->point[j].z);
			fprintf(fpam,"%6.1f %6.1f %6.1f\n",
				net->point[j].x - origin_shift[0],
				net->point[j].y - origin_shift[1],
				net->point[j].z - origin_shift[2]);
		}
	}
	printf("did points\n");
	fprintf(fpam,"\n@5\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		for (k=0;k<edge.npts;k++) {
			j = edge.pt[k];
			fprintf(fpam,"%6.2f\n",net->point[j].d);
		}
	}
	printf("did diameters\n");
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
// Look for edges with the same vertices
//-----------------------------------------------------------------------------------------------------
int amcheck_a(NETWORK *net)
{
	int ia, ib, npa, npb, kva[2], kvb[2], i, j;
	EDGE edgea, edgeb;
	POINT p0a, p1a, p0b, p1b, p, pa, pb;

	for (ia=0; ia<net->ne; ia++) {
		edgea = net->edgeList[ia];
		npa = edgea.npts;
		kva[0] = edgea.vert[0];
		kva[1] = edgea.vert[1];
		if (ia%100 == 0)
			printf("Edgea: %8d  kva: %8d %8d  npa: %3d\n", ia,kva[0],kva[1],npa);
		p0a = net->point[edgea.pt[1]];
		p1a = net->point[edgea.pt[npa-2]];
		for (ib=0; ib<net->ne; ib++) {
			if (ib == ia) continue;
			edgeb = net->edgeList[ib];
			npb = edgeb.npts;
			kvb[0] = edgeb.vert[0];
			kvb[1] = edgeb.vert[1];
			p0b = net->point[edgeb.pt[1]];
			p1b = net->point[edgeb.pt[npb-2]];

			// This produced no hits with A75.am
			if (kvb[0] == kva[0] && kvb[1] == kva[1]) {
				if (PTEQUIV(p0a,p0b) || PTEQUIV(p1a,p1b)) {
					printf("First two points repeated: %d %d\n",ia,ib);
					printf("edgea pts: \n");
					for (i=0;i<npa;i++) {
						j = edgea.pt[i];
						p = net->point[j];
						printf("%6d %8.1f %8.1f %8.1f\n",j,p.x,p.y,p.z);
					}
					printf("edgeb pts: \n");
					for (i=0;i<npb;i++) {
						j = edgeb.pt[i];
						p = net->point[j];
						printf("%6d %8.1f %8.1f %8.1f\n",j,p.x,p.y,p.z);
					}
					exit(1);
				}
			}
			if (kvb[0] == kva[1] && kvb[1] == kva[0]) {
				if (PTEQUIV(p0a,p1b) || PTEQUIV(p1a,p0b)) {
					printf("First two points repeated: %d %d\n",ia,ib);
					printf("edgea pts: \n");
					for (i=0;i<npa;i++) {
						j = edgea.pt[i];
						p = net->point[j];
						printf("%6d %8.1f %8.1f %8.1f\n",j,p.x,p.y,p.z);
					}
					printf("edgeb pts: \n");
					for (i=0;i<npb;i++) {
						j = edgeb.pt[i];
						p = net->point[j];
						printf("%6d %8.1f %8.1f %8.1f\n",j,p.x,p.y,p.z);
					}
					exit(1);
				}
			}
			// This produced no hits with A75.am
			for (i=1;i<npa-1;i++) {
				pa = net->point[edgea.pt[i]];
				for (j=1;j<npb-1;j++) {
					pb = net->point[edgeb.pt[j]];
					if (PTEQU(pa,pb)) {
						printf("    hit: %8d %3d %3d %8.1f %8.1f %8.1f\n",ib,i,j,pa.x,pa.y,pa.z);
						exit(1);
					}
				}
			}
		}
				//if (edge0.npts == edge1.npts) {
				//	printf("repeated edges: %d %d kv[0]: %d %d npts: %d\n",ie0,ie1,kv[0],edge1.vert[0],edge0.npts);
				//} else {
				//	printf("unequal npts: %d %d kv[0]: %d %d npts: %d %d\n",ie0,ie1,kv[0],edge1.vert[0],edge0.npts,edge1.npts);
				//}

//			}
	}		
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int amcheck_b(NETWORK *net)
{
	int iv, ie, kv, k, ncon, count[2], maxl;
	EDGE *ep;

	for (iv=0; iv<net->nv; iv++) {
		net->vertex[iv].nlinks = 0;
	}
	maxl = 0;
	for (ie=0; ie<net->ne; ie++) {
		ep = &net->edgeList[ie];
		for (k=0; k<2; k++) {
			kv = ep->vert[k];
			net->vertex[kv].link[net->vertex[kv].nlinks] = ie;
			net->vertex[kv].nlinks++;
			if (net->vertex[kv].nlinks > maxl) {
				maxl = net->vertex[kv].nlinks;
				printf("maxl: %d\n",maxl);
			}
			if (net->vertex[kv].nlinks == MAX_LINKS) {
				printf("ERROR: number of links exceeds MAX_LINKS\n");
				exit(1);
			}
		}
	}
	count[0] = count[1] = 0;
	for (ie=0; ie<net->ne; ie++) {
		ep = &net->edgeList[ie];
		ncon = 0;
		for (k=0; k<2; k++) {
			if (net->vertex[ep->vert[k]].nlinks > 1) ncon++;
		}
		if (ncon < 2) {
			if (ncon == 0) {
				printf("edge: %d has %d connections: %d %d\n",ie,ncon,ep->vert[0],ep->vert[1]);
				ep->used = false;
			}
			count[ncon]++;
		}
	}
	printf("edges with 0, 1 connection: %d %d\n",count[0],count[1]);
	printf("max number of nlinks: %d\n",maxl);
	kv = 0;
	printf("vertex %d nlinks: %d\n",kv,net->vertex[kv].nlinks);
	kv = 5016;
	printf("vertex %d nlinks: %d\n",kv,net->vertex[kv].nlinks);
	kv = 1;
	printf("vertex %d nlinks: %d\n",kv,net->vertex[kv].nlinks);
	kv = 244081;
	printf("vertex %d nlinks: %d\n",kv,net->vertex[kv].nlinks);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// A smoothed network NP1 is generated from NP0
// Edges and vertices are unchanged, but the number of points on an edge is reduced.
//-----------------------------------------------------------------------------------------------------
int amsmooth(NETWORK *net0, NETWORK *net1)
{
	int ie, iv, ip0, ip1, ne;
	EDGE edge0;

	printf("amsmooth\n");
//	net1->ne = net0->ne;
	net1->nv = net0->nv;
	net1->vertex = (VERTEX *)malloc(net1->nv*sizeof(VERTEX));
	net1->edgeList = (EDGE *)malloc(net0->ne*sizeof(EDGE));
	net1->point = (POINT *)malloc(net0->np*sizeof(POINT));
	for (iv=0; iv<net1->nv; iv++) {
		net1->vertex[iv].point = net0->vertex[iv].point;
	}
	net1->np = 0;
	ne = 0;
	for (ie=0; ie<net0->ne; ie++) {
		edge0 = net0->edgeList[ie];
		if (!edge0.used) continue;
//		printf("edge0: %d %d\n",ie,edge0.npts);
		net1->edgeList[ne].vert[0] = edge0.vert[0];
		net1->edgeList[ne].vert[1] = edge0.vert[1];
		net1->edgeList[ne].pt = (int *)malloc(net0->edgeList[ie].npts*sizeof(int));
		net1->edgeList[ne].used = true;
		int kfrom = edge0.pt[0];
		ip1 = 0;
		net1->point[net1->np] = net0->point[net0->edgeList[ie].pt[0]];
		net1->edgeList[ne].pt[0] = net1->np;
		net1->np++;
		ip1++;
		for (ip0=1; ip0<edge0.npts; ip0++) {
			int k2 = edge0.pt[ip0];
			double d = dist(net0,kfrom,k2);
			if (d > 0.5*(net0->point[kfrom].d/2 + net0->point[k2].d/2) || ip0 == edge0.npts-1) {
				net1->point[net1->np] = net0->point[net0->edgeList[ie].pt[ip0]];
				net1->point[net1->np].used = true;
				net1->edgeList[ne].pt[ip1] = net1->np;
				net1->np++;
				ip1++;
				kfrom = k2;
			}
		}
		net1->edgeList[ne].npts = ip1;
		ne++;
	}
	net1->ne = ne;
	printf("done: np: %d\n",net1->np);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// This code assumes that the selection of a region of interest (a cube) is carried out on the 3D image.
// This means that the cube centre (xc,yc,zc) and the width (diameter) are all specified in voxel
// coordinates.  The values are converted to um by multiplying by voxelsize in um.  This is necessary
// in order to specify the corresponding region in the Amira network (in which the distance unit is um).
// To enable comparison of the zoomed network file with the cropped 3D image, either the network must
// be scaled to use voxelsize as the distance unit (in which case direct comparison is possible in Amira),
// or the network file coordinates can be left in units of um (in which case the voxelsize must be specified
// when the image file is imported into Amira).  The second option has been adopted.
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *input_amfile;
	char drive[32], dir[128],filename[256], ext[32];
	char errfilename[256], output_amfile[256], result_file[256];
	char output_basename[256];
	float origin_shift[3];
	int cmgui_flag=1;
	NETWORK *NP0, *NP1;

	if (argc != 2) {
		printf("Usage: am_check input_amfile\n");
		fperr = fopen("am_check_error.log","w");
		fprintf(fperr,"Usage: am_check input_amfile\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(output_amfile,"checked.am");
	_splitpath(output_amfile,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_am_check.log",output_basename);
	sprintf(result_file,"%s_am_check.out",output_basename);
	fperr = fopen(errfilename,"w");

	fpout = fopen(result_file,"w");	
	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	NP1 = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;
//	amcheck_a(NP0);
	amcheck_b(NP0);
	return 0;

	amsmooth(NP0,NP1);

	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;
//	fprintf(fperr,"origin_shift: %f %f %f\n",origin_shift[0],origin_shift[1],origin_shift[2]);

	err = WriteAmiraFile(output_amfile,input_amfile,NP1,origin_shift);
	if (err != 0) return 4;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,NP1,origin_shift);
		if (err != 0) return 5;
	}
	return 0;
}