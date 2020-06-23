
// Modifications:
// 1. Start path at random position on a segment
// 2. Allow jump to nearby fibre

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
#include <random> // http://en.cppreference.com/w/cpp/numeric/random

#include "network.h"

int WriteCmguiData(char *basename, NETWORK *net, float origin_shift[]);

#define STR_LEN 128
#define MAX_LINKS 16
#define NB 10
#define NX NB
#define NY NB
#define NZ NB

struct fibre_str
{
	double L_actual, L_direct;
	double u[3];
	int kv[2];	// vertex indicies
	int pt[2];	// end point indicies - this is redundant, since pt = kv
	int nlinks[2];
	int link[2][MAX_LINKS];
	int near_fibre[2][2];		// [kend][0] is near fibre index (-1 means no near fibre), [kend][1] is end of near fibre, kend = (0,1) is end of this fibre
	double d[2];				// distance of near fibre node from this fibre node
};
typedef fibre_str FIBRE;	

int nfibres;
FIBRE *fibre;

float uvec[26][3];
POINT centre;

#define NTIMES 600
#define NDATAPTS 24

//int pathcnt[NPATHS];
//double path[NPATHS][NTIMES][3];
int *pathcnt;
POINT **path;

FILE *fpout, *fperr;

double shrink_factor;	// = 1.25;	// to correct for sample shrinkage
double max_len;
double min_len = 3;		// to compute fibre length distribution with a lower limit on length;

int npow = 4;					// probability of selecting a branch depends on power (npow) of cosine(turning angle)
int ntrials = 5000;				// number of random start locations
double start_radius = 100;		// radius of sphere within which paths start (um)
double deadend_radius = 0;		// radius of sphere within which deadends will be counted
bool save_paths;
int npaths;
double mean_speed=12;			// mean cell speed (um/min)
double CV = 0.10;				// coeff of variation of cell speed = std dev./mean
double deltat = 5;				// data point time interval (min)

#define EPSILON 0.001
#define PTEQU(a,b) (((a).x==(b).x)  && ((a).y==(b).y) && ((a).z==(b).z))
#define PTEQUIV(a,b) (fabs((a).x-(b).x) < EPSILON && fabs((a).y-(b).y) < EPSILON && fabs((a).z-(b).z) < EPSILON)

double PI = 4*atan(1.0);

double xmin, xmax, ymin, ymax, zmin, zmax;
double DX, DY, DZ;
int MAXBLOCK;
int *blocks;
int counter[NX][NY][NZ];
int ndead;
bool vertices_only = false;

bool use_len_limit, use_len_diam_limit;
float len_limit, len_diam_limit;
float ddiam, dlen;
#define NBOX 100

#define N1 MAXBLOCK
#define N2 2*MAXBLOCK
#define N3 2*MAXBLOCK*NX
#define N4 2*MAXBLOCK*NX*NY
#define B(j,kfe,ix,iy,iz) blocks[(iz)*N4 + (iy)*N3 + (ix)*N2 + (kfe)*N1 + (j)]	

int vcounter[NX][NY][NZ];
int *vertexList[NX][NY][NZ];


bool jumpy = true;		// jumping between fibres
bool v_jumpy = false;	// jumping vertex -> vertex
float jumpprobfactor = 1;
int njumpcalls;
double djumpsum;

using namespace std;

seed_seq seq;
mt19937 gen(seq);

uniform_real_distribution<double> uni_dist( 0.0, 1.0 ) ;

float pdist2(POINT p1, POINT p2)
{
	return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int SetupVertexLists(NETWORK *net)
{
	int ix, iy, iz, k;
	POINT p1;

	printf("SetupVertexLists\n");
	xmin = 1.0e10;
	ymin = 1.0e10;
	zmin = 1.0e10;
	xmax = -1.0e10;
	ymax = -1.0e10;
	zmax = -1.0e10;
	for (k=0; k<net->np; k++) {
		p1 = net->point[k];
		xmin = MIN(xmin,p1.x);
		xmax = MAX(xmax,p1.x);
		ymin = MIN(ymin,p1.y);
		ymax = MAX(ymax,p1.y);
		zmin = MIN(zmin,p1.z);
		zmax = MAX(zmax,p1.z);
	}
	DX = (xmax-xmin)/NX+1;
	DY = (ymax-ymin)/NY+1;
	DZ = (zmax-zmin)/NZ+1;

	fprintf(fpout,"\nDimension ranges and grid spacings\n");
	fprintf(fpout,"X: %f %f  NX: %d DX: %f\n",xmin,xmax,NX,DX);
	fprintf(fpout,"Y: %f %f  NY: %d DY: %f\n",ymin,ymax,NY,DY);
	fprintf(fpout,"Z: %f %f  NZ: %d DZ: %f\n",zmin,zmax,NZ,DZ);
	fflush(fpout);

	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				vcounter[ix][iy][iz] = 0;
			}
		}
	}
	for (k=0; k<net->nv; k++) {
		p1 = net->point[k];
		if (!p1.used) continue;
		ix = (int)((p1.x - xmin)/DX);
		iy = (int)((p1.y - ymin)/DY);
		iz = (int)((p1.z - zmin)/DZ);
		vcounter[ix][iy][iz]++;
	}
	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				vertexList[ix][iy][iz] = (int *)malloc(vcounter[ix][iy][iz]*sizeof(int));
			}
		}
	}

	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				vcounter[ix][iy][iz] = 0;
			}
		}
	}
	for (k=0; k<net->nv; k++) {
		p1 = net->point[k];
		if (!p1.used) continue;
		ix = (int)((p1.x - xmin)/DX);
		iy = (int)((p1.y - ymin)/DY);
		iz = (int)((p1.z - zmin)/DZ);
		*(vertexList[ix][iy][iz] + vcounter[ix][iy][iz]) = k;
		vcounter[ix][iy][iz]++;
		net->vertex[k].nlinks = 0;
	}
	for (k=0; k<net->ne; k++) {
		for (int iv=0; iv<2; iv++) {							// edge ends
			int kv = net->edgeList[k].vert[iv];					// end vertexes
			net->vertex[kv].link[net->vertex[kv].nlinks] = k;	// end vertex is connected to edge k
			net->vertex[kv].nlinks++;
		}
	}

	printf("did SetupVertexLists\n");

	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Uses fibre[].nlinks[], fibre[].link[][]
//-----------------------------------------------------------------------------------------------------
int SetupPointLists(NETWORK *net)
{
	int PLcounter[NX][NY][NZ];
	int *pointList[NX][NY][NZ];
	int i, j, k, j2, ix, iy, iz, npts;
	POINT p1, p2;
	EDGE *edge1, edge2;

	printf("SetupPointLists\n");

	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				PLcounter[ix][iy][iz] = 0;
			}
		}
	}
	for (k=0; k<net->np; k++) {
		p1 = net->point[k];
		if (!p1.used) continue;
		ix = (int)((p1.x - xmin)/DX);
		iy = (int)((p1.y - ymin)/DY);
		iz = (int)((p1.z - zmin)/DZ);
		PLcounter[ix][iy][iz]++;
	}
	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				pointList[ix][iy][iz] = (int *)malloc(PLcounter[ix][iy][iz]*sizeof(int));
			}
		}
	}
	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				PLcounter[ix][iy][iz] = 0;
			}
		}
	}
	for (k=0; k<net->np; k++) {
		p1 = net->point[k];
		if (!p1.used) continue;
		ix = (int)((p1.x - xmin)/DX);
		iy = (int)((p1.y - ymin)/DY);
		iz = (int)((p1.z - zmin)/DZ);
		*(pointList[ix][iy][iz] + PLcounter[ix][iy][iz]) = k;
		PLcounter[ix][iy][iz]++;
	} 

	// Check points in [5][4][5] for test.am
//	ix = 5;
//	iy = 4;
//	iz = 5;
//	fprintf(fpout,"points in block: %d %d %d\n",ix,iy,iz);
//	for (i=0; i<PLcounter[ix][iy][iz]; i++) {
//		j = *(pointList[ix][iy][iz] + i);
//		p1 = net->point[j];
//		fprintf(fpout,"%d  %f %f %f\n",j,p1.x,p1.y,p1.z);
//	}

	fprintf(fpout,"Find nearest points\n");
	POINT *plist = NULL;
	// Now find the nearest point to each point
	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			printf(".");
			for (iz=0; iz < NZ; iz++) {
				int n = PLcounter[ix][iy][iz];
				// Create the list of points in this block
				if (plist) free(plist);
				plist = (POINT *)malloc(n*sizeof(POINT));
				for (i=0; i<n; i++) {
					j = *(pointList[ix][iy][iz] + i);
					plist[i] = net->point[j];
				}
				for (int i1=0; i1<n; i1++) {
					int j1 = *(pointList[ix][iy][iz] + i1);
					p1 = plist[i1];
					int iedge1 = p1.iedge;
					int i2min = -1;
					float d2min = 1.0e10;
					for (int i2=0; i2<n; i2++) {
						if (i1 == i2) continue;
						j2 = *(pointList[ix][iy][iz] + i2);
						p2 = plist[i2];
						int iedge2 = p2.iedge;
						if (iedge1 == iedge2) continue;		// ignore other points on this edge

						// Also need to exclude edges that have a connection to this edge,
						// i.e. connected to the vertex at either end of iedge1
						bool connected = false;
						for (int iend = 0; iend<2; iend++) {
							int nlinks = fibre[iedge1].nlinks[iend];
							for (int ilink=0; ilink<nlinks; ilink++) {
								int ifib2 = fibre[iedge1].link[iend][ilink];
								if (ifib2 == iedge2) {
									connected = true;
									continue;
								}
							}
							if (connected) continue;
						}
						if (connected) continue;
						
						float d2 = pdist2(p1,p2);
						if (d2 < 1) printf("d2 < 1: iedge1: %d iedge2: %d\n",iedge1,iedge2);
						if (d2 < d2min) {
							d2min = d2;
							i2min = i2;
						}
					}
					if (i2min == -1) {	// no other edges in this block
						net->point[j1].nearestpt = -1;	// no edges nearby
						net->point[j1].d2near = 0;
					} else {
						// i2min is the index in plist[]
						net->point[j1].nearestpt = *(pointList[ix][iy][iz] + i2min);	// this is the pt index
						net->point[j1].d2near = d2min;
					}
				}
			}
		}
	}
	printf("\ndid set up nearestpt\n");
	int jcount[100] = {0};

	// Now look at each edge, and select the point with the smallest d2near
	for (k=0; k<net->ne; k++) {
		edge1 = &net->edgeList[k];
		if (!edge1->used) continue;
//		printf("edge: %d %p\n",k,edge1);
		npts = edge1->npts;
//		printf("npts: %d\n",npts);
		edge1->jumpable_pt = -1;
		edge1->dnearest = 0;
//		if (npts < 6) {
//			continue;
//		}
		float d2min = 1.0e10;
		int imin = -1;
//		for (i=3; i<=npts-3; i++) {
		for (i=1; i<=npts-2; i++) {
			int ip = edge1->pt[i];
//			printf("i: %d ip: %d\n",i,ip);
			p1 = net->point[ip];
			if (!p1.used) {
				printf("Error: in SetupPointLists: unused point used: %d\n",ip);
				exit(1);
			}
//			printf("d2near: %f\n",p1.d2near);
			if (p1.d2near == 0) continue;
			if (p1.d2near < d2min) {
				imin = i;
				d2min = p1.d2near;
//				printf("imin: %d d2min: %f\n",imin,d2min);
			}
		}
		if (imin == -1) continue;
		edge1->jumpable_pt = imin;	// index of edge1->pt[]
		int ip = edge1->pt[imin];
//		printf("imin: %d ip: %d\n",imin,ip);
		p1 = net->point[ip];
		float dj = sqrt(p1.d2near);
		float dmin = sqrt(d2min);
		edge1->dnearest = dmin;
//		if (edge1->dnearest == 1.0) {
//			printf("d2min: %f dnearest: %f\n",d2min,edge1->dnearest);
//			exit(1);
//		}
//		printf("d2min: %f p1.d2near: %f\n",d2min,p1.d2near);
//		printf("dmin: %f dj: %f\n",dmin,dj);
//		printf("dnearest: %f\n",edge1->dnearest);
		if (edge1->dnearest != dj) {
//			printf("dnearest (again): %f\n",edge1->dnearest);
			printf("edge: %d imin: %d d2min: %f dnearest: %f dj: %f\n",k,imin,d2min,edge1->dnearest,dj);
		}
//		if (imin > 0 && imin < 100) jcount[imin]++;
//		if (imin > 0 && imin < 100) {
//			int ii = npts - imin - 1;
//			jcount[ii]++;
//		}
//		fprintf(fpout,"edge1: %d  jumpable_pt: %d  npts: %d dnearest: %f\n",k,imin,npts,edge1.dnearest);

		// Check:
		if (imin < 0) {
			// no nearby edge
			continue;
		}
		p1 = net->point[edge1->pt[imin]];
		int ptnear = p1.nearestpt;
		int ienear = net->point[ptnear].iedge;
		edge2 = net->edgeList[ienear];
		// find jpt = index of ptnear on edge2
		int jpt = -1;
		npts = edge2.npts;
		for (i=0; i<npts; i++) {
			if (ienear == 0) printf("edge: %d i: %d pt[i]: %d\n",ienear,i,edge2.pt[i]);
			if (edge2.pt[i] == ptnear) {
				jpt = i;
				break;
			}
		}
		if (jpt < 0) {
			printf("Error: on edge: %d pt: %d ptnear: %d not found on edge2: %d\n",k,edge1->pt[imin],ptnear,ienear);
			printf("pt %d used: %d\n",ptnear,net->point[ptnear]);
			exit(1);
		}
//		exit(1);
	}

	// Now free pointList
	for (ix=0; ix < NX; ix++) {
		for (iy=0; iy < NY; iy++) {
			for (iz=0; iz < NZ; iz++) {
				free(pointList[ix][iy][iz]);
			}
		}
	}
	printf("Freed pointList\n");
	printf("did SetupPointLists\n");
		
	return 0;
}


//-----------------------------------------------------------------------------------------------------
// This assumes that edge dimensions (average diameter and length) have already been computed.
//-----------------------------------------------------------------------------------------------------
int CreateDistributions(NETWORK *net)
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double lsegadbox[NBOX];
	double ad, len, dave, ltot, dsum, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam;		// average point diameter
	double ave_seg_diam;	// average vessel diameter
	double ave_lseg_diam;	// length-weighted average vessel diameter
    int ie, ip, k, ka, kp, ndpts, nlpts, ndtot, nsegdtot, nptoowide, natoowide;
	EDGE edge;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lsegadbox[k] = 0;
		lvbox[k] = 0;
	}
	if (use_len_diam_limit) {
		printf("\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
		fprintf(fpout,"\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
	} else if (use_len_limit) {
		printf("\nUsing length lower limit = %6.1f um\n\n",len_limit);
		fprintf(fpout,"\nUsing length lower limit = %6.1f um\n\n",len_limit);
	}
	printf("\nDistributions\n");
	printf("shrinkage compensation factor: %6.2f\n\n",shrink_factor);
	fprintf(fpout,"\nDistributions\n");
	fprintf(fpout,"shrinkage compensation factor: %6.2f\n\n",shrink_factor);
	printf("Compute diameter distributions (length weighted)\n");
	fprintf(fpout,"Compute diameter distributions (length weighted)\n");

	// Diameters
//	ddiam = 0.5;
	ndtot = 0;
	nsegdtot = 0;
	lsegdtot = 0;
    nptoowide = 0;
    natoowide = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	ave_lseg_diam = 0;
	volume = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
//		printf("ie: %d\n",ie);
		if (!edge.used) continue;
		len = shrink_factor*edge.length_um;
		dave = edge.segavediam;
//		printf("ie,len,dave,npts: %d %f %f %d\n",ie,len,dave,edge.npts);
		if (use_len_limit && len < len_limit) continue;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = net->point[kp].d;
			ave_pt_diam += ad;
			if (ad < 0.001 || ad > 200) {
				printf("Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fpout,"Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				return 1;
			}
			ka = int(ad/ddiam + 0.5);
			if (ka >= NBOX) {
 //               printf("Vessel too wide (point): d: %f k: %d\n",ad,ka);
 //               fprintf(fpout,"Vessel too wide (point): d: %f k: %d\n",ad,ka);
				nptoowide++;
				continue;
			}
			adbox[ka]++;
			ndtot++;
		}
//		net->edgeList[ie].length_um = len;	// already computed 
		if (use_len_diam_limit && len/ad < len_diam_limit) continue;
		ave_seg_diam += dave;
		ave_lseg_diam += dave*len;
//		printf("ie: %6d dave,len,ave_lseg_diam: %8.1f %8.1f %10.1f\n",ie,dave,len,ave_lseg_diam);
//		fprintf(fpout,"ie: %6d dave,len,ave_lseg_diam: %8.1f %8.1f %10.1f\n",ie,dave,len,ave_lseg_diam);
		if (dave < 0.001 || dave > 200) {
			printf("Bad segment diameter: edge: %d ad: %f\n",ie,dave);
			fprintf(fpout,"Zero segment diameter: edge: %d ad: %f\n",ie,dave);
			return 1;
		}
		ka = int(dave/ddiam + 0.5);
		if (ka >= NBOX) {
//			printf("Vessel too wide (segment ave): d: %f k: %d\n",dave,ka);
//			fprintf(fpout,"Vessel too wide (segment ave): d: %f k: %d\n",dave,ka);
            natoowide++;
            continue;
		}
		segadbox[ka]++;
		nsegdtot++;
		lsegadbox[ka] += len;
		lsegdtot += len;
		volume += len*PI*(dave/2)*(dave/2);
	}
    if (natoowide > 0) {
        printf("Number of segment average diameters too big for NBOX: %d\n",natoowide);
        fprintf(fpout,"Number of segment average diameters too big for NBOX: %d\n",natoowide); 
    }
    if (nptoowide > 0) {
        printf("Number of point diameters too big for NBOX: %d\n",nptoowide);
        fprintf(fpout,"Number of point diameters too big for NBOX: %d\n",nptoowide);
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
	printf("\nCompute length distributions:\n");
	fprintf(fpout,"\nCompute length distributions:\n");
	// Lengths
//	dlen = 1;
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		len = shrink_factor*edge.length_um;
		k = int(len/dlen + 0.5);
		if (use_len_limit && len <= len_limit) continue;
		ad = edge.segavediam;
		if (use_len_diam_limit && len/ad < len_diam_limit) continue;
		if (k >= NBOX) {
			printf("Edge too long for boxes: len: %d  %6.1f  k: %d\n",ie,len,k);
			fprintf(fpout,"Edge too long for boxes: len: %d  %6.1f  k: %d\n",ie,len,k);
			continue;
		}
		lvbox[k]++;
		ave_len += len;
		ltot++;
	}

	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	ave_lseg_diam /= lsegdtot;
	fprintf(fpout,"Vertices: %d\n",net->nv);
	fprintf(fpout,"Points: %d\n",net->np);
	fprintf(fpout,"Vessels: %d\n",net->ne);
	fprintf(fpout,"ltot: %d\n",int(ltot));	// count of vessels used, excluding those dropped
	printf("Average pt diameter: %6.2f vessel diameter: %6.2f length-weighted: %6.2f\n",ave_pt_diam, ave_seg_diam, ave_lseg_diam);
	fprintf(fpout,"Ave_pt_diam: %6.2f\n",ave_pt_diam);
	fprintf(fpout,"Ave_vessel_diam: %6.2f\n",ave_seg_diam);
	fprintf(fpout,"Ave_wgt_vessel_diam: %6.2f\n",ave_lseg_diam);
	printf("Average vessel length: %6.1f\n",ave_len/ltot);
	fprintf(fpout,"Ave_vessel_len: %6.1f\n",ave_len/ltot);
	printf("Vessel_volume: %10.0f um3\n\n",volume);
	fprintf(fpout,"Vessel_volume: %10.0f um3\n\n",volume);

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

	for (k=NBOX-1; k>=0; k--) {
		if (lvbox[k] > 0) break;
	}
	nlpts = k+2;
	fprintf(fpout,"Vessel length distribution\n");
	fprintf(fpout,"   um    number  fraction\n");
	for (k=0; k<nlpts; k++) {
		fprintf(fpout,"%6.2f %8d %9.5f\n",k*dlen,lvbox[k],lvbox[k]/ltot);
		fflush(fpout);
	}
	fprintf(fpout,"Computed distributions\n");
	fflush(fpout);
	return 0;
} 

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float distance(NETWORK *net, int k1, int k2)
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
	fprintf(fpam,"# Created by conduit_analyse from: %s by joining up dead ends\n",amFileIn);
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
// Write Amira SpatialGraph file with some edges tagged as unused
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile_clean(char *amFileOut, char *amFileIn, NETWORK *net, float origin_shift[])
{
	int i, k, j, new_np, new_nv, new_ne, iv;
	int *new_iv, *old_iv;
	EDGE edge;

	printf("\nWriteAmiraFile_clean: %s\n",amFileOut);
	fprintf(fpout,"\nWriteAmiraFile_clean: %s\n",amFileOut);
	new_iv = (int *)malloc(net->nv*sizeof(int));
	old_iv = (int *)malloc(net->nv*sizeof(int));
	
	// Reset used flags on vertexes
	for (i=0; i<net->nv; i++)
		net->vertex[i].used = false;
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			if (vertices_only) net->edgeList[i].npts = 2;	// added
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				net->vertex[iv].used = true;
			}
		}
	}

	iv = 0;
	for (i=0; i<net->nv; i++) {
		if (net->vertex[i].used) {
			new_iv[i] = iv;
			old_iv[iv] = i;
			iv++;
		}
	}
	new_nv = iv;
	
	new_np = 0;
	new_ne = 0;
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			new_ne++;
			new_np += edge.npts;
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				if (!net->vertex[iv].used) {
					printf("Error: used edge: %d has unused vertex: %d %d\n",i,j,iv);
					exit(1);
				}
			}
		} else {
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				if (net->vertex[iv].used) {
					printf("Error: unused edge: %d has used vertex: %d %d\n",i,j,iv);
					exit(1);
				}
			}
		}
	}
//	new_np += new_nv;
	printf("new_ne: %d new_nv: %d new_np: %d\n",new_ne,new_nv,new_np);

	FILE *fpam = fopen(amFileOut,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
	fprintf(fpam,"# Created by conduit_analyse from: %s by removing unused edges\n",amFileIn);
	fprintf(fpam,"\n");
	fprintf(fpam,"define VERTEX %d\n",new_nv);
	fprintf(fpam,"define EDGE %d\n",new_ne);
	fprintf(fpam,"define POINT %d\n",new_np);
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
	for (iv=0;iv<new_nv;iv++) {
		i = old_iv[iv];
		fprintf(fpam,"%6.1f %6.1f %6.1f\n",
			net->vertex[i].point.x - origin_shift[0],
			net->vertex[i].point.y - origin_shift[1],
			net->vertex[i].point.z - origin_shift[2]);
	}
	printf("did vertices\n");
	fprintf(fpam,"\n@2\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			int iv0 = new_iv[edge.vert[0]];
			int iv1 = new_iv[edge.vert[1]];
	//		fprintf(fpam,"%d %d\n",edge.vert[0],edge.vert[1]);
			fprintf(fpam,"%d %d\n",iv0,iv1);
			if (iv0 == iv1) {
				printf("Error: identical vertex indexes: edge: %d vert: %d %d iv0,iv1: %d %d\n",i,edge.vert[0],edge.vert[1],iv0,iv1);
				exit(1);
			}
		}
	}
	printf("did edge vert\n");
	fprintf(fpam,"\n@3\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			fprintf(fpam,"%d\n",edge.npts);
		}
	}
	printf("did edge npts\n");
	fprintf(fpam,"\n@4\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			for (k=0;k<edge.npts;k++) {
				if (vertices_only)
					j = new_iv[edge.vert[k]];	// added
				else
					j = edge.pt[k];
				fprintf(fpam,"%6.1f %6.1f %6.1f\n",
					net->point[j].x - origin_shift[0],
					net->point[j].y - origin_shift[1],
					net->point[j].z - origin_shift[2]);
			}
		}
	}
	printf("did points\n");
	fprintf(fpam,"\n@5\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			for (k=0;k<edge.npts;k++) {
				if (vertices_only)
					j = new_iv[edge.vert[k]];	// added
				else
					j = edge.pt[k];
				fprintf(fpam,"%6.2f\n",net->point[j].d);
			}
		}
	}
	printf("did diameters\n");
	fclose(fpam);
	printf("Completed WriteAmiraFile_clean\n");
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Write Amira SpatialGraph file with some edges tagged as unused, no internal pts
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile_vertices(char *amFileOut, char *amFileIn, NETWORK *net, float origin_shift[])
{
	int i, k, j, new_np, new_nv, new_ne, iv;
	int *new_iv, *old_iv;
	EDGE edge;

	printf("\nWriteAmiraFile_vertices: %s\n",amFileOut);
	fprintf(fpout,"\nWriteAmiraFile_vertices: %s\n",amFileOut);
	new_iv = (int *)malloc(net->nv*sizeof(int));
	old_iv = (int *)malloc(net->nv*sizeof(int));
	// Reset used flags on vertexes
	for (i=0; i<net->nv; i++)
		net->vertex[i].used = false;
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			net->edgeList[i].npts = 2;	// added
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				net->vertex[iv].used = true;
			}
		}
	}

	iv = 0;
	for (i=0; i<net->nv; i++) {
		if (net->vertex[i].used) {
			new_iv[i] = iv;
			old_iv[iv] = i;
			iv++; 
		}
	}
	new_nv = iv;
	
	new_np = 0;
	new_ne = 0;
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			new_ne++;
			new_np += edge.npts;
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				if (!net->vertex[iv].used) {
					printf("Error: used edge: %d has unused vertex: %d %d\n",i,j,iv);
					exit(1);
				}
			}
		} else {
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				if (net->vertex[iv].used) {
					printf("Error: unused edge: %d has used vertex: %d %d\n",i,j,iv);
					exit(1);
				}
			}
		}
	}
//	new_np += new_nv;
	printf("new_ne: %d new_nv: %d new_np: %d\n",new_ne,new_nv,new_np);

	FILE *fpam = fopen(amFileOut,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
	fprintf(fpam,"# Created by conduit_analyse from: %s by removing unused edges\n",amFileIn);
	fprintf(fpam,"\n");
	fprintf(fpam,"define VERTEX %d\n",new_nv);
	fprintf(fpam,"define EDGE %d\n",new_ne);
	fprintf(fpam,"define POINT %d\n",new_np);
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
	for (iv=0;iv<new_nv;iv++) {
		i = old_iv[iv];
		fprintf(fpam,"%6.1f %6.1f %6.1f\n",
			net->vertex[i].point.x - origin_shift[0],
			net->vertex[i].point.y - origin_shift[1],
			net->vertex[i].point.z - origin_shift[2]);
	}
	printf("did vertices\n");
	fprintf(fpam,"\n@2\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			int iv0 = new_iv[edge.vert[0]];
			int iv1 = new_iv[edge.vert[1]];
	//		fprintf(fpam,"%d %d\n",edge.vert[0],edge.vert[1]);
			fprintf(fpam,"%d %d\n",iv0,iv1);
			if (iv0 == iv1) {
				printf("Error: identical vertex indexes: edge: %d vert: %d %d iv0,iv1: %d %d\n",i,edge.vert[0],edge.vert[1],iv0,iv1);
				exit(1);
			}
		}
	}
	printf("did edge vert\n");
	fprintf(fpam,"\n@3\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			fprintf(fpam,"%d\n",edge.npts);
		}
	}
	printf("did edge npts\n");
	fprintf(fpam,"\n@4\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			for (k=0;k<edge.npts;k++) {
				j = new_iv[edge.vert[k]];	// added
//				j = edge.pt[k];
				fprintf(fpam,"%6.1f %6.1f %6.1f\n",
					net->point[j].x - origin_shift[0],
					net->point[j].y - origin_shift[1],
					net->point[j].z - origin_shift[2]);
			}
		}
	}
	printf("did points\n");
	fprintf(fpam,"\n@5\n");
	for (i=0;i<net->ne;i++) {
		edge = net->edgeList[i];
		if (edge.used) {
			for (k=0;k<edge.npts;k++) {
				j = new_iv[edge.vert[k]];	// added
//				j = edge.pt[k];
				fprintf(fpam,"%6.2f\n",net->point[j].d);
			}
		}
	}
	printf("did diameters\n");
	fclose(fpam);
	printf("Completed WriteAmiraFile_clean\n");
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile(char *amFile, NETWORK *net)
{
	int i, j, k, kp, npts;
	int np_used, ne_used, nv_used;
	EDGE edge;
	char line[STR_LEN];

	fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
	fflush(fpout);
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
	net->point = (POINT *)malloc(net->np*sizeof(POINT));	// Note: this assumes that the largest pt index is < net->np
	bool *vused = (bool *)malloc(net->nv*sizeof(bool));		// to check usage of vertices
	printf("Allocated arrays: np: %d nv: %d ne: %d\n",net->np,net->nv,net->ne);
	fprintf(fpout,"Allocated arrays: np: %d nv: %d ne: %d\n",net->np,net->nv,net->ne);
	fflush(fpout);

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
					net->vertex[i].pt = kp;
					net->vertex[i].used = true;
				}
				kp++;
				printf("Got vertices\n");
			} else if (k == 2) {
				for (i=0;i<net->nv; i++)
					vused[i] = false;
				for (i=0;i<net->ne;i++) {
					if (fgets(line, STR_LEN, fpam) == NULL) {
						printf("ERROR reading section @2\n");
						return 1;
					}
					sscanf(line,"%d %d",&net->edgeList[i].vert[0],&net->edgeList[i].vert[1]);
					net->edgeList[i].used = true;
					vused[net->edgeList[i].vert[0]] = true;
					vused[net->edgeList[i].vert[1]] = true;
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
					for (int kk=0;kk<edge.npts;kk++) {
						if (fgets(line, STR_LEN, fpam) == NULL) { 
							printf("ERROR reading section @4\n");
							return 1;
						}
						if (kk > 0 && kk<edge.npts-1) {
							sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[kk] = kp;
							net->edgeList[i].pt_used[kk] = kp;
//							net->point[kp].iedge = i;
//							if (i == 0) printf("on edge: %d kk: %d pt: %d\n",i,kk,kp);
							if (kp == 6) printf("pt: %d is on edge: %d\n",kp,i);
							kp++;
						}
						if (kk > 0) {
							len = len + distance(net,net->edgeList[i].pt[kk-1],net->edgeList[i].pt[kk]);
						}
					}
					net->edgeList[i].length_um = len;
				}
				printf("Got edge points, made lengths, kp: %d\n",kp);
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
	for (j=0; j<net->np; j++) {
		net->point[j].used = false;
	}
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			net->point[j].iedge = i;
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
	// Count used vertexes
	nv_used = 0;
	for (i=0;i<net->nv;i++) {
		if (vused[i]) nv_used++;
	}
	printf("vertexes: nv: %d nv_used: %d\n",net->nv,nv_used);
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
int amcheck(NETWORK *net)
{
	int ia, ib, npa, npb, kva[2], kvb[2], i, j;
	EDGE edgea, edgeb;
	POINT p0a, p1a, p0b, p1b, p, pa, pb;

	for (ia=0; ia<net->ne; ia++) {
		edgea = net->edgeList[ia];
		npa = edgea.npts;
		kva[0] = edgea.vert[0];
		kva[1] = edgea.vert[1];
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
// A smoothed network net1 is generated from net0
// Edges and vertices are unchanged, but the number of points on an edge is reduced.
//-----------------------------------------------------------------------------------------------------
int amsmooth(NETWORK *net0, NETWORK *net1)
{
	int ie, iv, ip0, ip1;
	EDGE edge0;

	printf("amsmooth\n");
	net1->ne = net0->ne;
	net1->nv = net0->nv;
	net1->vertex = (VERTEX *)malloc(net1->nv*sizeof(VERTEX));
	net1->edgeList = (EDGE *)malloc(net1->ne*sizeof(EDGE));
	net1->point = (POINT *)malloc(net0->np*sizeof(POINT));
	for (iv=0; iv<net1->nv; iv++) {
		net1->vertex[iv].point = net0->vertex[iv].point;
	}
	net1->np = 0;
	for (ie=0; ie<net1->ne; ie++) {
		edge0 = net0->edgeList[ie];
//		printf("edge0: %d %d\n",ie,edge0.npts);
		net1->edgeList[ie].vert[0] = edge0.vert[0];
		net1->edgeList[ie].vert[1] = edge0.vert[1];
		net1->edgeList[ie].pt = (int *)malloc(net0->edgeList[ie].npts*sizeof(int));
		net1->edgeList[ie].used = true;
		int kfrom = edge0.pt[0];
		ip1 = 0;
		net1->point[net1->np] = net0->point[net0->edgeList[ie].pt[0]];
		net1->edgeList[ie].pt[0] = net1->np;
		net1->np++;
		ip1++;
		for (ip0=1; ip0<edge0.npts; ip0++) {
			int k2 = edge0.pt[ip0];
			float d = distance(net0,kfrom,k2);
			if (d > 0.5*(net0->point[kfrom].d/2 + net0->point[k2].d/2) || ip0 == edge0.npts-1) {
				net1->point[net1->np] = net0->point[net0->edgeList[ie].pt[ip0]];
				net1->point[net1->np].used = true;
				net1->edgeList[ie].pt[ip1] = net1->np;
				net1->np++;
				ip1++;
				kfrom = k2;
			}
		}
		net1->edgeList[ie].npts = ip1;
	}
	printf("done: np: %d\n",net1->np);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void make26directions()
{
	int idir, ix, iy, iz;
	float r;

	idir = 0;
	for (ix=-1; ix<2; ix++) {
		for (iy=-1; iy<2; iy++) {
			for (iz=-1; iz<2; iz++) {
				if (ix==0 && iy==0 && iz==0) continue;
				r = sqrt(float(ix*ix + iy*iy + iz*iz));
				uvec[idir][0] = ix/r;
				uvec[idir][1] = iy/r;
				uvec[idir][2] = iz/r;
				idir++;
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double dotproduct(double *u1, double *u2)
{
	double res=0;
	for (int i=0; i<3; i++) {
		res += u1[i]*u2[i];
	}
	return res;
}

//-----------------------------------------------------------------------------------------------------
// Get angle (degrees) between fibres kf1 and kf2 at vertex iv
//-----------------------------------------------------------------------------------------------------
double getAngle(int iv, int kf1, int kf2)
{
	int usign1, usign2;
	double *u1, *u2, L1, L2, dot, theta;

	L1 = fibre[kf1].L_direct;
	u1 = fibre[kf1].u;
	if (fibre[kf1].kv[0] == iv)
		usign1 = -1;
	else
		usign1 = 1;
	L2 = fibre[kf2].L_direct;
	u2 = fibre[kf2].u;
	if (fibre[kf2].kv[0] == iv)
		usign2 = 1;
	else
		usign2 = -1;
	dot = usign1*usign2*dotproduct(u1,u2);
	double rad;
	if (dot >= 1) {
		rad = 0;
	} else if (dot <= -1) {
		rad = 0.999*PI;
	} else {
		rad = acos(dot);
	}
	theta = 180.*rad/PI;	// turning angle
	theta = 180 - theta;	// fibre-fibre angle
	return theta;
}

//-----------------------------------------------------------------------------------------------------
// For block (ix,iy,iz)
// for j = 0,..,counter[ix][iy][iz]
// B(j,0,ix,iy,iz)  = index of fibre with an end in the block
// B(j,1,ix,iy,iz)  = fibre end in the block (0 = end 1, 1 = end 2)
// Note that a fibre may appear twice in the list for a block, if both ends are within the block
//-----------------------------------------------------------------------------------------------------
void setupBlockLists(NETWORK *net)
{
	POINT p1, p2;
	int ix, iy, iz, k;

	printf("setupBlockLists: nfibres: %d\n",nfibres);
	for (ix=0;ix<NX;ix++) {
		for (iy=0;iy<NY;iy++) {
			for (iz=0;iz<NZ;iz++) {
				counter[ix][iy][iz] = 0;
			}
		}
	}
	for (k=0; k<nfibres; k++) {
		int npts = net->edgeList[k].npts;
		p1 = net->point[net->edgeList[k].pt[0]];
		ix = (int)((p1.x - xmin)/DX);
		iy = (int)((p1.y - ymin)/DY);
		iz = (int)((p1.z - zmin)/DZ);
		B(counter[ix][iy][iz],0,ix,iy,iz) = k;
		B(counter[ix][iy][iz],1,ix,iy,iz) = 0;	// end 0
		counter[ix][iy][iz]++;
		//if (counter[ix][iy][iz] == MAXBLOCK) {
		//	printf("MAXBLOCK exceeded at: fibre k: %d  block: %d %d %d\n",k,ix,iy,iz);
		//	exit(5);
		//}
		p2 = net->point[net->edgeList[k].pt[npts-1]];
		ix = (int)((p2.x - xmin)/DX);
		iy = (int)((p2.y - ymin)/DY);
		iz = (int)((p2.z - zmin)/DZ);
		B(counter[ix][iy][iz],0,ix,iy,iz) = k;
		B(counter[ix][iy][iz],1,ix,iy,iz) = 1;	// end 1
		counter[ix][iy][iz]++;
	}
	printf("setupBlockLists\n");
//	fprintf(fpout,"setupBlockLists\n");
	int maxcount = 0;
	for (ix=0;ix<NX;ix++) {
		for (iy=0;iy<NY;iy++) {
			for (iz=0;iz<NZ;iz++) {
//				fprintf(fpout,"ix,iy,iz,counter: %3d %3d %3d %6d\n",ix,iy,iz,counter[ix][iy][iz]);
				maxcount = MAX(maxcount,counter[ix][iy][iz]);
			}
		}
	}
	printf("Max block count: %d\n",maxcount);
	fprintf(fpout,"Max block count: %d\n",maxcount);
	for (ix=0;ix<NX;ix++) {
		for (iy=0;iy<NY;iy++) {
			for (iz=0;iz<NZ;iz++) {
//				fprintf(fpout,"Block: %d %d %d  count: %d\n",ix,iy,iz,counter[ix][iy][iz]);
			}
		}
	}
	if (maxcount > MAXBLOCK) {
		printf("Error: MAXBLOCK exceeded\n");
		printf("MAXBLOCK: %d\n",MAXBLOCK);
		fprintf(fpout,"Error: MAXBLOCK exceeded\n");
		fprintf(fpout,"MAXBLOCK: %d\n",MAXBLOCK);
		ix = 7;
		iy = 8;
		iz = 8;
		POINT p;
//		fprintf(fpout,"block: %d %d %d  %8.2f %8.2f %8.2f\n",ix,iy,iz,xmin+ix*DX,ymin+iy*DX,zmin+iz*DX);
		for (k=0; k<counter[ix][iy][iz]; k++) {
			int npts = net->edgeList[k].npts;
			if (B(k,1,ix,iy,iz) == 0) {
				p = net->point[net->edgeList[k].pt[0]];
			} else {
				p = net->point[net->edgeList[k].pt[npts-1]];
			}
			fprintf(fpout,"fibre end: %6d %8.2f %8.2f %8.2f\n", B(k,0,ix,iy,iz),p.x,p.y,p.z);
		}
		fflush(fpout);

		exit(5);
	}
}

//-----------------------------------------------------------------------------------------------------
// Apply shrink_factor here only to the values that are output, and the length distribution.
// All node locations and fibre dimensions are left unscaled.
//-----------------------------------------------------------------------------------------------------
void MakeFibreList(NETWORK *net)
{
#define NBINS 50
#define NANGS 45
	int i, k, iv, idir, i_scaled, i_limit, nf_limit;
	int ix, iy, iz;
	int nlinks[2];
	POINT p1, p2;
	double L_actual_sum=0, L_direct_sum=0;
	double deltaL = 1.0;
	double v[3], dotsum[26];
	double Lmin, Lmax, Lsum, Lbin[NBINS+1], NBbin[MAX_LINKS], Lbin_scaled[NBINS+1], Lbin_limit[NBINS+1];
	double Angbin[NANGS], dang;
	int nvbin, nangbin;

	printf("MakeFibreList\n");
	fprintf(fpout,"MakeFibreList\n");
	dang = 180./NANGS;
	make26directions();
	nfibres = net->ne;
	
	Lmin = 1.0e10;
	Lmax = 0;
	Lsum = 0;
	for (i=0; i<=NBINS; i++) {
		Lbin[i] = 0;
		Lbin_scaled[i] = 0;
		Lbin_limit[i] = 0;
	}
	for (i=0; i<MAX_LINKS; i++)
		NBbin[i] = 0;
	
	fibre = (FIBRE *)malloc(nfibres*sizeof(FIBRE));
	nf_limit = 0;
	for (k=0; k<nfibres; k++) {
		fibre[k].kv[0] = net->edgeList[k].vert[0];
		fibre[k].kv[1] = net->edgeList[k].vert[1];
		fibre[k].pt[0] = net->edgeList[k].pt[0];
		fibre[k].pt[1] = net->edgeList[k].pt[net->edgeList[k].npts-1];
		fibre[k].L_actual = net->edgeList[k].length_um;
		L_actual_sum += fibre[k].L_actual;
		
		Lmin = MIN(Lmin,fibre[k].L_actual);
		Lmax = MAX(Lmax,fibre[k].L_actual);
		Lsum += fibre[k].L_actual;
		i = (int)(fibre[k].L_actual/deltaL + 0.5);
		i = MIN(i,NBINS);
		Lbin[i] += 1;
		i_scaled = (int)(shrink_factor*fibre[k].L_actual/deltaL + 0.5);
		i_scaled = MIN(i_scaled,NBINS);
		Lbin_scaled[i_scaled] += 1;
		if (shrink_factor*fibre[k].L_actual > min_len) {
			nf_limit++;
			i_limit = (int)(shrink_factor*fibre[k].L_actual/deltaL + 0.5);
			i_limit = MIN(i_limit,NBINS);
			Lbin_limit[i_limit] += 1;
		}
		
		int npts = net->edgeList[k].npts;
		fibre[k].L_direct = distance(net,net->edgeList[k].pt[0],net->edgeList[k].pt[npts-1]);
		L_direct_sum += fibre[k].L_direct;
		p1 = net->point[net->edgeList[k].pt[0]];
		p2 = net->point[net->edgeList[k].pt[npts-1]];
		fibre[k].u[0] = (p2.x - p1.x)/fibre[k].L_direct;
		fibre[k].u[1] = (p2.y - p1.y)/fibre[k].L_direct;
		fibre[k].u[2] = (p2.z - p1.z)/fibre[k].L_direct;
	}
	printf("\nShrinkage compensation factor: %6.2f\n\n",shrink_factor);
	printf("Average L_actual: %6.2f  scaled: %6.2f\n",L_actual_sum/nfibres,shrink_factor*L_actual_sum/nfibres);
	printf("Average L_direct: %6.2f  scaled: %6.2f\n",L_direct_sum/nfibres,shrink_factor*L_direct_sum/nfibres);
	/*
	printf("Min L_actual:     %6.2f  scaled: %6.2f\n",Lmin,shrink_factor*Lmin);
	printf("Max L_actual:     %6.2f  scaled: %6.2f\n",Lmax,shrink_factor*Lmax);
	printf("Total length: %10.0f  scaled: %10.0f\n",Lsum,shrink_factor*Lsum);
	fprintf(fpout,"\nShrinkage compensation factor: %6.2f\n\n",shrink_factor);
	fprintf(fpout,"Average L_actual: %6.2f  scaled: %6.2f\n",L_actual_sum/nfibres,shrink_factor*L_actual_sum/nfibres);
	fprintf(fpout,"Average L_direct: %6.2f  scaled: %6.2f\n",L_direct_sum/nfibres,shrink_factor*L_direct_sum/nfibres);
	fprintf(fpout,"Min L_actual:     %6.2f  scaled: %6.2f\n",Lmin,shrink_factor*Lmin);
	fprintf(fpout,"Max L_actual:     %6.2f  scaled: %6.2f\n",Lmax,shrink_factor*Lmax);
	fprintf(fpout,"Total length:   %8.0f  scaled: %8.0f\n",Lsum,shrink_factor*Lsum);
	fflush(fpout);
	fprintf(fpout,"Length distributions:\n");
	fprintf(fpout,"(Bin size: %3.1f)\n",deltaL);
	fprintf(fpout,"(the limited distribution is for fibres with length > %4.1f after scaling)\n",min_len);
	fprintf(fpout,"         unscaled   scaled     limited\n");
	for (i=1; i<=NBINS; i++) {
		fprintf(fpout,"%4.1f %10.4f %10.4f %10.4f\n",i*deltaL,Lbin[i]/nfibres,Lbin_scaled[i]/nfibres,Lbin_limit[i]/nf_limit);
	}
	*/

	// This code is to determine, for each fibre, for each end, the list of connected fibres.
	nvbin = 0;
	for (k=0; k<nfibres; k++) {
		for (int iend = 0; iend<2; iend++) {
			int kvt = fibre[k].kv[iend];
			nlinks[iend] = 0;
			VERTEX v = net->vertex[kvt];
			for (int i=0; i<v.nlinks; i++) {
				if (v.link[i] != k) {
					fibre[k].link[iend][nlinks[iend]] = v.link[i];
					nlinks[iend]++;
				}
			}
			fibre[k].nlinks[iend] = nlinks[iend];
			if (nlinks[iend] > 1) {
				nvbin++;
				NBbin[nlinks[iend]+1]++;
			}
//			printf("fibre: %d iend: %d nlinks: %d\n",k,iend,fibre[k].nlinks[iend]);
//			for (int i=0; i<nlinks[iend]; i++)
//				printf("i: %d link: %d\n",i,fibre[k].link[iend][i]);
		}
	}

			/*
			p1 = net->point[kvt];
			ix = (p1.x - xmin)/DX;
			iy = (p1.y - ymin)/DY;
			iz = (p1.z - zmin)/DZ;
			// Look at all vertexes in this block
			for (int j=0; j<vcounter[ix][iy][iz]; j++) {
				int kv = *(vertexList[ix][iy][iz] + j);
				VERTEX v = net->vertex[kv];
				for (int i=0; i<v.nlinks; i++) {
			
			
		}

		int npts = net->edgeList[k].npts;
		// Look at end 0 first
		p1 = net->point[net->edgeList[k].pt[0]];
		ix = (p1.x - xmin)/DX;
		iy = (p1.y - ymin)/DY;
		iz = (p1.z - zmin)/DZ;
		// Look at all fibre ends in this block
		for (int j=0; j<counter[ix][iy][iz]; j++) {
			int kk = B(j,0,ix,iy,iz);
			if (kk == k) continue;
			int kend = B(j,1,ix,iy,iz);
			if (fibre[kk].kv[kend] == kv[0]) {
				fibre[k].link[0][nlinks[0]] = kk;
				nlinks[0]++;
				if (nlinks[0] > MAX_LINKS) {
					fprintf(fpout,"Error: nlinks[0]: %d\n",nlinks[0]);
					fflush(fpout);
					exit(10);
				}
			}
		}
		// Now look at end 1
		p2 = net->point[net->edgeList[k].pt[npts-1]];
		ix = (p2.x - xmin)/DX;
		iy = (p2.y - ymin)/DY;
		iz = (p2.z - zmin)/DZ;
		// Look at all fibre ends in this block
		for (int j=0; j<counter[ix][iy][iz]; j++) {
			int kk = B(j,0,ix,iy,iz);
			if (kk == k) continue;
			int kend = B(j,1,ix,iy,iz);
			if (fibre[kk].kv[kend] == kv[1]) {
				fibre[k].link[1][nlinks[1]] = kk;
				nlinks[1]++;
				if (nlinks[1] > MAX_LINKS) {
					fprintf(fpout,"Error: nlinks[1]: %d\n",nlinks[1]);
					fflush(fpout);
					exit(10);
				}
			}
		}

		fibre[k].nlinks[0] = nlinks[0];
		fibre[k].nlinks[1] = nlinks[1];
		if (nlinks[0] > 1) {
			nvbin++;
			NBbin[nlinks[0]+1]++;
		}
		if (nlinks[1] > 1) {
			nvbin++;
			NBbin[nlinks[1]+1]++;
		}
	}
	*/

//	fprintf(fpout,"completed connected fibres\n");
//	fflush(fpout);

	if (centre.x==0 || centre.y==0 || centre.z==0) {
		// Get centre
		centre.x = 0;
		centre.y = 0;
		centre.z = 0;
		for (iv=0; iv<net->nv; iv++) {
			centre.x += net->vertex[iv].point.x;
			centre.y += net->vertex[iv].point.y;
			centre.z += net->vertex[iv].point.z;
		}
		centre.x /= net->nv;
		centre.y /= net->nv;
		centre.z /= net->nv;
	}

	// Evaluate angles at vertices, get centre
	for (i=0;i<NANGS;i++)
		Angbin[i] = 0;
	nangbin = 0;
	for (iv=0; iv<net->nv; iv++) {
		if (!net->vertex[iv].used) continue;
		if (net->vertex[iv].nlinks < 3) continue;
		for (int i1=0; i1<net->vertex[iv].nlinks; i1++) {
			int kf1 = net->vertex[iv].link[i1];
			for (int i2=0; i2<net->vertex[iv].nlinks; i2++) {
				if (i1 == i2) continue;
//				int kf2 = net->vertex[iv].fib[i2];	// Note float counting of pairs of fibres
				int kf2 = net->vertex[iv].link[i2];	// Note float counting of pairs of fibres
				double theta = getAngle(iv,kf1,kf2);
				i = (int)(theta/dang);
				if (i > NANGS) exit(1);
				Angbin[i]++;
				nangbin++;
			}
		}
	}
	fprintf(fpout,"\nFibre angle distribution: bin size = %f degrees\n",dang);
	for (i=0; i<NANGS; i++)
		fprintf(fpout,"%3d %6.1f %8.5f\n",i,dang*(i+0.5),Angbin[i]/nangbin);
	printf("made angle distribution\n");
	// Evaluate directional tendencies
	fprintf(fpout,"\nDirectional average dot-products\n");
	for (idir=0; idir<26; idir++) {
		dotsum[idir] = 0;
		for (i=0; i<3; i++)
			v[i] = uvec[idir][i];
		for (int k=0; k<nfibres; k++) {
			double dot = 0;
			for (i=0; i<3; i++) {
				dot += (v[i]*fibre[k].u[i]);
			}
			dotsum[idir] += abs(dot);
		}
		dotsum[idir] /= 2*nfibres;	// because each dot-product gets evaluated twice
		fprintf(fpout,"idir: %2d v: %6.3f %6.3f %6.3f dot ave: %8.4f\n",idir,v[0],v[1],v[2],dotsum[idir]); 
	}
	printf("\nNumber of vertices: %d\n",nvbin/2);
	printf("Length of fibre/vertex: %f\n",Lsum/net->nv);
	printf("Distribution of fibres/vertex:\n");
	fprintf(fpout,"\nNumber of vertices: %d\n",nvbin/2);
	fprintf(fpout,"Length of fibre/vertex: %f\n",Lsum/net->nv);
	fprintf(fpout,"Distribution of fibres/vertex:\n");
	for (i=0; i<MAX_LINKS; i++) {
		printf("%2d %d %8.5f\n",i,NBbin[i],float(NBbin[i])/nvbin);
		fprintf(fpout,"%2d %12.4e\n",i,float(NBbin[i])/nvbin);
	}

	// Now set up nearest vertex for each vertex
	for (int kv1=0; kv1<net->nv; kv1++) {
		VERTEX *v1 = &net->vertex[kv1];
		if (!v1->used) continue;
		POINT p1 = v1->point;
		ix = (int)((p1.x - xmin)/DX);
		iy = (int)((p1.y - ymin)/DY);
		iz = (int)((p1.z - zmin)/DZ);
		// Look at vertexes in [ix][iy][iz]
		double d2min = 1.0e10;
		int kvmin = -1;
		for (int i=0; i<vcounter[ix][iy][iz]; i++) {
			int kv2 = *(vertexList[ix][iy][iz] + i);
			if (kv2 == kv1) continue;
			VERTEX v2 = net->vertex[kv2];
			if (!v2.used) continue;
			bool connected = false;
			for (int ilink=0; ilink < v2.nlinks; ilink++) {
				int ie = v2.link[ilink];
				EDGE edge = net->edgeList[ie];
				if (edge.vert[0] == kv1 || edge.vert[1] == kv1) {
					connected = true;
					continue;
				}
			}
			if (connected) continue;
			POINT p2 = v2.point;
			double d2 = pdist2(p1,p2);
			if (d2 < d2min) {
				d2min = d2;
				kvmin = kv2;
			}
		}
		v1->kvnearest = kvmin;
		if (kvmin < 0)
			v1->dnearest = 0;
		else
			v1->dnearest = (float)sqrt(d2min);
	}
		
#if 0
	for (k=0; k<nfibres; k++) {
		if (k%10000 == 0) printf("fibre: %d\n",k);
		kv[0] = fibre[k].kv[0];
		kv[1] = fibre[k].kv[1];
		// Look at end 0 first
		int kend = 0;
		int k1 = fibre[k].pt[kend];
		p1 = net->point[k1];
		ix = (p1.x - xmin)/DX;
		iy = (p1.y - ymin)/DY;
		iz = (p1.z - zmin)/DZ;
		// end 0 of fibre[k] is in block[0][ix][iy][iz]
		// Look at all fibre ends in this block
		for (int j=0; j<counter[ix][iy][iz]; j++) {
			int kk = B(j,0,ix,iy,iz);
			if (kk == k) continue;
			// eliminate fibres connected to this fibre, this end
			int knear_end = B(j,1,ix,iy,iz);
			int k2 = fibre[kk].pt[1];
			double d = distance(net, k1, k2);
		}
	}
#endif
//	fprintf(fpout,"did MakeFibreList\n");
	return;
}

//-----------------------------------------------------------------------------------------------------
// Determine possible fibre-fibre jumps
//-----------------------------------------------------------------------------------------------------
int SetupJumps(NETWORK *net)
{
	int kv[2];
	for (int k=0; k<nfibres; k++) {
		if (k%10000 == 0) printf("fibre: %d\n",k);
		kv[0] = fibre[k].kv[0];
		kv[1] = fibre[k].kv[1];
		// Look at end 1 first
		int kend = 0;
		int k1 = fibre[k].pt[kend];
		POINT p1 = net->point[k1];
		int ix = (int)((p1.x - xmin)/DX);
		int iy = (int)((p1.y - ymin)/DY);
		int iz = (int)((p1.z - zmin)/DZ);
		// Look at all fibre ends in this block
		for (int j=0; j<counter[ix][iy][iz]; j++) {
			int kk = B(j,0,ix,iy,iz);
			if (kk == k) continue;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Generate a random integer in the range [n1,n2]
//-----------------------------------------------------------------------------------------------------
int random_int(int n1, int n2)
{
//	std::uniform_int_distribution<int> dist( n1, n2 ) ;
	uniform_int_distribution<int> dist( n1, n2 ) ;
	return dist(gen);
}

//-----------------------------------------------------------------------------------------------------
// jdir is the direction of motion: it is toward the end jdir
// Randomly select fibre ifib2 that is connected to end jdir of fibre ifib1
// Replace jdir with the index of the other end of the selected fibre
// Need to weight link probabilities on the basis of change of direction angle:
// more weight for links with less directional change.
//-----------------------------------------------------------------------------------------------------
int next_fibre(int ifib1, int *jdir)
{
//	printf("ifib1: %6d jdir: %d\n",ifib1,*jdir);
	int nlinks = fibre[ifib1].nlinks[*jdir];
	if (nlinks == 0) {
		ndead++;
//		printf("Dead end - reverse direction: %d\n",ndead);
		if (*jdir == 0)
			*jdir = 1;
		else
			*jdir = 0;
		return ifib1;
	}
//	printf("ifib1 links[0]: %d %d\n",fibre[ifib1].link[0][0],fibre[ifib1].link[0][1]);
//	printf("ifib1 links[1]: %d %d\n",fibre[ifib1].link[1][0],fibre[ifib1].link[1][1]);
	int kv1, kv2, ilink, ifib2;
	double theta, psum, p[MAX_LINKS];
	kv1 = fibre[ifib1].kv[*jdir];	// vertex
//	printf("kv: %d nlinks: %d\n",kv,nlinks);
	psum = 0;
	for (ilink=0; ilink<nlinks; ilink++) {
		ifib2 = fibre[ifib1].link[*jdir][ilink];
		theta = getAngle(kv1,ifib1,ifib2);
		theta = (180 - theta)*(PI/180.);	// turning angle
		if (theta > PI/2) {
			p[ilink] = 0.001;
		} else {
			p[ilink] = pow(cos(theta),npow);
		}
//		printf("ilink: %d theta: %f p: %f\n",ilink,theta,p[ilink]);
		psum += p[ilink];
	}
//	std::uniform_int_distribution<double> dist( 0, psum ) ;
//	uniform_real_distribution<double> dist( 0, psum ) ;
	double R = uni_dist(gen);
	R *= psum;
	psum = 0;
	for (ilink=0; ilink<nlinks; ilink++) {
		psum += p[ilink];
		if (psum >= R) {
			break;
		}
	}

//	ilink = random_int(0,nlinks-1);
	ifib2 = fibre[ifib1].link[*jdir][ilink];
//	printf("ilink: %d ifib2: %6d\n",ilink,ifib2);
//	printf("ifib2 links[0]: %d %d\n",fibre[ifib2].link[0][0],fibre[ifib2].link[0][1]);
//	printf("ifib2 links[1]: %d %d\n",fibre[ifib2].link[1][0],fibre[ifib2].link[1][1]);
//	printf("ifib1: %6d nlinks: %d ilink: %d ifib2: %6d\n",ifib1,nlinks,ilink,ifib2);
	if (fibre[ifib1].kv[*jdir] == fibre[ifib2].kv[0])
		*jdir = 1;
	else
		*jdir = 0;
	kv2 = fibre[ifib2].kv[*jdir];	// new vertex
	if (kv1 == kv2) {
		printf("vertex repeat: %d\n",kv1);
		exit(1);
	}
//	printf("new ifib jdir: %6d %d\n",ifib2,*jdir);
	return ifib2;
}

//-----------------------------------------------------------------------------------------------------
// Distance from end iend0 of fibre ifib0 to end iend1 (jdir) of fibre ifib1
// (because direction jdir on a fibre is the same as heading for end jdir)
// ifib0,iend0 is the path start point
//-----------------------------------------------------------------------------------------------------
double get_distance(NETWORK *net, int ifib0, int iend0, int ifib1, int iend1)
{
	int kv0 = fibre[ifib0].kv[iend0];
	int kv1 = fibre[ifib1].kv[iend1];
	POINT p0 = net->vertex[kv0].point;
	POINT p1 = net->vertex[kv1].point;
	float dx = p1.x - p0.x;
	float dy = p1.y - p0.y;
	float dz = p1.z - p0.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}


//-----------------------------------------------------------------------------------------------------
// The cell started on the current fibre (ifib1) at the end iend1.
// We first determine where it was after travelling distance dist
// then compute the straight-line distance of this point from the start vertex (iend0 of ifib0)
// This is returned
// Note: fibre[k] is actually edge edgelist[k]
//
// To accommodate jumping (the continuation of travel from jpt to the end of the fibre):
// For jump, need to start from ipt=jpt, the pt that was jumped to, on ifib1, and use jdir instead of iend0
// In usual case, ipt = 0 (if jdir = 1), or npts-1 (if jdir = 0)
//-----------------------------------------------------------------------------------------------------
double get_partial_distance(NETWORK *net, int ifib0, int iend0, int ifib1, int ipt, int jdir, double dist)
{
	double dd, dsum, dx, dy, dz, x, y, z, fac;
	int npts;
	EDGE *ep;

	ep = &net->edgeList[ifib1];
	npts = ep->npts;
	if (jdir == 1) {	// travel towards end 1
		dsum = 0;
		for (int kseg=ipt; kseg<npts-1; kseg++) {
			POINT p1 = net->point[ep->pt[kseg]];
			POINT p2 = net->point[ep->pt[kseg+1]];
			dx = p2.x - p1.x;
			dy = p2.y - p1.y;
			dz = p2.z - p1.z;
			dd = sqrt(dx*dx+dy*dy+dz*dz);
			if (dsum + dd > dist) {
				fac = (dist-dsum)/dd;
				x = p1.x + fac*dx;
				y = p1.y + fac*dy;
				z = p1.z + fac*dz;
				break;
			} else {
				dsum += dd;
			}
		}
	} else {	// travel towards end 0
		dsum = 0;
		for (int kseg=ipt; kseg>0; kseg--) {
			POINT p1 = net->point[ep->pt[kseg]];
			POINT p2 = net->point[ep->pt[kseg-1]];
			dx = p2.x - p1.x;
			dy = p2.y - p1.y;
			dz = p2.z - p1.z;
			dd = sqrt(dx*dx+dy*dy+dz*dz);
			if (dsum + dd > dist) {
				fac = (dist-dsum)/dd;
				x = p1.x + fac*dx;
				y = p1.y + fac*dy;
				z = p1.z + fac*dz;
				break;
			} else {
				dsum += dd;
			}
		}
	}
	int kv0 = fibre[ifib0].kv[iend0];
	POINT p0 = net->vertex[kv0].point;
	dx = x - p0.x;
	dy = y - p0.y;
	dz = z - p0.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
// The nearest pt on the fibre ifib is found, from the distance to the nearest fibre the jump probability
// is calculated, then a coin is tossed to determine if the jump happens.
// If the answer is yes, jpt is the jumped-to pt on jfib the jumped-to fibre (from pt ipt on ifib)
// (note that ipt and jpt are the indexes of the pts in edgeList[ifib] and edgeList[jfib], i.e. 0,..,npts-1,
// therefore the pt indexes in net->point[] are: edgeList[ifib].pt[ipt] and edgeList[jfib].pt[jpt])
// and the distance along the fibre ifib from the starting end (in the direction jdir) is calculated as jdist.
//-----------------------------------------------------------------------------------------------------
bool check_jump(NETWORK *net, int ifib, int jdir, int *jpt, int *jfib, double *dist)
{
	int i, ij, ip, npts, nearpt;
	POINT p1, p2;
	EDGE edge1, edge2;
	double djump, prob, R, d;

	njumpcalls++;

//	printf("check_jump: ifib: %d jdir: %d\n",ifib,jdir);
	edge1 = net->edgeList[ifib];
	ij = edge1.jumpable_pt;			// index of the jumpable pt on edge
	ip = edge1.pt[ij];
	djump = edge1.dnearest;			// jump distance
	djumpsum += djump;
	if (ij == 0) exit(1);
	if (ij < 0 || djump == 0) {
		printf("ij: %d djump: %f\n",ij,djump);
		return false;
	} else {
		prob = jumpprobfactor/(djump*djump);
	}
	if (djump < 1) printf("check_jump: djump < 1: %f\n",djump);
	R = uni_dist(gen);
	if (R > prob) {
		return false;
	}
	// Jumping
	npts = edge1.npts;
//	printf("Jumping: ifib: %d ij: %d djump: %f\n",ifib,ij,djump);
	// Distance from the start vertex depends on jdir, ij
	if (jdir == 1) {	// travel from end 0 towards end 1
//		printf("jdir: %d\n",jdir);
		d = 0;
		p1 = net->point[edge1.pt[0]];
		for (i=1; i<=ij; i++) {
			p2 = net->point[edge1.pt[i]];
			d += sqrt(pdist2(p1,p2));
			p1 = p2;
		}
	} else {
//		printf("jdir: %d\n",jdir);
		d = 0;
		p1 = net->point[edge1.pt[npts-1]];
		for (i=npts-2; i>=ij; i--) {
			p2 = net->point[edge1.pt[i]];
			d += sqrt(pdist2(p1,p2));
			p1 = p2;
		}
	}
	*dist = d;
//	printf("d: %f\n",d);
	nearpt = net->point[ip].nearestpt;	// this is the index in net->point[]
//	printf("nearpt: %d\n",nearpt);
	*jfib = net->point[nearpt].iedge;	// the pt nearpt on fibre jfib is nearest to the jumpable pt on ifib
//	printf("jfib: %d  nearpt: %d\n",*jfib,nearpt);
	// need the index jpt in the pt list for jfib
	edge2 = net->edgeList[*jfib];
	npts = edge2.npts;
	*jpt = -1;
	for (i=0; i<npts; i++) {
		if (edge2.pt[i] == nearpt) {
			*jpt = i;
			break;
		}
	}
//	printf("ifib: %d ij: %d ip: %d djump: %f prob: %f R: %f\n",ifib,ij,ip,djump,prob,R);
//	fprintf(fpout,"ifib: %d ij: %d jfib: %d jpt: %d\n",ifib,ij,*jfib,*jpt);
//	printf("jfib: %d  nearpt: %d  jpt: %d\n",*jfib,nearpt,*jpt);
	return true;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool check_vjump(NETWORK *net, int ifib1, int jdir1, int *ifib, int *jdir, int *kv, float *djump)
{
	njumpcalls++;
	if (jumpprobfactor == 0) return false;
	int kv1 = net->edgeList[ifib1].vert[jdir1];
	VERTEX v1 = net->vertex[kv1];
	int kv2 = v1.kvnearest;
	if (kv2 < 0) return false;
	double d = v1.dnearest;
	djumpsum += d;
	double prob = jumpprobfactor/(d*d);
	double R = uni_dist(gen);
	if (R > prob) {
		return false;
	}
//	printf("kv1: %d kv2: %d d: %f prob: %f R: %f\n",kv1,kv2,d,prob,R);
	VERTEX v2 = net->vertex[kv2];
	int nlinks = v2.nlinks;
	if (nlinks == 1) {
		*ifib = v2.link[0];
	} else {
		// directional choice - choose the connected edge that is most in the same direction as the jump vj[] to the vertex v2
		double vj[3], vf[3];
		vj[0] = v2.point.x - v1.point.x;
		vj[1] = v2.point.y - v1.point.y;
		vj[2] = v2.point.z - v1.point.z;
		double dotmax = -1.0e10;
		int imax;
		int kve;
		for (int i=0; i<v2.nlinks; i++) {
			int fib = v2.link[i];
			EDGE e = net->edgeList[fib];
			// kve is the vertex at the other end
			if (e.vert[0] == kv2)
				kve = e.vert[1];
			else
				kve = e.vert[0];
			VERTEX ve = net->vertex[kve];
			vf[0] = ve.point.x - v2.point.x;
			vf[1] = ve.point.y - v2.point.y;
			vf[2] = ve.point.z - v2.point.z;
			double d2 = 0;
			for (int j=0; j<3; j++) 
				d2 += vf[j]*vf[j];
			d = sqrt(d2);
			double dot = 0;
			for (int j=0; j<3; j++)
				dot += vj[j]*vf[j]/d;
			if (dot > dotmax) {
				dotmax = dot;
				imax = i;
			}
		}
		*ifib = v2.link[imax];
			
		/*	random choice
		prob = 1.0/nlinks;
		R = uni_dist(gen);
		double psum = 0;
		for (int i=0; i<v2.nlinks; i++) {
			psum += prob;
			if (R < psum) {
				*ifib = v2.link[i];
				continue;
			}
		}
		*/
	}
	if (net->edgeList[*ifib].vert[0] == kv2) {
		*jdir = 1;
	} else {
		*jdir = 0;
	}
	*kv = kv2;

	return true;
}


//-----------------------------------------------------------------------------------------------------
// First set up a list of all fibres with an end inside the starting sphere, then restrict attention
// to these.
//-----------------------------------------------------------------------------------------------------
void traverse(NETWORK *net, int nsteps, double *tpt, double *res_d2sum, double *Cm, int *nsims)
{
	int i, istep, itpt;
	bool reached;
	int *nvisits=NULL;
	double *res_d2tsum;
	double *Lsum=NULL;
	double Ltotal, R, stepsum;
	double t, dt, speed, std_speed;
	bool dbug;

	int nsphere_fibres;
	int sphere_fibre[100000];

	printf("traverse: jumpprobfactor: %f\n",jumpprobfactor);
	fprintf(fpout,"traverse: jumpy: %d v_jumpy: %d\n",jumpy,v_jumpy);
	fprintf(fpout,"traverse: jumpprobfactor: %f\n",jumpprobfactor);
	mean_speed /= shrink_factor;	// the mean speed is scaled to account for the fact that the network description
									// is based on the unscaled measurements - not compensated for shrinkage
									// This is equivalent to scaling fibre lengths by multiplying by shrink_factor
	std_speed = CV*mean_speed;
	res_d2tsum = new double[NDATAPTS];
	for (i=0; i<NDATAPTS; i++) {
		res_d2sum[i] = 0;
		res_d2tsum[i] = 0;
	}

	nvisits = new int[nfibres];
	Lsum = new double[nfibres];
	nsphere_fibres = 0;
	for (i=0; i<nfibres; i++) {
		if (!net->edgeList[i].used) continue;
		//if (i == 0)
		//	Lsum[0] = fibre[0].L_actual;
		//else
		//	Lsum[i] = Lsum[i-1] + fibre[i].L_actual;
		nvisits[i] = 0;
		for (int j=0; j<2; j++) {
			int kp = fibre[i].pt[j];
			float dx = centre.x - net->point[kp].x;
			float dy = centre.y - net->point[kp].y;
			float dz = centre.z - net->point[kp].z;
			if (sqrt(dx*dx+dy*dy+dz*dz) < start_radius) {
				sphere_fibre[nsphere_fibres] = i;
				nsphere_fibres++;
				break;
			}
		}
	}
	for (i=0; i<nsphere_fibres; i++) {
		if (i == 0) 
			Lsum[0] = fibre[sphere_fibre[0]].L_actual;
		else
			Lsum[i] = Lsum[i-1] + fibre[sphere_fibre[i]].L_actual;
	}
	Ltotal = Lsum[nsphere_fibres-1];
	printf("nsphere_fibres: %d\n",nsphere_fibres);

    normal_distribution<double> norm_dist( 0.0, 1.0 ) ;
	
	stepsum = 0;
	njumpcalls = 0;
	djumpsum = 0;
	double d;
	int np = 0;
	for(;;) {
		dbug = (np <= -1);
		R = uni_dist(gen);
		R *= Ltotal;
		// Randomly choose initial direction on the fibre
		int jdir = random_int(0,1);

		// Randomly select a fibre i
		for (i=0; i<nsphere_fibres; i++) {
			if (R <= Lsum[i]) {
				break;
			}
		}
		int ifib = sphere_fibre[i];
		int ifib0 = ifib;

		int iend0;
		if (jdir == 0)
			iend0 = 1;
		else
			iend0 = 0;
		// Starting on fibre ifib0, at the end 0 if jdir=1, at the end 1 if jdir=0
		// i.e. moving towards end jdir of ifib (=ifib0)

		// Check distance of start point from the centre
		int kv0 = fibre[ifib0].kv[iend0];
		POINT p0 = net->vertex[kv0].point;

		float dx = centre.x - p0.x;
		float dy = centre.y - p0.y;
		float dz = centre.z - p0.z;
		d = sqrt(dx*dx+dy*dy+dz*dz);
//		printf("ifib0,iend0,kv0: %d %d %d\n",ifib0,iend0,kv0);
//		printf("centre,p0: %f %f %f   %f %f %f\n",centre.x,centre.y,centre.z,p0.x,p0.y,p0.z);
//		printf("k,d,start_radius: %d %f %f\n\n",k,d,start_radius);
		if (d > start_radius) continue;
		if (save_paths) {
			if (np < npaths) {
				pathcnt[np] = 1;
				path[np][0].x = 0;
				path[np][0].y = 0;
				path[np][0].z = 0;
			}
		}

		if (dbug) fprintf(fpout,"np: %d\n",np);

		speed = mean_speed;
		R = norm_dist(gen);
		speed += R*std_speed;
		bool jump = false;			// jump or not
		int ipt, jpt, jfib;	// when jumpable pt is reached, where jump to jpt on jfib is possible (jpt indices on the edges)
		double jdist;	// distance from the starting end of the current fibre
		reached = false;
		itpt = 0;
		t = 0;
		if (v_jumpy) {
			for (istep=0; istep<nsteps; istep++) {
				if (!net->edgeList[ifib].used) {
					printf("Error: in traverse: unused edge: ifib: %d\n",ifib);
					exit(1);
				}
				dt = fibre[ifib].L_actual/speed;	// time to reach next vertex
				if (t + dt > tpt[itpt]) {
					dt = tpt[itpt] - t;
					// we need to find where the cell was dt after leaving the other end (not jdir) on fibre ifib
					t += dt;
					// ipt is the starting pt index
					if (jdir == 0) {
						ipt = net->edgeList[ifib].npts-1;
					} else {
						ipt = 0;
					}
					d = get_partial_distance(net,ifib0,iend0,ifib,ipt,jdir,dt*speed);	// distance of the time point location from p0
					d *= shrink_factor;
					reached = true;
				} else {
					t += dt;
				}
				if (reached) {	// reached the next time point
					res_d2sum[itpt] += d*d;
					res_d2tsum[itpt] += d*d/t;
					itpt++;
					if (itpt == NDATAPTS) break;
					reached = false;
				}
				// check to see if we want to jump to a nearby vertex
				int ifib1 = ifib;
				int jdir1 = jdir;
				int kv;
				float djump;
				jump = check_vjump(net,ifib1,jdir1,&ifib,&jdir,&kv, &djump);
//				printf("istep: %d t: %f ifib1: %d ifib: %d jump: %d\n",istep,t,ifib1,ifib,jump);
				// returns new fibre (ifib), direction (jdir), vertex (kv), jump distance (djump)
				if (jump) {
					dt = djump/speed;
					t += dt;
					if (t > tpt[itpt]) {	// assume the time point was reached when the next fibre was reached (small error at worst)
						POINT p = net->point[kv];
						//	d = dist from p to start pt p0
						d = sqrt(pdist2(p,p0));
						printf("time point reached: kv: %d d: %f\n",kv,d);
						d *= shrink_factor;
						res_d2sum[itpt] += d*d;
						res_d2tsum[itpt] += d*d/t;
						itpt++;
						if (itpt == NDATAPTS) break;
					}
				} else {
					ifib = ifib1;
					jdir = jdir1;
					ifib = next_fibre(ifib1,&jdir);
				}
				nvisits[ifib]++;
				if (save_paths && np < npaths && istep+1 < NTIMES) {
					int kv1 = fibre[ifib1].kv[jdir1];
					POINT p1 = net->vertex[kv1].point;
					dx = p1.x - p0.x;
					dy = p1.y - p0.y;
					dz = p1.z - p0.z;
					pathcnt[np]++;
					path[np][istep+1].x = dx;
					path[np][istep+1].y = dy;
					path[np][istep+1].z = dz;
				}
			}
		} else {
			for (istep=0; istep<nsteps; istep++) {
				EDGE *edge1 = &net->edgeList[ifib];
				if (dbug) {
					int ip;
					POINT p;
					if (jdir == 1) {	// we are at end 0
						ip = edge1->pt[0];
						p = net->point[ip];
					} else {			// we are at end 1
						ip = edge1->pt[edge1->npts-1];
						p = net->point[ip];
					}
					d = sqrt(pdist2(p,p0));
					fprintf(fpout,"start istep: %d ifib: %d jdir: %d d: %f\n",istep,ifib,jdir,d);
				}
				if (!net->edgeList[ifib].used) {
					printf("Error: in traverse: unused edge: ifib: %d\n",ifib);
					exit(1);
				}
				jump = false;
				int ij = edge1->jumpable_pt;			// index of the jumpable pt on edge
				if (jumpy && (ij >= 0) ) {
					jump = check_jump(net, ifib, jdir, &jpt, &jfib, &jdist);	// jdist = distance travelled on the fibre to the jump pt
				}
				if (jump) {
					dt = jdist/speed;					// time to reach jump pt
					if (dbug) {
						int ip1 = edge1->pt[ij];
						POINT p1 = net->point[ip1];
						EDGE *edge2 = &net->edgeList[jfib];
						int ip2 = edge2->pt[jpt];
						POINT p2 = net->point[ip2];
						d = sqrt(pdist2(p1,p2));
						fprintf(fpout,"jump: from ifib: %d ij: %d --> jfib: %d jpt: %d distance: %f\n",ifib,ij,jfib,jpt,d);
					}
				} else {
					dt = fibre[ifib].L_actual/speed;	// time to reach next vertex
				}
	//			printf("istep,dt: %d %f\n",istep,dt);
				if (t + dt > tpt[itpt]) {
					dt = tpt[itpt] - t;
					// we need to find where the cell was dt after leaving the other end (not jdir) on fibre ifib
					t += dt;
					// ipt is the starting pt index
					if (jdir == 0) {
						ipt = net->edgeList[ifib].npts-1;
					} else {
						ipt = 0;
					}
					d = get_partial_distance(net,ifib0,iend0,ifib,ipt,jdir,dt*speed);	// distance of the time point location from p0
					d *= shrink_factor;
					reached = true;
				} else {
					t += dt;
				}
				if (reached) {	// reached the next time point
					if (dbug) {
						if (jump) 
							fprintf(fpout,"itpt: %d d: %f (before jump pt)\n",itpt,d);
						else
							fprintf(fpout,"itpt: %d d: %f (before vertex)\n",itpt,d);
					}

					res_d2sum[itpt] += d*d;
					res_d2tsum[itpt] += d*d/t;
					itpt++;
					if (itpt == NDATAPTS) break;
					reached = false;
				}
				// either arrived at the vertex, or, if jump, at the jump pt
				// if jump, need to progress to the vertex on the next fibre, otherwise continue as usual
				if (jump) {
//					EDGE *edge1 = &net->edgeList[ifib];
//					int ij = edge1->jumpable_pt;			// index of the jumpable pt on edge ( = ipt)
					double djump = sqrt(edge1->dnearest);	// jump distance
//					printf("jfib: %d jpt: %d djump: %f\n",jfib,jpt,djump);
					// new fibre
					ifib = jfib;
					// note: jpt is the index of pt in jfib list.  
					// check if time point has been reached
					dt = djump/speed;
					t += dt;
					if (t > tpt[itpt]) {	// assume the time point was reached when the next fibre was reached (small error at worst)
						// Let the point[] index be kp
						int kp = net->edgeList[ifib].pt[jpt];
						POINT p = net->point[kp];
						//	d = dist from p to start pt p0
						d = sqrt(pdist2(p,p0));
//						printf("time point reached: kp: %d d: %f\n",kp,d);
						d *= shrink_factor;
						if (dbug) fprintf(fpout,"itpt: %d d: %f (in jump)\n",itpt,d);
						res_d2sum[itpt] += d*d;
						res_d2tsum[itpt] += d*d/t;
						itpt++;
						if (itpt == NDATAPTS) break;
					}
					EDGE *edge2 = &net->edgeList[ifib];
					int npts = edge2->npts;
					// motion continues along fibre ifib (edge2) in direction jdir, to the vertex
					// need to compute remaining distance to the vertex
					// don't allow another jump
					// note: if jpt = 0 or npts-1, we are already at a vertex, must set jdir appropriately
					// to correspond to motion towards that vertex
					if (jpt == 0) {
						jdir = 0;
						if (dbug) fprintf(fpout,"at end 0\n");
						continue;	// really should goto where path is recorded
					} else if (jpt == npts-1) {
						jdir = 1;
						if (dbug) fprintf(fpout,"at end 1\n");
						continue;	// really should goto where path is recorded
					} else {
//						printf("random jdir\n");
						jdir = random_int(0,1);
					}
//					printf("set jdir: %d\n",jdir);
					double dleft = 0;
					POINT p1, p2;
					if (jdir == 1) {	// travel from jpt towards end 1
						p1 = net->point[edge2->pt[jpt]];
						for (i=jpt+1; i<=npts-1; i++) {
							p2 = net->point[edge2->pt[i]];
							dleft += sqrt(pdist2(p1,p2));
							p1 = p2;
						}
					} else {			// travel from jpt towards end 0
						p1 = net->point[edge2->pt[jpt]];
						for (i=jpt-1; i>=0; i--) {
							p2 = net->point[edge2->pt[i]];
							dleft += sqrt(pdist2(p1,p2));
							p1 = p2;
						}
					}
					dt = dleft/speed;
					if (dbug) fprintf(fpout,"jdir: %d jpt: %d npts: %d dleft: %f dt: %f\n",jdir,jpt,npts,dleft,dt);
					if (t + dt > tpt[itpt]) {
//						printf("time point reached: dt: %f t+dt: %f tpt[itpt]: %f\n",dt,t+dt,tpt[itpt]);
						dt = tpt[itpt] - t;
//						printf("get_partial_distance with reduced dt: %f jpt: %d npts: %d\n",dt,jpt,npts);
						// we need to find where the cell was dt after leaving the other end (not jdir) on fibre ifib
						t += dt;
						d = get_partial_distance(net,ifib0,iend0,ifib,jpt,jdir,dt*speed);	// distance of the time point location from p0
//						printf("time point reached: itpt: %d  d: %f\n",itpt,d);
						d *= shrink_factor;
						reached = true;
					} else {
						t += dt;
					}
					if (reached) {	// reached the next time point
						if (dbug) fprintf(fpout,"itpt: %d d: %f (after jump pt)\n",itpt,d);
						res_d2sum[itpt] += d*d;
						res_d2tsum[itpt] += d*d/t;
						itpt++;
						if (itpt == NDATAPTS) break;
						reached = false;
					}
				}
				if (save_paths && np < npaths && istep+1 < NTIMES) {
					int kv1 = fibre[ifib].kv[jdir];
					POINT p1 = net->vertex[kv1].point;
					dx = p1.x - p0.x;
					dy = p1.y - p0.y;
					dz = p1.z - p0.z;
					pathcnt[np]++;
					path[np][istep+1].x = dx;
					path[np][istep+1].y = dy;
					path[np][istep+1].z = dz;
				}
				int ifib1 = ifib;
				int iend1= jdir;
				ifib = next_fibre(ifib,&jdir);

				if (ifib == ifib1 && jdir == iend1) {
					printf("No movement! %d %d\n",ifib,jdir);
					exit(1);
				}
				nvisits[ifib]++;
			}
		}
		np++;
		if (np%1000 == 0) printf("Did path: %d  istep: %d\n",np, istep);
		stepsum += istep;
		if (np == ntrials) break;
	}
	if (istep >= nsteps) {
		printf("Used up all steps!  istep: %d\n",istep);
		exit(1);
	}
	printf("npaths: %d\n",np);
	for (itpt=0; itpt<NDATAPTS; itpt++) {
		Cm[itpt] = res_d2tsum[itpt]/(6*np);		// this is with fixed time - actual t is used, not tmax
	}
	*nsims = np;
	printf("Number of jump calls: %d\n",njumpcalls);
	if (njumpcalls > 0) 
		printf("Average djump: %f\n",djumpsum/njumpcalls);
	printf("Average number of steps: %d\n",(int)(stepsum/np));
}

//-----------------------------------------------------------------------------------------------------
// Count deadends in the network.  
// If deadend_radius = 0 
//    count the whole network
// else
//    count in a sphere of radius deadend_radius centred at (x0,y0,z0)
//-----------------------------------------------------------------------------------------------------
int NumberOfDeadends(NETWORK *net)
{
	float len;
	int k1, k2;
	int nshort = 0;
	int ndeadends = 0;

	if (deadend_radius > 0) {
		fprintf(fpout,"\nCounting deadends in a sphere: radius: %6.1f centre: %6.1f %6.1f %6.1f\n",deadend_radius, centre.x, centre.y, centre.z);
	} else {
		fprintf(fpout,"\nCounting deadends in the whole network\n");
	}
	for (int i=0; i<nfibres; i++) {
		for (int k=0; k<2; k++) {
			if (fibre[i].nlinks[k] == 0) {	// deadend
				if (deadend_radius > 0) {
					int kp = fibre[i].pt[k];
					float dx = centre.x - net->point[kp].x;
					float dy = centre.y - net->point[kp].y;
					float dz = centre.z - net->point[kp].z;
					if (sqrt(dx*dx+dy*dy+dz*dz) < deadend_radius) {
						ndeadends++;
						k1 = fibre[i].pt[0];
						k2 = fibre[i].pt[1];
						len = distance(net, k1, k2);
						if (len <= min_len) nshort++;
					}
				} else {
					ndeadends++;
					k1 = fibre[i].pt[0];
					k2 = fibre[i].pt[1];
					len = distance(net, k1, k2);
					if (len <= min_len) nshort++;
				}
			}
		}
	}
	printf("Number of dead ends: %d\n",ndeadends);
	fprintf(fpout,"Number of dead ends: %d\n",ndeadends);
	printf("Number of short dead ends: %d\n",nshort);
	fprintf(fpout,"Number of short dead ends: %d\n",nshort);
	return ndeadends;
}

//-----------------------------------------------------------------------------------------------------
// Try to find a vertex to connect a dead end to.
//-----------------------------------------------------------------------------------------------------
void Connect(NETWORK *NP1)
{
	int i, j, j0, j1, k, kdead;
	int deadend[2], nd, nfibres_temp, kv0, kv1, kmax, kv, pt0, pt1, ntdead, njoins, npts, nshort, nfree;
	double v[3], w[3], dv, dw, cosa, value, valuemax;
	POINT p0, p1, q;
	double a = 0.5;
	bool joined;
	bool do_join = true;
	bool noneg = false;		// noneg = true if connect link angle > 90 not allowed (cosa < 0)
	bool usecosa = true;	// take the deviation angle into account

	struct join_str
	{
		int kv[2];
		int pt[2];
		bool duplicate;
	};
	typedef join_str JOIN;	
	JOIN *join;

	printf("\nConnect\n");
	double dwmin, dwmax, dwmin_max = 0, dmax = 0;

	// Determine nlinks for each vertex
	for (i=0; i<NP1->nv; i++) {
		NP1->vertex[i].nlinks = 0;
	}
	for (i=0; i<nfibres; i++) {
		kv = NP1->edgeList[i].vert[0];
		NP1->vertex[kv].nlinks++;
		kv = NP1->edgeList[i].vert[1];
		NP1->vertex[kv].nlinks++;
	}
	join = (JOIN *)malloc(nfibres*sizeof(JOIN));	// guaranteed big enough
	ntdead = 0;
	njoins = 0;
	nshort = 0;
	nfree = 0;
	nfibres_temp = nfibres;
	for (i=0; i<nfibres_temp; i++) {
		nd = 0;
		for (int kend=0; kend<2; kend++) {
			if (fibre[i].nlinks[kend] == 0) {	// no connected fibres at this end kend (= 0,1)
				//k1 = fibre[i].pt[0];			// index of end point of one vertex of the fibre
				//k2 = fibre[i].pt[1];			// index of end point of other vertex of the fibre
				//len = distance(NP1, k1, k2);	// distance between vertices
				//if (len <= min_len) {		// Here short deadends were dropped, but this makes a bad .am file
				//	NP1->edgeList[i].used = false;
				//	nshort++;
				//	nd = 0;
				//	break;
				//}
				deadend[nd] = kend;
				nd++;
			}
		}
		if (nd == 0) continue;	// fibre has no deadend
		if (nd == 2) {
			printf("Fibre unconnected at both ends! %d\n",i);
			printf("Error exit\n");
			fprintf(fpout,"Fibre unconnected at both ends! %d\n",i);
			fprintf(fpout,"Error exit\n");
			fclose(fpout);
			exit(1);
		}
		if (!do_join) continue;

		// The fibre has a deadend at end kend - note: we can only get nd = 1
		joined = false;	
		dwmin = 9999;
		for (kdead=0; kdead<nd; kdead++) {		// Probably nd always = 1, kdead = 0, unless there are fibres disconnected at both ends
			ntdead++;
			kv1 = fibre[i].kv[deadend[kdead]];		// vertex index of the deadend vertex
			pt1 = fibre[i].pt[deadend[kdead]];		// this is the pt index of the deadend vertex
			if (deadend[kdead] == 0) {
				kv0 = fibre[i].kv[1];				// vertex index of the other vertex, used to get p0 and then v[]
//				pt0 = fibre[i].pt[1];				// this is the pt index of the other vertex (NOT USED)
			} else {
				kv0 = fibre[i].kv[0];
//				pt0 = fibre[i].pt[0];				// ditto (NOT USED)
			}
			//p0 = NP1->vertex[kv0].point;
			int npts = NP1->edgeList[i].npts;
			if (ntdead%100 == 0)
				printf("fibre: %d npts: %d nlinks: %d %d  kv1: %d  ntdead: %d\n",i,npts,fibre[i].nlinks[0],fibre[i].nlinks[1],kv1,ntdead);
			// To estimate v[], unit vector in the direction of the deadend fibre, now use pts separated by 2 pts, if possible
			if (deadend[kdead] == 0) {
				if (npts >= 4)
					pt0 = NP1->edgeList[i].pt[3];
				else
					pt0 = NP1->edgeList[i].pt[npts-1];
			} else {
				if (npts >= 4)
					pt0 = NP1->edgeList[i].pt[npts-4];
				else
					pt0 = NP1->edgeList[i].pt[0];
			}
			p0 = NP1->point[pt0];
			p1 = NP1->vertex[kv1].point;
			// Checked that this is the same as using pt0, pt1.  No need for kv0
			v[0] = p1.x - p0.x;
			v[1] = p1.y - p0.y;
			v[2] = p1.z - p0.z;
	//		printf("deadend: fibre: %d verticies: %d %d v: %f %f %f\n",i,kv0,kv1,v[0],v[1],v[2]);
			dv = 0;
			for (k=0; k<3; k++) {
				dv += v[k]*v[k];
			}
			dv = sqrt(dv);
			for (k=0; k<3; k++) {
				v[k] /= dv;
			}
	//		printf("normalised v: %f %f %f\n",v[0],v[1],v[2]);
			kmax = -1;
			valuemax = -9999;
			double cosamax;
			for (kv=0; kv<NP1->nv; kv++) {
				if (kv == kv1 || kv == kv0) continue;	// don't try to connect to the parent fibre
				if (NP1->vertex[kv].nlinks == 0) continue;	// this is an unconnected vertex
				q = NP1->vertex[kv].point;
				w[0] = q.x - p1.x;
				w[1] = q.y - p1.y;
				w[2] = q.z - p1.z;
				dw = 0;
				for (k=0; k<3; k++) {
					dw += w[k]*w[k];
				}
				dw = sqrt(dw);	// distance to the vertex
				if (dw < dwmin) dwmin = dw;
				if (dw > max_len) continue;		// no added fibres longer than max_len
				cosa = 0;
				if (usecosa) {
					for (k=0; k<3; k++) {
						cosa += v[k]*w[k];
					}
					if (noneg && cosa < 0) continue;
					cosa /= dw;
					value = (a + cosa)/dw;
				} else {
					value = 1/dw;
				}
				// where theta = angle between v[] and w[], dw = length of w[]
				double valuethreshold = -1000;	// was 0.1;
				if (value > valuemax && value > valuethreshold) {
					kmax = kv;
					valuemax = value;
					dwmax = dw;
					cosamax = cosa;
	//				printf("\ncosa,d,valuemax: %f %f %f\n",cosa,dw,valuemax);
				}
			}
			if (kmax < 0) {
				printf("No suitable vertex for end vertex: %d on fibre: %d\n",kv1,i);
				fprintf(fpout,"No suitable vertex for end vertex: %d on fibre: %d\n",kv1,i);
			} else {
//				fprintf(fpout,"adding fibre to deadend: dw: %4.1f cosa: %5.3f valuemax: %f\n",dwmax,cosamax,valuemax);
				joined = true;
				dmax = MAX(dmax,dwmax);
				join[njoins].kv[0] = kv1;
				join[njoins].kv[1] = kmax;
				join[njoins].pt[0] = pt1;
//				join[njoins].pt[1] = ?;	// what is the point index corresponding to vertex index kmax?
				// need to look at any fibre connected to the vertex.  Choose #0
//				int kf = NP1->vertex[kmax].fib[0];
				int kf = NP1->vertex[kmax].link[0];
				if (fibre[kf].kv[0] == kmax) {
					join[njoins].pt[1] = fibre[kf].pt[0];
				} else {
					join[njoins].pt[1] = fibre[kf].pt[1];
				}

				join[njoins].duplicate = false;
				njoins++;
				if (ntdead%100 == 0) 
					printf("Nearest vertex to: %d is: %d\n",kv1,kmax);
			}
		}

		if (!joined) {
//			NP1->edgeList[i].used = false;	// this is to drop edges that can't be connected // TRY REMOVING
			nfree++;	// what is nfree?  Number of deadends that can't be connected
		}
		if (dwmin > dwmin_max) dwmin_max = dwmin;
	}
//	printf("Total short dead ends: %d\n",nshort);
//	fprintf(fpout,"Total short dead ends: %d\n",nshort);
	printf("Total unconnectable dead ends: %d\n",nfree);
	fprintf(fpout,"Total unconnectable dead ends: %d\n",nfree);
	printf("Total candidate joins: %d\n",njoins);
	fprintf(fpout,"Total candidate joins: %d\n",njoins);
	printf("longest join: dmax: %f\n",dmax);
	fprintf(fpout,"longest join: dmax: %f\n",dmax);
	//printf("dwmin_max: %f\n",dwmin_max);
	//fprintf(fpout,"dwmin_max: %f\n",dwmin_max);
//	fprintf(fpout,"\nCandidate joins:\n");
//	for (i=0; i<njoins; i++) {
//		fprintf(fpout,"%8d %8d %8d\n",i,join[i].kv[0],join[i].kv[1]);
//	}

	// Now clear out duplicate joins
	printf("Clearing out duplicates\n");
	for (j0=0; j0<njoins; j0++) {
		if (join[j0].duplicate) continue;
		kv0 = join[j0].kv[0];
		kv1 = join[j0].kv[1];
		for (j1=0; j1<njoins; j1++) {
			if (join[j1].duplicate) continue;
			if (kv0 == join[j1].kv[1] && kv1 == join[j1].kv[0]) {
//				fprintf(fpout,"Duplicate: %8d %8d %8d %8d\n",j0,j1,kv0,kv1);
				join[j1].duplicate = true;
			}
		}
	}
	printf("Duplicates have been tagged\n");
	//printf("Write out valid joins\n");
	//fprintf(fpout,"\nExtra fibres:\n");
	int nvalid = 0;
	for (j0=0; j0<njoins; j0++) {
		if (join[j0].duplicate) continue;
		nvalid++;
	//	fprintf(fpout,"%8d %8d\n",join[j0].kv[0],join[j0].kv[1]);
	}
	printf("\nNumber of valid extra fibres: %d\n",nvalid);
	fprintf(fpout,"\nNumber of valid extra fibres: %d\n",nvalid);
	fflush(fpout);

	// Now add the extra fibres to NP1
	// The number of edges ne is increased, and the edgelist is supplemented
	// Need to allocate more space for edgelist
	int ne = NP1->ne; //- nfree;	// - nshort
	EDGE *temp_edgeList;
	temp_edgeList = (EDGE *)malloc(ne*sizeof(EDGE));
	j = 0;
	for (i=0; i<NP1->ne; i++) {
		if (!NP1->edgeList[i].used) {
			free(NP1->edgeList[i].pt);
			continue;
		}
		npts = NP1->edgeList[i].npts;
		//if (npts == 0) {
		//	printf("(a) npts = 0: i,j: %d %d\n",i,j);
		//	exit(1);
		//}
		temp_edgeList[j].npts = npts;
		temp_edgeList[j].pt = (int *)malloc(npts*sizeof(int));
		for (k=0; k<npts; k++) {
			temp_edgeList[j].pt[k] = NP1->edgeList[i].pt[k];
		}
		temp_edgeList[j].vert[0] = NP1->edgeList[i].vert[0];
		temp_edgeList[j].vert[1] = NP1->edgeList[i].vert[1];
		temp_edgeList[j].segavediam = NP1->edgeList[i].segavediam;
		temp_edgeList[j].length_um = NP1->edgeList[i].length_um;
		j++;
		free(NP1->edgeList[i].pt);
	}
	free(NP1->edgeList);
	printf("Copied used edges to temp_edgelist\n");
	fprintf(fpout,"Copied used edges to temp_edgelist\n");
	NP1->edgeList = (EDGE *)malloc((ne + nvalid)*sizeof(EDGE));
	printf("Allocated new edgelist: %d\n",ne + nvalid);
	fprintf(fpout,"Allocated new edgelist: %d\n",ne + nvalid);
	for (i=0; i<ne; i++) {
		npts = temp_edgeList[i].npts;
		NP1->edgeList[i].npts = npts;
		NP1->edgeList[i].pt = (int *)malloc(npts*sizeof(int));
		for (k=0; k<npts; k++) {
			NP1->edgeList[i].pt[k] = temp_edgeList[i].pt[k];
		}
		NP1->edgeList[i].vert[0] = temp_edgeList[i].vert[0];
		NP1->edgeList[i].vert[1] = temp_edgeList[i].vert[1];
		NP1->edgeList[i].segavediam = temp_edgeList[i].segavediam;
		NP1->edgeList[i].length_um = temp_edgeList[i].length_um;
		NP1->edgeList[i].used = true;
	}
	free(temp_edgeList);
	printf("Copied temp_edgelist to edgelist\n");
	fprintf(fpout,"Copied temp_edgelist to edgelist\n");

	// Add new edges: edge.vert[0],edge.vert[1] and edge.npts
	printf("Add valid edges\n");
	fprintf(fpout,"Add valid edges\n");
	fflush(fpout);
	k = 0;
	for (j=0; j<njoins; j++) {
		if (join[j].duplicate) continue;
		NP1->edgeList[ne+k].vert[0] = join[j].kv[0];
		NP1->edgeList[ne+k].vert[1] = join[j].kv[1];
		NP1->edgeList[ne+k].npts = 2;
		NP1->edgeList[ne+k].pt = (int *)malloc(2*sizeof(int));
		NP1->edgeList[ne+k].pt[0] = join[j].pt[0];
		NP1->edgeList[ne+k].pt[1] = join[j].pt[1];
//		fprintf(fpout,"edge: %d %d vert: %d %d pt: %d %d\n",k,ne+k,join[j].kv[0],join[j].kv[1],join[j].pt[0],join[j].pt[1]);
//		fflush(fpout);
		k++;
	}
	NP1->ne = ne + nvalid;
	free(join);
	printf("Added valid joins to edgelist: %d\n",nvalid);
	fprintf(fpout,"Added valid joins to edgelist: %d\n",nvalid);
	return;

// Check deadends
	printf("Checking for deadends\n");
	fprintf(fpout,"Checking for deadends\n");
	int *nvertex_links = (int *)malloc(NP1->nv*sizeof(int));
	for (i=0; i<NP1->ne; i++) {
		for (int k = 0; k<2; k++) {
			int kv = NP1->edgeList[i].vert[k];
			nvertex_links[kv]++;
		}
	}
	int ndeadends = 0;
	for (i=0; i<NP1->nv; i++) {
		if (nvertex_links[i] == 1) {
			ndeadends++;
			printf("vertex %d is connected to only one edge\n",i);
			fprintf(fpout,"vertex %d is connected to only one edge\n",i);
		}
	}
	free(nvertex_links);
	printf("ndeadends: %d\n",ndeadends);
	fprintf(fpout,"ndeadends: %d\n",ndeadends);
}

//-----------------------------------------------------------------------------------------------------
// Check how vertex and point are related
//-----------------------------------------------------------------------------------------------------
void check_vertex(NETWORK *net)
{
	printf("check_vertex\n");
	for (int kv=0; kv<100; kv++) {
		POINT pv = net->vertex[kv].point;
		printf("vertex kv: %d at: %f %f %f\n",kv,pv.x,pv.y,pv.z);
//		int ipt = net->vertex[kv].pt;
//		POINT p = net->point[ipt];
//		printf("pt: %d %f %f %f\n",ipt,p.x,p.y,p.z);
		int nlinks = net->vertex[kv].nlinks;
//		int nf = net->vertex[kv].nf;
		printf("nlinks: %d\n",nlinks);
		for (int i=0; i<nlinks; i++) {
//			int ifib = net->vertex[kv].fib[i];
			int ilink = net->vertex[kv].link[i];
			FIBRE fib = fibre[ilink];
			printf("fibre: %d  kv[2]: %d %d pt[2]: %d %d\n",ilink,fib.kv[0],fib.kv[1],fib.pt[0],fib.pt[1]);
		}
	}
	exit(1);
}

//-----------------------------------------------------------------------------------------------------
// Check jump distances
//-----------------------------------------------------------------------------------------------------
void check_dist(NETWORK *net)
{
	int i, j, ie, ip, iv;
	float dd = 1;
	float nbin[20];
	int n = 20;
	printf("min edge distances\n");
	for (j=0; j<n; j++)
		nbin[j] = 0;
	int ntot = 0;
	float dsum = 0;
	for (ie=0; ie<net->ne; ie++) {
		float d = net->edgeList[ie].dnearest;
		if (d == 0) continue;
//		i = net->edgeList[ie].jumpable_pt;
//		ip = net->edgeList[ie].pt[i];
//		float dj = sqrt(net->point[ip].d2near);
//		if (dj != d) printf("ie: %d i: %d ip: %d d: %f dj: %f\n",ie,i,ip,d,dj);
		ntot++;
		dsum += d;
		j = (int)(d/dd);
		if (j > n-1) j = n-1;
		nbin[j]++;
	}
	for (j=0; j<n; j++) {
		printf("j: %d fraction: %f\n",j,nbin[j]/ntot);
	}
	printf("average d: %f\n",dsum/ntot);
	printf("min point distances\n");
	for (j=0; j<n; j++)
		nbin[j] = 0;
	ntot = 0;
	for (ip=0; ip<net->np; ip++) {
		float d = net->point[ip].d2near;
		if (d == 0) continue;
		ntot++;
		d = sqrt(d);
		j = (int)(d/dd);
		if (j > n-1) j = n-1;
		nbin[j]++;
	}
	for (j=0; j<n; j++) {
		printf("j: %d fraction: %f\n",j,nbin[j]/ntot);
	}
	printf("min vertex distances\n");
	dsum = 0;
	for (j=0; j<n; j++)
		nbin[j] = 0;
	ntot = 0;
	for (iv=0; iv<net->nv; iv++) {
		float d = net->vertex[iv].dnearest;
		if (d <= 0) continue;
		ntot++;
		dsum += d;
		j = (int)(d/dd);
		if (j > n-1) j = n-1;
		nbin[j]++;
	}
	for (j=0; j<n; j++) {
		printf("j: %d fraction: %f\n",j,nbin[j]/ntot);
	}
	printf("average d: %f\n",dsum/ntot);

	fprintf(fpout,"\nVertexes with 2 links:\n");
	n=0;
	for (int iv=0; iv<net->nv; iv++) {
		VERTEX v = net->vertex[iv];
		if (v.nlinks == 2) {
			n++;
//			fprintf(fpout,"iv: %d\n",iv);
		}
	}
	fprintf(fpout,"n: %d nv: %d\n",n,net->nv);

	// check fibres 16,17
	for (int ifib=16; ifib<=17; ifib++) {
		fprintf(fpout,"\n");
		EDGE e = net->edgeList[ifib];
		int npts = e.npts;
		fprintf(fpout,"ifib: %d npts: %d ij: %d\n",ifib,npts,e.jumpable_pt);
		for (i=0; i<npts; i++)
			fprintf(fpout,"i: %d pt: %d\n",i,e.pt[i]);

		VERTEX v = net->vertex[e.vert[0]];
		int nlinks = v.nlinks;
		fprintf(fpout,"vert[0]: %d nlinks: %d\n",e.vert[0],nlinks);
		for (i=0; i<nlinks; i++)
			fprintf(fpout,"i: %d link: %d\n",i,v.link[i]);
		v = net->vertex[e.vert[1]];
		nlinks = v.nlinks;
		fprintf(fpout,"vert[1]: %d nlinks: %d\n",e.vert[1],nlinks);
		for (i=0; i<nlinks; i++)
			fprintf(fpout,"i: %d link: %d\n",i,v.link[i]);
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool inlist(int ie, int n, int *list)
{
	for (int i=0; i<n; i++) {
		if (list[i] == ie) 
			return true;
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void Fragments(NETWORK *net)
{
	int list[10000], badlist[10000];
	int n, ie0, ie, iv, nlinks, nt, i, j, nadded;
	EDGE *e;
	VERTEX *v;
	bool dbug;

	printf("Fragments: ne: %d\n",net->ne);
	nt = 10;
	int nbad = 0;
	for (ie0=0; ie0<net->ne; ie0++) {
//	for (ie0=0; ie0<1; ie0++) {
		dbug = false; //(ie0 == -34282);
		if (dbug) printf("ie0: %d\n",ie0);
		e = &net->edgeList[ie0];
		if (!e->used) continue;
		n = 0;
		list[n] = ie0;
		n++;

		for (int k=0; k<30; k++) {
			if (dbug) printf("k: %d n: %d\n",k,n);
			int nold = n;
			if (n > 200) break;
			nadded = 0;
			for (i=0; i<nold; i++) {
				ie = list[i];
				e = &net->edgeList[ie];
				if (dbug) printf("    i: %d ie: %d\n",i,ie);
				for (int iend=0; iend<2; iend++) {
					iv = e->vert[iend];
					v = &net->vertex[iv];
					nlinks = v->nlinks;
					if (dbug) printf("    iend: %d iv: %d\n",iend,iv);
					for (j=0; j<nlinks; j++) {
						int ie1 = v->link[j];
						if (dbug) printf("      j: %d ie1: %d\n",j,ie1);
						if (!inlist(ie1,n,list)) {
							list[n] = ie1;
							n++;
							if (dbug) printf("add: %d\n",ie1);
							nadded++;
						}
					}
				}
			}
			if (dbug) printf("        nadded: %d\n",nadded);
			if (nadded == 0) break;
		}
		if (n < 190) {
			if (nadded != 0) printf("nadded != 0: nadded: %d\n",nadded);
			badlist[nbad] = ie0;
			nbad++;
			if (n > 10) printf("ie0: %d n: %d\n",ie0,n);
			if (dbug) {
			for (i=0; i<n; i++)
				printf("ie0: %d i: %d list[i]: %d\n",ie0,i,list[i]);
			}
		}
	}
	printf("nbad: %d ne: %d\n",nbad,net->ne);
	fprintf(fpout,"\nFragments: nbad: %d out of ne: %d\n",nbad,net->ne);
	// Now mark bad edges and all points and vertexes as unused
	for (i=0; i<nbad; i++) {
		ie = badlist[i];
		e = &net->edgeList[ie];
		e->used = false;
		for (j=0; j<2; j++) {
			iv = e->vert[j];
			v = &net->vertex[iv];
			v->used = false;
		}
		fprintf(fpout,"dropped edge: %d vertexes: %d %d\n",ie,e->vert[0],e->vert[1]);
		for (j=0; j<e->npts; j++) {
			int ip = e->pt[j];
			net->point[ip].used = false;
		}
	}

	// Redo association of edges and vertexes
	int kv;
	for (kv=0; kv<net->nv; kv++)
		net->vertex[kv].nlinks = 0;
	for (int k=0; k<net->ne; k++) {
		if (!net->edgeList[k].used) continue;
		for (int iv=0; iv<2; iv++) {							// edge ends
			kv = net->edgeList[k].vert[iv];						// end vertexes
			net->vertex[kv].link[net->vertex[kv].nlinks] = k;	// end vertex is connected to edge k
			net->vertex[kv].nlinks++;
		}
	}

}

//-----------------------------------------------------------------------------------------------------
// OLD PURPOSE
// This code assumes that the selection of a region of interest (a cube) is carried out on the 3D image.
// This means that the cube centre (xc,yc,zc) and the width (diameter) are all specified in voxel
// coordinates.  The values are converted to um by multiplying by voxelsize in um.  This is necessary
// in order to specify the corresponding region in the Amira network (in which the distance unit is um).
// To enable comparison of the zoomed network file with the cropped 3D image, either the network must
// be scaled to use voxelsize as the distance unit (in which case direct comparison is possible in Amira),
// or the network file coordinates can be left in units of um (in which case the voxelsize must be specified
// when the image file is imported into Amira).  The second option has been adopted.
// NOW
// (1) Read the network from the .am file
// (2) Determine length and direction (unit vector u(k)) of each segment S(k), k=1,NS
// (3) For a regular grid of angles: P(i,j) = (i*da, j*da), i=1,..,N, j=1,..,N, da = pi/N
//     generate unit vectors v(i,j) then compute Sum(over k) {u(k).v(i,j)}
//
// Q:
// What is a segment?
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *input_amfile, *output_file;
	char output_amfile[256];
	char drive[32], dir[128],filename[256], ext[32];
//	char errfilename[256];
//	char output_basename[256];
	int mode, do_connect, ijumpy;
	float origin_shift[3];
	int limit_mode;
	float limit_value;
	int ivert;

	int itpt;
	double *tpt;
	double d2ave, Cm1;
	double res_d2sum[NDATAPTS];
	double Cm[NDATAPTS];
	int cmgui_flag=1;
	NETWORK *NP0;

	printf("Three functions: (1) To strip unconnected edges out of the network\n");
	printf("                 (2) To join dead ends in the network (if possible), \n                   and trim unconnected dead ends\n");
	printf("                 (3) To estimate the coefficient of motility Cm\n");
	printf("\n");
	printf("These three functions are carried out by separate program executions.  \nThe network should be cleaned up and healed (joined, trimmed) before Cm estimation.\n");
	printf("For the healing pass, set max_len \n(unscaled upper limit on length of an added connection).\n");
	printf("The shrinkage compensation factor must be specified\n(Note that the dimensions in the am file are left unscaled)\n");
	printf("\n");

	printf("argc: %d\n", argc);
	if (argc != 4 && argc != 9 && argc != 20) {
		printf("To strip out unconnected edges:\n");
		printf("Usage: conduit_analyse input_amfile output_file ivert\n");
		printf("       ivert       = flag for am file with vertices only\n");
		printf("       (the cleaned network file will be created as: \n");
		printf("            clean.am if ivert=0, clean_vert.am if ivert=1)\n");
		printf("\nTo perform joining and trimming of dead ends:\n");
		printf("Usage: conduit_analyse input_amfile output_file sfactor max_len limit_mode limit_value ddiam dlen\n");
		printf("       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		printf("       max_len     = upper limit on (unscaled) connecting fibre length (um)\n");
		printf("       limit_mode  = restriction for computing fibre statistics:\n");
		printf("                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		printf("       limit_value = value of limit for the chosen mode\n");
		printf("       ddiam       = bin size for diameter distribution\n");
		printf("       dlen        = bin size for length distribution\n");
		printf("\nTo simulate cell paths and estimate Cm:\n");
		printf("Usage: conduit_analyse input_amfile output_file sfactor npow ntrials x0 y0 z0 radius deadend_radius speed CV npaths ijumpy jpfactor limit_mode limit_value ddiam dlen\n");
		printf("       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		printf("       npow        = power of cos(theta) in weighting function for probability\n");
		printf("                     of taking a branch (theta = turning angle)\n");
		printf("       ntrials     = number of cell paths simulated\n");
		printf("       x0,y0,z0    = centre of sphere within which the cell paths start (um)\n");
		printf("       radius      = radius of sphere within which the cell paths start (um)\n");
		printf("       deadend_radius = radius of sphere within which deadends are counted (um)\n");
		printf("       speed       = mean cell speed (um/min)\n");
		printf("       CV          = coefficient of variation of speed (std dev/mean)\n");
		printf("       npaths      = number of cell paths to save\n");
		printf("       ijumpy      = controls jumping:\n");
		printf("                     0=no jumps, 1=jump along edge, 2=jump at vertex\n");
		printf("       jpfactor    = jump probability factor\n");
		printf("       limit_mode  = restriction for computing fibre statistics:\n");
		printf("                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		printf("       limit_value = value of limit for the chosen mode\n");
		printf("       ddiam       = bin size for diameter distribution\n");
		printf("       dlen        = bin size for length distribution\n");
		printf("\n");
		printf("Writing conduit_analyse_error.log\n");
		fpout = fopen("conduit_analyse_error.log","w");
		fprintf(fpout,"To strip out unconnected edges:\n");
		fprintf(fpout,"Usage: conduit_analyse input_amfile output_file ivert\n");
		fprintf(fpout,"       ivert       = flag for am file with vertices only\n");
		fprintf(fpout,"       (the cleaned network file will be created as: \n");
		fprintf(fpout,"            clean.am if ivert=0, clean_vert.am if ivert=1)\n");
		fprintf(fpout,"\nTo perform joining and trimming of dead ends:\n");
		fprintf(fpout,"Usage: conduit_analyse input_amfile output_file sfactor max_len limit_mode limit_value ddiam dlen\n");
		fprintf(fpout,"       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		fprintf(fpout,"       max_len     = upper limit on (unscaled) connecting fibre length (um)\n");
		fprintf(fpout,"       limit_mode  = restriction for computing fibre statistics:\n");
		fprintf(fpout,"                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		fprintf(fpout,"       limit_value = value of limit for the chosen mode\n");
		fprintf(fpout,"       ddiam       = bin size for diameter distribution\n");
		fprintf(fpout,"       dlen        = bin size for length distribution\n");
		fprintf(fpout,"\nTo simulate cell paths and estimate Cm:\n");
		fprintf(fpout,"Usage: conduit_analyse input_amfile output_file sfactor npow ntrials x0 y0 z0 radius speed CV npaths ijumpy jpfactor limit_mode limit_value ddiam dlen\n");
		fprintf(fpout,"       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		fprintf(fpout,"       npow        = power of cos(theta) in weighting function for probability\n");
		fprintf(fpout,"                     of taking a branch (theta = turning angle)\n");
		fprintf(fpout,"       ntrials     = number of cell paths simulated\n");
		fprintf(fpout,"       x0,y0,z0    = centre of sphere within which the cell paths start (um)\n");
		fprintf(fpout,"       radius      = radius of sphere within which the cell paths start (um)\n");
		fprintf(fpout,"       speed       = mean cell speed (um/min)\n",CV);
		fprintf(fpout,"       CV          = coefficient of variation of speed (std dev/mean)\n");
		fprintf(fpout,"       npaths      = number of cell paths to save\n");
		fprintf(fpout,"       ijumpy      = controls jumping:\n");
		fprintf(fpout,"                     0=no jumps, 1=jump along edge, 2=jump at vertex\n");
		fprintf(fpout,"       jpfactor    = jump probability factor\n");
		fprintf(fpout,"       limit_mode  = restriction for computing fibre statistics:\n");
		fprintf(fpout,"                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		fprintf(fpout,"       limit_value = value of limit for the chosen mode\n");
		fprintf(fpout,"       ddiam       = bin size for diameter distribution\n");
		fprintf(fpout,"       dlen        = bin size for length distribution\n");
		fprintf(fpout,"\n");
		fprintf(fpout,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fpout,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fpout);
		return 1;	// Wrong command line
	} 
	input_amfile = argv[1];
	output_file = argv[2];
	fpout = fopen(output_file,"w");	

	_splitpath(output_file,drive,dir,filename,ext);
	fprintf(fpout,"output_file: %s\n",output_file);
//	fprintf(fpout,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);

	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
//	NP1 = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;
	printf("Read Amira file\n");

	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;

	if (argc == 4) {
		// testing
//		WriteAmiraFile_vertices("vert.am",input_amfile,NP0,origin_shift);
//		return 0;
//		vertices_only = true;
		sscanf(argv[3],"%d",&ivert);
		vertices_only = (ivert != 0);
		mode = 1;
		SetupVertexLists(NP0);
		// To look for unconnected fragments and create a network with all edges used.
		Fragments(NP0);
		if (vertices_only)
			strcpy(filename,"clean_vert.am");
		else
			strcpy(filename,"clean.am");
		sprintf(output_amfile,"%s%s%s",drive,dir,filename);
		printf("output_amfile: %s\n",output_amfile);
		fprintf(fpout,"output_amfile: %s\n",output_amfile);
		WriteAmiraFile_clean(output_amfile,input_amfile,NP0,origin_shift);
		return 0;
	} else if (argc == 9) {
		mode = 2;
		sscanf(argv[3],"%lf",&shrink_factor);
		sscanf(argv[4],"%lf",&max_len);
		sscanf(argv[5],"%d",&limit_mode);
		sscanf(argv[6],"%f",&limit_value);
		sscanf(argv[7],"%f",&ddiam);
		sscanf(argv[8],"%f",&dlen);
		do_connect = true;
		deadend_radius = 0;
		//char add_str[13];
		//if (max_len < 10)
		//	sprintf(a/dd_str,"_jt_max%3.1f.am",max_len);
		//else
		//	sprintf(add_str,"_jt_max%4.1f.am",max_len);
		//strcpy(output_amfile,input_amfile);
		//int len;
		//for(len=strlen(output_amfile);len >= 0 && output_amfile[len] != '.';len--);
		//output_amfile[len] = '\0';
		//strcat(output_amfile,add_str);
		//printf("\noutput_amfile: %s\n",output_amfile);
		//fprintf(fpout,"\noutput_amfile: %s\n",output_amfile);
		//return 0;

	} else if (argc == 20) {
		mode = 3;
		sscanf(argv[3],"%lf",&shrink_factor);
		sscanf(argv[4],"%d",&npow);
		sscanf(argv[5],"%d",&ntrials);
		sscanf(argv[6],"%f",&centre.x);
		sscanf(argv[7],"%f",&centre.y);
		sscanf(argv[8],"%f",&centre.z);
		sscanf(argv[9],"%lf",&start_radius);
		sscanf(argv[10],"%lf",&deadend_radius);
		sscanf(argv[11],"%lf",&mean_speed);
		sscanf(argv[12],"%lf",&CV);
		sscanf(argv[13],"%d",&npaths);
		sscanf(argv[14],"%d",&ijumpy);
		sscanf(argv[15],"%f",&jumpprobfactor);
		sscanf(argv[16],"%d",&limit_mode);
		sscanf(argv[17],"%f",&limit_value);
		sscanf(argv[18],"%f",&ddiam);
		sscanf(argv[19],"%f",&dlen);
		if (ijumpy == 0) {
			jumpy = false;
			v_jumpy = false;
		} else if (ijumpy == 1) {
			jumpy = true;
			v_jumpy = false;
			fprintf(fpout,"\nedge jumping\n");
		} else if (ijumpy == 2) {
			jumpy = false;
			v_jumpy = true;
			fprintf(fpout,"\nvertex jumping\n");
		} else {
			printf("Error: ijumpy: %d must be 0,1,2\n",ijumpy);
			exit(1);
		}
		do_connect = false;
		save_paths = (npaths > 0);
		if (save_paths) {
			// Need to allocate space for paths
			pathcnt = (int *)malloc(npaths*sizeof(int));
			path = (POINT **)malloc(npaths*sizeof(POINT *));
			for (int i=0; i<npaths; i++) {
				path[i] = (POINT *)malloc(NTIMES*sizeof(POINT));
			}
		}
	}

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

	if (ddiam > 0) {
		err = CreateDistributions(NP0);
		if (err != 0) return 4;
	}
	SetupVertexLists(NP0);
	MakeFibreList(NP0);
//	check_vertex(NP0); 
	ndead = NumberOfDeadends(NP0);
	if (mode == 2 || (mode == 3 && deadend_radius == 0)) {
		fprintf(fpout,"Number of dead ends in the whole network: %d\n",ndead);
	} else {
		fprintf(fpout,"Number of dead ends in a sphere of radius: %f  %d\n",deadend_radius,ndead);
	}
	if (mode ==3 && jumpy)
		SetupPointLists(NP0);
//	check_dist(NP0);

	if (do_connect) {
		printf("do_connect\n");
		fprintf(fpout,"do_connect\n");
		char add_str[13];
		if (max_len < 10)
			sprintf(add_str,"_jt_max%3.1f.am",max_len);
		else
			sprintf(add_str,"_jt_max%4.1f.am",max_len);
		strcpy(output_amfile,input_amfile);
		int len;
		for(len=(int)strlen(output_amfile);len >= 0 && output_amfile[len] != '.';len--);
		filename[len] = '\0';
		strcat(output_amfile,add_str);
		Connect(NP0);
		err = WriteAmiraFile(output_amfile,input_amfile,NP0,origin_shift);
		if (err != 0) return 3;

		fprintf(fpout,"\nRecompute distributions with connected network\n");
		for (int i=0;i<NP0->ne;i++) {
			EDGE edge = NP0->edgeList[i];
			double dave = 0;
			for (int k=0;k<edge.npts;k++) {
				int j = edge.pt[k];
				dave += NP0->point[j].d;
				NP0->edgeList[i].segavediam = (float)(dave/edge.npts);
			}
		}
		if (ddiam > 0) {
			err = CreateDistributions(NP0);
			if (err != 0) return 4;
		}
		return 0;
	}

	// tpt[] holds the equi-spaced time points for the saved paths
	tpt = new double[NDATAPTS];
	for (itpt=0; itpt<NDATAPTS; itpt++) {
		tpt[itpt] = (itpt+1)*deltat;
	}

	fprintf(fpout,"\nSimulating cell paths\n");
	fprintf(fpout,"ntrials: %d deltat: %6.2f\n",ntrials,deltat);
	fprintf(fpout,"Starting sphere centre: %6.1f %6.1f %6.1f\n",centre.x,centre.y,centre.z);
	fprintf(fpout,"Starting sphere radius: %6.1f\n",start_radius);
	fprintf(fpout,"npow: %d\n",npow);
	fprintf(fpout,"mean speed: %6.1f CV: %6.2f\n",mean_speed,CV); 

	int nruns = 1;
	for (int irun=0; irun<nruns; irun++) {
		fprintf(fpout,"\n");
	ndead = 0;
	int nsims;		// this is the actual number of paths simulated, possibly less than the requested number: ntrials
	traverse(NP0,1000,tpt,res_d2sum,Cm,&nsims);
	printf("# of deadends encountered: %d\n",ndead);
	printf("\nCm estimation results:\n");
	fprintf(fpout,"# of paths simulated: %d\n",nsims);
	fprintf(fpout,"# of deadends encountered: %d\n",ndead);
	fprintf(fpout,"\nCm estimation results:\n");
	for (itpt=0; itpt<NDATAPTS; itpt++) {
//		d2ave = res_d2sum[itpt]/ntrials;
		d2ave = res_d2sum[itpt]/nsims;
		Cm1 = d2ave/(6*tpt[itpt]);		// this is with fixed time - actual t is used, not tmax
		printf("t: %6.1f d2: %9.1f Cm: %6.2f\n",tpt[itpt],d2ave,Cm1);
		fprintf(fpout,"t: %6.1f d2: %9.1f Cm: %6.2f\n",tpt[itpt],d2ave,Cm1);		// dropped Cm[itpt] - no difference now
	}
	}
	if (save_paths) {
//		fpout = fopen("path.dat","w");
//		printf("opened path.dat\n");
		fprintf(fpout,"\nCell paths\n");
		fprintf(fpout,"%d\n",npaths);
		for (int k=0; k<npaths; k++) {
			fprintf(fpout,"%d\n",pathcnt[k]);
			for (int i=0; i<pathcnt[k]; i++) {
				fprintf(fpout,"%8.1f %8.1f %8.1f\n",path[k][i].x,path[k][i].y,path[k][i].z);
			}
		}
	}

	fclose(fpout);
	return 0;
	/*
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,NP0,origin_shift);
		if (err != 0) return 5;
	}
	return 0;
	*/
}
