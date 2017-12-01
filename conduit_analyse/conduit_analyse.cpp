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
#define MAX_LINKS 8
#define NB 10
#define NX NB
#define NY NB
#define NZ NB

struct fibre_str
{
	double L_actual, L_direct;
	double u[3];
	int kv[2];	// vertex indicies
	int pt[2];	// end point indicies
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

//#define NPATHS 200
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

bool use_len_limit, use_len_diam_limit;
float len_limit, len_diam_limit;
float ddiam, dlen;
#define NBOX 100

#define N1 MAXBLOCK
#define N2 2*MAXBLOCK
#define N3 2*MAXBLOCK*NX
#define N4 2*MAXBLOCK*NX*NY
#define B(j,kfe,ix,iy,iz) blocks[(iz)*N4 + (iy)*N3 + (ix)*N2 + (kfe)*N1 + (j)]	

int jumpy = 0;

using namespace std;

//std::seed_seq seq;
//std::mt19937 gen(seq);
seed_seq seq;
mt19937 gen(seq);


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
//		printf("ie,len,dave: %d %f %f\n",ie,len,dave);
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
	}
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
	printf("Allocated arrays: np: %d nv: %d ne: %d\n",net->np,net->nv,net->ne);

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
					for (int kk=0;kk<edge.npts;kk++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
						if (kk > 0 && kk<edge.npts-1) {
							sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[kk] = kp;
							net->edgeList[i].pt_used[kk] = kp;
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
// A smoothed network NP1 is generated from NP0
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
	float rad;
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
		ix = (p1.x - xmin)/DX;
		iy = (p1.y - ymin)/DY;
		iz = (p1.z - zmin)/DZ;
		B(counter[ix][iy][iz],0,ix,iy,iz) = k;
		B(counter[ix][iy][iz],1,ix,iy,iz) = 0;	// end 1
		counter[ix][iy][iz]++;
		p2 = net->point[net->edgeList[k].pt[npts-1]];
		ix = (p2.x - xmin)/DX;
		iy = (p2.y - ymin)/DY;
		iz = (p2.z - zmin)/DZ;
		B(counter[ix][iy][iz],0,ix,iy,iz) = k;
		B(counter[ix][iy][iz],1,ix,iy,iz) = 1;	// end 2
		counter[ix][iy][iz]++;
	}
	int maxcount = 0;
	for (ix=0;ix<NX;ix++) {
		for (iy=0;iy<NY;iy++) {
			for (iz=0;iz<NZ;iz++) {
				maxcount = MAX(maxcount,counter[ix][iy][iz]);
			}
		}
	}
	printf("Max block count: %d\n",maxcount);
	//for (ix=0;ix<NX;ix++) {
	//	for (iy=0;iy<NY;iy++) {
	//		for (iz=0;iz<NZ;iz++) {
	//			printf("Block: %d %d %d  count: %d\n",ix,iy,iz,counter[ix][iy][iz]);
	//		}
	//	}
	//}
	if (maxcount > MAXBLOCK) {
		printf("Error: MAXBLOCK exceeded\n");
		printf("MAXBLOCK: %d\n",MAXBLOCK);
		exit(1);
	}
}

//-----------------------------------------------------------------------------------------------------
// Apply shrink_factor here only to the values that are output, and the length distribution.
// All node locations and fibre dimensions are left unscaled.
//-----------------------------------------------------------------------------------------------------
void makeFibreList(NETWORK *net)
{
#define NBINS 50
#define NANGS 45
	int i, k, iv, idir, i_scaled, i_limit, nf_limit;
	int ix, iy, iz;
	int kv[2], nlinks[2];
	POINT p1, p2;
	float L_actual_sum=0, L_direct_sum=0;
	float deltaL = 1.0;
	float v[3], dotsum[26];
	float Lmin, Lmax, Lsum, Lbin[NBINS+1], NBbin[10], Lbin_scaled[NBINS+1], Lbin_limit[NBINS+1];
	float Angbin[NANGS], dang;
	int nvbin, nangbin;

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
	for (i=0; i<10; i++)
		NBbin[i] = 0;
	
	fibre = (FIBRE *)malloc(nfibres*sizeof(FIBRE));
	xmin = ymin = zmin = 1.0e10;
	xmax = ymax = zmax = 0;
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
		i = fibre[k].L_actual/deltaL + 0.5;
		i = MIN(i,NBINS);
		Lbin[i] += 1;
		i_scaled = shrink_factor*fibre[k].L_actual/deltaL + 0.5;
		i_scaled = MIN(i_scaled,NBINS);
		Lbin_scaled[i_scaled] += 1;
		if (shrink_factor*fibre[k].L_actual > min_len) {
			nf_limit++;
			i_limit = shrink_factor*fibre[k].L_actual/deltaL + 0.5;
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
		xmin = MIN(xmin,p1.x);
		xmin = MIN(xmin,p2.x);
		xmax = MAX(xmax,p1.x);
		xmax = MAX(xmax,p2.x);
		ymin = MIN(ymin,p1.y);
		ymin = MIN(ymin,p2.y);
		ymax = MAX(ymax,p1.y);
		ymax = MAX(ymax,p2.y);
		zmin = MIN(zmin,p1.z);
		zmin = MIN(zmin,p2.z);
		zmax = MAX(zmax,p1.z);
		zmax = MAX(zmax,p2.z);
	}
	DX = (xmax-xmin)/NX+1;
	DY = (ymax-ymin)/NY+1;
	DZ = (zmax-zmin)/NZ+1;

	printf("X: %f %f  DX: %f\n",xmin,xmax,DX);
	printf("Y: %f %f  DY: %f\n",ymin,ymax,DY);
	printf("Z: %f %f  DZ: %f\n",zmin,zmax,DZ);

	MAXBLOCK = (50*nfibres)/(NX*NY*NZ*2);
	printf("MAXBLOCK: %d\n",MAXBLOCK);
	blocks = (int *)malloc(NX*NY*NZ*2*MAXBLOCK*sizeof(int));
	setupBlockLists(net);

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
	fflush(fpout);

	// This code is to determine, for each fibre, for each end, the list of connected fibres.
	nvbin = 0;
	for (k=0; k<nfibres; k++) {
		if (k%10000 == 0) printf("fibre: %d\n",k);
		kv[0] = fibre[k].kv[0];
		nlinks[0] = 0;
		kv[1] = fibre[k].kv[1];
		nlinks[1] = 0;
		int npts = net->edgeList[k].npts;
		// Look at end 1 first
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
					printf("Error: nlinks[0]: %d\n",nlinks[0]);
					exit(1);
				}
			}
		}
		// Now look at end 2
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
					printf("Error: nlinks[1]: %d\n",nlinks[1]);
					exit(1);
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

/*
	nvbin = 0;
	for (k=0; k<nfibres; k++) {
		if (k%1000 == 0) printf("fibre: %d\n",k);
		kv[0] = fibre[k].kv[0];
		nlinks[0] = 0;
		kv[1] = fibre[k].kv[1];
		nlinks[1] = 0;
		for (int kk=0; kk<nfibres; kk++) {
			if (kk == k) continue;
			if (fibre[kk].kv[0] == kv[0]) {
				fibre[k].link[0][nlinks[0]] = kk;
				nlinks[0]++;
				if (nlinks[0] > MAX_LINKS) {
					printf("Error: nlinks[0]: %d\n",nlinks[0]);
					exit(1);
				}
			}
			if (fibre[kk].kv[1] == kv[0]) {
				fibre[k].link[0][nlinks[0]] = kk;
				nlinks[0]++;
				if (nlinks[0] > MAX_LINKS) {
					printf("Error: nlinks[0]: %d\n",nlinks[0]);
					exit(1);
				}
			}
			if (fibre[kk].kv[0] == kv[1]) {
				fibre[k].link[1][nlinks[1]] = kk;
				nlinks[1]++;
				if (nlinks[1] > MAX_LINKS) {
					printf("Error: nlinks[1]: %d\n",nlinks[1]);
					exit(1);
				}
			}
			if (fibre[kk].kv[1] == kv[1]) {
				fibre[k].link[1][nlinks[1]] = kk;
				nlinks[1]++;
				if (nlinks[1] > MAX_LINKS) {
					printf("Error: nlinks[1]: %d\n",nlinks[1]);
					exit(1);
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
//		fprintf(fpout,"fibre: %6d nlinks1,nlinks2: %2d %2d\n",k,nlinks1,nlinks2);
	}

*/
	printf("set up nlinks\n");

	if (centre.x==0 || centre.y==0 || centre.z==0) {
		// Get centre
		centre.x = 0;
		centre.y = 0;
		centre.z = 0;
		for (iv=0; iv<net->nv; iv++) {
			net->vertex[iv].nf = 0;
			centre.x += net->vertex[iv].point.x;
			centre.y += net->vertex[iv].point.y;
			centre.z += net->vertex[iv].point.z;
		}
		centre.x /= net->nv;
		centre.y /= net->nv;
		centre.z /= net->nv;
	}

	// Evaluate angles at vertices, get centre
	for (k=0; k<nfibres; k++) {
		kv[0] = fibre[k].kv[0];
		kv[1] = fibre[k].kv[1];
		net->vertex[kv[0]].fib[net->vertex[kv[0]].nf] = k;
		net->vertex[kv[0]].nf++;
		net->vertex[kv[1]].fib[net->vertex[kv[1]].nf] = k;
		net->vertex[kv[1]].nf++;
	}
	printf("set up vertex.nf\n");
	for (i=0;i<NANGS;i++)
		Angbin[i] = 0;
	nangbin = 0;
	for (iv=0; iv<net->nv; iv++) {
		if (net->vertex[iv].nf < 3) continue;
		for (int i1=0; i1<net->vertex[iv].nf; i1++) {
			int kf1 = net->vertex[iv].fib[i1];
			for (int i2=0; i2<net->vertex[iv].nf; i2++) {
				if (i1 == i2) continue;
				int kf2 = net->vertex[iv].fib[i2];	// Note float counting of pairs of fibres
				double theta = getAngle(iv,kf1,kf2);
				i = theta/dang;
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
			float dot = 0;
			for (i=0; i<3; i++) {
				dot += (v[i]*fibre[k].u[i]);
			}
			dotsum[idir] += abs(dot);
		}
		dotsum[idir] /= 2*nfibres;	// because each dot-product gets evaluated twice
		fprintf(fpout,"idir: %2d v: %6.3f %6.3f %6.3f dot ave: %8.4f\n",idir,v[0],v[1],v[2],dotsum[idir]); 
	}
//	return;
	printf("\nNumber of vertices: %d\n",nvbin/2);
	printf("Length of fibre/vertex: %f\n",Lsum/net->nv);
	printf("Distribution of fibres/vertex:\n");
	fprintf(fpout,"\nNumber of vertices: %d\n",nvbin/2);
	fprintf(fpout,"Length of fibre/vertex: %f\n",Lsum/net->nv);
	fprintf(fpout,"Distribution of fibres/vertex:\n");
	for (i=0; i<10; i++) {
		printf("%2d %8.5f\n",i,float(NBbin[i])/nvbin);
		fprintf(fpout,"%2d %8.5f\n",i,float(NBbin[i])/nvbin);
	}

	// Now set up near fibres
	if (jumpy == 0) return;
	for (k=0; k<nfibres; k++) {
		if (k%10000 == 0) printf("fibre: %d\n",k);
		kv[0] = fibre[k].kv[0];
		kv[1] = fibre[k].kv[1];
		// Look at end 1 first
		int kend = 0;
		int k1 = fibre[k].pt[kend];
		p1 = net->point[k1];
		ix = (p1.x - xmin)/DX;
		iy = (p1.y - ymin)/DY;
		iz = (p1.z - zmin)/DZ;
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
	return;
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
	uniform_real_distribution<double> dist( 0, psum ) ;
	double R = dist(gen);
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
// The cell started at ifib1, iend1, we first determine where it was after travelling distance dist
// Note: fibre[k] is actually edge edgelist[k]
//-----------------------------------------------------------------------------------------------------
double get_partial_distance(NETWORK *net, int ifib0, int iend0, int ifib1, int iend1, double dist)
{
	double dd, dsum, dx, dy, dz, x, y, z, fac;
	EDGE *ep;

	ep = &net->edgeList[ifib1];
	dsum = 0;
	for (int kseg=0; kseg<ep->npts-1; kseg++) {
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
	int kv0 = fibre[ifib0].kv[iend0];
	POINT p0 = net->vertex[kv0].point;
	dx = x - p0.x;
	dy = y - p0.y;
	dz = z - p0.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
// Travel randomly on the network.
// nsteps = max number of steps
// tmax = maximum time (min)
//-----------------------------------------------------------------------------------------------------
void traverse(NETWORK *net, int nsteps, double *tpt, double *res_d2sum, double *Cm)
{
	int i, k, iR, istep, itpt, iseg, reached;
	int *nvisits=NULL;
	double *res_d2tsum;
	double *Lsum=NULL;
	double Ltotal, R, sum;
	double t, dt, speed, std_speed;

	printf("traverse\n");
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
	for (i=0; i<nfibres; i++) {
		if (i == 0)
			Lsum[0] = fibre[0].L_actual;
		else
			Lsum[i] = Lsum[i-1] + fibre[i].L_actual;
		nvisits[i] = 0;
	}
	Ltotal = Lsum[nfibres-1];

//    std::normal_distribution<double> norm_dist( 0.0, 1.0 ) ;
    normal_distribution<double> norm_dist( 0.0, 1.0 ) ;
	
//	std::uniform_real_distribution<double> uni_dist( 0.0, 1.0 ) ;
	uniform_real_distribution<double> uni_dist( 0.0, 1.0 ) ;
	double d, d2sum=0;
	int np = 0;
	for (k=0; k<1000000; k++) {
		R = uni_dist(gen);
		R *= Ltotal;

		// Randomly choose initial direction on the fibre
		int jdir = random_int(0,1);

		// Randomly select a fibre i
		for (i=0; i<nfibres; i++) {
			if (R <= Lsum[i]) {
				break;
			}
		}
		int ifib = i;
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

		speed = mean_speed;
		R = norm_dist(gen);
		speed += R*std_speed;
		reached = 0;
		itpt = 0;
		t = 0;
		for (istep=0; istep<nsteps; istep++) {
			dt = fibre[ifib].L_actual/speed;
//			printf("istep,dt: %d %f\n",istep,dt);
			if (t + dt > tpt[itpt]) {
//				if (itpt == NDATAPTS-1) {	// This is the last point, we can locate it more accurately on the fibre ifib
				dt = tpt[itpt] - t;
				// we need to find where the cell was dt after leaving the other end (not jdir) on fibre ifib
				t += dt;
				d = get_partial_distance(net,ifib0,iend0,ifib,1-jdir,dt*speed);
				reached = 1;
			} else {
				t += dt;
				d = get_distance(net,ifib0,iend0,ifib,jdir);	// d = distance from path start point p0
			}
			d *= shrink_factor;
			d2sum += d*d;
			if (reached) {	// reached the next time point
				res_d2sum[itpt] += d*d;
				res_d2tsum[itpt] += d*d/t;
//				printf("istep: %d t: %f d: %f d2sum: %f\n",istep,t,d,d2sum);
				itpt++;
				if (itpt == NDATAPTS) break;
				reached = 0;
			}
//			} else {
//				t += dt;
//			}
			int kv1 = fibre[ifib].kv[jdir];
			POINT p1 = net->vertex[kv1].point;
			dx = p1.x - p0.x;
			dy = p1.y - p0.y;
			dz = p1.z - p0.z;
			if (save_paths && np < npaths && istep+1 < NTIMES) {
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
		np++;
		printf("Did path: %d  istep: %d\n",np, istep);
		if (np == ntrials) break;
	}
	if (istep >= nsteps) {
		printf("Used up all steps!  istep: %d\n",istep);
		exit(1);
	}
	printf("Number of reversals: %d\n",ndead);
	fprintf(fpout,"Number of reversals: %d\n",ndead);
	printf("npaths: %d\n",np);
//	*Cm = d2sum/(6*npaths*tmax);
	for (itpt=0; itpt<NDATAPTS; itpt++) {
		Cm[itpt] = res_d2tsum[itpt]/(6*np);		// this is with fixed time - actual t is used, not tmax
	}
	//for (i=0; i<nfibres; i++) {
	//	if (nvisits[i] == 0) printf("No visit to fibre: %d\n",i);
	//}
}

//-----------------------------------------------------------------------------------------------------
// Count deadends in the network
//-----------------------------------------------------------------------------------------------------
int NumberOfDeadends(NETWORK *net)
{
	float len;
	int k1, k2;
	int nshort = 0;
	int ndeadends = 0;
	for (int i=0; i<nfibres; i++) {
		for (int k=0; k<2; k++) {
			if (fibre[i].nlinks[k] == 0) {
				ndeadends++;
				k1 = fibre[i].pt[0];
				k2 = fibre[i].pt[1];
				len = distance(net, k1, k2);
				if (len <= min_len) nshort++;
			}
		}
	}
	printf("Number of short dead ends: %d\n",nshort);
	fprintf(fpout,"Number of short dead ends: %d\n",nshort);
	return ndeadends;
}

//-----------------------------------------------------------------------------------------------------
// Try to find a vertex to connect a dead end to.
//-----------------------------------------------------------------------------------------------------
void connect(NETWORK *NP1)
{
	int i, j, j0, j1, k, kdead, k1, k2;
	int deadend[2], nd, nfibres_temp, kv0, kv1, kmax, kv, pt0, pt1, ptmax, ntdead, njoins, npts, nshort, nfree;
	double v[3], w[3], dv, dw, cosa, value, valuemax, len;
	POINT p0, p1, q;
	double a = 0.5;
	bool joined;
	bool do_join = true;

	struct join_str
	{
		int kv[2];
		int pt[2];
		bool duplicate;
	};
	typedef join_str JOIN;	
	JOIN *join;

	printf("\nconnect\n");
	join = (JOIN *)malloc(nfibres*sizeof(JOIN));	// guaranteed big enough
	ntdead = 0;
	njoins = 0;
	nshort = 0;
	nfree = 0;
	nfibres_temp = nfibres;
	for (i=0; i<nfibres_temp; i++) {
		nd = 0;
		for (k=0; k<2; k++) {
			if (fibre[i].nlinks[k] == 0) {
				k1 = fibre[i].pt[0];
				k2 = fibre[i].pt[1];
				len = distance(NP1, k1, k2);
				//if (len <= min_len) {	// Here short deadends were dropped, now they are left to make possible connections
				//	NP1->edgeList[i].used = false;
				//	nshort++;
				//	nd = 0;
				//	break;
				//}
				deadend[nd] = k;
				nd++;
			}
		}
		if (nd == 0) continue;

		if (!do_join) continue;

		joined = false;
		for (kdead=0; kdead<nd; kdead++) {
			ntdead++;
			kv1 = fibre[i].kv[deadend[kdead]];
			pt1 = fibre[i].pt[deadend[kdead]];
			if (deadend[kdead] == 0) {
				kv0 = fibre[i].kv[1];
				pt0 = fibre[i].pt[1];
			} else {
				kv0 = fibre[i].kv[0];
				pt0 = fibre[i].pt[0];
			}
			if (ntdead%100 == 0)
				printf("fibre: %d nlinks: %d %d  kv0, kv1: %d %d  ntdead: %d\n",i,fibre[i].nlinks[0],fibre[i].nlinks[1],kv0,kv1,ntdead);
			p0 = NP1->vertex[kv0].point;
			p1 = NP1->vertex[kv1].point;
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
			valuemax = 0;
			for (kv=0; kv<NP1->nv; kv++) {
				if (kv == kv1) continue;
				q = NP1->vertex[kv].point;
				w[0] = q.x - p1.x;
				w[1] = q.y - p1.y;
				w[2] = q.z - p1.z;
				dw = 0;
				for (k=0; k<3; k++) {
					dw += w[k]*w[k];
				}
				dw = sqrt(dw);	// distance to the vertex
				if (dw > max_len) continue;
				cosa = 0;
				for (k=0; k<3; k++) {
					cosa += v[k]*w[k];
				}
				if (cosa < 0) continue;
				cosa /= dw;
				value = (a + cosa)/dw;
				if (value > valuemax) {
					kmax = kv;
					valuemax = value;
	//				printf("\ncosa,d,valuemax: %f %f %f\n",cosa,dw,valuemax);
				}
			}
			if (kmax < 0) {
//				printf("No suitable vertex for: %d\n",kv1);
			} else {
				joined = true;
				join[njoins].kv[0] = kv1;
				join[njoins].kv[1] = kmax;
				join[njoins].pt[0] = pt1;
//				join[njoins].pt[1] = ?;	// what is the point index corresponding to vertex index kmax?
				// need to look at any fibre connected to the vertex.  Choose #0
				int kf = NP1->vertex[kmax].fib[0];
				if (fibre[kf].kv[0] == kmax) {
					join[njoins].pt[1] = fibre[kf].pt[0];
				} else {
					join[njoins].pt[1] = fibre[kf].pt[1];
				}

				join[njoins].duplicate = false;
				njoins++;
//				fprintf(fpout,"%8d %8d\n",kv1,kmax);
				if (ntdead%100 == 0) 
					printf("Nearest vertex to: %d is: %d\n",kv1,kmax);
			}
		}
		if (!joined) {
			NP1->edgeList[i].used = false;
			nfree++;
		}
	}
	printf("Total short dead ends: %d\n",nshort);
	fprintf(fpout,"Total short dead ends: %d\n",nshort);
	printf("Total unconnectable dead ends: %d\n",nfree);
	fprintf(fpout,"Total unconnectable dead ends: %d\n",nfree);
	printf("Total candidate joins: %d\n",njoins);
	fprintf(fpout,"Total candidate joins: %d\n",njoins);

	// Now clear out duplicate joins
	printf("Clearing out duplicates\n");
	for (j0=0; j0<njoins; j0++) {
		if (join[j0].duplicate) continue;
		kv0 = join[j0].kv[0];
		kv1 = join[j0].kv[1];
		for (j1=0; j1<njoins; j1++) {
			if (join[j1].duplicate) continue;
			if (kv0 == join[j1].kv[1] && kv1 == join[j1].kv[0]) {
//				printf("Duplicate: %d %d\n",j0,j1);
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
	int ne = NP1->ne - nfree;	// - nshort
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
	char drive[32], dir[128],filename[256], ext[32];
	char errfilename[256], output_amfile[256];
	char output_basename[256];
	int do_connect;
	float origin_shift[3];
	int limit_mode;
	float limit_value;

	int itpt;
	double *tpt;
	double d2ave, Cm1;
	double res_d2sum[NDATAPTS];
	double Cm[NDATAPTS];
	int cmgui_flag=1;
	NETWORK *NP0, *NP1;

	printf("Two functions: (1) To join dead ends in the network (if possible), \n                   and trim unconnected dead ends\n");
	printf("               (2) To estimate the coefficient of motility Cm\n");
	printf("\n");
	printf("These two functions are carried out by separate program executions.  \nThe network should be healed (joined, trimmed) before Cm estimation.\n");
	printf("For the healing pass, set max_len \n(unscaled upper limit on length of an added connection).\n");
	printf("The shrinkage compensation factor must be specified\n(Note that the dimensions in the am file are left unscaled)\n");
	printf("\n");
	if (argc != 9 && argc != 16) {
		printf("To perform joining and trimming of dead ends:\n");
		printf("Usage: conduit_analyse input_amfile output_file sfactor max_len limit_mode limit_value ddiam dlen\n");
		printf("       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		printf("       max_len     = upper limit on (unscaled) connecting fibre length (um)\n");
		printf("       limit_mode  = restriction for computing fibre statistics:\n                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		printf("       limit_value = value of limit for the chosen mode\n");
		printf("       ddiam       = bin size for diameter distribution\n");
		printf("       dlen        = bin size for length distribution\n");
		printf("\nTo simulate cell paths and estimate Cm:\n");
		printf("Usage: conduit_analyse input_amfile output_file sfactor npow ntrials x0 y0 z0 radius speed npaths limit_mode limit_value ddiam dlen\n");
		printf("       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		printf("       npow        = power of cos(theta) in weighting function for prob. of \n                     taking a branch (theta = turning angle)\n");
		printf("       ntrials     = number of cell paths simulated\n");
		printf("       x0,y0,z0    = centre of sphere within which the cell paths start (um)\n");
		printf("       radius      = radius of sphere within which the cell paths start (um)\n");
		printf("       speed       = mean cell speed (um/min) (Note that CV = %6.2f)\n",CV);
		printf("       npaths      = number of cell paths to save\n");
		printf("       limit_mode  = restriction for computing fibre statistics:\n                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		printf("       limit_value = value of limit for the chosen mode\n");
		printf("       ddiam       = bin size for diameter distribution\n");
		printf("       dlen        = bin size for length distribution\n");
		printf("\n");
		fpout = fopen("conduit_analyse_error.log","w");
		fprintf(fpout,"To perform joining and trimming of dead ends:\n");
		fprintf(fpout,"Usage: conduit_analyse input_amfile output_file sfactor max_len\n");
		fprintf(fpout,"       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		fprintf(fpout,"       max_len     = upper limit on (unscaled) connecting fibre length (um)\n");
		fprintf(fpout,"       limit_mode  = restriction for computing fibre statistics:\n                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
		fprintf(fpout,"       limit_value = value of limit for the chosen mode\n");
		fprintf(fpout,"       ddiam       = bin size for diameter distribution\n");
		fprintf(fpout,"       dlen        = bin size for length distribution\n");
		fprintf(fpout,"\nTo simulate cell paths and estimate Cm:\n");
		fprintf(fpout,"Usage: conduit_analyse input_amfile output_file sfactor npow ntrials x0 y0 z0 radius speed npaths\n");
		fprintf(fpout,"       sfactor     = shrinkage compensation factor e.g. 1.25\n");
		fprintf(fpout,"       npow        = power of cos(theta) in weighting function for prob. of \n                     taking a branch (theta = turning angle)\n");
		fprintf(fpout,"       ntrials     = number of cell paths simulated\n");
		fprintf(fpout,"       x0,y0,z0    = centre of sphere within which the cell paths start (um)\n");
		fprintf(fpout,"       radius      = radius of sphere within which the cell paths start (um)\n");
		fprintf(fpout,"       speed       = mean cell speed (um/min) (Note that CV = %6.2f)\n",CV);
		fprintf(fpout,"       npaths      = number of cell paths to save\n");
		fprintf(fpout,"       limit_mode  = restriction for computing fibre statistics:\n                     0 = no limit, 1 = len limit, 2 =  len/diam limit\n");
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
	sscanf(argv[3],"%lf",&shrink_factor);
	if (argc == 9) {
		sscanf(argv[4],"%lf",&max_len);
		sscanf(argv[5],"%d",&limit_mode);
		sscanf(argv[6],"%f",&limit_value);
		sscanf(argv[7],"%f",&ddiam);
		sscanf(argv[8],"%f",&dlen);
		do_connect = true;

		//char add_str[13];
		//if (max_len < 10)
		//	sprintf(add_str,"_jt_max%3.1f.am",max_len);
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

	} else if (argc == 16) {
		sscanf(argv[4],"%d",&npow);
		sscanf(argv[5],"%d",&ntrials);
		sscanf(argv[6],"%f",&centre.x);
		sscanf(argv[7],"%f",&centre.y);
		sscanf(argv[8],"%f",&centre.z);
		sscanf(argv[9],"%lf",&start_radius);
		sscanf(argv[10],"%lf",&mean_speed);
		sscanf(argv[11],"%d",&npaths);
		sscanf(argv[12],"%d",&limit_mode);
		sscanf(argv[13],"%f",&limit_value);
		sscanf(argv[14],"%f",&ddiam);
		sscanf(argv[15],"%f",&dlen);
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

	/*
	_splitpath(input_amfile,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(result_file,"%s_fd.out",output_basename);
	printf("output_basename: %s\n",output_basename);
	*/
	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	NP1 = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;
	printf("Read Amira file\n");

//	for (int ie=0; ie<NP0->ne; ie++) {
//		printf("fibre: %d npts: %d\n",ie,NP0->edgeList[ie].npts);
//	}
//	exit(0);

	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;
//	printf("ne: %d\n",NP0->ne);
	err = CreateDistributions(NP0);
	if (err != 0) return 4;

	makeFibreList(NP0);
	printf("did makeFibreList\n");
	ndead = NumberOfDeadends(NP0);
	fprintf(fpout,"Number of fibres, deadends in the network: %d  %d\n",nfibres,ndead);

	if (do_connect) {
		char add_str[13];
		if (max_len < 10)
			sprintf(add_str,"_jt_max%3.1f.am",max_len);
		else
			sprintf(add_str,"_jt_max%4.1f.am",max_len);
		strcpy(output_amfile,input_amfile);
		int len;
		for(len=strlen(output_amfile);len >= 0 && output_amfile[len] != '.';len--);
		output_amfile[len] = '\0';
		strcat(output_amfile,add_str);
		connect(NP0);
		printf("\noutput_amfile: %s\n",output_amfile);
		fprintf(fpout,"\noutput_amfile: %s\n",output_amfile);
		err = WriteAmiraFile(output_amfile,input_amfile,NP0,origin_shift);
		if (err != 0) return 3;
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

	ndead = 0;
	traverse(NP0,1000,tpt,res_d2sum,Cm);
	printf("# of deadends encountered: %d\n",ndead);
	printf("\nCm estimation results:\n");
	fprintf(fpout,"# of deadends encountered: %d\n",ndead);
	fprintf(fpout,"\nCm estimation results:\n");
	for (itpt=0; itpt<NDATAPTS; itpt++) {
		d2ave = res_d2sum[itpt]/ntrials;
		Cm1 = d2ave/(6*tpt[itpt]);		// this is with fixed time - actual t is used, not tmax
		printf("t: %6.1f d2: %9.1f Cm: %6.2f\n",tpt[itpt],d2ave,Cm1);
		fprintf(fpout,"t: %6.1f d2: %9.1f Cm: %6.2f\n",tpt[itpt],d2ave,Cm1);		// dropped Cm[itpt] - no difference now
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
	err = WriteAmiraFile(output_amfile,input_amfile,NP1,origin_shift);
	if (err != 0) return 4;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,NP1,origin_shift);
		if (err != 0) return 5;
	}
	return 0;
	*/
}
