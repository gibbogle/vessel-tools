// prune.cpp
// This program has an error that causes dead-ends to be messed up (wrapped back on themselves). FIXED

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

#include "prune.h"

struct pair_str
{
	int i1, i2;
};
typedef pair_str PAIR;

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


int WriteCmguiData(char *basename);

int nv, ne, np;
int nv_used, ne_used, np_used;
EDGE *edgeList;
VERTEX *vertex;
POINT *point;

FILE *fperr, *fpout;
FILE *exelem, *exnode;
char output_basename[128];
double ratio_limit;

#define NBOX 400
#define STR_LEN 128
//#define RATIO_LIMIT 4	// should be an input parameter
#define PI 3.14159

#define JOIN_LOOSE_ENDS false

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float dist(int k1, int k2)
{
	float dx = point[k2].x - point[k1].x;
	float dy = point[k2].y - point[k1].y;
	float dz = point[k2].z - point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

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
//-----------------------------------------------------------------------------------------------------
bool samevalue(float v1, float v2)
{
	float tol = 0.01;

	return (fabs(v1-v2) < tol);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool samepoint(POINT *p1, POINT *p2)
{
	if (!samevalue(p1->x,p2->x)) return false;
	if (!samevalue(p1->y,p2->y)) return false;
	if (!samevalue(p1->z,p2->z)) return false;
	return true;
}

//-----------------------------------------------------------------------------------------------------
// Check consistency between edge vert[] and pt[]
//-----------------------------------------------------------------------------------------------------
int checkEdgeEndPts()
{
	EDGE *edge;
	POINT *pv0, *pv1, *pt0, *pt1;
	int ie, npts, k, kp0, kp1, kv0, kv1, kvp0, kvp1, err;

	err = 0;
	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		npts = edge->npts;
		pt0 = &point[edge->pt[0]];
		pt1 = &point[edge->pt[npts-1]];
		kv0 = edge->vert[0];
		kv1 = edge->vert[1];
		pv0 = &vertex[kv0].point;
		pv1 = &vertex[kv1].point;
//		if ((kp0==kvp0 && kp1==kvp1) || (kp0==kvp1 && kp1==kvp0)) continue;
		if ((samepoint(pt0,pv0) && samepoint(pt1,pv1)) || (samepoint(pt0,pv1) && samepoint(pt1,pv0))) continue;
		err = 1;
		fprintf(fpout,"checkEdgeEndPts: edge: %d vert: %d %d -> %d %d pt: %d %d\n",ie,kv0,kv1,kvp0,kvp1,kp0,kp1);
		fflush(fpout);
	}
	return err;
}

//-----------------------------------------------------------------------------------------------------
// The average diameter of a vessel (edge) is now estimated by dividing the volume by the length.
// All distances in the .am file are in um.
//-----------------------------------------------------------------------------------------------------
int CreateDistributions(double delta_diam, double delta_len)
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double lsegadbox[NBOX];
//	double ddiam, dlen;
	double ad, len, dlen, ltot, lsum, dsum, dvol, r2, r2prev, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam, ave_seg_diam;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot;
	double lenlimit = 4.0;
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
//	ddiam = 0.5;
	ndtot = 0;
	nsegdtot = 0;
	lsegdtot = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	volume = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
//		printf("ie: %d npts: %d\n",ie,edge.npts);
		fprintf(fperr,"ie: %d npts: %d\n",ie,edge.npts);
		fflush(fperr);
		bool dbug = false;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		dvol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = point[kp].d;
//			ad = avediameter[kp];
			ave_pt_diam += ad;
			if (dbug) {
				printf("%d  %d  %f  %f\n",ip,kp,ad,delta_diam);
				fprintf(fperr,"%d  %d  %f  %f\n",ip,kp,ad,delta_diam);
			}
			fflush(fperr);
//			dsum += ad;
			if (ad < 0.001) {
				printf("Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fperr,"Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				return 1;
			}
			ka = int(ad/delta_diam + 0.5);
			if (ka >= NBOX) {
				printf("Vessel too wide (point): d: %f k: %d\n",ad,ka);
				fprintf(fperr,"Vessel too wide (point): d: %f k: %d\n",ad,ka);
				continue;
			}
			adbox[ka]++;
			ndtot++;
			if (ip > 0) {
				dlen = dist(kp,kpprev);
				r2 = ad*ad/4;
				dvol += PI*dlen*(r2 + r2prev)/2;
				lsum += dlen;
			}
			kpprev = kp;
			r2prev = r2;
		}
		edgeList[ie].length_um = lsum;
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
//		ka = int(ad/ddiam + 0.5);
		ka = int(ad/delta_diam );
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
			d95 = (k-1)*delta_diam;
			break;
		}
	}
	printf("Compute length distributions: lower limit = %6.1f um\n",lenlimit);
	fprintf(fperr,"Compute length distributions: lower limit = %6.1f um\n",lenlimit);
	// Lengths
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/delta_len);
		if (len <= lenlimit) continue;
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
	fprintf(fpout,"Total vertices: %d  points: %d\n",nv,np);
	fprintf(fpout,"Vessels: %d\n",ne);
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
	fprintf(fpout,"'um'    number  fraction    length  fraction\n");
	for (k=0; k<ndpts; k++) {
		fprintf(fpout,"'%6.2f -%6.2f' %8d %9.5f  %8.0f %9.5f\n",k*delta_diam,(k+1)*delta_diam,segadbox[k],segadbox[k]/float(nsegdtot),
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
	fprintf(fpout,"\nVessel length distribution\n");
	fprintf(fpout,"'um'    number  fraction\n");
	for (k=0; k<nlpts; k++) {
		fprintf(fpout,"'%6.2f -%6.2f' %8d %9.5f\n",k*delta_len,(k+1)*delta_len,lvbox[k],lvbox[k]/ltot);
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// The average diameter of a segment (edge) is now estimated by dividing the volume by the length.
//-----------------------------------------------------------------------------------------------------
int oldCreateDistributions(void)
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double ad, len_um, ddiam, dlen, dtot, ltot, lsum, dsum, dvol, r2, r2prev;
	double segdtot;
	double ave_len_um, volume, d95;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts;
	double lenlimit = 3.0;
	EDGE edge;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lvbox[k] = 0;
	}
	printf("Compute diameter distributions\n");
	fprintf(fperr,"Compute diameter distributions\n");
	// Diameters
	ddiam = 0.5;
	dtot = 0;
	segdtot = 0;
	volume = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		dvol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = point[kp].d;
//			dsum += ad;
			if (ad < 0.001) {
				fprintf(fperr,"Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				return 1;
			}
			ka = ad/ddiam;
			if (ka >= NBOX) {
				printf("Vessel too wide (point): d: %f k: %d\n",ad,ka);
				continue;
			}
			adbox[ka]++;
			dtot++;
			if (ip > 0) {
				dlen = dist(kp,kpprev);
				r2 = ad*ad/4;
				dvol += PI*dlen*(r2 + r2prev)/2;
				lsum += dlen;
			}
			kpprev = kp;
			r2prev = r2;
		}
		edgeList[ie].length_um = lsum;
		volume += dvol;
//		ad = dsum/edge.npts;
		ad = 2*sqrt(dvol/(PI*lsum));
		if (ad < 0.001) {
			fprintf(fperr,"Zero edge diameter: edge: %d ad: %f\n",ie,ad);
			return 1;
		}
		ka = ad/ddiam;
		if (ka >= NBOX) {
			printf("Vessel too wide (segment ave): d: %f k: %d\n",ad,ka);
			continue;
		}
		segadbox[ka]++;
		segdtot++;
	}
	// Determine d95, the diameter that >95% of points exceed.
	dsum = 0;
	for (k=0; k<NBOX; k++) {
		dsum += adbox[k]/dtot;
		if (dsum > 0.05) {
			d95 = (k-1)*ddiam;
			break;
		}
	}
//	printf("Compute length distributions: lower limit = d95: %f\n",d95);
//	printf("Compute length distributions: no lower limit\n");
	printf("Compute length distributions: lower limit = %6.1f um\n",lenlimit);
	// Lengths
	dlen = 1;
	ltot = 0;
	ave_len_um = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		len_um = edge.length_um;
//		if (len_um <= d95) continue;
		k = len_um/dlen + 0.5;
		if (k*dlen <= lenlimit) continue;
		if (k >= NBOX) {
			printf("Edge too long: len_um: %d  %f  k: %d\n",ie,len_um,k);
			continue;
		}
		lvbox[k]++;
		ave_len_um += len_um;
		ltot++;
	}
	fprintf(fpout,"Total vertices: %d  points: %d\n",nv,np);
	fprintf(fpout,"Edges: %d\n",ne);
	printf("Average length: vox: %f\n",ave_len_um/ltot);
	fprintf(fpout,"Average length: vox: %f\n",ave_len_um/ltot);
	fprintf(fpout,"Volume: %10.0f\n\n",volume);

	ndpts = 0;
	for (k=0; k<NBOX; k++) {
		if (adbox[k]/dtot >= 0.0005) {
			ndpts = k+1;
		}
	}
	fprintf(fpout,"   um   ave_diam  seg_ave\n");
	for (k=0; k<ndpts; k++) {
		fprintf(fpout,"%6.2f %9.4f %9.4f\n",k*ddiam,adbox[k]/dtot,segadbox[k]/segdtot);
	}

	nlpts = 0;
	for (k=0; k<NBOX; k++) {
		if (lvbox[k]/ltot >= 0.0005) {
			nlpts = k+1;
		}
	}
	fprintf(fpout,"   um   vox_len\n");
	for (k=0; k<nlpts; k++) {
		fprintf(fpout,"%6.2f %9.4f\n",k*dlen,lvbox[k]/ltot);
//		printf("%3d %8.3f %8.3f %8.3f %8.3f\n",k,adbox[k]/dtot,mdbox[k]/dtot,lvbox[k]/ltot,ljbox[k]/ltot);
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
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile(char *amFile)
{
	int i, j, k, kp, npts, nee, npp;
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
			sscanf(line+13,"%d",&nv);
			k++;
		}
		if (strncmp(line,"define EDGE",11) == 0) {
			sscanf(line+11,"%d",&ne);
			nee = 2*ne;
			k++;
		}
		if (strncmp(line,"define POINT",12) == 0) {
			sscanf(line+12,"%d",&np);
			npp = 2*np;
			k++;
		}
	}

	vertex = (VERTEX *)malloc(nv*sizeof(VERTEX));
	edgeList = (EDGE *)malloc(nee*sizeof(EDGE));	// 2* for added joining edges
	point = (POINT *)malloc(npp*sizeof(POINT));		// 2* for added joining edges

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
					edgeList[i].length_um = len;
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
						if (j == 13) {
							printf("j=13: d: %f\n",point[j].d);
						}
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
// Check that all edges are completed, i.e. no vertex appears within an edge.
//-----------------------------------------------------------------------------------------------------
int adjoinEdges(void)
{
	int *nvrefs;
	PAIR *pair;
	int ie, iv, kv0, kv1, ivmax, nloose, i1, i2, n1, n2, k1, k2, i, ntwo, nvunused, neunused;
	int temp[1000];	// should be big enough
	EDGE edge1, edge2;

	printf("adjoinEdges\n");
	nvrefs = (int *)malloc(nv*sizeof(int));
	pair = (PAIR *)malloc(nv*sizeof(PAIR));
	for (;;) {
		for (iv=0; iv<nv; iv++) {
			nvrefs[iv] = 0;
			pair[iv].i1 = 0;
			pair[iv].i2 = 0;
		}
		neunused = 0;
		ivmax = 0;
		for (ie=0; ie<ne; ie++) {
			edge1 = edgeList[ie];
			if (!edge1.used) {
				neunused++;
//				fprintf(fpout,"Unused edge: %d  %d\n",neunused,ie);
				continue;
			}
			n1 = edge1.npts;
			kv0 = edge1.pt[0];
			kv1 = edge1.pt[n1-1];
			if (nvrefs[kv0] == 0)
				pair[kv0].i1 = ie;
			else if (nvrefs[kv0] == 1)
				pair[kv0].i2 = ie;
			nvrefs[kv0]++;
			if (nvrefs[kv1] == 0)
				pair[kv1].i1 = ie;
			else if (nvrefs[kv1] == 1)
				pair[kv1].i2 = ie;
			nvrefs[kv1]++;
			if (kv0 > ivmax) ivmax = kv0;
			if (kv1 > ivmax) ivmax = kv1;
		}
		ntwo = 0;
		nloose = 0;
		nvunused = 0;
		for (iv=0; iv<=ivmax; iv++) {
			if (nvrefs[iv] == 0) {
	//			printf("Unused: %d\n",iv);
				nvunused++;
			} else if (nvrefs[iv] == 1) {
				nloose++;
			} else if (nvrefs[iv] == 2) {
				i1 = pair[iv].i1;
				i2 = pair[iv].i2;
				edge1 = edgeList[i1];
				edge2 = edgeList[i2];
				if (!edge1.used || !edge2.used) continue;	// already processed, adjoined
				ntwo++;
				// The two edges are edge1 and edge2
				// They will be combined into one, and replace edge1, edge2 will be deleted
				n1 = edge1.npts;
				kv0 = edge1.pt[0];
				kv1 = edge1.pt[n1-1];
				if (kv0 == iv)
					k1 = 0;
				else
					k1 = n1-1;
				n2 = edge2.npts;
				kv0 = edge2.pt[0];
				kv1 = edge2.pt[n2-1];
	//			fprintf(fpout,"Two occurrences: %d  %d %d\n",iv,i1,i2);
	//			if ((kv0 == edge1.pt[0] && kv1 == edge1.pt[n1-1]) ||
	//				(kv1 == edge1.pt[0] && kv0 == edge1.pt[n1-1])) {	// this is a loop
	//				printf("Duplicated edge: %d\n",i2);
	//				edgeList[i2].used = false;
	//				continue;
	//			}
				if (kv0 == iv)
					k2 = 0;
				else
					k2 = n2-1;
				if (k1 == 0 && k2 == 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[n1-i-1];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[i];
				} else if (k1 == 0 && k2 != 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[n1-i-1];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[n2-i-1];
				} else if (k1 != 0 && k2 == 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[i];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[i];
				} else if (k1 != 0 && k2 != 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[i];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[n2-i-1];
				}
				edge1.npts = n1 + n2 - 1;
				edge1.npts_used = edge1.npts;
				free(edge1.pt);
				free(edge1.pt_used);
				edge1.pt = (int *)malloc(edge1.npts*sizeof(int));
				edge1.pt_used = (int *)malloc(edge1.npts*sizeof(int));
				for (i=0; i<edge1.npts; i++) {
					edge1.pt[i] = temp[i];
					edge1.pt_used[i] = temp[i];
				}
				edge1.vert[0] = temp[0];
				edge1.vert[1] = temp[edge1.npts-1];
				edge1.used = true;
				edgeList[i1] = edge1;
				edgeList[i2].used = false;
//				fprintf(fpout,"Adjoined edge: %d to edge: %d\n",i2,i1);
			}
		}
		printf("Number of two-healings: %d\n",ntwo);
		printf("Number of unused edges: %d\n",neunused);
		printf("Number of unused vertices: %d\n",nvunused);
		if (ntwo == 0) break;
	}
	free(nvrefs);
	free(pair);
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
// Set up the othogonal unit vectors as an axis system with origin at P0.
// Vx is towards P1, Vy is normal to Vx (and a global axis), Vz is orthogonal to Vx and Vy
//-----------------------------------------------------------------------------------------------------
int makeAxes(double *e0, double p0[], double p1[], double Vx[], double Vy[], double Vz[])
{
	int i;
	double a[3], d;

	for (i=0; i<3; i++) {
		Vx[i] = p1[i] - p0[i];
	}
	d = sqrt(dotproduct(Vx,Vx));
	for (i=0; i<3; i++) {
		Vx[i] /= d;
	}
	a[0] = 1; a[1] = 0; a[2] = 0;	// try X axis
	d = sqrt(fabs(dotproduct(Vx,a)));
	if (d > 0.9) {
		a[0] = 0; a[1] = 1; a[2] = 0;	// try Y axis
	}
	crossProduct(Vy,Vx,a);
	d = sqrt(dotproduct(Vy,Vy));
	for (i=0; i<3; i++) {
		Vy[i] /= d;
	}
	crossProduct(Vz,Vx,Vy);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int getP(double d, double theta, double phi, double p0[], double Vx[], double Vy[], double Vz[], double p[])
{
	int i;
	double L, x, y, z;

	L = d/(2*sin(theta)*cos(phi));
	x = sin(theta)*cos(phi);
	y = sin(theta)*sin(phi);
	z = cos(theta);
	for (i=0; i<3; i++) {
		p[i] = p0[i] + L*(x*Vx[i] + y*Vy[i] + z*Vz[i]);
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// The quantity we seek to minimise is the sum of the squares of the three angles, alpha1, alpha2, alpha3,
// between the edge segments.
// alpha1 is the angle between e0 and P0-P
// alpha2 is the angle between P0-P and P-P1
// alpha3 is the angle between P-P1 and -e1
//-----------------------------------------------------------------------------------------------------
double Q(double *e0, double *e1, double p0[], double p1[], double Vx[], double Vy[], double Vz[],
	double d, double theta, double phi)
{
	double p[3], v0[3], v1[3];
	double alpha1, alpha2, alpha3;
	double d0, d1, cosa;
	int i;

	getP(d,theta,phi,p0,Vx,Vy,Vz,p);

	for (i=0; i<3; i++) {
		v0[i] = p[i] - p0[i];
		v1[i] = p1[i] - p[i];
	}
	d0 = sqrt(dotproduct(v0,v0));
	d1 = sqrt(dotproduct(v1,v1));
	cosa = dotproduct(e0,v0)/d0;
	alpha1 = acos(cosa);
	cosa = dotproduct(v0,v1)/(d0*d1);
	alpha2 = acos(cosa);
	cosa = -dotproduct(v1,e1)/d1;
	alpha3 = acos(cosa);
	return alpha1*alpha1 + alpha2*alpha2 + alpha3*alpha3;
}

//-----------------------------------------------------------------------------------------------------
// Point kp0 is the end of an edge that has a final segment with orientation of the unit vector e0,
// point kp1 is the end of an edge that has a final segment with orientation of the unit vector e1.
// We seek a point P equidistant from P0 and P1 such that the sum of squares of angles along the
// connection P0-P-P1 is minimised.  The three angles are:
//   the angle e0 makes with P0-P = alpha1
//   the angle P0-P makes with P-P1 = alpha2
//   the angle P-P1 makes with -e1 = alpha3
//   Q = alpha1^2 + alpha2^2 + alpha3^2 is to be minimised.
// 
// At the end P0 we start by defining three orthogonal unit vectors that become the axes of
// a local coordinate system, represented conveniently in spherical coordinates.
// Define Vx as the unit vector corresponding to P0-P1.
// Define Vy as any unit vector perpendicular to Vx.
// Define Vz as the unit vector perpendicular to Vx and Vy, Vz = Vx x Vy (cross-product).
//
// Then a general unit vector is given by:
//   Vu = x.Vx + y.Vy + z.Vz
// subject to:
//   x^2 + y^2 + z^2 = 1 (which implies that -1<=x<=1, -1<=y<=1, -1<=z<=1)
// A general unit vector in this coordinate system can be parametrised by theta and phi,
// where theta is the angle with the Z axis (Vx) 0 <= theta <= PI
// and phi is the angle about the Z axis, starting from 0 at the X axis, 0 <= phi <= 2PI
// where:
//   x = sin(theta).cos(phi)
//   y = sin(theta).sin(phi)
//   z = cos(theta)
// For phi = 0, a line from P0 in the direction Vu, with L = d/(2.sin(theta)), 
// ends at a point P that is equidistant from P0 and P1.
// If the line segment length is L, for any phi, x = d/2 ==> L = d/(2.sin(theta).cos(phi))
// Then the parameters (theta, phi) give the point P, offset by the corresponding vector from P0, as:
//   P = P0 + L.(x.Vx + y.Vy + z.Vz)
// From e0, P-P0, P1-P, and e1 the angles alpha1, alpha2, and alpha3 can be computed, hence Q.
// We need to determine the (theta, phi) pair that minimises Q.
// Let t = theta, p = phi.
// Q = Q(t,p)
// The minimum occurs where dQ/dt = 0 and dQ/dp = 0.
// Let F1(t,p) = dQ/dt,
//     F2(t,p) = dQ/dp.
// We now have a two-dimensional root-finding problem, where the evaluation of F1 and F2 is numerical.
// Using Newton-Raphson, the derivatives dFi/dt and dFi/dp (i = 1,2) are numerical.
// dF1/dt = d/dt(dQ/dt) = (Q(t+dt,p) - 2.Q(t,p) + Q(t-dt,p))/(dt^2)
// dF1/dp = d/dp(dQ/dt) = ((Q(t+dt,p+dp) - Q(t,p+dp))/dt - (Q(t+dt,p) - Q(t,p))/dt)/dp
// dF2/dp = d/dp(dQ/dp) = (Q(t,p+dp) - 2.Q(t,p) + Q(t,p-dp))/(dp^2)
// dF2/dt = d/dt(dQ/dp) = ((Q(t+dt,p+dp) - Q(t+dt,p))/dp - (Q(t,p+dp) - Q(t,p))/dp)/dt
//
// We solve:
//   |J11 J12| |dt| = |-F1|
//   |J21 J22| |dp|   |-F2|
//-----------------------------------------------------------------------------------------------------
int joining_edge(double *e0, double *e1, int kp0, int kp1)
{
	double p0[3], p1[3], p[3];		// locations of P0 and P1
	double Vx[3], Vy[3], Vz[3];		// orthogonal unit vectors
	double d;	// distance from P0 -> P1
	double theta;	// theta
	double phi;	// phi
	double L;	// distance from P0 -> P, and from P -> P1
	double J11, J12, J21, J22;		// dF1/dt, dF1/dp, dF2/dt, dF2/dp
	double Qzz, Qzm, Qmz, Qzp, Qpz, Qpp;	// z = 0, m = -1, p = +1
	double F1, F2;
	double del, dt, dp, len;
	int i, k, kp;
	double epsilon = 0.00001;

	p0[0] = point[kp0].x; p0[1] = point[kp0].y; p0[2] = point[kp0].z; 
	p1[0] = point[kp1].x; p1[1] = point[kp1].y; p1[2] = point[kp1].z; 
	makeAxes(e0,p0,p1,Vx,Vy,Vz);

	d = 0;
	for (i=0; i<3; i++) {
		del = p0[i]-p1[i];
		d += del*del;
	}
	d = sqrt(d);
	theta = PI/2;	// initial guess of theta 
	phi = 0;		// and phi: point midway between P0 and P1
//	printf("P0: %8.4f %8.4f %8.4f P1: %8.4f %8.4f %8.4f\n",p0[0],p0[1],p0[2],p1[0],p1[1],p1[2]);
//	printf("d = %8.4f\n",d);
	del = 0.01;

	k = 0;
	for (;;) {
		k++;
		Qzz = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta,phi);
		Qmz = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta-del,phi);
		Qzm = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta,phi-del);
		Qpp = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta+del,phi+del);
		Qpz = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta+del,phi);
		Qzp = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta,phi+del);
		F1 = (Qpz - Qzz)/del;
		F2 = (Qzp - Qzz)/del;
		/*
		printf("Qzz: %10.6f\n",Qzz);
		printf("Qmz: %10.6f\n",Qmz);
		printf("Qzm: %10.6f\n",Qzm);
		printf("Qpp: %10.6f\n",Qpp);
		printf("Qpz: %10.6f\n",Qpz);
		printf("Qzp: %10.6f\n",Qzp);
		
		printf("Qzz: %10.6f F1: %10.6f  F2: %10.6f  fabs: %10.6f %10.6f\n",Qzz,F1,F2,fabs(F1),fabs(F2));
		*/
		if (fabs(F1) < epsilon && fabs(F2) < epsilon) break;
		J11 = (Qpz - 2*Qzz + Qmz)/(del*del);
		J22 = (Qzp - 2*Qzz + Qzm)/(del*del);
		J12 = (Qpp - Qzp - Qpz + Qzz)/(del*del);
		J21 = (Qpp - Qpz - Qzp + Qzz)/(del*del);
		/*
		printf("J11: %10.6f\n",J11);
		printf("J12: %10.6f\n",J12);
		printf("J21: %10.6f\n",J21);
		printf("J22: %10.6f\n",J22);
		*/
		dt = (J12*F2 - J22*F1)/(J11*J22 - J12*J21);
		dp = (-F1 - J11*dt)/J12;
		theta += dt;
		phi += dp;
		L = d/(2*sin(theta)*cos(phi));
		if (k == 20) {
			printf("Too many iterations\n");
			exit(1);
		}
		if (L > 0.75*d || (fabs(dt) < epsilon && fabs(dp) < epsilon)) break;
	}

	getP(d,theta,phi,p0,Vx,Vy,Vz,p);

	kp = np;
	point[kp].x = p[0];
	point[kp].y = p[1];
	point[kp].z = p[2];
	point[kp].d = (point[kp0].d + point[kp1].d)/2;
	point[kp].used = true;
	np++;
	edgeList[ne].npts = 3;
	edgeList[ne].pt = (int *)malloc(edgeList[ne].npts*sizeof(int));
	edgeList[ne].pt_used = (int *)malloc(edgeList[ne].npts*sizeof(int));
	edgeList[ne].npts_used = 3;
	edgeList[ne].pt[0] = kp0;
	edgeList[ne].pt[1] = kp;
	edgeList[ne].pt[2] = kp1;
	len = dist(kp0,kp) + dist(kp,kp1);
	if (len > 300) {
		printf("Too long! %d %d %d  %f\n",kp0,kp,kp1,len);
		printf("P0, P, P1\n");
		for (i=0; i<3; i++) {
			printf("%6.3f  %6.3f  %6.3f\n",p0[i],p[i],p1[i]);
		}
		exit(1);
	}
	edgeList[ne].pt_used[0] = kp0;
	edgeList[ne].pt_used[1] = kp;
	edgeList[ne].pt_used[2] = kp1;
	edgeList[ne].vert[0] = kp0;
	edgeList[ne].vert[1] = kp1;
	edgeList[ne].used = true;
	ne++;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Try to join loose ends.  There are several criteria making up the joining score for a loose end pair:
// (1) closeness - close ends get a high score
// (2) angle - small angle made a by straight joining line at each end gets a high score
// (3) diameter match - closely matched diameters get a high score (NOT USED)
//
// Currently a single segment is used to join the two vertices.  Need to use intermediate points
// to create a curve.
//-----------------------------------------------------------------------------------------------------
int joinLooseEnds(LOOSEND end[], int nloose)
{
	int i0, i1, k, kp0, kp1, njoined, imax;
	double *dir0, *dir1;
	double dave0, dave1, d, vec[3], ang0, ang1, s, smax;
	EDGE edge0, edge1;
#define MAX_END_SEPARATION 100.
#define THRESHOLD_SCORE 1./MAX_END_SEPARATION
#define MAX_DIAMETER 14

	njoined = 0;
	for (i0=0; i0<nloose; i0++) {
		if (end[i0].joined >= 0) continue;
		edge0 = edgeList[end[i0].iedge];
		kp0 = edge0.vert[end[i0].iend];
		dir0 = end[i0].dir;
		dave0 = end[i0].dave;
		if (dave0 > MAX_DIAMETER) continue;
		imax = -1;
		smax = 0;
		for (i1=i0; i1<nloose; i1++) {
			if (i0 == i1) continue;
			if (end[i1].joined >= 0) continue;

			edge1 = edgeList[end[i1].iedge];
			kp1 = edge1.vert[end[i1].iend];
			dir1 = end[i1].dir;
			dave1 = end[i1].dave;
			if (dave1 > MAX_DIAMETER) continue;
			d = dist(kp0,kp1);
			if (d > MAX_END_SEPARATION) continue;
			vec[0] = point[kp1].x - point[kp0].x;
			vec[1] = point[kp1].y - point[kp0].y;
			vec[2] = point[kp1].z - point[kp0].z;
			ang0 = angle(dir0,vec);
			if (ang0 < -1.001*PI || ang0 > 1.001*PI) {
				printf("Bad ang0: %f\n",ang0);
				return 1;
			}
			for (k=0; k<3; k++)
				vec[k] = -vec[k];
			ang1 = angle(dir1,vec);
			if (ang1 < -1.001*PI || ang1 > 1.001*PI) {
				printf("Bad ang1: %f\n",ang1);
				return 1;
			}
			s = fang(ang0)*fang(ang1)*fdist(d);
			if (s < THRESHOLD_SCORE) continue;
			if (s > smax) {
				imax = i1;
				smax = s;
			}
		}
		if (imax >= 0) {
			njoined += 2;
			end[i0].joined = imax;
			end[imax].joined = i0;
//			printf("Joined: %6d  %6d %6d  %6d %6d  %8.4f\n",njoined,i0,imax,kp0,kp1,smax);
			edge1 = edgeList[end[imax].iedge];
			kp1 = edge1.vert[end[imax].iend];
			dir1 = end[imax].dir;

			joining_edge(dir0,dir1,kp0,kp1);
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// To remove obviously unwanted dead-end twigs.
// A twig is the segment from a terminus to a vertex (junction).
// The first criterion is that short fat twigs should be removed, since being fat such a twig
// is very unlikely to have lost its connection because of insufficient label intensity, and
// is much more likely to be an artifact of the thinning algorithm.
// One problem is that we do not want to remove the end of a main artery or vein, which
// may branch soon after the terminus.
//
// A secondary criterion, for thin vessels, is that there is no apparent matching twig nearby,
// i.e. another thin twig within a feasible distance, and possibly with feasible orientation.
//
// I don't believe this is finding all the loose ends.
//-----------------------------------------------------------------------------------------------------
int pruner(int iter)
{
	int i, j, k, ii, kv0, kv1, n0, n1, j1, j2, k1, k2, nloose, npruned, kp0, err;
	float len, dave, dsum;
	EDGE edge, edge0;
	LOOSEND end[20000];

	printf("pruner\n");
	nloose = 0;
	npruned = 0;
	// Step through edges, looking for loose ends.
	for (i=0; i<ne; i++) {
		if (!edgeList[i].used) continue;
		edge = edgeList[i];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
//		bool hit = false;
//		if (kv0 == 33239 || kv1 == 33239) {
//			hit = true;
//			printf("Hit point 33239 on edge: %d\n",i);
//		}
		n0 = 0;
		n1 = 0;
		for (ii=0; ii<ne; ii++) {
			if (!edgeList[ii].used) continue;
			if (i == ii) continue;
			if (edgeList[ii].vert[0] == kv0 || edgeList[ii].vert[1] == kv0) n0=1;
			if (edgeList[ii].vert[0] == kv1 || edgeList[ii].vert[1] == kv1) n1=1;
		}
		if (n0+n1 == 0) {
			printf("Error: pruner: edge: %d is unconnected\n",i);
			fprintf(fperr,"Error: pruner: edge: %d is unconnected\n",i);
			edgeList[i].used = false;
//			return 1;
		} else if (n0+n1 == 1) {	// This is a loose end
			// Need the length and average diameter
			len = 0;
			dsum = 0;
			for (k=0; k<edge.npts; k++) {
				j = edge.pt[k];
				dsum += point[j].d;
				if (k > 0) {
					len += dist(j,edge.pt[k-1]);
				}
			}
			dave = dsum/edge.npts;
			end[nloose].iedge = i;
			if (n0 == 0) {
				end[nloose].iend = 0;
				j2 = 0;
				j1 = 1;
			} else {
				end[nloose].iend = 1;
				j2 = edge.npts-1;
				j1 = j2 - 1;
			}
			end[nloose].dave = dave;
			end[nloose].len = len;
			k1 = edge.pt[j1];
			k2 = edge.pt[j2];
			len = dist(k1,k2);
			end[nloose].dir[0] = (point[k2].x - point[k1].x)/len;
			end[nloose].dir[1] = (point[k2].y - point[k1].y)/len;
			end[nloose].dir[2] = (point[k2].z - point[k1].z)/len;
			end[nloose].joined = -1;
			edge0 = edgeList[end[nloose].iedge];
			kp0 = edge0.vert[end[nloose].iend];
//			fprintf(fpout,"%6d  %6d %2d %2d %2d %6d %6d %6.3f %6.3f %6.3f  %6d\n",nloose,i,end[nloose].iend,j1,j2,k1,k2,
//				end[nloose].dir[0],end[nloose].dir[1],end[nloose].dir[2],kp0);
			nloose++;
		}
	}

	if (iter == 0 && JOIN_LOOSE_ENDS) {
		err = joinLooseEnds(end,nloose);
		if (err != 0) return err;
	}

	npruned = 0;
	for (i=0; i<nloose; i++) {
		if (end[i].joined < 0) {
			if (end[i].len/end[i].dave < ratio_limit) {	// prune this twig 
				edgeList[end[i].iedge].used = false;
				npruned++;
			}
		}
	}
	printf("Number of loose ends: %d, number pruned: %d\n",nloose,npruned);
	fprintf(fperr,"Number of loose ends: %d, number pruned: %d\n",nloose,npruned);
	return 0;
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
// A loop is broken if its edges are still used.
// Each vertex must be checked to see if it has external connections.  There are three cases:
// (1) One vertex is connected externally (this is a dead-end).
// (2) Two vertices are connected externally (side loop)
// (3) Three vertices are connected externally (3-way junction loop)
//-----------------------------------------------------------------------------------------------------
int fixloop(int j1, int j2, int j3)
{
	int e[3];
	int vert[3];
	bool used[3];
	int nconn[3];
	int i, j, kv0, kv1, ncmin, nc2, imax;
	double d, dmax;

	e[0] = j1;
	e[1] = j2;
	e[2] = j3;
	for (i=0; i<3; i++) {
		used[i] = edgeList[e[i]].used;
		nconn[i] = 0;
	}
	vert[0] = edgeList[e[0]].vert[0];
	vert[1] = edgeList[e[0]].vert[1];
	if (edgeList[e[1]].vert[0] == vert[0] || edgeList[e[1]].vert[0] == vert[1]) {
		vert[2] = edgeList[e[1]].vert[1];
	} else {
		vert[2] = edgeList[e[1]].vert[0];
	}
	for (i=0; i<ne; i++) {
		kv0 = edgeList[i].vert[0];
		kv1 = edgeList[i].vert[1];
		for (j=0; j<3; j++) {
			if (vert[j] == kv0 || vert[j] == kv1) {
				nconn[j] += 1;
			}
		}
	}
//	printf("fixloop: %6d %6d %6d  used: %2d %2d %2d  nconn: %2d %2d %2d\n",
//		j1,j2,j3,used[0],used[1],used[2],nconn[0],nconn[1],nconn[2]);

	nc2 = 0;
	ncmin = 999;
	for (i=0; i<3; i++) {
		if (nconn[i] < ncmin) ncmin = nconn[i];
		if (nconn[i] == 2) nc2++;
	}
	if (nc2 > 0) {
		for (i=0; i<3; i++) {
			kv0 = edgeList[e[i]].vert[0];
			kv1 = edgeList[e[i]].vert[1];
			for (j=0; j<3; j++) {
				if (nconn[j] == 2 && (vert[j] == kv0 || vert[j] == kv1)) {
					edgeList[e[i]].used = false;
					break;
				}
			}
		}
	} else if (nc2 == 0) {	// DO NOT arbitrarily choose edge j1, choose longest
		dmax = 0;
		for (i=0; i<3; i++) {
			kv0 = edgeList[e[i]].vert[0];
			kv1 = edgeList[e[i]].vert[1];
			d = dist(kv0,kv1);
			if (d > dmax) {
				imax = i;
				dmax = d;
			}
		}
		edgeList[e[imax]].used = false;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Search for loops
// deloop works fine.  Should remove longest side when all three have nconn = 3
//-----------------------------------------------------------------------------------------------------
int deloop(void)
{
#define NE2MAX 100000
	int i, ii, kv0, kv1, kkv0, kkv1;
	int npairs, ne2, j1, j2, j3, nloops, k1, k3;
	EDGE edge, eedge;
	PAIR pair[NE2MAX];
	int e2[NE2MAX];
	bool dup;

	printf("deloop\n");

// Does any non-vertex point occur in more than one edge?  No.
	
	ne2 = 0;
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		if (!edge.used) continue;
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		if (edge.npts == 2) {
			if (ne2 == NE2MAX) {
				printf("Array dimension e2 exceeded\n");
				return 1;
			}
//			printf("2-edge: %4d %6d  %6d %6d\n",ne2,i,kv0,kv1);
			e2[ne2] = i;
			ne2++;
			dup = false;
			for (ii=i+1; ii<ne; ii++) {
				eedge = edgeList[ii];
				if (!eedge.used) continue;
				kkv0 = eedge.vert[0];
				kkv1 = eedge.vert[1];
				// This finds a number of edges
				if ((kkv0 == kv0 && kkv1 == kv1) || (kkv0 == kv1 && kkv1 == kv0)) {
//					printf("Duplicate edges: %6d  %6d  %4d %4d\n",i,ii,edge.npts,eedge.npts);
//					fprintf(fpout,"%6d  %6d  %4d %4d\n",i,ii,edge.npts,eedge.npts);
					dup = true;
					if (eedge.npts == 2) {
						point[kkv0].d = 1.414*point[kkv0].d;
						point[kkv1].d = 1.414*point[kkv1].d;
					}
					break;
				}
			}
			if (dup) {
				ne2--;
				edgeList[i].used = false;
			}
		}
	}
	printf("Number of 2-edges (edges with two points): %d\n",ne2);

	npairs = 0;
	for (j1=0; j1<ne2; j1++) {
		edge = edgeList[e2[j1]];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		for (j2=j1+1; j2<ne2; j2++) {
			eedge = edgeList[e2[j2]];
			if (eedge.vert[0] == kv1) {			// -->  -->
				pair[npairs].i1 = e2[j1];
				pair[npairs].i2 = e2[j2];
				npairs++;
			} else if (eedge.vert[1] == kv1) {	// -->  <--
				pair[npairs].i1 = e2[j1];
				pair[npairs].i2 = -e2[j2];
				npairs++;
			} else if (eedge.vert[0] == kv0) {	// <--  -->
				pair[npairs].i1 = -e2[j1];
				pair[npairs].i2 = e2[j2];
				npairs++;
			} else if (eedge.vert[1] == kv0) {	// <--  <--
				pair[npairs].i1 = -e2[j1];
				pair[npairs].i2 = -e2[j2];
				npairs++;
			}
		}
	}
	printf("Number of 2-edge pairs (connections): %d\n",npairs);
	
	bool hit;
	nloops = 0;
	for (i=0; i<npairs; i++) {
		for (ii=i+1; ii<npairs; ii++) {
			hit = false;
			if (abs(pair[i].i2) == abs(pair[ii].i1)) {
				j1 = pair[i].i1;
				j2 = abs(pair[i].i2);
				j3 = pair[ii].i2;
				hit = true;
			} else if (abs(pair[i].i1) == abs(pair[ii].i2)) {
				j1 = pair[i].i2;
				j2 = abs(pair[i].i1);
				j3 = pair[ii].i1;
				hit = true;
			} else if (abs(pair[i].i2) == abs(pair[ii].i2)) {
				j1 = pair[i].i1;
				j2 = abs(pair[i].i2);
				j3 = pair[ii].i1;
				hit = true;
			} else if (abs(pair[i].i1) == abs(pair[ii].i1)) {
				j1 = pair[i].i2;
				j2 = abs(pair[i].i1);
				j3 = pair[ii].i2;
				hit = true;
			}
			if (hit) {
				if (j1 > 0)
					k1 = edgeList[j1].vert[0];
				else
					k1 = edgeList[-j1].vert[1];
				if (j3 > 0)
					k3 = edgeList[j3].vert[1];
				else
					k3 = edgeList[-j3].vert[0];
				if (k1 == k3) {
					nloops++;
					fixloop(abs(j1),abs(j2),abs(j3));
				}
			}
		}
	}
	printf("Number of loops: %d\n",nloops);
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
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int n_prune_cycles, err;
	char *input_amfile;
	char drive[32], dir[1024],filename[256], ext[32];
	char errfilename[1024], output_amfile[1024], outfilename[1024], result_file[1024];
	int prune_flag, cmgui_flag;
	double ddiam, dlen;

	if (argc != 9) {
		printf("Usage: prune input_amfile output_amfile ratio_limit n_prune_cycles prune_flag cmgui_flag delta_diam delta_len\n");
		fperr = fopen("prune_error.log","w");
		fprintf(fperr,"Usage: prune input_amfile output_amfile ratio_limit n_prune_cycles prune_flag cmgui_flag delta_diam delta_len\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(outfilename,argv[2]);
	sscanf(argv[3],"%lf",&ratio_limit);
	sscanf(argv[4],"%d",&n_prune_cycles);
	sscanf(argv[5],"%d",&prune_flag);
	sscanf(argv[6],"%d",&cmgui_flag);
	sscanf(argv[7],"%lf",&ddiam);
	sscanf(argv[8],"%lf",&dlen);
	if (prune_flag == 0) n_prune_cycles = 0;
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_prune.log",output_basename);
	sprintf(output_amfile,"%s.am",output_basename);
	sprintf(result_file,"%s.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	fpout = fopen(result_file,"w");	
	err = ReadAmiraFile(input_amfile);
	if (err != 0) return 2;
	if (n_prune_cycles > 0) {
		err = adjoinEdges();
		if (err != 0) return 3;
		err = deloop();
		if (err != 0) return 4;
		err = adjoinEdges();
		if (err != 0) return 3;
		for (int k=0; k<n_prune_cycles; k++) {
			err = pruner(k);
			if (err != 0) return 5;
			err = adjoinEdges();
			if (err != 0) return 3;
			err = checkEdgeEndPts();
			if (err != 0) return 6;
		}
		err = squeezer();	// must squeeze, or SpatialGraph and CMGUI files are not consistent
		if (err != 0) return 7;
	}

	err = WriteAmiraFile(output_amfile,input_amfile);
	if (err != 0) return 8;
//	err = oldCreateDistributions();
	err = CreateDistributions(ddiam, dlen);
	if (err != 0) return 9;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename);
		if (err != 0) return 10;
	}
	return 0;
}