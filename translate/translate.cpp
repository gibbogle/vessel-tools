// trans.cpp
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

//#include "translate.h"
#include "network.h"

//int WriteCmguiData(char *basename);
//float dist(int k1, int k2);

//int nv, ne, np;
//int nv_used, ne_used, np_used;
//EDGE *edgeList;
//VERTEX *vertex;
//POINT *point;

FILE *fperr, *fpout;
FILE *exelem, *exnode;
char output_basename[128];

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

/*
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
		if (!edge.used || edge.npts == 2) continue;
//		printf("ie: %d npts: %d\n",ie,edge.npts);
//		fprintf(fperr,"ie: %d npts: %d\n",ie,edge.npts);
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
		if (!edge.used || edge.npts == 2) continue;
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
	printf("Total vertices: %d  points: %d\n",nv,np);
	fprintf(fpout,"Total vertices: %d  points: %d\n",nv,np);
	printf("Vessels: %d  used: %d\n",ne,int(ltot));
	fprintf(fpout,"Vessels: %d  used: %d\n",ne,int(ltot));
	printf("Average pt diameter: %6.2f vessel diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	fprintf(fpout,"Average pt diameter: %6.2f vessel diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	printf("Average vessel length: %6.1f\n",ave_len/ltot);
	fprintf(fpout,"Average vessel length: %6.1f\n",ave_len/ltot);
	printf("Total vessel length: %10.0f\n",ave_len);
	fprintf(fpout,"Total vessel length: %10.0f\n",ave_len);
	printf("Volume: %10.0f\n\n",volume);
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
int oldCreateDistributions()
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double ad, len_vox, ddiam, dlen, ltot, lsum, dsum, dvol, r2, r2prev;
	double ave_len_vox, volume, d95;
	double ave_pt_diam, ave_seg_diam;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot;
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
	ndtot = 0;
	nsegdtot = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	volume = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
//		printf("ie: %d npts: %d\n",ie,edge.npts);
//		fprintf(fperr,"ie: %d npts: %d\n",ie,edge.npts);
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
			ka = int(ad/ddiam);
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
		ka = int(ad/ddiam);
		if (ka >= NBOX) {
			printf("Vessel too wide (segment ave): d: %f k: %d\n",ad,ka);
			fprintf(fperr,"Vessel too wide (segment ave): d: %f k: %d\n",ad,ka);
			continue;
		}
		segadbox[ka]++;
		nsegdtot++;
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
	ave_len_vox = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		len_vox = edge.length_um;
		k = int(len_vox/dlen + 0.5);
		if (k*dlen <= lenlimit) continue; 
		if (k >= NBOX) {
			printf("Edge too long: len_vox: %d  %f  k: %d\n",ie,len_vox,k);
			fprintf(fperr,"Edge too long: len_vox: %d  %f  k: %d\n",ie,len_vox,k);
			continue;
		}
		lvbox[k]++;
		ave_len_vox += len_vox;
		ltot++;
	}
	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	fprintf(fpout,"Total vertices: %d  points: %d\n",nv,np);
	fprintf(fpout,"Segments: %d\n",ne);
	printf("Average pt diameter: %6.2f segment diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	fprintf(fpout,"Average pt diameter: %6.2f segment diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	printf("Average length: %6.1f\n",ave_len_vox/ltot);
	fprintf(fpout,"Average length: %6.1f\n",ave_len_vox/ltot);
	fprintf(fpout,"Volume: %10.0f\n\n",volume);

	ndpts = 0;
	for (k=0; k<NBOX; k++) {
		if (adbox[k]/float(ndtot) >= 0.0005) {
			ndpts = k+1;
		}
	}
	fprintf(fpout,"   um   ave_diam  seg_ave\n");
	for (k=0; k<ndpts; k++) {
		fprintf(fpout,"%6.2f %9.4f %9.4f\n",k*ddiam,adbox[k]/float(ndtot),segadbox[k]/float(nsegdtot));
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
	}
	return 0;
}
*/


/*
//-----------------------------------------------------------------------------------------------------
// Write Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile1(char *amFileOut, char *amFileIn)
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
	fprintf(fpam,"# Created by translater.exe from: %s\n",amFileIn);
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
*/

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
		printf("%s\n",line);
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
//					printf("edge: %d  npts: %d\n",i,edge.npts);
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

/*
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile1(char *amFile)
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

	vertex = (VERTEX *)malloc(2*nv*sizeof(VERTEX));
	edgeList = (EDGE *)malloc(nee*sizeof(EDGE));	// 2* for added joining edges
	point = (POINT *)malloc(npp*sizeof(POINT));	// 2* for added joining edges

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
					vertex[i].point.d = 0;
					vertex[i].point_index = -1;
				}
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
				}
				printf("Got edge npts, total: %d\n",npts);
			} else if (k == 4) {
				kp = 0;
				for (i=0;i<ne;i++) {
					edge = edgeList[i];
					float len = 0;
//					printf("edge: %d  npts: %d\n",i,edge.npts);
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
						sscanf(line,"%f %f %f",&x,&y,&z);
						if (k == 0 || k == edge.npts-1) {
							get_vertex_ptindex(x,y,z,kp,&kpv);
							point[kpv].x = x;
							point[kpv].y = y;
							point[kpv].z = z;
							edgeList[i].pt[k] = kpv;
							edgeList[i].pt_used[k] = kpv;
						} else {
							point[kp].x = x;
							point[kp].y = y;
							point[kp].z = z;
							edgeList[i].pt[k] = kp;
							edgeList[i].pt_used[k] = kp;
						}
//						printf("point: %d  %f %f %f\n",kp,x,y,z);
						kp++;
						if (k > 0) {
							len = len + dist(edgeList[i].pt[k-1],edgeList[i].pt[k]);
						}
					}
					edgeList[i].length_um = len;
				}
				printf("got point positions\n");
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
						if (point[j].d == 0) {
							printf("Error: ReadAmiraFile: zero diameter: i: %d npts: %d k: %d j: %d\n",i,edge.npts,k,j);
							return 1;
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
	int err = 0;
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			point[j].used = true;
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
			if (point[j1].x == point[j2].x && point[j1].y == point[j2].y && point[j1].z == point[j2].z) {
				printf("edge: %d vertices: %d %d same pos: %f %f %f\n",i,j1,j2,point[j1].x,point[j1].y,point[j1].z);
				fprintf(fperr,"edge: %d vertices: %d %d same pos: %f %f %f\n",i,j1,j2,point[j1].x,point[j1].y,point[j1].z);
//				err = 1;
			}
		}
	}
	fclose(fpam);
	return err;
}

//-----------------------------------------------------------------------------------------------------
// Reorder element and node identifiers
//-----------------------------------------------------------------------------------------------------
int squeezer(void)
{
	int i, j, k;
	int ne_x, nv_x, np_x, knew, i_x, kv;
	EDGE edge;
	EDGE *edgeList_x;
	VERTEX *vertex_x;
	POINT *point_x;
	int *oldpt;
	POINT p;

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
			if (knew == -1) {
				oldpt[np_x] = edgeList[i].vert[j];
				vertex_x[np_x].point = point[oldpt[np_x]];
				point_x[np_x] = point[oldpt[np_x]];
				knew = np_x;
//				printf("Vertex: %6d %6d %3d %6d\n",nv_x, i,j,edgeList[i].vert[j]);
//				fprintf(fpout,"Vertex: %6d %6d %3d %6d\n",nv_x, i,j,edgeList[i].vert[j]);
				np_x++;
			}
			edge.vert[j] = knew;
		}
		edgeList_x[ne_x] = edge;
		ne_x++;
	}
	nv_x = np_x;

	//fprintf(fpout, "Original vertices: %d\n",nv);
	//for (kv=0; kv<nv; kv++) {
	//	p = point[kv];
	//	fprintf(fpout,"%6d %6d %8.1f %8.1f %8.1f %8.2f\n",kv,oldpt[kv],p.x,p.y,p.z,p.d);
	//}

	// Now add the edge interior points to the list.  These occur once only.
	i_x = 0;
	for (i=0; i<ne; i++) {
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
	printf("Did squeezer\n");
	return 0;
}
*/

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
	char errfilename[256], output_amfile[256], outfilename[256];
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
//	strcpy(output_amfile,"zzz.am");
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
	if (err != 0) return 7;
//	err = CreateDistributions(ddiam,dlen);
//	if (err != 0) return 8;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,net,origin_shift);
		if (err != 0) return 9;
	}
	err = EdgeDimensions(net->edgeList,net->point,net->ne);
	if (err != 0) return 8;
	err = CreateDistributions(net);
	if (err != 0) return 9;
	return 0;
}