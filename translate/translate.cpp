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

#include "translate.h"

int WriteCmguiData(char *basename);
float dist(int k1, int k2);

int nv, ne, np;
int nv_used, ne_used, np_used;
EDGE *edgeList;
VERTEX *vertex;
APOINT *point;

FILE *fperr, *fpout;
FILE *exelem, *exnode;
char output_basename[128];
double ratio_limit;

#define NBOX 400
#define STR_LEN 128
//#define RATIO_LIMIT 4	// should be an input parameter
#define PI 3.14159

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
// The average diameter of a vessel (edge) is now estimated by dividing the volume by the length.
// All distances in the .am file are in um.
//-----------------------------------------------------------------------------------------------------
int CreateDistributions()
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
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
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
	point = (APOINT *)malloc(npp*sizeof(APOINT));		// 2* for added joining edges

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
//			fclose(fpam);
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
				err = 1;
			}
			// Check for same end points with different indices
			int j1 = edge.pt[0];
			int j2 = edge.pt[1];
			if (point[j1].x == point[j2].x && point[j1].y == point[j2].y && point[j1].z == point[j2].z) {
				printf("edge: %d vertices: %d %d same pos: %f %f %f\n",i,j1,j2,point[j1].x,point[j1].y,point[j1].z);
				fprintf(fperr,"edge: %d vertices: %d %d same pos: %f %f %f\n",i,j1,j2,point[j1].x,point[j1].y,point[j1].z);
				err = 1;
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
	int ne_x, nv_x, np_x, knew, i_x;
	EDGE edge;
	EDGE *edgeList_x;
	VERTEX *vertex_x;
	APOINT *point_x;
	int *oldpt;

	printf("squeezer\n");
	vertex_x = (VERTEX *)malloc(nv*sizeof(VERTEX));
	edgeList_x = (EDGE *)malloc(2*ne*sizeof(EDGE));		// 2* for added joining edges
	point_x = (APOINT *)malloc(2*np*sizeof(APOINT));
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
		if (ne_x == 1737) {
			printf("new edge 1737 was: %d\n",i);
			fprintf(fperr,"new edge 1737 was: %d\n",i);
		}
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
int checkUnconnected(void)
{
	int i, kv0, kv1, n0, n1, ii;
	EDGE edge;

	// Step through edges, looking for unconnected edges.
	for (i=0; i<ne; i++) {
		if (!edgeList[i].used) continue;
		edge = edgeList[i];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		n0 = 0;
		n1 = 0;
		for (ii=0; ii<ne; ii++) {
			if (!edgeList[ii].used) continue;
			if (i == ii) continue;
			if (edgeList[ii].vert[0] == kv0 || edgeList[ii].vert[1] == kv0) n0=1;
			if (edgeList[ii].vert[0] == kv1 || edgeList[ii].vert[1] == kv1) n1=1;
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
	char errfilename[256], output_amfile[256], outfilename[256], result_file[256];
	int cmgui_flag;


	if (argc != 4) {
		printf("Usage: translate input_amfile output_file cmgui_flag\n");
		fperr = fopen("translate_error.log","w");
		fprintf(fperr,"Usage: translate input_amfile output_file cmgui_flag\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(outfilename,argv[2]);
	strcpy(output_amfile,"zzz.am");
	sscanf(argv[3],"%d",&cmgui_flag);
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_translate.log",output_basename);
//	sprintf(result_file,"%s_prune.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

//	fpout = fopen(result_file,"w");	
	fpout = fopen(outfilename,"w");	
	err = ReadAmiraFile(input_amfile);
	if (err != 0) return 2;

//	err = squeezer();	// must squeeze, or SpatialGraph and CMGUI files are not consistent
//	if (err != 0) return 6;

	err = WriteAmiraFile(output_amfile,input_amfile);
	if (err != 0) return 7;
	err = CreateDistributions();
	if (err != 0) return 8;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename);
		if (err != 0) return 9;
	}
	return 0;
}