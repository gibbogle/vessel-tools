// Locate very short edges and remove them by merging the two end junction nodes.
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
int ne_net[NEMAX];

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

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
bool EqualPoints(POINT p1, POINT p2)
{
	if (p1.x != p2.x) return false;
	if (p1.y != p2.y) return false;
	if (p1.z != p2.z) return false;
	return true;
}

//-----------------------------------------------------------------------------------
// mode = 'N' for numbers
//        'D' for diameters
//        'P' for positions
//-----------------------------------------------------------------------------------
void showedge(NETWORK *net, int ie, char mode)
{
	int k, kp;

	printf("edge: %d npts: %d\n",ie,net->edgeList[ie].npts);
	for (k=0; k<net->edgeList[ie].npts; k++) {
		kp = net->edgeList[ie].pt[k];
		if (mode == 'N')
			printf("%7d ",kp);
		else if (mode == 'D')
			printf("%7.1f ",net->point[kp].d);
		else 
			printf("%7.1f %7.1f %7.1f   ",net->point[kp].x,net->point[kp].y,net->point[kp].z);
		printf("\n");
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
// Find vertex nodes with only two edges connected
//-----------------------------------------------------------------------------------------------------
int CheckNetwork1(NETWORK *net, char *str)
{
	int iv, i, ne, nv, n2;

	fprintf(fpout,"CheckNetwork: %s\n",str);
	nv = net->nv;
	ne = net->ne;
	for (iv=0; iv<nv; iv++) {
		net->vertex[iv].nlinks = 0;
	}
	// Set up net edges and count vertex edges
	for (i=0; i<ne; i++) {
		for (int k=0; k<2;k++) {
			iv = net->edgeList[i].vert[k];
			net->vertex[iv].nlinks++;
		}
	}
	n2 = 0;
	for (iv=0; iv<nv; iv++) {
		if (net->vertex[iv].nlinks == 2) {
			n2++;
//			fprintf(fpout,"vertex: %6d nlinks: %2d\n",iv,net->vertex[iv].nlinks);
		}
	}
	fprintf(fpout,"Number of 2-link vertices: %d\n",n2);
	return 0;
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
				fprintf(fpout,"CheckNetwork: repeated (0): ie,npts,v0,ip,kp0,kp: %d %d %d %d %d %d\n",ie,npts,0,ip,kp0,kp);
				fprintf(fpout,"removed\n");
//				fprintf(fpout,"netID: %d\n",edge.netID);
				repeated++;
				net->edgeList[ie].used = false;
			}
		}
		for (ip=0; ip<npts-1; ip++) {
			kp = edge.pt[ip];
			p = net->point[kp];
			if (EqualPoints(p,p1)) {
				fprintf(fpout,"CheckNetwork: repeated (1): ie,npts,ip,v1,kp,kp1: %d %d %d %d %d %d\n",ie,npts,ip,npts-1,kp,kp1);
				fprintf(fpout,"removed\n");
//				fprintf(fpout,"netID: %d\n",edge.netID);
				repeated++;
				net->edgeList[ie].used = false;
			}
		}
	}
//	if (repeated > 0) return repeated;
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
				} else {
					net->edgeList[ie1].used = false;
				}
			}
		}
	}
	return 0;
}	

/*
//-----------------------------------------------------------------------------------------------------
// Compute vessel lengths and average diameters.
//-----------------------------------------------------------------------------------------------------
int EdgeDimensions(EDGE *edges, POINT *points, int ne)
{
	int ie, ip, kp, kprev;
	float dx, dy, dz, len, dlen, dlen2, ad, r2, r2prev, dsum, lsum, vol, diam;
	EDGE edge;

	printf("EdgeDimensions:\n");
	fprintf(fpout,"EdgeDimensions:\n");
	for (ie=0; ie<ne; ie++) {
		edge = edges[ie];
		if (!edge.used) continue;
		kprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		vol = 0;
//		printf("\nedge, npts: %4d %4d\n",ie,edge.npts);
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = points[kp].d;
			r2 = ad*ad/4;
//			printf("ip,kp: %4d %4d x,y,z,d: %6.1f %6.1f %6.1f %6.1f\n",ip,kp,points[kp].x,points[kp].y,points[kp].z,points[kp].d);
			if (ad < 0.001 || ad > 100) {
				printf("Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fperr,"Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				return 1;
			}
			if (ip > 0) {
				dx = points[kprev].x - points[kp].x;
				dy = points[kprev].y - points[kp].y;
				dz = points[kprev].z - points[kp].z;
				dlen2 = dx*dx+dy*dy+dz*dz;
//				printf("dx,dy,dz,dlen2: %6.1f %6.1f %6.1f %6.1f\n",dx,dy,dz,dlen2);
				if (dlen2 == 0) {
					printf("EdgeDimensions: error: dlen = 0: ie,npts,ip,kp: %d %d %d %d point: %6.1f %6.1f %6.1f\n",
						 ie,edge.npts,ip,kp,points[kp].x,points[kp].y,points[kp].z);
					fprintf(fpout,"EdgeDimensions: error: dlen = 0: ie,npts,ip,kp: %d %d %d %d point: %6.1f %6.1f %6.1f\n",
						 ie,edge.npts,ip,kp,points[kp].x,points[kp].y,points[kp].z);
					return 2;
				}
				dlen = sqrt(dlen2);
				vol += PI*dlen*(r2 + r2prev)/2;
				lsum += dlen;
			}
			kprev = kp;
			r2prev = r2;
		}
		len = lsum;
		diam = 2*sqrt(vol/(PI*len));
		edges[ie].length_um = len;	
		edges[ie].segavediam = diam;	
//		printf("edge: %4d len,diam: %6.1f %6.1f\n",ie,len,diam);
	}
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
	double ad, len, dave, ddiam, dlen, ltot, dsum, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam;		// average point diameter
	double ave_seg_diam;	// average vessel diameter
	double ave_lseg_diam;	// length-weighted average vessel diameter
	int ie, ip, k, ka, kp, ndpts, nlpts, ndtot, nsegdtot;
	EDGE edge;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lsegadbox[k] = 0;
		lvbox[k] = 0;
	}
	if (use_len_diam) {
		printf("\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
		fprintf(fpout,"\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
	} else if (use_len) {
		printf("\nUsing length lower limit = %6.1f um\n\n",len_limit);
		fprintf(fpout,"\nUsing length lower limit = %6.1f um\n\n",len_limit);
	}
	printf("Compute diameter distributions (length weighted)\n");
	fprintf(fperr,"Compute diameter distributions (length weighted)\n");
	// Diameters
	ddiam = 0.5;
	ndtot = 0;
	nsegdtot = 0;
	lsegdtot = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	ave_lseg_diam = 0;
	volume = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		dave = edge.segavediam;
		if (use_len && len < len_limit) continue;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = net->point[kp].d;
			ave_pt_diam += ad;
			if (ad < 0.001 || ad > 200) {
				printf("Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fperr,"Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
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
		}
		net->edgeList[ie].length_um = len;	// already computed
		if (use_len_diam && len/ad < len_diam_limit) continue;
		ave_seg_diam += dave;
		ave_lseg_diam += dave*len;
//		printf("ie: %6d dave,len,ave_lseg_diam: %8.1f %8.1f %10.1f\n",ie,dave,len,ave_lseg_diam);
//		fprintf(fpout,"ie: %6d dave,len,ave_lseg_diam: %8.1f %8.1f %10.1f\n",ie,dave,len,ave_lseg_diam);
		if (dave < 0.001 || dave > 200) {
			printf("Bad segment diameter: edge: %d ad: %f\n",ie,dave);
			fprintf(fperr,"Zero segment diameter: edge: %d ad: %f\n",ie,dave);
			return 1;
		}
		ka = int(dave/ddiam + 0.5);
		if (ka >= NBOX) {
			printf("Vessel too wide (segment ave): d: %f k: %d\n",dave,ka);
			fprintf(fperr,"Vessel too wide (segment ave): d: %f k: %d\n",dave,ka);
			continue;
		}
		segadbox[ka]++;
		nsegdtot++;
		lsegadbox[ka] += len;
		lsegdtot += len;
		volume += len*PI*(dave/2)*(dave/2);
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
	fprintf(fperr,"\nCompute length distributions:\n");
	// Lengths
	dlen = 1;
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/dlen + 0.5);
		if (use_len && k*dlen <= len_limit) continue;
		ad = edge.segavediam;
		if (use_len_diam && len/ad < len_diam_limit) continue;
		if (k >= NBOX) {
			printf("Edge too long for boxes: len: %d  %6.1f  k: %d\n",ie,len,k);
			fprintf(fperr,"Edge too long for boxes: len: %d  %6.1f  k: %d\n",ie,len,k);
			continue;
		}
		lvbox[k]++;
		ave_len += len;
		ltot++;
	}
	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	ave_lseg_diam /= lsegdtot;
	fprintf(fpout,"Total vertices: %d  points: %d\n",net->nv,net->np);
	fprintf(fpout,"Vessels: %d\n",net->ne);
	printf("\nAverage pt diameter: %6.2f vessel diameter: %6.2f length-weighted: %6.2f\n",
		ave_pt_diam, ave_seg_diam, ave_lseg_diam);
	fprintf(fpout,"\nAverage pt diameter: %6.2f vessel diameter: %6.2f length-weighted vessel diameter: %6.2f\n",
		ave_pt_diam, ave_seg_diam, ave_lseg_diam);
	printf("Average vessel length: %6.1f\n",ave_len/ltot);
	fprintf(fpout,"Average vessel length: %6.1f\n",ave_len/ltot);
	printf("Volume: %10.0f um3\n\n",volume);
	fprintf(fpout,"Volume: %10.0f um3\n\n",volume);

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
*/

/*
//-----------------------------------------------------------------------------------------------------
// The average diameter of a vessel (edge) is now estimated by dividing the volume by the length.
// All distances in the .am file are in um.
//-----------------------------------------------------------------------------------------------------
int CreateDistributions1(NETWORK *net)
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double lsegadbox[NBOX];
	double ad, len, ddiam, dlen, ltot, lsum, dsum, dvol, r2, r2prev, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam;		// average point diameter
	double ave_seg_diam;	// average vessel diameter
	double ave_lseg_diam;	// length-weighted average vessel diameter
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot;
	EDGE edge;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lsegadbox[k] = 0;
		lvbox[k] = 0;
	}
	if (use_len_diam) {
		printf("\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
		fprintf(fpout,"\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
	} else {
		printf("\nUsing length lower limit = %6.1f um\n\n",len_limit);
		fprintf(fpout,"\nUsing length lower limit = %6.1f um\n\n",len_limit);
	}
	printf("Compute diameter distributions (length weighted)\n");
	fprintf(fperr,"Compute diameter distributions (length weighted)\n");
	// Diameters
	ddiam = 0.5;
	ndtot = 0;
	nsegdtot = 0;
	lsegdtot = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	volume = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
//		printf("ie: %d npts: %d\n",ie,edge.npts);
//		fprintf(fperr,"ie: %d npts: %d\n",ie,edge.npts);
//		fflush(fperr);
		bool dbug = false;
		len = edge.length_um;
		if (!use_len_diam && len < len_limit) continue;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		dvol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = net->point[kp].d;
			ave_pt_diam += ad;
			if (dbug) {
				printf("%d  %d  %f  %f\n",ip,kp,ad,ddiam);
				fprintf(fperr,"%d  %d  %f  %f\n",ip,kp,ad,ddiam);
			}
			fflush(fperr);
			if (ad < 0.001 || ad > 200) {
				printf("Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fperr,"Bad point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
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
				dlen = dist(net,kp,kpprev);
				r2 = ad*ad/4;
				dvol += PI*dlen*(r2 + r2prev)/2;
				lsum += dlen;
			}
			kpprev = kp;
			r2prev = r2;
		}
		len = lsum;
//		net->edgeList[ie].length_um = len;	// already computed
		volume += dvol;
		if (dbug) {
			printf("lsum: %f\n",len);
			fprintf(fperr,"lsum: %f\n",len);
			fflush(fperr);
			if (lsum == 0) return 1;
		}
		ad = 2*sqrt(dvol/(PI*len));	// segment diameter
		net->edgeList[ie].segavediam = ad;
		if (use_len_diam && len/ad < len_diam_limit) continue;
//		printf("ie: %d dvol,lsum,ad: %8.0f %6.1f %6.1f\n",ie,dvol,lsum,ad);
		fprintf(fpout,"ie: %d dvol,lsum,ad,total: %8.0f %6.1f %6.1f %8.1f\n",ie,dvol,lsum,ad,ave_seg_diam);
		if (dvol <= 0) {
			printf("Bad diameter: %8.2f\n",ad);
			return 1;
		}
		ave_seg_diam += ad;
		ave_lseg_diam += ad*len;
		if (ad < 0.001 || ad > 200) {
			printf("Bad segment diameter: edge: %d ad: %f\n",ie,ad);
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
		lsegadbox[ka] += len;
		lsegdtot += len;
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
	fprintf(fperr,"\nCompute length distributions:\n");
	// Lengths
	dlen = 1;
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/dlen + 0.5);
		if (!use_len_diam && k*dlen <= len_limit) continue;
		ad = edge.segavediam;
		if (use_len_diam && len/ad < len_diam_limit) continue;
		if (k >= NBOX) {
			printf("Edge too long for boxes: len: %d  %6.1f  k: %d\n",ie,len,k);
			fprintf(fperr,"Edge too long for boxes: len: %d  %6.1f  k: %d\n",ie,len,k);
			continue;
		}
		lvbox[k]++;
		ave_len += len;
		ltot++;
	}
	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	fprintf(fpout,"ave_seg_diam,nsegdtot: %12.3e %8d\n",ave_seg_diam,nsegdtot);
	ave_lseg_diam /= lsegdtot;
	fprintf(fpout,"Total vertices: %d  points: %d\n",net->nv,net->np);
	fprintf(fpout,"Vessels: %d\n",net->ne);
	printf("\nAverage pt diameter: %6.2f vessel diameter: %6.2f length-weighted: %6.2f\n",
		ave_pt_diam, ave_seg_diam, ave_lseg_diam);
	fprintf(fpout,"\nAverage pt diameter: %6.2f vessel diameter: %6.2f length-weighted vessel diameter: %6.2f\n",
		ave_pt_diam, ave_seg_diam, ave_lseg_diam);
	printf("Volume: %10.0f\n\n",volume); 
	printf("Average vessel length: %6.1f\n",ave_len/ltot);
	fprintf(fpout,"Average vessel length: %6.1f\n",ave_len/ltot);
	printf("Volume: %10.0f\n\n",volume);
	printf("Volume: %10.0f\n\n",volume);
	printf("Volume: %10.0f\n\n",volume);
	printf("Volume: %10.0f\n\n",volume);
	printf("Volume: %10.0f\n\n",volume);
	fprintf(fpout,"Volume: %10.0f\n\n",volume);

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
*/

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int CreateReducedNet(NETWORK *net0, NETWORK *net1)
{
	int i, j, ie, iv, kp, kv, ne, nv, npts;
	float d, dave;
	EDGE edge;
	EDGE *elist;
	int *ivlist;

	// Create a list of used edges, and a list of their vertex indices
	elist = (EDGE *)malloc(net0->ne*sizeof(EDGE));
	ivlist = (int *)malloc(2*net0->ne*sizeof(int));
	nv = 0;
	ne = 0;
	for (i=0; i<net0->ne; i++) {
		edge = net0->edgeList[i];
		if (edge.used) {
			elist[ne] = edge;
			for (j=0; j<2; j++) {
				iv = edge.vert[j];
				if (!net0->vertex[iv].used) {
					printf("Error: CreateReducedNet: vertex not used: %d %d %d\n",i,j,iv);
					return 1;
				}
				kv = ivlistAdd(iv,ivlist,&nv);
				elist[ne].vert[j] = kv;
			}
			ne++;
		}
	}
	// Now from this list of edges we need to create the new network net1.
	// Note that the vertex indices in ivlist[] are for the net0 vertex list.
	net1->vertex = (VERTEX *)malloc(nv*sizeof(VERTEX));
	net1->edgeList = (EDGE *)malloc(ne*sizeof(EDGE));
	net1->point = (POINT *)malloc(net0->np*sizeof(POINT));	// keeping full original list of points
	// Set up net1 vertices
	for (iv=0; iv<nv; iv++) {
		net1->vertex[iv] = net0->vertex[ivlist[iv]];
		net1->vertex[iv].nlinks = 0;
	}
	// Set up net1 edges and count vertex edges
	for (i=0; i<ne; i++) {
		net1->edgeList[i] = elist[i];
		for (int k=0; k<2;k++) {
			iv = net1->edgeList[i].vert[k];
			net1->vertex[iv].nlinks++;
		}
	}
	// Set up vertex edge lists
	for (iv=0; iv<nv; iv++) {
		net1->vertex[iv].edge = (int *)malloc(net1->vertex[iv].nlinks*sizeof(int));
		net1->vertex[iv].nlinks = 0;
	}
	for (i=0; i<ne; i++) {
		for (int k=0; k<2;k++) {
			iv = net1->edgeList[i].vert[k];
			net1->vertex[iv].edge[net1->vertex[iv].nlinks] = i;
			net1->vertex[iv].nlinks++;
		}
	}

	// Copy net0 points to net1, initially set all points unused
	for (i=0; i<net0->np; i++) {
		net1->point[i] = net0->point[i];
		net1->point[i].used = false;
	}
	// Now flag those points that are used.
	for (ie=0; ie<ne; ie++) {
		dave = 0;
		npts = net1->edgeList[ie].npts;
		for (kp=0; kp<npts; kp++) {
			i = net1->edgeList[ie].pt[kp];
			d = net1->point[i].d;
			if (d < 0.001 || d > 200) {
				printf("Error: CreateReducedNet: ie,kp,i: %d %d %dd = 0\n",ie,kp,i);
				return 1;
			}
			dave += d;
			net1->point[i].used = true;
		}
//		fprintf(fpout,"edge diameter: %d %d %6.1f\n",ie,npts,dave/npts);
	}

	net1->ne = ne;
	net1->nv = nv;
	net1->np = net0->np;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Find vertex nodes with only two edges connected
//-----------------------------------------------------------------------------------------------------
int CompleteNetwork(NETWORK *net)
{
	int ie, iv, i, j, kp0, kp1, kp2, ne, nv, n2;
	int npts, iv0, iv1;
	float len, dsum;
	EDGE edge;
	POINT pv, pn;

	fprintf(fpout,"CompleteNetwork\n");
	nv = net->nv;
	ne = net->ne;
	for (iv=0; iv<nv; iv++) {
		net->vertex[iv].nlinks = 0;
	}
	// Just for safety, search for redundant edges, i.e. iv0 -> iv1 and iv1 -> iv0
	// Mark as unused any redundant edge
	for (ie=0; ie<ne; ie++) {
		iv0 = net->edgeList[ie].vert[0];
		iv1 = net->edgeList[ie].vert[1];
		for (i=0; i<ne; i++) {
			if (i == ie) continue;
			if (iv1 == net->edgeList[i].vert[0] && iv0 == net->edgeList[ie].vert[1]) {
				printf("redundant edge: %d  vert: %d %d\n",ie,iv0,iv1);
				fprintf(fpout,"redundant edge: %d  vert: %d %d\n",ie,iv0,iv1);
				net->edgeList[i].used = false;
			}
		}
	}
	// Set up net edges and count vertex edges
	for (ie=0; ie<ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		for (int k=0; k<2;k++) {
			iv = edge.vert[k];
			net->vertex[iv].nlinks++;
		}
		npts = edge.npts;
		if (npts == 0) {
			printf("Error: npts = 0: edge: %d\n",ie);
			return 1;
		}
		iv0 = edge.vert[0];
		iv1 = edge.vert[1];
		kp0 = edge.pt[0];
		kp1 = edge.pt[npts-1];
		// This is to check that for an edge, vert[0] corresponds to pt[0], vert[1] to pt[npts-1]
		pv = net->vertex[iv0].point;
		pn = net->point[kp0];
		if (!EqualPoints(pv,pn)) {
			fprintf(fpout,"Not the same points (0): ie: %d vert: %d %d kp: %d %d\n",ie,iv0,iv1,kp0,kp1);
		}
		pv = net->vertex[iv1].point;
		pn = net->point[kp1];
		if (!EqualPoints(pv,pn)) {
			fprintf(fpout,"Not the same points (1): ie: %d vert: %d %d kp: %d %d\n",ie,iv0,iv1,kp0,kp1);
		}
	}
	n2 = 0;
	for (iv=0; iv<nv; iv++) {
		if (net->vertex[iv].nlinks == 2) {
			n2++;
//			fprintf(fpout,"vertex: %6d nlinks: %2d\n",iv,net->vertex[iv].nlinks);
		}
	}
	fprintf(fpout,"Number of 2-link vertices: %d\n",n2);
	
	// Set up vertex edge lists
	for (iv=0; iv<nv; iv++) {
//		net->vertex[iv].edge = (int *)malloc(net->vertex[iv].nlinks*sizeof(int));
		net->vertex[iv].edge = (int *)malloc(20*sizeof(int));	// 20 should be enough
		net->vertex[iv].nlinks = 0;
		
	}
	for (i=0; i<ne; i++) {
		if (!net->edgeList[i].used) continue;
		for (int k=0; k<2;k++) {
			iv = net->edgeList[i].vert[k];
			net->vertex[iv].edge[net->vertex[iv].nlinks] = i;
			net->vertex[iv].nlinks++;
		}
	}
	for (i=0; i<ne; i++) {
		if (!net->edgeList[i].used) continue;
		npts = net->edgeList[i].npts;
		iv = net->edgeList[i].vert[0];
		if (net->vertex[iv].nlinks > 3) {
			pv = net->vertex[iv].point;
//			fprintf(fpout,"Vertex: %6d pos: %6.1f %6.1f %6.1f  nlinks: %4d\n",iv,pv.x,pv.y,pv.z,net->vertex[iv].nlinks);
		}
		len = 0;
		dsum = 0;
		for (j=0; j<npts; j++) {
			kp1 = net->edgeList[i].pt[j];
			if (j < npts-1) {
				kp2 = net->edgeList[i].pt[j+1];
				len += dist(net,kp1,kp2);
			}
			dsum += net->point[kp1].d;
			if ( net->point[kp1].d < 0.001 || net->point[kp1].d > 200) {
				printf("Error: CompleteNetwork: bad d: %f\n",net->point[kp1].d);
				return 1;
			}
		}
		net->edgeList[i].length_um = len;
		net->edgeList[i].segavediam = dsum/npts;
		if (i == 0) printf("length of edge 0: %f\n",len);
		if (len == 0) {
			printf("Error: CompleteNetwork: len = 0: edge: %d\n",ie);
			return 1;
		}
//		fprintf(fpout,"edge: %6d npts: %3d len: %6.1f diam: %6.1f\n",i,npts,len,dsum/npts);
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// The reduction of vertices iv0 and iv1 to a single vertex is not valid if the connected edges have
// a common vertex.
//-----------------------------------------------------------------------------------------------------
void CheckReduction(NETWORK *net, int ie, bool *ok)
{
	int iv0, iv1, j0, j1, ie0, ie1, nlinks0, nlinks1, ivend0, ivend1;
	EDGE edge0, edge1;

	iv0 =  net->edgeList[ie].vert[0];
	iv1 =  net->edgeList[ie].vert[1];
	nlinks0 = net->vertex[iv0].nlinks;
	nlinks1 = net->vertex[iv1].nlinks;
	for (j0=0; j0<nlinks0; j0++) {
		ie0 = net->vertex[iv0].edge[j0];
		if (ie0 == ie) continue;
		edge0 = net->edgeList[ie0];
		if (edge0.vert[0] == iv0) {
			ivend0 = edge0.vert[1];
		} else {
			ivend0 = edge0.vert[0];
		}
		for (j1=0; j1<nlinks1; j1++) {
			ie1 = net->vertex[iv1].edge[j1];
			if (ie1 == ie) continue;
			edge1 = net->edgeList[ie1];
			if (edge1.vert[0] == iv1) {
				ivend1 = edge1.vert[1];
			} else {
				ivend1 = edge1.vert[0];
			}
			if (ivend0 == ivend1) {
//				printf("edges: %d %d  vert: %d %d  %d %d\n",ie0,ie1,edge0.vert[0],edge0.vert[1],edge1.vert[0],edge1.vert[1]);
				*ok = false;
				return;
			}
		}
	}
	*ok = true;
}
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int ReduceNetwork(NETWORK *NP0, NETWORK *NP1, float len_diam_limit)
{
	int ie, ie0, j, kp0, kp1, kp, ndropped, err;
	int npts, iv0, iv1, nlinks0, nlinks1, npts1;
	int jj, elist[20];
	float len, diam;
	bool ok;
	EDGE edge0, edge;
	POINT pmid, p0, p1, p;

	err = CompleteNetwork(NP0);	// ensure that edge lengths and diameters and vertex edge lists are computed
	if (err != 0) return err;

	fprintf(fpout,"ReduceNetwork\n");

	ndropped = 0;
	for (ie0=0; ie0<NP0->ne; ie0++) {
		edge0 = NP0->edgeList[ie0];
		if (!edge0.used) {
			printf("edge not used: %d\n",ie0);
			continue;
		}
		diam = edge0.segavediam;
		len = edge0.length_um;
		if (len == 0) {
			printf("Error: ReduceNetwork: len = 0: edge: %d\n",ie0);
			return 1;
		}
		if (len/diam > len_diam_limit) continue;
		npts = edge0.npts;
		iv0 = edge0.vert[0];
		if (!NP0->vertex[iv0].used) {
			printf("vertex iv0 not used: %d\n",iv0);
			fprintf(fpout,"vertex iv0 not used: %d\n",iv0);
			return 1;
		}
		iv1 = edge0.vert[1];
		if (!NP0->vertex[iv1].used) {
			printf("vertex iv1 not used: %d\n",iv1);
			fprintf(fpout,"vertex iv1 not used: %d\n",iv1);
			return 1;
		}

		CheckReduction(NP0,ie0,&ok);
		if (!ok) {
			continue;
		}

		nlinks0 = NP0->vertex[iv0].nlinks;
		nlinks1 = NP0->vertex[iv1].nlinks;
		if (nlinks1 == 1 || nlinks1 > 3 || nlinks0 > 3) continue;

		kp0 = edge0.pt[0];
		kp1 = edge0.pt[npts-1];
		pmid.x = (NP0->vertex[iv0].point.x + NP0->vertex[iv1].point.x)/2;
		pmid.y = (NP0->vertex[iv0].point.y + NP0->vertex[iv1].point.y)/2;
		pmid.z = (NP0->vertex[iv0].point.z + NP0->vertex[iv1].point.z)/2;
		pmid.d = (NP0->point[kp0].d + NP0->point[kp1].d)/2;
		NP0->vertex[iv0].point = pmid;
		NP0->point[kp0] = pmid;
		// THIS IS NOT ENOUGH!!!  Need to move the end points of all connected edges.
		// This is Amira stupidity - a separate point for every edge connected at the same vertex.

		printf("Dropping edge: %d vertex: %d %d  kp0,kp1: %d %d\n",ie0,iv0,iv1,kp0,kp1);
//		printf("npts: %d len: %6.1f diam: %6.1f len/diam: %6.3f\n",npts,len,diam,len/diam);
//		fprintf(fpout,"Dropping edge: %d vertex: %d %d  kp0,kp1: %d %d\n",ie0,iv0,iv1,kp0,kp1);
//		fprintf(fpout,"npts: %d len: %6.1f diam: %6.1f len/diam: %6.3f\n",npts,len,diam,len/diam);
//		fprintf(fpout,"pmid: %6.1f %6.1f %6.1f\n",pmid.x,pmid.y,pmid.z);
		// Need to drop the ie0 edge from the edge list for iv0
		jj = 0;
		for (j=0; j<nlinks0; j++) {
			ie = NP0->vertex[iv0].edge[j];
			if (ie == ie0) continue;
			elist[jj] = ie;
			jj++;
		}
		nlinks0 = jj;
		NP0->vertex[iv0].nlinks = nlinks0;
		for (j=0; j<nlinks0; j++) {
			NP0->vertex[iv0].edge[j] = elist[j];
			// Can move the end points here:
			ie = NP0->vertex[iv0].edge[j];
			edge = NP0->edgeList[ie];
			if (edge.vert[0] == iv0) {
				kp = edge.pt[0];
				NP0->point[kp] = pmid;
			} else {
				kp = edge.pt[edge.npts-1];
				NP0->point[kp] = pmid;
			}
		}

		for (j=0; j<nlinks1; j++) {
			ie = NP0->vertex[iv1].edge[j];
			if (ie == ie0 || !NP0->edgeList[ie].used) continue;
			NP0->vertex[iv0].edge[NP0->vertex[iv0].nlinks] = ie;
			NP0->vertex[iv0].nlinks++;
			if (NP0->edgeList[ie].vert[0] == iv1) {
				NP0->edgeList[ie].vert[0] = iv0;
				NP0->edgeList[ie].pt[0] = kp0;
//				fprintf(fpout,"connected edge: %d at vert[0], pt: %d -> kp0: %d\n",ie,0,kp0);
				kp = NP0->edgeList[ie].pt[1];
//				printf("v0: ie0,ie,kp: %d %d %d\n",ie0,ie,kp);
//				fprintf(fpout,"v0: ie0,ie,kp: %d %d %d\n",ie0,ie,kp);
				p0 = NP0->point[kp0];
				p = NP0->point[kp];
				if (EqualPoints(p0,p)) {
					printf("Error: same points (0): ie0,ie,kp0: %d %d %d\n",ie0,ie,kp0);
					fprintf(fpout,"Error: same points (0): ie0,ie,kp0: %d %d %d\n",ie0,ie,kp0);
					err = 2;
				}
			} else {
				npts1 = NP0->edgeList[ie].npts;
				NP0->edgeList[ie].vert[1] = iv0;
				NP0->edgeList[ie].pt[npts1-1] = kp0;
//				fprintf(fpout,"connected edge: %d at vert[1], pt: %d -> kp0: %d\n",ie,npts1-1,kp0);
				kp = NP0->edgeList[ie].pt[npts1-2];
//				printf("v1: ie0,ie,kp: %d %d %d\n",ie0,ie,kp);
//				fprintf(fpout,"v1: ie0,ie,kp: %d %d %d\n",ie0,ie,kp);
				p = NP0->point[kp];
				p1 = NP0->point[kp0];
				if (EqualPoints(p,p1)) {
					printf("Error: same points (1): ie0,ie,kp0: %d %d %d\n",ie0,ie,kp0);
					fprintf(fpout,"Error: same points (1): ie0,ie,kp0: %d %d %d\n",ie0,ie,kp0);
					err = 2;
				}
			}
		}
		NP0->vertex[iv1].used = false;
		NP0->edgeList[ie0].used = false;
		ndropped++;
	}
	printf("\nNumber of edges reduced: %d\n",ndropped);
	fprintf(fpout,"\nNumber of edges reduced: %d\n",ndropped);
	err = CreateReducedNet(NP0,NP1);
	return err;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	char *input_amfile;
	char drive[32], dir[2048],filename[256], ext[32];
	char errfilename[2048], output_amfile[2048], result_file[2048];
	char output_basename[2048];
	int cmgui_flag;
	float origin_shift[3];
	NETWORK *NP0, *NP1;

	if (argc != 6) {
		printf("Usage: reduce input_amfile output_amfile len_diam_limit len_limit cmgui_flag\n");
		fperr = fopen("reduce_error.log","w");
		fprintf(fperr,"Usage: reduce input_amfile output_amfile len_diam_limit len_limit cmgui_flag\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(output_amfile,argv[2]);
	sscanf(argv[3],"%f",&len_diam_limit);
	sscanf(argv[4],"%f",&len_limit);
	sscanf(argv[5],"%d",&cmgui_flag);

	_splitpath(output_amfile,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_reduce.log",output_basename);
	sprintf(result_file,"%s_reduce.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	fpout = fopen(result_file,"w");	
	printf("Input .am file: %s\n",input_amfile);
	printf("Minimum vessel length/diameter: %6.1f\n",len_diam_limit);
	printf("Minimum vessel length (for distributions): %6.1f\n",len_limit);
	printf("cmgui_flag: %d\n",cmgui_flag);
	fprintf(fpout,"Input .am file: %s\n",input_amfile);
	fprintf(fpout,"Minimum vessel length/diameter: %6.3f\n",len_diam_limit);
	fprintf(fpout,"Minimum vessel length (for distributions): %6.1f\n",len_limit);
	fprintf(fpout,"cmgui_flag: %d\n",cmgui_flag);

	if (len_limit < 0) {		// this signals that length/diameter limit should be used for the distributions
		use_len_diam_limit = true;
		use_len_limit = false;
	} else {					// use this as the length limit for the distributions
		use_len_diam_limit = false;
		use_len_limit = true;
	}
	ddiam = 0.5;
	dlen = 1.0;

	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;

	err = CheckNetwork(NP0,"Original network");
	if (err != 0) return 3;

	NP1 = (NETWORK *)malloc(sizeof(NETWORK));

	err = ReduceNetwork(NP0,NP1,len_diam_limit);
	if (err != 0) return 4;

	printf("\nNP1: ne, nv, np: %d %d %d\n\n",NP1->ne,NP1->nv,NP1->np);
	err = CheckNetwork(NP1,"Reduced network");
	if (err != 0) return 5;

	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;

	err = WriteAmiraFile(output_amfile,input_amfile,NP1,origin_shift);
	if (err != 0) return 6;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,NP1,origin_shift);
		if (err != 0) return 7;
	}
	err = EdgeDimensions(NP1->edgeList,NP1->point,NP1->ne);
	if (err != 0) return 8;
	err = CreateDistributions(NP1);
	if (err != 0) return 9;
	return 0;
}