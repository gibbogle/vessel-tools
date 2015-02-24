#include <cstdio>
#include <math.h>
#include "network.h"

extern FILE *fperr, *fpout;
extern bool use_len_limit, use_len_diam_limit;
extern float len_limit, len_diam_limit;
extern float ddiam, dlen;

//-----------------------------------------------------------------------------------------------------
// Compute vessel lengths and average diameters.
//-----------------------------------------------------------------------------------------------------
int EdgeDimensions(EDGE *edges, POINT *points, int ne)
{
	int ie, ip, kp, kprev;
	float dx, dy, dz, len, deltalen, deltalen2, ad, r2, r2prev, dsum, lsum, vol, diam;
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
				deltalen2 = dx*dx+dy*dy+dz*dz;
//				printf("dx,dy,dz,dlen2: %6.1f %6.1f %6.1f %6.1f\n",dx,dy,dz,dlen2);
				if (deltalen2 == 0) {
					printf("EdgeDimensions: error: dlen = 0: ie,npts,ip,kp: %d %d %d %d point: %6.1f %6.1f %6.1f\n",
						 ie,edge.npts,ip,kp,points[kp].x,points[kp].y,points[kp].z);
					fprintf(fpout,"EdgeDimensions: error: dlen = 0: ie,npts,ip,kp: %d %d %d %d point: %6.1f %6.1f %6.1f\n",
						 ie,edge.npts,ip,kp,points[kp].x,points[kp].y,points[kp].z);
//					return 2;
				}
				deltalen = sqrt(deltalen2);
				vol += PI*deltalen*(r2 + r2prev)/2;
				lsum += deltalen;
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
	double ad, len, dave, ltot, dsum, lsegdtot;
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
	if (use_len_diam_limit) {
		printf("\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
		fprintf(fpout,"\nUsing length/diameter lower limit = %6.1f\n\n",len_diam_limit);
	} else if (use_len_limit) {
		printf("\nUsing length lower limit = %6.1f um\n\n",len_limit);
		fprintf(fpout,"\nUsing length lower limit = %6.1f um\n\n",len_limit);
	}
	printf("Compute diameter distributions (length weighted)\n");
	fprintf(fperr,"Compute diameter distributions (length weighted)\n");
	// Diameters
//	ddiam = 0.5;
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
		if (use_len_limit && len < len_limit) continue;
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
		if (use_len_diam_limit && len/ad < len_diam_limit) continue;
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
//	dlen = 1;
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/dlen + 0.5);
		if (use_len_limit && k*dlen <= len_limit) continue;
		ad = edge.segavediam;
		if (use_len_diam_limit && len/ad < len_diam_limit) continue;
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