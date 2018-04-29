/*
 * To extract topology from a skeletonized .tif  
 */

#include <cstdio>
#include <vector>
/*
#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
*/
#if (defined (LINUX) || defined (__linux__))
#include <libgen.h>

#endif
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

#include "network2.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//#define MAXNBRS 30
//#define MAXEDGES 1000000
#define PI 3.14159
#define NBOX 800

bool FRC_fix = false;
bool METHOD2 = false;

struct xyz_str
{
	double x, y, z;
};
typedef xyz_str XYZ;

#define V3Dskel(a,b,c) pskel[(c)*imsize_xy+(b)*width+(a)]
#define V3D(a,b,c)     p[(c)*imsize_xy+(b)*width+(a)]
#define Vindex(a,b,c)  pindex[(c)*imsize_xy+(b)*width+(a)]
typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im, imskel;
long long width, height, depth, imsize_xy;
unsigned char *pskel, *p;
int *pindex;
VOXEL *Vlist;
EDGE *edgeList;
NETWORK *net;
VERTEX *vertexList;
float *avediameter;
int nlit, ne, np, nv, ne_used;
int lendist[100];
int diamdist[100];
float ddiam = 0.5;
bool fixed_diam_flag = false;
bool use_object;
float FIXED_DIAMETER = 0;
float max_diameter;
bool use_max_diameter;
float vsize[3];	// in um
float min_end_len = 5;
float len_limit = 2;
int nzerodiameters;
double calib_param = 0;

bool use_pt_diameters = true;	// this is how it was done in topo

FILE *fperr, *fpout;

int WriteAmiraFile(char *outFile, char *vessFile, char *skelFile, float *origin_shift);
int getPointDiameters();
void freeEdgeList();

//-----------------------------------------------------------------------------------------------------
// For Linux to create output_basename without extension
//-----------------------------------------------------------------------------------------------------
char *remove(char* mystr) {
    char *retstr;
    char *lastdot;
    if (mystr == NULL)
         return NULL;
    if ((retstr = (char *)malloc (strlen (mystr) + 1)) == NULL)
        return NULL;
    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, '.');
    if (lastdot != NULL)
        *lastdot = '\0';
    return retstr;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
//float zdist(int k1, int k2)
//{
//	float dx = vsize[0]*(voxel[k2].pos[0] - voxel[k1].pos[0]);
//	float dy = vsize[0]*(voxel[k2].pos[1] - voxel[k1].pos[1]);
//	float dz = vsize[0]*(voxel[k2].pos[2] - voxel[k1].pos[2]);
//	return sqrt(dx*dx+dy*dy+dz*dz);
//}

//-----------------------------------------------------------------------------------------------------
// Distance between voxels in um (because vsize is voxel size)
//-----------------------------------------------------------------------------------------------------
float dist_um(int p1[], int p2[])
{
	float d2=0;
	d2 += vsize[0]*vsize[0]*(p1[0]-p2[0])*(p1[0]-p2[0]);
	d2 += vsize[1]*vsize[1]*(p1[1]-p2[1])*(p1[1]-p2[1]);
	d2 += vsize[2]*vsize[2]*(p1[2]-p2[2])*(p1[2]-p2[2]);
	return sqrt(d2);
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
// Rotate v0 about N through angle to get v
//-----------------------------------------------------------------------------------------------------
int Rotate(double *v0, double *N, double angle, double *v)
{
	double cosa, sina, d;
	XYZ q1, q2, u;

	cosa = cos(angle);
	sina = sin(angle);

	q1.x = v0[0]; q1.y = v0[1]; q1.z = v0[2];
	u.x = N[0]; u.y = N[1]; u.z = N[2];

   d = sqrt(u.y*u.y + u.z*u.z);

   /* Step 2 */
   if (d != 0) {
      q2.x = q1.x;
      q2.y = q1.y * u.z / d - q1.z * u.y / d;
      q2.z = q1.y * u.y / d + q1.z * u.z / d;
   } else {
      q2 = q1;
   }

   /* Step 3 */
   q1.x = q2.x * d - q2.z * u.x;
   q1.y = q2.y;
   q1.z = q2.x * u.x + q2.z * d;

   /* Step 4 */
   q2.x = q1.x * cosa - q1.y * sina;
   q2.y = q1.x * sina + q1.y * cosa;
   q2.z = q1.z;

   /* Inverse of step 3 */
   q1.x =   q2.x * d + q2.z * u.x;
   q1.y =   q2.y;
   q1.z = - q2.x * u.x + q2.z * d;

   /* Inverse of step 2 */
   if (d != 0) {
      q2.x =   q1.x;
      q2.y =   q1.y * u.z / d + q1.z * u.y / d;
      q2.z = - q1.y * u.y / d + q1.z * u.z / d;
   } else {
      q2 = q1;
   }

   /* Inverse of step 1 */
   v[0] = q2.x;
   v[1] = q2.y;
   v[2] = q2.z;
   return 0;
}

//-----------------------------------------------------------------------------------------------------
// The average diameter of a vessel (edge) is now estimated by dividing the volume by the length.
// The scale value is no longer required since the distance calculations now yield values in um 
// - to account for dx=dy!=dz 
//-----------------------------------------------------------------------------------------------------
int CreateDistributions_topo()
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double lsegadbox[NBOX];
	double ad, len, ddiam, dlen, ltot, lsum, dsum, dvol, vol, r2, r2prev, lsegdtot;
	double tot_len, volume, d95;
	double ave_pt_diam, ave_seg_diam;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot, nptstot, nptsusedtot;
	EDGE edge;
	bool dbug;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lsegadbox[k] = 0;
		lvbox[k] = 0;
	}
	printf("Compute diameter distributions\n");
	fprintf(fperr,"Compute diameter distributions\n");
	fprintf(fpout,"Compute diameter distributions\n");
	fflush(fpout);
	nptstot = 0;
	nptsusedtot = 0;
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
//		fprintf(fperr,"ie: %d npts: %d npts_used: %d\n",ie,edge.npts,edge.npts_used);
		fflush(fperr);
		nptstot += edge.npts;
//		nptsusedtot += edge.npts_used;
		dbug = false;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		vol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = avediameter[kp];
			ave_pt_diam += ad;
			if (dbug) {
				printf("%d  %d  %f  %f\n",ip,kp,ad,ddiam);
				fprintf(fperr,"%d  %d  %f  %f\n",ip,kp,ad,ddiam);
			}
			fflush(fperr);
			if (ad < 0.001) {
				printf("Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				fprintf(fperr,"Zero point diameter: edge: %d point: %d ad: %f\n",ie,ip,ad);
				return 1;
			}
			ka = int(ad/ddiam + 0.5);
			if (ka >= NBOX) {
				printf("Vessel too wide (point): d: %f ka: %d\n",ad,ka);
				fprintf(fperr,"Vessel too wide (point): d: %f ie,ip,kp,ka: %d %d %d %d  npts: %d\n",ad,ie,ip,kp,ka,edge.npts);
				fflush(fperr);
				ka = NBOX-1;
				ad  = ddiam*(NBOX-1);
//				continue;
			}
			adbox[ka]++;
			ndtot++;
			r2 = ad*ad/4;
			if (ip > 0) {
//				dlen = dist_um(kp,kpprev);
				dlen = dist_um(Vlist[kp].initial_pos,Vlist[kpprev].initial_pos);
				dvol = PI*dlen*(r2 + r2prev)/2;
				vol += dvol;
				lsum += dlen;
//				fprintf(fperr,"ie: %d ip: %d kp: %d ad: %6.2f dlen: %6.2f r2: %6.2f r2prev: %6.2f dvol: %6.2f\n",
//					ie,ip,kp,ad,dlen,r2,r2prev,dvol);
			}
			kpprev = kp;
			r2prev = r2;
		}
		edgeList[ie].length_um = lsum;
		volume += vol;
		if (dbug) {
			printf("lsum: %f\n",lsum);
			fprintf(fperr,"lsum: %f\n",lsum);
			fflush(fperr);
			if (lsum == 0) return 1;
		}
		ad = 2*sqrt(vol/(PI*lsum));	// segment diameter
		edgeList[ie].segavediam = ad;
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
//		fprintf(fperr,"edge: %4d len,diam,vol,ka: %6.1f %6.1f %6.1f %d\n",ie,lsum,ad,vol,ka);
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
	printf("Compute length distributions: lower limit = %6.1f um\n",len_limit);
	fprintf(fperr,"Compute length distributions: lower limit = %6.1f um\n",len_limit);
	fflush(fperr);
	fprintf(fpout,"Compute length distributions: lower limit = %6.1f um\n",len_limit);
	fflush(fpout);
	// Lengths
	dlen = 1;
	ltot = 0;
	tot_len = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/dlen + 0.5);	// was k = int(len/dlen)
		if (len <= len_limit) continue;
		if (k >= NBOX) {
			printf("Edge too long: ie: %d  k: %d len: %f\n",ie,k,len);
			fprintf(fperr,"Edge too long: ie: %d  k: %d len: %f\n",ie,k,len);
			fflush(fperr);
			continue;
		}
		lvbox[k]++;
		tot_len += len;
		ltot++;
	}
	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	fprintf(fpout,"Total vertices: %d  points: %d\n",nv,np);
	fprintf(fpout,"Vessels: %d ltot: %d\n",ne,int(ltot));
	printf("Average pt diameter: %6.2f vessel diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	fprintf(fpout,"Average pt diameter: %6.2f vessel diameter: %6.2f\n",ave_pt_diam, ave_seg_diam);
	printf("Average vessel length: %6.1f\n",tot_len/ltot);
	fprintf(fpout,"Average vessel length: %6.1f\n",tot_len/ltot);
	printf("Total vessel length: %6.1f\n",tot_len);
	fprintf(fpout,"Total vessel length: %6.1f\n",tot_len);
	printf("Total vessel volume: %10.0f\n\n",volume);
	fprintf(fpout,"Total vessel volume: %10.0f\n\n",volume);
	fflush(fpout);
	for (k=NBOX-1; k>=0; k--) {
		if (segadbox[k] > 0) break;
	}
	ndpts = k+2;
	fprintf(fpout,"Vessel diameter distribution\n");
	fprintf(fpout,"   um    number  fraction    length  fraction\n");
	for (k=0; k<ndpts; k++) {
		fprintf(fpout,"%6.2f %8d %9.5f  %8.0f %9.5f\n",k*ddiam,segadbox[k],segadbox[k]/float(nsegdtot),
			lsegadbox[k],lsegadbox[k]/lsegdtot);
		fflush(fpout);
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
	fflush(fpout);
	printf("nptstot: %d\n",nptstot);
	printf("nptsusedtot: %d\n",nptsusedtot);
	fprintf(fpout,"did CreateDistributions_topo\n");
	fflush(fpout);
	return 0;
}


//-----------------------------------------------------------------------------------------------------
// The average diameter of a vessel (edge) is now estimated by dividing the volume by the length.
// The scale value is no longer required since the distance calculations now yield values in um 
// - to account for dx=dy!=dz 
//-----------------------------------------------------------------------------------------------------
int CreateDistributions()
{
	int adbox[NBOX], lvbox[NBOX];
	int segadbox[NBOX];
	double lsegadbox[NBOX];
	double ad, len, ddiam, dlen, lsum, dsum, dvol, vol, r2, r2prev, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam, ave_seg_diam;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot, nptstot, nptsusedtot, ltot;
	EDGE edge;
	bool dbug;

	for (k=0;k<NBOX;k++) {
		adbox[k] = 0;
		segadbox[k] = 0;
		lsegadbox[k] = 0;
		lvbox[k] = 0;
	}
	if (use_object) {
	printf("Compute diameter distributions\n");
	fprintf(fpout,"Compute diameter distributions\n");
	if (use_max_diameter) fprintf(fpout,"\nUSING maximum diameter: %f\n\n",max_diameter);
	nptstot = 0;
	nptsusedtot = 0;
	// Diameters
	ddiam = 0.5;
	ndtot = 0;
	nsegdtot = 0;
	lsegdtot = 0;
	ave_pt_diam = 0;
	ave_seg_diam = 0;
	volume = 0;
	lsum = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		fflush(fperr);
		nptstot += edge.npts;
		/*
//		nptsusedtot += edge.npts_used;
		dbug = false;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		vol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
			ad = avediameter[kp];
			ave_pt_diam += ad;
			if (dbug) {
				printf("%d  %d  %f  %f\n",ip,kp,ad,ddiam);
				fprintf(fperr,"%d  %d  %f  %f\n",ip,kp,ad,ddiam);
			}
			fflush(fperr);
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
			r2 = ad*ad/4;
			if (ip > 0) {
				dlen = dist_um(kp,kpprev);
				dvol = PI*dlen*(r2 + r2prev)/2;
				vol += dvol;
				lsum += dlen;
//				fprintf(fperr,"ie: %d ip: %d kp: %d ad: %6.2f dlen: %6.2f r2: %6.2f r2prev: %6.2f dvol: %6.2f\n",
//					ie,ip,kp,ad,dlen,r2,r2prev,dvol);
			}
			kpprev = kp;
			r2prev = r2;
		}
		edgeList[ie].length_um = lsum;
		volume += vol;
		if (dbug) {
			printf("lsum: %f\n",lsum);
			fprintf(fperr,"lsum: %f\n",lsum);
			fflush(fperr);
			if (lsum == 0) return 1;
		}
		ad = 2*sqrt(vol/(PI*lsum));	// segment diameter
		*/
		ad = edge.segavediam;
		if (!edge.used) continue;
		lsum += edge.length_um;
//		volume += edge.length_um*PI*(ad/2)*(ad/2);
		volume += edge.volume;
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
//		fprintf(fperr,"edge: %4d len,diam,vol,ka: %6.1f %6.1f %6.1f %d\n",ie,lsum,ad,vol,ka);
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
	}

	printf("Compute length distributions: lower limit = %6.1f um\n",len_limit);
	fprintf(fpout,"Compute length distributions: lower limit = %6.1f um\n",len_limit);
	// Lengths
	dlen = 1;
	ltot = 0;
	ave_len = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		ad = edge.segavediam;
		if (!edge.used) continue;
//		if (!edge.used) continue;
		len = edge.length_um;
		k = int(len/dlen + 0.5);	// was k = int(len/dlen)
		if (len <= len_limit) continue;
		if (k >= NBOX) {
			printf("Edge too long: ie: %d  k: %d len: %f  k: %d\n",ie,k,len);
			fprintf(fperr,"Edge too long: ie: %d  k: %d len: %f  k: %d\n",ie,k,len);
			continue;
		}
		lvbox[k]++;
		ave_len += len;
		ltot++;
	}

	if (use_object) {
	ave_pt_diam /= ndtot;
	ave_seg_diam /= nsegdtot;
	printf("Average vessel diameter: %6.2f um\n",ave_seg_diam);
	fprintf(fpout,"Average vessel diameter: %6.2f um\n",ave_seg_diam);
	printf("Total vessel volume: %10.0f um3\n",volume);
	fprintf(fpout,"Total vessel volume: %10.0f um3\n",volume);
	printf("Total vessel length: %10.0f um\n\n",lsum);
	fprintf(fpout,"Total vessel length: %10.0f um\n\n",lsum);
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
	}

//	fprintf(fpout,"Total vertices: %d  points: %d\n",nv,np);
	fprintf(fpout,"\nNumber of vessels: %d  number exceeding min length: %d\n",ne,ltot);
	printf("Average vessel length: %6.2f\n",ave_len/ltot);
	fprintf(fpout,"Average vessel length: %6.2f\n",ave_len/ltot);
	for (k=NBOX-1; k>=0; k--) {
		if (lvbox[k] > 0) break;
	}
	nlpts = k+2;
	fprintf(fpout,"Vessel length distribution\n");
	fprintf(fpout,"   um    number  fraction\n");
	for (k=0; k<nlpts; k++) {
		fprintf(fpout,"%6.2f %8d %9.5f\n",k*dlen,lvbox[k],lvbox[k]/float(ltot));
	}
	printf("nptstot: %d\n",nptstot);
//	printf("nptsusedtot: %d\n",nptsusedtot);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Get distance from p1[] in direction of unit vector v[] to the first black voxel.
// For FRC, try reducing the distance by 1/2 to get a better estimate: fix = true
//-----------------------------------------------------------------------------------------------------
double GetRadius2(double p1[3], double v[3])
{
	int i, k, ixyz[3], ixyzmax[3];
	double r, dr, xyz, vox_um, rfix, ir;
	bool out;

//	printf("GetRadius2\n");
	vox_um = vsize[0];
	dr = vox_um/10;
	ixyzmax[0] = width-1;
	ixyzmax[1] = height-1;
	ixyzmax[2] = depth-1;
	i = 0;
	for(;;) {
		r = i*dr;
//		printf("r: %d %6.1f\n",i,r);
		out = false;
		for (k=0; k<3; k++) {
			xyz = p1[k] + r*v[k];
			ixyz[k] = xyz/vsize[k];
			if (ixyz[k] < 0) out = true;
			if (ixyz[k] > ixyzmax[k]) out = true;
		}
//		printf("ixyz: %d %d %d\n",ixyz[0],ixyz[1],ixyz[2]);
		if (out) break;
		if (V3D(ixyz[0],ixyz[1],ixyz[2]) == 0) {
			if (r == 0) {
//				printf("GetRadius2 (a): r2=0: i,x,y,z,V: %d  %d %d %d  %d\n",i,ixyz[0],ixyz[1],ixyz[2],V3D(ixyz[0],ixyz[1],ixyz[2]));
//				fprintf(fperr,"GetRadius2 (a): r2=0: i,x,y,z,V: %d  %d %d %d  %d\n",i,ixyz[0],ixyz[1],ixyz[2],V3D(ixyz[0],ixyz[1],ixyz[2]));
				return 0;	// was 1.0;
			}
			break;
		}
		i++;
	}
	if (FRC_fix) 
		rfix = MAX(0.0,r-0.5*vox_um);
	else
		rfix = r;
	return rfix*rfix;		// was r*r
}

//-----------------------------------------------------------------------------------------------------
// At the point of interest, p1, the centreline connects p0 and p2.
// We want to estimate the mean square radius of the vessel at p1.
// We do this by shooting lines equispaced through 360 degrees from p1,
// in the plane normal to p0-p2, to find distance to the first black voxel.
// This method should be used on a simplified edge, i.e. with sufficient
// separation between p0 and p2 to provide a reasonable estimate of the 
// centreline direction vector N[].
//-----------------------------------------------------------------------------------------------------
int EstimateDiameter(double p0[3], double p1[3], double p2[3], double *r2ave, double *r2min, bool *zero)
{
	int i;
	double N[3], sum, dp, r2, r2sum;
	double v0[3], v[3];
	double angle;
	int nrays;
	double RMIN = 0.5;	// was 0.5

//	printf("EstimateDiameter: %f %f %f\n",p1[0],p1[1],p1[2]);
	sum = 0;
	for (i=0; i<3; i++) {
		N[i] = p2[i] - p0[i];	// axis of vessel, normal to the plane
		sum += N[i]*N[i];
	}
	sum = sqrt(sum);
	for (i=0; i<3; i++) {
		N[i] = N[i]/sum;		// unit normal
	}
	// N is normal to the plane that approximates the slice through the vessel at p1
	// Need to determine a a base vector through p1 in this plane
	v[0] = 1; v[1] = 0; v[2] = 0;	// try this
	dp = dotproduct(v,N);
	dp = abs(dp);
	if (dp > 0.9) {					// too close to N, use this
		v[0] = 0; v[1] = 1; v[2] = 0;
	}
	// Define v0 as the vector normal to both v and N, i.e. lies in the plane
	crossProduct(v0,v,N);
	dp = dotproduct(v0,v0);
	dp = sqrt(dp);
	for (i=0; i<3; i++) {
		v0[i] = v0[i]/dp;		// unit vector in the plane
	}
	r2sum = 0;
	*r2min = 1.0e10;
	nrays = 32;
	*zero = false;
	for (i=0;i<nrays;i++) {
		angle = (2*PI*i)/nrays;
		Rotate(v0,N,angle,v);
		r2 = GetRadius2(p1,v);		// problem if voxels are on the boundary of the image
//		if (r2 == 0) *zero = true;
		r2sum += r2;
		if (r2 < *r2min) *r2min = r2;
	}
	*r2ave = r2sum/(2*nrays);
	*r2ave = MAX(RMIN*RMIN,*r2ave);
	if (*zero) nzerodiameters++;
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// This method uses the preceding and following pts on an edge to estimate the centreline direction
// vector.  To improve the estimate of the cutting plane perpendicular to the vessel, the pts need 
// to be well separated, i.e. this must be applied to the simplified network.
// It can be used only for the interior pts, therefore not when npts = 2.
// NOT USED NOW
//-----------------------------------------------------------------------------------------------------
double getDiameter(int kp0, int kp1, int kp2)
{
	int i;
	double p1[3], p2[3], p0[3];
	double r2_ave, r2_min, diam;
	double dlim = 50.0;
	bool zero;

	for (i=0; i<3; i++) {
//		p0[i] = vsize[i]*(Vlist[kp0].pos[i] + 0.5);		// Note: centres of voxel cubes
//		p1[i] = vsize[i]*(Vlist[kp1].pos[i] + 0.5);
//		p2[i] = vsize[i]*(Vlist[kp2].pos[i] + 0.5);
		p0[i] = vsize[i]*(Vlist[kp0].initial_pos[i] + 0.5);		// Note: these pos are guaranteed to fall within the object image
		p1[i] = vsize[i]*(Vlist[kp1].initial_pos[i] + 0.5);
		p2[i] = vsize[i]*(Vlist[kp2].initial_pos[i] + 0.5);
	}
	// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
	EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min,&zero);
	if (zero) {
		int *pos;
		printf("getDiameter: EstimateDiameter gave a zero for: %d %d %d\n",kp0,kp1,kp2);
		fprintf(fpout,"getDiameter: EstimateDiameter gave a zero for: %d %d %d\n",kp0,kp1,kp2);
		pos = Vlist[kp0].initial_pos;
		fprintf(fpout,"kp0: %d %d %d  V3D: %d\n", pos[0],pos[1],pos[2],V3D(pos[0],pos[1],pos[2]));
		pos = Vlist[kp1].initial_pos;
		fprintf(fpout,"kp1: %d %d %d  V3D: %d\n", pos[0],pos[1],pos[2],V3D(pos[0],pos[1],pos[2]));
		pos = Vlist[kp2].initial_pos;
		fprintf(fpout,"kp2: %d %d %d  V3D: %d\n", pos[0],pos[1],pos[2],V3D(pos[0],pos[1],pos[2]));
		exit(99);
	}
	diam = 2*sqrt(r2_ave);
	if (calib_param != 0) {
		diam *= calib_param;
	}
	return diam;
}

//-----------------------------------------------------------------------------------------------------
// Estimates diameters for all interior points.
// Note that an edge with npts=2 (no intermediate points) is not processed.
// This is how it was done in topo
// Note: if this is used, then the createDistributions_topo must be used.
//-----------------------------------------------------------------------------------------------------
int getDiameters_topo(void)
{
	int ie, npts, k, kp0, kp1, kp2;
	EDGE *edge;

	printf("getDiameters_topo\n");
	fprintf(fpout,"getDiameters_topo: ne: %d  np: %d\n",ne,np);
	fflush(fpout);
	for (k=1; k<=np; k++) {
		avediameter[k] = 0;
//		Vlist[k].diameter = 0;
	}
	ne_used = 0;
	for (ie=0;ie<ne;ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
		ne_used++;
		npts = edge->npts;
		for (k=1;k<npts-1;k++) {
			kp0 = edge->pt[k-1];
			kp1 = edge->pt[k];
			kp2 = edge->pt[k+1];
			if (npts < 3) {
				printf("npts < 3: ie: %d %d\n",ie,npts);
				return 1;
			}
			avediameter[kp1] = getDiameter(kp0,kp1,kp2);
		}
	}
	fprintf(fpout,"getDiameters_topo done\n");
	fflush(fpout);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Estimates average diameter for all segments (edges).
// Note that lit voxel numbering now starts at 1
// The segment length has already been calculated, here it is checked.
//-----------------------------------------------------------------------------------------------------
int getDiameters(void)
{
	int ie, npts, k, kp0, kp1, kp2, i;
	double p1[3], p2[3], p0[3];
	double r2_ave, r2_min, r2_prev, diam, vol;
	double d, dave, dlen, dlenx, len, dvol, area;
	EDGE *edge;
	VOXEL *pv;
	bool zero;

	printf("getDiameters\n");
	fprintf(fpout,"getDiameters\n");
	nzerodiameters = 0;

	if (!use_object) {
		for (ie=0;ie<ne;ie++) {
			edge = &edgeList[ie];
			npts = edge->npts;
			edge->segavediam = FIXED_DIAMETER;
			for (k=0; k<npts; k++) {
				kp1 = edge->pt[k];
				pv = &Vlist[kp1];
				pv->diameter = FIXED_DIAMETER;
			}
		}
		return 0;
	}

	for (i=0; i<100; i++)
		diamdist[i] = 0;
	if (!METHOD2) {
		for (k=1; k<=np; k++) {
			Vlist[k].diameter = 0;
		}
	}
	ne_used = 0;
	for (ie=0;ie<ne;ie++) {
		edge = &edgeList[ie];
		npts = edge->npts;
//		printf("ie: %d npts: %d\n",ie,npts);
		if (npts == 2) {
			kp0 = edge->pt[0];
			kp1 = edge->pt[1];
//			printf("kp0% %d kp1: %d diams: %f %f\n",kp0,kp1,Vlist[kp0].diameter,Vlist[kp1].diameter);
			len = dist_um(Vlist[kp0].initial_pos,Vlist[kp1].initial_pos);
			if (METHOD2) {
				d = 0.5*(Vlist[kp0].diameter + Vlist[kp1].diameter);
			} else {
				for (i=0; i<3; i++) {
					p0[i] = vsize[i]*(Vlist[kp0].initial_pos[i] + 0.5);		// Note: these pos are not guaranteed to fall within the object image
					p1[i] = vsize[i]*(Vlist[kp1].initial_pos[i] + 0.5);
					p2[i] = 2*p1[i] - p0[i];
				}
				EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min,&zero);
				if (zero) {
					printf("getDiameters: got a zero for ie: %d  npts=2: %d %d\n",ie,kp0,kp2);
					fprintf(fpout,"getDiameters: got a zero for ie: %d  npts=2: %d %d\n",ie,kp0,kp2);
					exit(1);
				}
				d = 2*sqrt(r2_ave);
			}
			area = PI*pow(d/2,2);
		} else {	// sum volumes and lengths of pieces, average area = vol/len
			dave = 0;
			vol = 0;
			len = 0;
			for (k=1;k<npts-1;k++) {
				kp0 = edge->pt[k-1];
				kp1 = edge->pt[k];
				kp2 = edge->pt[k+1];
				dlen = dist_um(Vlist[kp0].initial_pos,Vlist[kp1].initial_pos);
//				printf("kp0,kp1: %d %d diams: %f %f\n",kp0,kp1,Vlist[kp0].diameter,Vlist[kp1].diameter);
				if (METHOD2) {
					r2_ave = 0.5*(pow(Vlist[kp0].diameter/2,2) + pow(Vlist[kp1].diameter/2,2));
				} else {
					// Basic check for insideness
					pv = &Vlist[kp1];
					if (V3D(pv->pos[0],pv->pos[1],pv->pos[2]) == 0) {
						printf("In getDiameters: edge: %d voxel: %d is not inside the object\n",ie,kp1);
						fprintf(fpout,"In getDiameters: edge: %d voxel: %d is not inside the object\n",ie,kp1);
						exit(1);
					}
					for (i=0; i<3; i++) {
						p0[i] = vsize[i]*(Vlist[kp0].initial_pos[i] + 0.5);		// Note: these pos are guaranteed to fall within the object image
						p1[i] = vsize[i]*(Vlist[kp1].initial_pos[i] + 0.5);
						p2[i] = vsize[i]*(Vlist[kp2].initial_pos[i] + 0.5);
					}
					// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
					EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min,&zero);
					if (zero) {
						printf("getDiameters: got a zero for ie: %d  %d %d %d\n",ie,kp0,kp1,kp2);
						fprintf(fpout,"getDiameters: got a zero for ie: %d  %d %d %d\n",ie,kp0,kp1,kp2);
						exit(1);
					}
				}

				if (k == 1) {
					dvol = PI*dlen*r2_ave;
				} else {
					dvol = PI*dlen*(r2_ave + r2_prev)/2;
				}
				if (k == npts-2) {	// need to add the last piece
					dlenx = dist_um(Vlist[kp1].initial_pos,Vlist[kp2].initial_pos);
					dlen += dlenx;
					if (METHOD2) {
						r2_ave = 0.5*(pow(Vlist[kp1].diameter/2,2) + pow(Vlist[kp2].diameter/2,2));
					}
					dvol += PI*dlenx*r2_ave;
				}
				vol += dvol;
				len += dlen;
				r2_prev = r2_ave;
			}
//			d = dave/(npts-2);
			if (fabs(len-edge->length_um) > 0.001) {
				printf("getDiameters: length check: %d %f %f\n",ie,len,edge->length_um);
				exit(1);
			}
			area = vol/len;
			d = 2*sqrt(area/PI);
		}
		edge->segavediam = d;
		if (d == 0) {
			printf("Error: d = 0: ie: %d\n",ie);
			exit(1);
		}
		if (use_max_diameter && d > max_diameter) {
			edge->used = false;
		} else {
			ne_used++;
			edge->used = true;
			edge->volume = area*len;
			if (!METHOD2) {
				for (k=0; k<npts; k++) {
					kp1 = edge->pt[k];
					pv = &Vlist[kp1];
					pv->diameter = d;
				}
			}
		}
		//pv = &Vlist[edge->pt[0]];
		//pv->diameter = MAX(d,pv->diameter);
		//pv = &Vlist[edge->pt[npts-1]];
		//pv->diameter = MAX(d,pv->diameter);
		//int n = d/ddiam + 0.5;	// needs to be length weighted
		//n = MAX(n,1);
		//n = MIN(n,100);
		//diamdist[n-1]++;	
	}
	//printf("n  diam count\n");
	//fprintf(fpout,"n  diam count\n");
	//for (i=0; i<100; i++) {
	//	printf("%d %6.1f %d\n",i+1,(i+1)*ddiam,diamdist[i]);
	//	fprintf(fpout,"%d %6.1f %d\n",i+1,(i+1)*ddiam,diamdist[i]);
	//}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int	removeShortEnds()
{
#define nblist 1000
	int k0, k, k1, nbrs, n, ib, i, nends;
	int blist[nblist];
	double endlen;
	VOXEL *pv, *pv0;

	printf("removeShortEnds\n");
	nends = 0;
	for (k0=1; k0<=nlit; k0++) {
		pv = &Vlist[k0];
		if (pv->nbrs != 1) continue;
//		printf("\nk0: %d nbrs: %d\n",k0,pv->nbrs);
		n = 0;
		blist[n] = k0;
		n = 1;
		// Check the length of this loose end
		k = pv->nbr[0];
		k1 = k0;
		endlen = dist_um(pv->pos,Vlist[k].pos);
//		printf("endlen: %f\n",endlen);
		blist[n] = k;
		n++;
		for(;;) {
			pv = &Vlist[k];
			nbrs = pv->nbrs;
//			printf("k: %d nbrs: %d\n",k,nbrs);
			if (nbrs > 2) break;	// reached a vertex
			if (nbrs < 2) {
				printf("Error: k0: %d k: %d nbrs: %d\n",k0,k,nbrs);
				for (i=0; i<n; i++) {
					printf("%d %d\n",i,blist[i]);
				}
				exit(1);
			}
			// Find the next neighbour on the segment
			if (pv->nbr[0] == k1) {
				k1 = k;
				k = pv->nbr[1];
			} else {
				k1 = k;
				k = pv->nbr[0];
			}
			endlen += dist_um(pv->pos,Vlist[k].pos);
//			printf("endlen: %f\n",endlen);
			blist[n] = k;
			n++;
		}
		if (endlen < min_end_len) {
			// Remove this short loose end
			// First remove the voxels in the list, except for the last one
			for (i=0; i<n-1; i++) {
				k = blist[i];
				pv = &Vlist[k];
				pv->nbrs = 0;
				Vindex(pv->pos[0],pv->pos[1],pv->pos[2]) = 0;
//				printf("removed: %d\n",k);
			}
			pv0 = &Vlist[blist[n-1]];	// This is the vertex pointer
			k = blist[n-2];				// This is the last voxel on the loose end before the vertex
//			printf("vertex is: %d last voxel: %d\n",blist[n-1],k);
//			printf("nbr list: ");
//			for (i=0; i<pv0->nbrs; i++)
//				printf("  %d",pv0->nbr[i]);
//			printf("\n");
			for (ib=0; ib<pv0->nbrs; ib++) {
				if (pv0->nbr[ib] == k) break;
			}
			// The loose end is on branch ib
			// Now remove this branch from pv0->nbr[]
//			printf("nbr to remove: ib: %d\n",ib);
			for (i=ib; i<pv0->nbrs-1; i++) {
				pv0->nbr[i] = pv0->nbr[i+1];
			}
			pv0->nbrs--;
//			printf("new nbr list: ");
//			for (i=0; i<pv0->nbrs; i++)
//				printf("  %d",pv0->nbr[i]);
//			printf("\n");
			nends++;
		}
	}
	return nends;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int	removeLoops()
{
#define nblist 1000
	int k0, k, k1, kb, nbrs, nbrs0, n, ib, i, nloops;
	int blist[nblist];
	int nbr[30];
	VOXEL *pv, *pv0;

	printf("removeLoops\n");
	nloops = 0;
	for (k0=1; k0<=nlit; k0++) {
		pv0 = &Vlist[k0];
		nbrs0 = pv0->nbrs;
		if (nbrs0 < 3) continue;
//		printf("\nk0: %d nbrs: %d\n",k0,nbrs0);
		for (ib=0; ib<nbrs0; ib++) {
			n = 0;
			k1 = k0;
			k = pv0->nbr[ib];
			blist[n] = k;
			n++;
			for(;;) {
				pv = &Vlist[k];
				nbrs = pv->nbrs;
	//			printf("k: %d nbrs: %d\n",k,nbrs);
				if (nbrs > 2) break;	// reached a vertex
				if (nbrs < 2) break;	// a loose end
				// Find the next neighbour on the segment
				if (pv->nbr[0] == k1) {
					k1 = k;
					k = pv->nbr[1];
				} else {
					k1 = k;
					k = pv->nbr[0];
				}
				blist[n] = k;
				n++;
			}
			if (nbrs < 2) continue;
			if (k == k0) {  // Loops back to starting vertex
				n--;	// ignore last blist, which is k0
//				printf("nbrs0: %d nbrs: %d\n",nbrs0,nbrs);
				// Need to remove the voxels at blist[0] and blist[n-1] from the nbr[] for k0,
				// then remove the voxels 0 .. n-1
//				printf("n: %d ends of loop: %d %d\n",n,blist[0],blist[n-1]);
				// First fix the nbr list at k0
				kb = 0;
				for (ib=0; ib<nbrs0; ib++) {
					if (pv0->nbr[ib] != blist[0] && pv0->nbr[ib] != blist[n-1]) {
						nbr[kb] = pv0->nbr[ib];
						kb++;
					}
				}
				pv0->nbrs = kb;
				for (ib=0; ib<kb; ib++) {
					pv0->nbr[ib] = nbr[ib];
				}

				// Remove the voxels in the blist
				for (i=0; i<n; i++) {
					k = blist[i];
					pv = &Vlist[k];
					pv->nbrs = 0;
					Vindex(pv->pos[0],pv->pos[1],pv->pos[2]) = 0;
	//				printf("removed: %d\n",k);
				}
//				printf("k0: %d new nbrs0: %d\n",k0,pv0->nbrs);

	//			printf("new nbr list: ");
	//			for (i=0; i<pv0->nbrs; i++)
	//				printf("  %d",pv0->nbr[i]);
	//			printf("\n");
				nloops++;
				break;	// remove at most one loop per vertex
			}
		}
	}
	return nloops;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int createVlist()
{
#define nblist 1000
	int k, k0, kp, ke, kl, x, y, z, ib, ib1, ib2, i, j, len, n, nb, nloops, nends;
	int blist[nblist];
	int nbr_temp[30];
	VOXEL *pv, *pv0;
	bool repeat, dbug;

	printf("createVlist\n");
	k = 0;
	for (x=0;x<width; x++) {
		for (y=0;y<height; y++) {
			for (z=0;z<depth; z++) {
				if (V3Dskel(x,y,z) != 0) {
					k++;
					Vindex(x,y,z) = k;
					Vlist[k].pos[0] = x;
					Vlist[k].pos[1] = y;
					Vlist[k].pos[2] = z;
					Vlist[k].initial_pos[0] = x;
					Vlist[k].initial_pos[1] = y;
					Vlist[k].initial_pos[2] = z;
					Vlist[k].diameter = 0;
				}
			}
		}
	}
	printf("Create nbr lists\n");
	// Find neighbour voxels
	for (k=1; k<=nlit; k++) {
//		if (k%1000 == 0) printf("k: %d\n",k);
		int xmin, xmax, ymin, ymax, zmin, zmax;
		x = Vlist[k].pos[0]; 
		y = Vlist[k].pos[1]; 
		z = Vlist[k].pos[2]; 
		xmin = MAX(0,x-1);
		xmax = MIN(width-1,x+1);
		ymin = MAX(0,y-1);
		ymax = MIN(height-1,y+1);
		zmin = MAX(0,z-1);
		zmax = MIN(depth-1,z+1);
		int nbrs = 0;
		for (int xx=xmin; xx<=xmax; xx++) {
			for (int yy=ymin; yy<=ymax; yy++) {
				for (int zz=zmin; zz<=zmax; zz++) {
					int kk = Vindex(xx,yy,zz);
					if (kk == k) continue;
					if (kk > 0) {
						nbrs++;
						if (nbrs > MAXNBRS) {
							printf("Error: too many neighbours\n");
							fprintf(fperr,"Error: too many neighbours\n");
							return 1;
						}
						Vlist[k].nbr[nbrs-1] = kk;
					}
				}
			}
		}
		Vlist[k].nbrs = nbrs;
	}

//	imskel = nullptr;
	imskel = ImageType::New();

	for (i=0; i<4; i++) {
		nends = removeShortEnds();	
		printf("Removed short ends: %d\n",nends);
		fprintf(fpout,"Removed short ends: %d\n",nends);
		nloops = removeLoops();	
		printf("Removed loops: %d\n",nloops);
		fprintf(fpout,"Removed loops: %d\n",nloops);
		if (nends == 0 && nloops == 0) break;
	}

	avediameter = (float *)malloc((nlit+1)*sizeof(float));
	fprintf(fpout,"Allocated avediameter: %d\n",nlit+1);

	return 0;
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int createVlist1()
{
#define nblist 1000
	int k, k0, kp, ke, kl, x, y, z, ib, ib1, ib2, i, j, len, n, nb, nloops, nends;
	int blist[nblist];
	int nbr_temp[30];
	VOXEL *pv, *pv0;
	bool repeat, dbug;
	int ndbug = 6;
	int kdbug[6] = {5346350, 5345611, 5344788, 5347081, 5347867, 5347082};

	k = 0;
	for (x=0;x<width; x++) {
		for (y=0;y<height; y++) {
			for (z=0;z<depth; z++) {
				if (V3Dskel(x,y,z) != 0) {
					k++;
					Vindex(x,y,z) = k;
					Vlist[k].pos[0] = x;
					Vlist[k].pos[1] = y;
					Vlist[k].pos[2] = z;
					Vlist[k].initial_pos[0] = x;
					Vlist[k].initial_pos[1] = y;
					Vlist[k].initial_pos[2] = z;
					Vlist[k].diameter = 0;
				}
			}
		}
	}
	// Find neighbour voxels
	for (k=1; k<=nlit; k++) {
		int xmin, xmax, ymin, ymax, zmin, zmax;
		x = Vlist[k].pos[0]; 
		y = Vlist[k].pos[1]; 
		z = Vlist[k].pos[2]; 
		xmin = MAX(0,x-1);
		xmax = MIN(width-1,x+1);
		ymin = MAX(0,y-1);
		ymax = MIN(height-1,y+1);
		zmin = MAX(0,z-1);
		zmax = MIN(depth-1,z+1);
		int nbrs = 0;
		for (int xx=xmin; xx<=xmax; xx++) {
			for (int yy=ymin; yy<=ymax; yy++) {
				for (int zz=zmin; zz<=zmax; zz++) {
					int kk = Vindex(xx,yy,zz);
					if (kk == k) continue;
					if (kk > 0) {
						nbrs++;
						if (nbrs > MAXNBRS) {
							printf("Error: too many neighbours\n");
							fprintf(fperr,"Error: too many neighbours\n");
							return 1;
						}
						Vlist[k].nbr[nbrs-1] = kk;
					}
				}
			}
		}
		Vlist[k].nbrs = nbrs;
	}

//	imskel = nullptr;
	imskel = ImageType::New();

	for (i=0; i<ndbug; i++) {
		k = kdbug[i];
		printf("voxel: %d nbrs: %d\n",k,Vlist[k].nbrs);
		for (j=0; j<Vlist[k].nbrs; j++) {
			printf("   %d %d",j,Vlist[k].nbr[j]);
		}
		printf("\n");
	}

	printf("Start cleanup\n");
	fprintf(fpout,"Start cleanup\n");
	// Do basic cleanup - remove simple loops and short ends
	nloops = 0;
	nends = 0;
	for (k0=1; k0<=nlit; k0++) {
		pv0 = &Vlist[k0];
		if (k0 == 5346350) printf("createVlist: k0: %d nbrs: %d\n",k0,pv0->nbrs);
		if (pv0->nbrs < 3) continue;
		dbug = (k0 == 5346350);
		if (dbug) {

			for (i=0; i<ndbug; i++) {
				k = kdbug[i];
				printf("voxel: %d nbrs: %d\n",k,Vlist[k].nbrs);
				for (j=0; j<Vlist[k].nbrs; j++) {
					printf("   %d %d",j,Vlist[k].nbr[j]);
				}
				printf("\n");
			}

			printf("k0: %d  nbrs: %d  ", k0,pv0->nbrs);
			for (i=0; i<pv0->nbrs; i++)
				printf("%d",pv0->nbr[i]);
			printf("\n");
		}
		do {
			repeat = false;
			for (ib=0; ib<pv0->nbrs; ib++) {	// look at every branch ib
				kp = k0;
				k = pv0->nbr[ib];
				len = 0;
				for (;;) {	// traverse the branch
					pv = &Vlist[k];
					if (dbug) {
						printf("ib: %d k: %d  nbrs: %d\n",ib,k,pv->nbrs);
					}
					if (len >= nblist) {
						printf("createVlist: voxel k0: %d nbrs: %d ib: %d branch length len: %d exceeds nblist: %d\n",k0,pv0->nbrs,ib,len,nblist); 
						exit(1);
					}
					blist[len] = k;	// keep track of the branch voxels
					len++;
					if (pv->nbrs > 2) {		// another vertex, end of segment
						if (k < k0) {
							// segment end voxel is an already processed voxel, i.e. this is a duplicate
						} else if (k == k0) {
							// a loop - this is to be removed IF none of the loop voxels has nbrs > 2
							// Note: when a loop is removed, because pv0->nbr[] has changed we must exit this vertex
							repeat = true;
							for (i=0; i<len-1; i++) {
								if (i >= nblist) exit(1);
								kl = blist[i];
								pv = &Vlist[kl];
								if (pv->nbrs > 2) repeat = false;
							}
							if (!repeat) break;
							nloops++;
							for (i=0; i<len-1; i++) {
								if (i >= nblist) exit(1);
								kl = blist[i];
								if (dbug) printf("remove kl: %d\n",kl);
								pv = &Vlist[kl];
								pv->nbrs = 0;
								Vindex(pv->pos[0],pv->pos[1],pv->pos[2]) = 0;
							}
							// Now remove the two branches of the loop from pv0->nbr[]
							// One branch is ib1 = ib, the other, ib2, has pv0->nbr[ib2] = kp
							ib1 = ib;
							for (ib=0; ib<pv0->nbrs; ib++) {
								if (ib >= 30) exit(1);
								nbr_temp[ib] = pv0->nbr[ib];
								if (pv0->nbr[ib] == kp) {
									ib2 = ib;
								}
							}
							n = 0;
							for (ib=0; ib<pv0->nbrs; ib++) {
								if (ib != ib1 && ib != ib2) {
									pv0->nbr[n] = pv0->nbr[ib];
									n++;
								}
							}
							pv0->nbrs = pv0->nbrs-2;
						} else {
							// We have completed an internal segment, from k0 to k
						}
						break;
					} else if (pv->nbrs == 2) {
						// Find the next neighbour on the segment
						if (pv->nbr[0] == kp) {
							kp = k;
							k = pv->nbr[1];
						} else {
							kp = k;
							k = pv->nbr[0];
						}
						if (dbug) printf("next nbr: %d\n",k);
					} else if (pv->nbrs == 1) {
						// We have completed an end segment, from k0 to k
						if (len < min_end_len) {
							if (dbug) printf("Removing short loose end\n");
							// a short end, to be removed
							// Note: when an end is removed, because pv0->nbr[] has changed we must exit this vertex
							nends++;
							repeat = true;
							for (i=0; i<len; i++) {
								if (i >= nblist) exit(1);
								ke = blist[i];
								if (dbug) printf("remove ke: %d, set nbrs = 0\n",ke);
								pv = &Vlist[ke];
								pv->nbrs = 0;
								Vindex(pv->pos[0],pv->pos[1],pv->pos[2]) = 0;
							}
							// Now remove this branch from pv0->nbr[]
							for (i=ib; i<pv0->nbrs-1; i++) {
								pv0->nbr[i] = pv0->nbr[i+1];
							}
							pv0->nbrs--;
							if (dbug) {
								printf("removed branch from pv0, now pv0->nbrs = %d\n",pv0->nbrs);
								for (int i=0; i<pv0->nbrs; i++) {
									printf("i: %d %d\n",i,pv0->nbr[i]);
								}
							}
						}
						break;
					} else if (pv->nbrs == 0) {
						printf("Error: nbrs = 0: voxel: %d\n",k);
						fprintf(fperr,"Error: nbrs = 0: voxel: %d\n",k);
						return 1;
					}
				}
//				if (repeat) continue;
				if (repeat) break;
			}
		} while (repeat);
	}
	printf("Removed simple loops: %d and short ends: %d\n",nloops,nends);
	fprintf(fpout,"Removed simple loops: %d and short ends: %d\n",nloops,nends);
	avediameter = (float *)malloc((nlit+1)*sizeof(float));
	fprintf(fpout,"Allocated avediameter: %d\n",nlit+1);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int traceSegments()
{
	int k0, ib, k, kp, nsegvoxels, nsegs, nloops, nends, nshow=0;
	int segvoxel[1000];	//[200];
	float len;
	VOXEL *pv, *pv0, *pv1;
	EDGE *pe;
	bool duplicate, loop, flag, dbug;

	printf("\ntraceSegments\n");
	fprintf(fpout,"\ntraceSegments\n");
	fflush(fpout);
	net->edgeList = (EDGE *)malloc(MAXEDGES*sizeof(EDGE));
	edgeList = net->edgeList;

	for (int i=0; i<100; i++)
		lendist[i] = 0;
	flag = true;
	nsegs = 0;
	nloops = 0;
	nends = 0;
	for (k0=1; k0<=nlit; k0++) {
		dbug = (k0 == -1);
		pv0 = &Vlist[k0];
		if (dbug) fprintf(fpout,"k0: %d nbrs: %d\n",k0,pv0->nbrs);
		if (pv0->nbrs < 3) continue;	// not a 3-vertex
		if (nsegs < nshow) {
			printf("\nstart vertex: %d  nbrs: %d\n",k0,pv0->nbrs);
		}
//		fprintf(fpout,"k0: %d nsegs: %d nbrs: %d\n",k0,nsegs,pv0->nbrs);
//		fflush(fpout);
		for (ib=0; ib<pv0->nbrs; ib++) {	// branch i
			kp = k0;
			pv1 = pv0;
			k = pv0->nbr[ib];
			nsegvoxels = 0;
			segvoxel[nsegvoxels] = k0;
			nsegvoxels++;
			duplicate = false;
			loop = false;
			len = 0;
//			if (nsegs == 0) printf("seg 0: k: %d len: %f\n",k0,len);
			for (;;) {
				if (dbug) printf("branch: %d next: %d\n",ib,k);
				segvoxel[nsegvoxels] = k;
				nsegvoxels++;
				pv = &Vlist[k];
				if (dbug) printf("nbrs: %d\n",pv->nbrs);
				if (pv->nbrs > 2) {		// another vertex, end of segment
					if (k < k0) {
						// segment end voxel is an already processed voxel, i.e. this is a duplicate
						duplicate = true;
					} else if (k == k0) {
						// a loop - this is a problem!
						loop = true;
						nloops++;
					} else {
						// We have completed an internal segment, from k0 to k (save it later)
						float dlen = dist_um(pv->pos,pv1->pos);
						if (dlen == 0) {
							printf("zero dlen: branch: %d kp: %d  k: %d\n",ib,kp,k);
							printf("pos: %d %d %d   %d %d %d\n",pv1->pos[0],pv1->pos[1],pv1->pos[2],pv->pos[0],pv->pos[1],pv->pos[2]);
							fprintf(fperr,"zero dlen: branch: %d kp: %d  k: %d\n",ib,kp,k);
							fprintf(fperr,"pos: %d %d %d   %d %d %d\n",pv1->pos[0],pv1->pos[1],pv1->pos[2],pv->pos[0],pv->pos[1],pv->pos[2]);
							return 1;
						}
						len += dlen;
//						if (nsegs == 0) printf("seg 0: k: %d len: %f\n",k,len);
						int klen = MIN(99,int(len + 0.5));
						lendist[klen]++;

						if (nsegs >= MAXEDGES) {
							printf("Error: MAXEDGES exceeded: %d\n",nsegs);
							fprintf(fperr,"Error: MAXEDGES exceeded: %d\n",nsegs);
							return 1;
						}
						pe = &edgeList[nsegs];
						pe->npts = nsegvoxels;
						pe->pt = (int *)malloc(pe->npts*sizeof(int));
						for (int i=0; i<pe->npts; i++)
							pe->pt[i] = segvoxel[i];
						pe->length_um = len;
						pe->used = true;
//						Vlist[k0].diameter = 0;
//						Vlist[k].diameter = 0;

						nsegs++;
					}
					break;
				} else if (pv->nbrs == 2) {
					float dlen = dist_um(pv->pos,pv1->pos);
					if (dlen == 0) {
						printf("zero dlen: branch: %d k0: %d  k: %d\n",ib,k0,k);
						printf("pos: %d %d %d   %d %d %d\n",pv1->pos[0],pv1->pos[1],pv1->pos[2],pv->pos[0],pv->pos[1],pv->pos[2]);
						fprintf(fperr,"zero dlen: branch: %d k0: %d  k: %d\n",ib,k0,k);
						fprintf(fperr,"pos: %d %d %d   %d %d %d\n",pv1->pos[0],pv1->pos[1],pv1->pos[2],pv->pos[0],pv->pos[1],pv->pos[2]);
						return 1;
					}
					len += dlen;
//					if (nsegs == 0) printf("seg 0: k: %d len: %f\n",k,len);
					pv1 = pv;
					// Find the next neighbour on the segment
					if (pv->nbr[0] == kp) {
						kp = k;
						k = pv->nbr[1];
					} else {
						kp = k;
						k = pv->nbr[0];
					}
				} else if (pv->nbrs == 1) {
					// We have completed an end segment, from k0 to k (save it later)
					float dlen = dist_um(pv->pos,pv1->pos);
					if (dlen == 0) {
						printf("zero dlen: branch: %d k0: %d  k: %d\n",ib,k0,k);
						printf("pos: %d %d %d   %d %d %d\n",pv1->pos[0],pv1->pos[1],pv1->pos[2],pv->pos[0],pv->pos[1],pv->pos[2]);
						fprintf(fperr,"zero dlen: branch: %d k0: %d  k: %d\n",ib,k0,k);
						fprintf(fperr,"pos: %d %d %d   %d %d %d\n",pv1->pos[0],pv1->pos[1],pv1->pos[2],pv->pos[0],pv->pos[1],pv->pos[2]);
						return 1;
					}
					len += dlen;
//					if (nsegs == 0) printf("seg 0: k: %d len: %f\n",k,len);
//					fprintf(fpout,"end: %d  len: %f\n",nends,len);
					nends++;
					if (len > min_end_len) {	// drop short ends
						int klen = MIN(99,int(len + 0.5));
						lendist[klen]++;

						if (nsegs >= MAXEDGES) {
							printf("Error: MAXEDGES exceeded: %d\n",nsegs);
							fprintf(fperr,"Error: MAXEDGES exceeded: %d\n",nsegs);
							return 1;
						}
						pe = &edgeList[nsegs];
						pe->npts = nsegvoxels;
						pe->pt = (int *)malloc(pe->npts*sizeof(int));
						for (int i=0; i<pe->npts; i++)
							pe->pt[i] = segvoxel[i];
						pe->length_um = len;
						pe->used = true;
//						Vlist[k0].diameter = 0;
//						Vlist[k].diameter = 0;

						nsegs++;
					} else {	// remove these voxels
						for (int i=1; i<nsegvoxels; i++) {
							Vlist[segvoxel[i]].nbrs == 0;
						}
					}
					break;
				} else if (pv->nbrs == 0) {
					printf("Error: nbrs = 0: voxel: %d\n",k);
					printf("k0: %d ib: %d k: %d\n",k0,ib,k);
					fprintf(fperr,"Error: nbrs = 0: voxel: %d\n",k);
					fprintf(fperr,"k0: %d ib: %d k: %d\n",k0,ib,k);
					return 1;
				}
			}
			if (nsegs < nshow) {
				if (duplicate) {
					printf("branch: %d  duplicate: ",ib);
					fprintf(fpout,"branch: %d  duplicate: ",ib);
				} else if (loop) {
					printf("branch: %d  loop\n",ib);
					fprintf(fpout,"branch: %d  loop\n",ib);
				} else {
					printf("branch: %d  segment: %d: ",ib,nsegs);
					fprintf(fpout,"branch: %d  segment: %d: ",ib,nsegs);
				}
				printf(" %d -> %d  nvoxels: %d: ",k0,k,nsegvoxels);
				fprintf(fpout," %d -> %d  nvoxels: %d: ",k0,k,nsegvoxels);
				if (nsegvoxels > 2) {
					for (int i=0; i<nsegvoxels; i++) {
						printf("%d  ",segvoxel[i]);
						fprintf(fpout,"%d  ",segvoxel[i]);
					}
				}
				printf("\n");
				fprintf(fpout,"\n");
			}
		}
		fflush(fpout);
	}
	printf("nsegs: %d\n",nsegs);
	printf("nloops: %d\n",nloops);
	fprintf(fpout,"nsegs: %d\n",nsegs);
	fprintf(fpout,"nloops: %d\n",nloops);
	printf("did traceSegments\n");
	fflush(fpout);
//	fprintf(fpout,"length distribution\n");
//	for (int i=0; i<100; i++) {
//		printf("%d %f\n",i,lendist[i]/float(nsegs));
//		fprintf(fpout,"%d %f\n",i,lendist[i]/float(nsegs));
//	}
	ne = nsegs;
	net->ne = nsegs;
	return 0;
}

/*
//-----------------------------------------------------------------------------------------------------
// Count # of vertices in neighbourhood of each vertex
//-----------------------------------------------------------------------------------------------------
void vertexDensity1()
{
	int kv, kv0, x, y, z, n, nvt, nc;
	int *count;
	int d = 2;
	VOXEL *pv, *pv0;

	nc = (2*d+1)*(2*d+1)*(2*d+1);
	count = (int *)malloc(nc*sizeof(int));
	for (int i=0;i<nc;i++) 
		count[i] = 0;
	nvt = 0;
	for (kv0=1; kv0<=nlit; kv0++) {
		pv0 = &Vlist[kv0];
		n = 0;
		if (pv0->nbrs < 3) continue;
		nvt++;
		for (x = pv0->pos[0]-d; x <= pv0->pos[0]+d; x++) {
			for (y = pv0->pos[1]-d; y <= pv0->pos[1]+d; y++) {
				for (z = pv0->pos[2]-d; z <= pv0->pos[2]+d; z++) {
					kv = Vindex(x,y,z);
					if (kv == 0) continue;
					pv = &Vlist[kv];
					if (pv->nbrs > 2) n++;
				}
			}
		}
		count[n]++;
	}
	for (n=1; n<nc;n++) {
		printf("n: %d count[n]: %d\n",n,count[n]);
		fprintf(fpout,"n: %d count[n]: %d\n",n,count[n]);
	}
	printf("nvt: %d\n",nvt);
	fprintf(fpout,"nvt: %d\n",nvt);
}
*/

//-----------------------------------------------------------------------------------------------------
// Returns true if kvd is in kvlist[]
//-----------------------------------------------------------------------------------------------------
bool in_list(int kvd, int kvlist[], int n)
{
	for (int i=0; i<n; i++) {
		if (kvd == kvlist[i]) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------
// In the 3x3x3 block at (ix,iy,iz), join the n vertices together
//-----------------------------------------------------------------------------------------------------
int coalesce(int kvlist[], int n, int kmid)
{
//	int newpos[3];
	int kvhat;
	VOXEL *pv0, *pv, *pvd;
	bool dbug = false;

//	printf("\ncoalesce: %d %d %d  n: %d\n",ix,iy,iz,n);
//	dbug = (ix==198) && (iy == 96) && (iz==228);
	dbug = (kmid == -1);
	if (dbug) printf("kvlist: %d %d %d %d\n",kvlist[0],kvlist[1],kvlist[2],kvlist[3]);

//	kvhat = kvlist[0];
	kvhat = kmid;
	pv0 = &Vlist[kvhat];
	// ------- Try leaving kvhat where it was, to avoid stepping on another point (should not be necessary) -------
	//if (Vindex(newpos[0],newpos[1],newpos[2]) == 0) {
	//	pv0->pos[0] = newpos[0];
	//	pv0->pos[1] = newpos[1];
	//	pv0->pos[2] = newpos[2];
	//	Vindex(newpos[0],newpos[1],newpos[2]) = 255;	// Note that the image is changed
	//}

	// Now look at all neighbours of the vertices in kvlist[]
	int hat_nbrs = 0;
	int hat_nbr[50];	// this will store the new list of neighbours for vertex kvhat
	for (int i=0; i<n; i++) {
		int kv = kvlist[i];
		pv = &Vlist[kv];
		if (dbug) printf("interior vertex #: %d %d  nbrs: %d\n",i,kv,pv->nbrs);
		for (int j=0; j<pv->nbrs; j++) {
			int kvd = pv->nbr[j];	// this is a neighbour of the vertex in the cube
			if (dbug) printf("kv: %d #: %d kvd: %d\n",kv,j,kvd);
			// We need to ignore any neighbour that is in kvlist[]
			if (in_list(kvd,kvlist,n)) continue;
			// We also need to ignore any neighbour that has already been added to hat_nbr[]!
			if (in_list(kvd,hat_nbr,hat_nbrs)) continue;
			hat_nbr[hat_nbrs] = kvd;	// add the neighbour to the new list of kvhat neighbours
			pvd = &Vlist[kvd];
			hat_nbrs++;
			// Now we need to adjust the neighbour list for kvd to replace all references to 
			// members of kvlist by a single entry for kvhat
			int nbrs = 0;
			int nbr[30];
			nbr[nbrs] = kvhat;
			nbrs++;
			if (pvd->nbrs > MAXNBRS) {
				printf("coalesce error: too many neighbours\n");
				fprintf(fperr,"coalesce error: too many neighbours\n");
				return 1;
			}
			for (int k=0; k<pvd->nbrs; k++) {
				int kvdd = pvd->nbr[k];
				if (in_list(kvdd,kvlist,n)) continue;
				nbr[nbrs] = kvdd;
				nbrs++;
			}
			// Now copy the new neighbour list for the outside neighbour
			if (dbug) printf("new nbrs for kvd: %d  %d\n",kvd,nbrs);
			pvd->nbrs = nbrs;
			for (int k=0; k<nbrs; k++) {
				pvd->nbr[k] = nbr[k];
				if (dbug) printf("%d %d\n",k,nbr[k]);
			}

		}
		if (kv != kvhat) {
			// Set the # of neighbours for the surplus vertices to 0.  This flags the vertex as no longer used. Also set Vindex(x,y,z) = 0
			pv->nbrs = 0;
			Vindex(pv->pos[0],pv->pos[1],pv->pos[2]) = 0;	// Note that the image is changed
		}
	}
	// Now copy the new neighbour list for kvhat
//	printf("coalesce: (%d %d %d) n: %d  kvhat: %d nbrs: %d\n",ix,iy,iz,n,kvhat,hat_nbrs);
//	fprintf(fpout,"coalesce: n: %d  kvhat: %d nbrs: %d\n",n,kvhat,hat_nbrs);
	pv0->nbrs = hat_nbrs;
	if (hat_nbrs > MAXNBRS) {
		printf("Error: # of neighbours exceeds limit of: %d\n",MAXNBRS);
		fprintf(fperr,"Error: # of neighbours exceeds limit of: %d\n",MAXNBRS);
		return 1;
	}
	if (dbug) printf("kvhat pos: %d %d %d\n",pv0->pos[0],pv0->pos[1],pv0->pos[2]);
	for (int j=0; j<pv0->nbrs; j++) {
		pv0->nbr[j] = hat_nbr[j];
		if (dbug) {
			pv = &Vlist[hat_nbr[j]];
			printf("nbr: %d %d pos: %d %d %d\n",j,hat_nbr[j],pv->pos[0],pv->pos[1],pv->pos[2]);
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Find voxel in kvlist nearest to the CoM of the n voxels
//-----------------------------------------------------------------------------------------------------
void selectMidVoxel(int kvlist[], int n, int *kmid)
{
	float sum[3];
	sum[0] = sum[1] = sum[2] = 0;
	for (int i=0; i<n; i++) {
		sum[0] += Vlist[kvlist[i]].pos[0];
		sum[1] += Vlist[kvlist[i]].pos[1];
		sum[2] += Vlist[kvlist[i]].pos[2];
	}
	// Find nearest k
	int kmin;
	float d2min = 1.0e10;
	for (int i=0; i<n; i++) {
		float d2 = 0;
		float del;
		for (int j=0; j<3; j++) {
			del = Vlist[kvlist[i]].pos[j] - sum[j]/n;
			d2 += del*del;
		}
		if (d2 < d2min) {
			d2 = d2min;
			kmin = i;
		}
	}
	*kmid = kvlist[kmin];
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void group(int kvlist1[],int n1, int kvlist2[], int *n2)
{
	int nn[20], imax, nmax, value[20];
	int list[20][20];
	VOXEL *pv;

	nmax = 0;
	for (int i=0; i<n1; i++) {
		value[i] = 2<<i;
		nn[i] = 0;
		list[i][nn[i]] = i;
		nn[i]++;
		int k0 = kvlist1[i];
		for (int j=0; j<n1; j++) {
			int k = kvlist1[j];
			pv = &Vlist[k];
			for (int ib=0; ib<pv->nbrs; ib++) {
				if (pv->nbr[ib] == k0) {
					list[i][nn[i]] = j;
					nn[i]++;
					value[i] += 2<<j;
					break;
				}
			}
		}
		if (nn[i] > nmax) {
			nmax = nn[i];
			imax = i;
		}
	}
	*n2 = nn[imax];
	for (int i=0; i<nn[imax]; i++) {
		kvlist2[i] = kvlist1[list[imax][i]];
	}
}

//-----------------------------------------------------------------------------------------------------
// Count # of vertices in every 3x3x3 cube, coalesce vertices
//-----------------------------------------------------------------------------------------------------
int vertexDensity2(int noffsets)
{
	int kv, kv0, x, y, z, n1, n2, nc, offset, ndropped, kmid, err;
	int nx, ny, nz, ix, iy, iz, xx, yy, zz;
	int nvt;
	int nvt_max = 10000;
	int *count;
	int kvlist1[27], kvlist2[27];
	int offsetx[4] = {0, 1, 0, 0};
	int offsety[4] = {0, 0, 1, 0};
	int offsetz[4] = {0, 0, 0, 1};
	VOXEL *pv, *pv0;
	bool random_vertex = false;

	if (random_vertex) {	// randomly select a vertices as foci for coalescence - doesn't work well, creates short segments
		printf("\nvertexDensity2: random vertex selection: nvt_max: %d\n",nvt_max);
		nc = 3*3*3;
		count = (int *)malloc(nc*sizeof(int));
		for (int i=0;i<nc;i++) count[i] = 0;
		ndropped = 0;
		nvt = 0;
		for (;;) {
			if (nvt > nvt_max) break;
			double R = (double)rand() / (double)(RAND_MAX+1);
			kv0 = R*nlit;
			pv = &Vlist[kv0];
			if (pv->nbrs != 3) continue;
			nvt++;
			ix = pv->pos[0];
			iy = pv->pos[1];
			iz = pv->pos[2];
			n1 = 0;
			for (x=MAX(ix-1,0); x <= MIN(width-1,ix+1); x++) {
				for (y=MAX(iy-1,0); y <= MIN(height-1,iy+1); y++) {
					for (z=MAX(iz-1,0); z <= MIN(depth-1,iz+1); z++) {
						kv = Vindex(x,y,z);
						if (kv == 0) continue;
						pv = &Vlist[kv];
						if (pv->nbrs > 0) {		// try adding all lit voxels to the list.  Need to use nbrs to avoid vertices dropped in coalesce
							kvlist1[n1] = kv;
							n1++;
						}
					}
				}
			}
			count[n1]++;
			if (n1 > 1) {
				// Need to select those that are connected, and the collapse voxel location
				group(kvlist1,n1,kvlist2,&n2);
				selectMidVoxel(kvlist2,n2,&kmid);
				err = coalesce(kvlist2,n2,kmid);
				if (err != 0) return 1;
				ndropped += n2-1;
			}
		}
		free(count);
		printf("Dropped points: %d\n",ndropped);
	} else {
		printf("\nvertexDensity2: noffsets: %d\n",noffsets);
		fprintf(fpout,"\nvertexDensity2: noffsets: %d\n",noffsets);
		nx = width/3-1;
		ny = height/3-1;
		nz = depth/3-1;
		for (offset=0; offset<noffsets; offset++) {
			printf("\noffset: %d\n",offset);
			fprintf(fpout,"\noffset: %d\n",offset);
			nc = 3*3*3;
			count = (int *)malloc(nc*sizeof(int));
			for (int i=0;i<nc;i++) count[i] = 0;
			ndropped = 0;
			for (ix=0;ix<nx;ix++) {
				for (iy=0;iy<ny;iy++) {
					for (iz=0;iz<nz;iz++) {
						n1 = 0;
						for (x=3*ix; x<3*(ix+1); x++) {
							for (y=3*iy; y<3*(iy+1); y++) {
								for (z=3*iz; z<3*(iz+1); z++) {
									xx = x + offsetx[offset];
									yy = y + offsety[offset];
									zz = z + offsetz[offset];
									kv = Vindex(xx,yy,zz);
									if (kv == 0) continue;
									pv = &Vlist[kv];
									if (pv->nbrs > 0) {		// try adding all lit voxels to the list.  Need to use nbrs to avoid vertices dropped in coalesce
										kvlist1[n1] = kv;
										n1++;
									}
								}
							}
						}
						count[n1]++;
						//if (n1 == 8) {
						//	fprintf(fpout,"big count: %d  ix,iy,iz: %d %d %d\n",n1,ix,iy,iz);
						//}
						// At this point, we know that there are n lit voxels in the 3x3x3 cube, with indices stored in kvlist[]
						if (n1 > 1) {
							// Need to select those that are connected, and the collapse voxel location
							group(kvlist1,n1,kvlist2,&n2);
							selectMidVoxel(kvlist2,n2,&kmid);
							err = coalesce(kvlist2,n2,kmid);
							if (err != 0) return 1;
							ndropped += n2-1;
						}
					}
				}
			}
		}
		free(count);
		printf("Dropped points: %d\n",ndropped);
	}	// end of offset loop
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Process vertices in NxNxN cube around every vertex with nbrs > 2 (N = 2*nr+1)
// This does not work as expected.  nbrs can continue to grow.
//-----------------------------------------------------------------------------------------------------
int vertexDensity3()
{
	int kv, kv0, nbrs, x, y, z, n1, n2, nc, offset, ndropped, kmid, err;
	int nx, ny, nz, ix, iy, iz, xx, yy, zz, x0, y0, z0;
	int xmin, xmax, ymin, ymax, zmin, zmax;
	int nr=1;
	int kvlist1[27], kvlist2[27];
	VOXEL *pv, *pv0;

	printf("\nvertexDensity3\n");
	ndropped = 0;
	for (kv0=1; kv0<=nlit; kv0++) {
		pv0 = &Vlist[kv0];
		nbrs = pv0->nbrs;
		if (nbrs < 3) continue;
		x0 = pv0->pos[0];
		y0 = pv0->pos[1];
		z0 = pv0->pos[2];
		xmin = MAX(0,x0-nr);
		xmax = MIN(width-1,x0+nr);
		ymin = MAX(0,y0-nr);
		ymax = MIN(height-1,y0+nr);
		zmin = MAX(0,z0-nr);
		zmax = MIN(depth-1,z0+nr);
		n1 = 0;
		for (xx=xmin; xx<=xmax; xx++) {
			for (yy=ymin; yy<=ymax; yy++) {
				for (zz=zmin; zz<=zmax; zz++) {
					kv = Vindex(xx,yy,zz);
					if (kv == 0) continue;
					pv = &Vlist[kv];
					if (pv->nbrs > 0) {		// try adding all lit voxels to the list.  Need to use nbrs to avoid vertices dropped in coalesce
						kvlist1[n1] = kv;
						n1++;
					}
				}
			}
		}
		if (n1 > 1) {
			// Need to select those that are connected, and the collapse voxel location
			group(kvlist1,n1,kvlist2,&n2);
			selectMidVoxel(kvlist2,n2,&kmid);
			nbrs = Vlist[121710].nbrs;
			printf("kv0: %d n1: %d kmid: %d nbrs of 121710: %d\n",kv0,n1,kmid,nbrs);
			err = coalesce(kvlist2,n2,kmid);
			if (err != 0) return 1;
			ndropped += n2-1;
		}
	}
	printf("Dropped points: %d\n",ndropped);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// The nj used points on an edge have index in pe->pt[] of jj[]
//-----------------------------------------------------------------------------------------------------
void saveNetwork()
{
	int ie, npts, deln, ndel, i, j, nj;
	int jj[50];
	float diam;
	EDGE *pe;
	VOXEL *pv;

	for (ie=0; ie<ne; ie++) {
		pe = &edgeList[ie];
		npts = pe->npts;
		diam = pe->segavediam;
		npts = ie+2;
		if (npts < 5) {
			deln = 1;
		} else if (npts < 10) {
			deln = 2;
		} else {
			deln = 3;
		}
		ndel = int((npts-1)/deln) + 1;
		printf("npts: %d deln: %d ndel: %d: ",npts,deln,ndel);
		nj = 0;
		for (i=0; i<ndel; i++) {
			j = i*deln;
//			printf("%d ",j);
			jj[nj] = j;
			nj++;
		}
		if (j < npts-1) {		
//			printf("%d ",npts-1);
			jj[nj] = npts-1;
			nj++;
		}
//		printf("\n");
	}
}

//-----------------------------------------------------------------------------------------------------
// Drop voxel k from the nbr list for pv
//-----------------------------------------------------------------------------------------------------
void dropnbr(int kdrop, VOXEL *pv)
{
	for (int i=0; i<pv->nbrs; i++) {
		if (pv->nbr[i] == kdrop) {
			for (int j=i; j<pv->nbrs-1;j++) {
				pv->nbr[j] = pv->nbr[j+1];
			}
			break;
		}
	}
	pv->nbrs--;
}

//-----------------------------------------------------------------------------------------------------
// Try to eliminate small loops
// First, tiny triangles: o   o
//                          o
//-----------------------------------------------------------------------------------------------------
int deloop()
{
	int k, k0, i, i1, i2, k1, k2, j, kj1, kj2, kdrop, nbrs, x0, y0, z0, x, y, z;
	int xmin, xmax, ymin, ymax, zmin, zmax;
	int nnear, nearlist[20], n3dropped, n4dropped;
	int nr = 3;		// the searched cube has side = 2*nr+1
	float d1, d2;
	bool triangle, loop4;
	VOXEL *pv, *pv0, *pv1, *pv2, *pvdrop;

	// Tiny triangles
	n3dropped = 0;
	printf("deloop: removing tiny triangles\n");
	for (k0=1; k0<=nlit; k0++) {
		pv0 = &Vlist[k0];
		nbrs = pv0->nbrs;
		if (nbrs < 3) continue;
		x0 = pv0->pos[0];
		y0 = pv0->pos[1];
		z0 = pv0->pos[2];
		nnear = 0;
		triangle = false;
		xmin = MAX(0,x0-nr);
		xmax = MIN(width-1,x0+nr);
		ymin = MAX(0,y0-nr);
		ymax = MIN(height-1,y0+nr);
		zmin = MAX(0,z0-nr);
		zmax = MIN(depth-1,z0+nr);
		// Look at lit voxels within (2*nr+1) cube 
		for (x=xmin; x<=xmax; x++) {
			for (y=ymin; y<=ymax; y++) {
				for (z=zmin; z<=zmax; z++) {
					k = Vindex(x,y,z);
					if (k != 0) {
						pv = &Vlist[k];
						if (pv->nbrs > 0) {
							// does this voxel have k0 as a nbr?
							for (i=0; i<pv->nbrs; i++) {
								if (pv->nbr[i] == k0) {
									nearlist[nnear] = k;
									nnear++;
									break;
								}
							}
						}
					}
				}
			}
		}
		// nearlist[] is the list of voxels that have k0 as a nbr
		// now need to see if any pair (k1,k2) in nearlist[] are neighbours
		for (i1=0; i1<nnear; i1++) {
			k1 = nearlist[i1];
			pv1 = &Vlist[k1];
			for (j=0; j<pv1->nbrs; j++) {
				kj1 = pv1->nbr[j];
				for (i2=i1+1; i2<nnear; i2++) {
					k2 = nearlist[i2];
					if (k2 == kj1) {
						triangle = true;
						break;
					}
				}
				if (triangle) break;
			}
			if (triangle) break;
		}
		if (triangle) {
			// k0 - k1 - k2 makes a triangle loop.  Drop the longest link of k0-k1 and k0-k2
	//		printf("k0 - k1 - k2 makes a triangle loop: %d %d %d\n",k0,k1,k2);
			pv2 = &Vlist[k2];
			d1 = dist_um(pv0->pos,pv1->pos);
			d2 = dist_um(pv0->pos,pv2->pos);
			if (d1 > d2) {
				kdrop = k1;
				pvdrop = pv1;
			} else {
				kdrop = k2;
				pvdrop = pv2;
			}
			// drop nbr kdrop from k0 nbr[]
			dropnbr(kdrop,pv0);
			// drop nbr k0 from kdrop nbr[]
			dropnbr(k0,pvdrop);
			n3dropped++;
		}
		continue;

// Now look for bigger loops
// To be safe, need to recreate nearlist[]
		nnear = 0;
		// Look at lit voxels within (2*nr+1) cube 
		for (x=xmin; x<=xmax; x++) {
			for (y=ymin; y<=ymax; y++) {
				for (z=zmin; z<=zmax; z++) {
					k = Vindex(x,y,z);
					if (k != 0) {
						pv = &Vlist[k];
						if (pv->nbrs > 0) {
							// does this voxel have k0 as a nbr?
							for (i=0; i<pv->nbrs; i++) {
								if (pv->nbr[i] == k0) {
									nearlist[nnear] = k;
									nnear++;
									break;
								}
							}
						}
					}
				}
			}
		}
// Check if two of the nearlist k1, k2 have a nbr in common (not k0)
		for (i1=0; i1<nnear; i1++) {
			k1 = nearlist[i1];
			pv1 = &Vlist[k1];
			for (int j1=0; j1<pv1->nbrs; j1++) {
				kj1 = pv1->nbr[j1];
				for (i2=i1+1; i2<nnear; i2++) {
					k2 = nearlist[i2];
					pv2 = &Vlist[k2];
					for (int j2=0; j2<pv2->nbrs; j2++) {
						kj2 = pv2->nbr[j2];
						if ((kj2 == kj1) && kj2 != k0) {
							loop4 = true;
							break;
						}
					}
					if (loop4) break;
				}
				if (loop4) break;
			}
			if (loop4) break;
		}
		if (loop4) {
			if (Vlist[k1].nbrs == 2) {	
			} else if (Vlist[k2].nbrs == 2) { 
			} else if (Vlist[kj1].nbrs == 2) { 
			} else {	// need to choose which voxel to drop - arbitrarily drop kj1
			}
		}
	}
	printf("dropped %d links to break small triangles\n",n3dropped);

	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Create mapping between Vlist and vertexList
//-----------------------------------------------------------------------------------------------------
int createVertexList()
{
	int k, ie, i, npts, iv, iv0, iv1, kp;
//	float *dsum;
//	int *nsum;
	VOXEL *pv;
	EDGE *pe;


	// Count vertices and points
	nv = 0;
	np = 0;
	for (k=1; k<=nlit; k++) {
		pv = &Vlist[k];
		if (pv->nbrs == 2) np++;
		if (pv->nbrs == 1 || pv->nbrs > 2) nv++;
	}
	printf("counted vertices: %d\n",nv);
	vertexList = (VERTEX *)malloc(nv*sizeof(VERTEX));		// this will hold voxel number corresponding to vertex number
	// Assign first pt sequence numbers to the vertices
	kp = 0;
	for (k=1; k<=nlit; k++) {
		pv = &Vlist[k];
		if (pv->nbrs == 1 || pv->nbrs > 2) {
			pv->vertex_num = kp;
			vertexList[kp].voxel_index = k;
			kp++;
		}
	}
	printf("set vertex voxel indexes\n");

	// Now assign pt numbers to the other lit voxels  NOT NEEDED
	//for (k=1; k<=nlit; k++) {
	//	pv = &Vlist[k];
	//	if (pv->nbrs == 2) {
	//		pv->kp = kp;
	//		kp++;
	//	}
	//}
	//printf("kp: %d  np: %d\n",kp,np);

	// Set edge vert[] to their kp sequence numbers
	for (ie=0;ie<ne;ie++) {
		pe = &edgeList[ie];
		npts = pe->npts;
		iv0 = Vlist[pe->pt[0]].vertex_num;
		pe->vert[0] = iv0;
		int iv1 = Vlist[pe->pt[npts-1]].vertex_num;
		pe->vert[1] = iv1;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Need to adjust diameters near junctions, because the diameter estimation method does not work well
// where vessels come together.
// Method:
// For each edge:
//    Find the points that are less than alpha (e.g. 1.5) diameters from each end, set edge.avediameter[k] = 0
// For each edge:
//    Set zero-diameter points to diameter of nearest point on the edge with non-zero diameter (if it exists)
//    except at end points - here set to max of this and current value
// For each edge:
//    If 
//-----------------------------------------------------------------------------------------------------
int FixDiameters_topo()
{
	EDGE *edge;
	int ie, k, kp0, kp1, kp, npts, n0, n1, ipass, nzero, iv, err;
	double d0, d1, diam, d2ave, dave, vsum, dmax, dmin;
	double len, len0, len1, lsum;
	bool done;
	double alpha = 0.3;
	bool FIX_JUNCTIONS = true;

	printf("FixDiameters\n");
	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
		npts = edge->npts;
		kp0 = edge->pt[0];
		kp1 = edge->pt[npts-1];
		for (k=1; k<npts-1; k++) {
			kp = edge->pt[k];
			diam = avediameter[kp];
//			d0 = dist_um(kp,kp0);
//			d1 = dist_um(kp,kp1);
			d0 = dist_um(Vlist[kp].initial_pos,Vlist[kp0].initial_pos);
			d1 = dist_um(Vlist[kp].initial_pos,Vlist[kp1].initial_pos);
			if (d0 < alpha*diam || d1 < alpha*diam) {
				avediameter[kp] = 0;
			}
		}
	}
	printf("did stage 1\n");
	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
		npts = edge->npts;
		kp0 = edge->pt[0];
		kp1 = edge->pt[npts-1];
		d0 = 0;
		n0 = 0;
		for (k=1; k<npts-1; k++) {
			kp = edge->pt[k];
			d0 = avediameter[kp];
			if (d0 > 0) break;
			n0 = k;
		}
		d1 = 0;
		n1 = npts-1;
		for (k=npts-2; k>0; k--) {
			kp = edge->pt[k];
			d1 = avediameter[kp];
			if (d1 > 0) break;
			n1 = k;
		}
		if (d0 == 0) {
			if (d1 != 0) {
				printf("d1 != 0\n");
				return 1;
			}
			continue;
		}
		for (k=1; k<=n0; k++) {
			kp = edge->pt[k];
			avediameter[kp] = d0;
		}
		avediameter[kp0] = MAX(avediameter[kp0],d0);
		for (k=n1; k<npts-1; k++) {
			kp = edge->pt[k];
			avediameter[kp] = d1;
		}
		avediameter[kp1] = MAX(avediameter[kp1],d1);
	}
	printf("did stage 2\n");
	// There may still remain some edges with unset diameters
	ipass = 0;
	for (;;) {
		printf("pass: %d\n",ipass);
		ipass++;
		done = true;
		for (ie=0; ie<ne; ie++) {
			edge = &edgeList[ie];
			if (!edge->used) continue;
			npts = edge->npts;
			kp0 = edge->pt[0];
			kp1 = edge->pt[npts-1];
			d0 = avediameter[kp0];
			d1 = avediameter[kp1];
			if (d0 == 0 && d1 > 0) {
				done = false;
				for (k=0; k<npts; k++) {
					kp = edge->pt[k];
					avediameter[kp] = d1;
				}
			}
			if (d1 == 0 && d0 > 0) {
				done = false;
				for (k=0; k<npts; k++) {
					kp = edge->pt[k];
					avediameter[kp] = d0;
				}
			}
		}
		if (done) break;
	}
	printf("did stage 3\n");
	nzero = 0;
	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
		npts = edge->npts;
		err = 0;
		for (k=0; k<npts; k++) {
			kp = edge->pt[k];
			if (avediameter[kp] == 0) err = 1;
		}
		if (err != 0) {
			nzero++;
			n0 = 0;
			for (k=0; k<npts; k++) {
				kp = edge->pt[k];
				if (avediameter[kp] == 0) {
					n0 = k-1;
					break;
				}
				d0 = avediameter[kp];
			}
			n1 = npts-1;
			for (k=npts-1; k>=0; k--) {
				kp = edge->pt[k];
				if (avediameter[kp] == 0) {
					n1 = k+1;
					break;
				}
				d1 = avediameter[kp];
			}
			for (k=n0+1; k<n1; k++) {
				kp = edge->pt[k];
				avediameter[kp] = d0 + (d1-d0)*(float)(k-n0)/(n1-n0);
			}
//			exit(1);
		}
	}
	// set uniform (average) diameter along each edge
	// the total vessel volume is conserved, and end pt diameters are not changed
	// interior point diameters are all set to the same value = dave = sqrt(d2ave)
	// where: (d0*d0+d2ave)*len0/2 + d2ave*(lsum - len0 - len1) + (d1*d1+d2ave)*len1/2 = vsum
	// d2ave = (vsum - d0*d0*len0/2 - d1*d1*len1/2)/(lsum - len0/2 - len1/2)
	// d0 = start pt diameter, len0 = start segment length,
	// d1 = end pt diameter, len1 = end segment length,
	bool uniform_diameter = true;
	if (uniform_diameter) {	
		printf("Creating uniform edge diameters\n");
		for (ie=0; ie<ne; ie++) {
			edge = &edgeList[ie];
			if (!edge->used) continue;
			npts = edge->npts;
			if (npts < 4) continue;
			lsum = 0;
			vsum = 0;
			for (k=0; k<npts-1; k++) {
				kp0 = edge->pt[k];
				kp1 = edge->pt[k+1];
				d0 = avediameter[kp0];
				d1 = avediameter[kp1];
//				len = zdist(kp0,kp1);
				len = dist_um(Vlist[kp0].initial_pos,Vlist[kp1].initial_pos);
				lsum += len;
				vsum += len*(d0*d0 + d1*d1)/2;
				if (k == 0) len0 = len;
				if (k == npts-2) len1 = len;
			}
			kp0 = edge->pt[0];
			d0 = avediameter[kp0];
			kp1 = edge->pt[npts-1];
			d1 = avediameter[kp1];
			d2ave = (vsum - d0*d0*len0/2 - d1*d1*len1/2)/(lsum - len0/2 - len1/2);
			dave = sqrt(d2ave);
			for (k=1; k<npts-1; k++) {
				kp = edge->pt[k];
				avediameter[kp] = dave;
			}
		}
	}
	if (!FIX_JUNCTIONS) return 0;
	// Finally, need to fix junction node diameters.
	// At each junction, find the connected edges. Set the point avediameter to the maximum 
	// or average of the computed diameters of the connected edges.
#if 0
	// Need to identify the point corresponding to the vertex!
	// This is edge->pt[0] if iv == edge->vert[0], edge->pt[npts-1] if iv = edge->vert[1]
	for (iv=0; iv<nv; iv++) {
//		if (!vertex[iv].used || vertex[iv].nlinks_used == 0) continue;
		dmin = 1.0e10;
		dmax = 0;
		dave = 0;
		n1 = 0;
		for (k=0; k<vertexList[iv].nlinks; k++) {
			ie = vertexList[iv].edge[k];
			edge = &edgeList[ie];
			if (!edge->used) continue;
			n1++;
			npts = edge->npts;
			if (iv == edge->vert[0]) {
				kp = edgeList[ie].pt[1];
			} else if (iv == edge->vert[1]) {
				kp = edgeList[ie].pt[npts-2];
			} else {
				printf("Error: FixDiameters: iv: %d vert: %d %d\n",iv,edge->vert[0],edge->vert[1]);
				return 1;
			}
			if (avediameter[kp] < dmin) dmin = avediameter[kp];
			if (avediameter[kp] > dmax) dmax = avediameter[kp];
			dave += avediameter[kp];
		}
		dave /= n1;
//		if (dave == 0) continue;
		for (k=0; k<vertexList[iv].nlinks; k++) {
			ie = vertexList[iv].edge[k];
			edge = &edgeList[ie];
			if (!edge->used) continue;
			npts = edge->npts;
			if (iv == edge->vert[0]) {
				kp = edgeList[ie].pt[0];
			} else {
				kp = edgeList[ie].pt[npts-1];
			}
//			printf("iv: %6d  %3d  %6d %6d %6.1f %6.1f\n",iv,n1,k,kp,avediameter[kp],dmax);
			//if (junction_max || dmax/dmin > max_ratio) {
			//	avediameter[kp] = dmax;
			//} else {
			//	avediameter[kp] = dave;
			//}
			avediameter[kp] = dave;
		}
	}
#endif

// Instead, look at edges
	EDGE *pe;
	for (ie=0; ie<ne; ie++) {
		pe = &edgeList[ie];
		npts = pe->npts;
		for (iv=0; iv<2; iv++) {
			if (iv == 0)
				kp0 = pe->pt[0];
			else
				kp0 = pe->pt[npts-1];
			if (avediameter[kp0] > 0) continue;
			int nbrs = Vlist[kp0].nbrs;
			dave = 0;
			for (int i=0; i<nbrs; i++) {
				kp1 = Vlist[kp0].nbr[i];
				if (Vlist[kp1].nbrs != 2) {
					double p0[3],p2[3],r2_ave,r2_min,d;
					bool zero;
					for (int i=0; i<3; i++) {
						p0[i] = vsize[i]*(Vlist[kp0].pos[i] + 0.5);		// Note: centres of voxel cubes
						p2[i] = vsize[i]*(Vlist[kp1].pos[i] + 0.5);
					}
					// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
					EstimateDiameter(p0,p0,p2,&r2_ave,&r2_min,&zero);
					d = 2*sqrt(r2_ave);
					dave += d;
				} else {
					dave += avediameter[kp1];
				}
			}
			avediameter[kp0] = dave/nbrs;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Looking for edges to compare with topo to see why volume is so different
//-----------------------------------------------------------------------------------------------------
void showEdges()
{
	EDGE *pe;
	int iv0, iv1, *pos0, *pos1;
	double ad;
	for (int ie=0;ie<ne;ie++) {
		pe = &edgeList[ie];
		ad = pe->segavediam;
		iv0 = pe->vert[0];
		iv1 = pe->vert[1];
		pos0 = Vlist[vertexList[iv0].voxel_index].pos;
		pos1 = Vlist[vertexList[iv1].voxel_index].pos;
		if (pos0[0] == 33 && (pos0[1] == 41 || pos0[1] == 42) && ad > 14) {
			fprintf(fpout,"ie: %d vert: %d %d ad: %f\n",ie,iv0,iv1,ad);
			fprintf(fpout,"%d %d %d   %d %d %d\n",pos0[0],pos0[1],pos0[2],pos1[0],pos1[1],pos1[2]);
		}
		if (pos1[0] == 33 && (pos1[1] == 41 || pos1[1] == 42) && ad > 14) {
			fprintf(fpout,"ie: %d vert: %d %d ad: %f\n",ie,iv0,iv1,ad);
			fprintf(fpout,"%d %d %d   %d %d %d\n",pos0[0],pos0[1],pos0[2],pos1[0],pos1[1],pos1[2]);
		}
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	float voxelsize_x, voxelsize_y, voxelsize_z;
	char *vessFile, *skelFile;
	int count, err;
	char *outfilename;
	char drive[32], dir[256],filename[256], ext[32];
	char amfilename[256];
	char output_basename[256], errfilename[256];
	char *temp_name;
	int noffsets = 1;
	int cmgui_flag, diam_flag;
	float origin_shift[3] = {0, 0, 0};

	if (argc != 12) {
		printf("Usage: topo2 skel_tiff object_tiff output_file voxelsize_x voxelsize_y voxelsize_z fixed_diam_flag fixed_diam len_limit min_end_len max_diam\n");
		fperr = fopen("topo2_error.log","w");
		fprintf(fperr,"Usage: topo2 skel_tiff object_tiff output_file voxelsize_x voxelsize_y voxelsize_z fixed_diam_flag fixed_diam len_limit min_end_len max_diam\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	outfilename = argv[3];

#if (defined (_WIN32) || defined (_WIN64))
    // windows code
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
#elif (defined (LINUX) || defined (__linux__))
    // linux code
	strcpy(output_basename, basename(outfilename));
	temp_name = remove(output_basename);
	strcpy(output_basename,temp_name);
#endif
	sprintf(errfilename,"%s%s",output_basename,"_topo2.log");
	fperr = fopen(errfilename,"w");

	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
	fprintf(fperr,"Basename: %s\n",output_basename);
	sprintf(amfilename,"%s.am",output_basename);

	skelFile = argv[1];
	vessFile = argv[2];
	sscanf(argv[4],"%f",&voxelsize_x);
	sscanf(argv[5],"%f",&voxelsize_y);
	sscanf(argv[6],"%f",&voxelsize_z);
	sscanf(argv[7],"%d",&diam_flag);
	sscanf(argv[8],"%f",&FIXED_DIAMETER);
	sscanf(argv[9],"%f",&len_limit);
	sscanf(argv[10],"%f",&min_end_len);
	sscanf(argv[11],"%f",&max_diameter);

	fixed_diam_flag = (diam_flag == 1);
	use_max_diameter = (max_diameter > 0);
	vsize[0] = voxelsize_x;
	vsize[1] = voxelsize_y;
	vsize[2] = voxelsize_z;
	cmgui_flag = 1;

	fpout = fopen(outfilename,"w");
	printf("Input skeleton file: %s\n",skelFile);
	fprintf(fpout,"Input skeleton file: %s\n",skelFile);

	printf("Input object file: %s\n",vessFile);
	fprintf(fpout,"Input object file: %s\n",vessFile);
	if (vessFile[0] == '-') {
		fixed_diam_flag = true;
		if (FIXED_DIAMETER == 0) FIXED_DIAMETER = 1;
		printf("No object file is supplied, no diameters will be estimated\n");
		fprintf(fpout,"No object file is supplied, no diameters will be estimated\n");
		use_object = false;
	} else {
		use_object = true;
	}
	
	printf("Output basename: %s\n",output_basename);
	fprintf(fpout,"Output basename: %s\n",output_basename);
	printf("Voxel size: x,y,z: %f %f %f\n",voxelsize_x, voxelsize_y,voxelsize_z);
	fprintf(fpout,"Voxel size: x,y,z: %f %f %f\n",voxelsize_x, voxelsize_y,voxelsize_z);
	fprintf(fpout,"Distribution length limit: %6.1f  Minimum dead-end length: %6.1f\n",len_limit,min_end_len);
	printf("Fixed diameter flag: %d\n",fixed_diam_flag);		
	fprintf(fperr,"Fixed diameter flag: %d\n",fixed_diam_flag);
	if (fixed_diam_flag) {
		printf("Using a fixed vessel diameter: %f\n",FIXED_DIAMETER);
		fprintf(fpout,"Using a fixed vessel diameter: %f\n",FIXED_DIAMETER);
	} else if (FRC_fix) {
		printf("\nNOTE: This version of the program is applicable to FRC networks.\n");
		printf("A small correction is made to the estimate of vessel radius (-0.5 um)\n");
		printf("Accurate diameter estimation is impossible with fewer than 5 pixels per diameter\n\n");
		fprintf(fpout,"\nNOTE: This version of the program is applicable to FRC networks.\n");
		fprintf(fpout,"A small correction is made to the estimate of vessel radius (-0.5 um)\n");
		fprintf(fpout,"Accurate diameter estimation is impossible with fewer than 5 pixels per diameter\n\n");
	}
	fprintf(fpout,"Amira file: %s\n",amfilename);
	fflush(fpout);

	net = (NETWORK *)malloc(sizeof(NETWORK));

	typedef itk::ImageFileReader<ImageType> FileReaderType;

	// Read skeleton file
	{
		FileReaderType::Pointer readerskel = FileReaderType::New();
		readerskel->SetFileName(skelFile);
		try
		{
			printf("Reading input skeleton file: %s\n",skelFile);
			readerskel->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cout << e << std::endl;
			fprintf(fperr,"Read error on skeleton file\n");
			fclose(fperr);
			return 2;	// Read error on input file
		}

		imskel = readerskel->GetOutput();
	}

	width = imskel->GetLargestPossibleRegion().GetSize()[0];
	height = imskel->GetLargestPossibleRegion().GetSize()[1];
	depth = imskel->GetLargestPossibleRegion().GetSize()[2];
	imsize_xy = width*height;

	printf("Image dimensions: width, height, depth: %lld %lld %lld\n",width,height,depth);
	pskel = (unsigned char *)(imskel->GetBufferPointer());
	fflush(fpout);

	nlit = 0;
	for (long long i=0; i<width*height*depth; i++) {
		if (pskel[i] > 0) nlit++;
	}
	// nlit = number of lit voxels in the skeleton image
	np = nlit;
	net->np = nlit;

	pindex = (int *)malloc(width*height*depth*sizeof(int));		// maps location (i,j,k) to a skeleton lit voxel index (k>0) or 0 (unlit)
	net->point = (VOXEL *)malloc((nlit+1)*sizeof(VOXEL));
	Vlist = net->point;

	printf("Number of lit skeleton voxels: %d\n",nlit);

	printf("Creating the lit voxel index list\n");
	fflush(fpout);
	err = createVlist();
	if (err != 0) {
		printf("Error in createVlist\n");
		fprintf(fperr,"Error in createVlist\n");
		return 5;
	} else {
		printf("Created Vlist\n");
	}


	if (use_object) {
		// Read original object file (binary)
		FileReaderType::Pointer reader = FileReaderType::New();
		reader->SetFileName(vessFile);
		try
		{
			printf("Reading input object file: %s\n",vessFile);
			reader->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cout << e << std::endl;
			fprintf(fperr,"Read error on object file\n");
			fclose(fperr);
			return 3;	// Read error on input file
		}

		im = reader->GetOutput();

		int w = im->GetLargestPossibleRegion().GetSize()[0];
		int h = im->GetLargestPossibleRegion().GetSize()[1];
		int d = im->GetLargestPossibleRegion().GetSize()[2];

		if (w != width || h != height || d != depth) {
			printf("Error: skeleton and object files differ in size\n");
			fprintf(fperr,"Error: skeleton and object files differ in size\n");
			fclose(fperr);
			return 4;
		}

		p = (unsigned char *)(im->GetBufferPointer());
		int nlit_obj = 0;
		for (long long i=0; i<width*height*depth; i++) {
			if (p[i] > 0) nlit_obj++;
		}
		// nlit_obj = number of lit voxels in the object image
		printf("Number of lit object voxels: %d\n",nlit_obj);
		fprintf(fpout,"Number of lit object voxels: %d\n",nlit_obj);
		fflush(fpout);
	}

	if (METHOD2) {
		err = deloop();
		if (err != 0) {
			printf("Error in deloop (method2)\n");
			fprintf(fperr,"Error in deloop (method2)\n");
			return 7;
		} else {
			printf("Eliminated small triangular loops\n");
		}
		err = traceSegments();
		if (err != 0) {
			printf("Error in traceSegments (method2)\n");
			fprintf(fperr,"Error in traceSegments (method2)\n");
			return 8;
		} else {
			printf("Traced segments (method2): ne: %d\n",ne);
		}
		printf("Did traceSegments\n");
		err = getPointDiameters();
		if (err != 0) {
			printf("Error in getPointDiameters\n");
			fprintf(fperr,"Error in getPointDiameters\n");
			return 9;
		} else {
			printf("Estimated point diameters\n");
		}
		freeEdgeList();
	}

//	fprintf(fpout,"DROP vertexDensity2 for TESTING\n");
	err = vertexDensity2(noffsets);
	if (err != 0) {
		printf("Error in vertexDensity2\n");
		fprintf(fperr,"Error in vertexDensity2\n");
		return 6;
	} else {
		printf("Coalesced lit voxels\n");
	}
	fflush(fpout);

//	fprintf(fpout,"DROP deloop for TESTING\n");
	err = deloop();
	if (err != 0) {
		printf("Error in deloop\n");
		fprintf(fperr,"Error in deloop\n");
		return 7;
	} else {
		printf("Eliminated small triangular loops\n");
	}
	fflush(fpout);

	err = traceSegments();
	if (err != 0) {
		printf("Error in traceSegments\n");
		fprintf(fperr,"Error in traceSegments\n");
		return 8;
	} else {
		printf("Traced segments\n");
	}

	fprintf(fpout,"Traced segments\n");
//	fprintf(fpout,"Early exit\n");
	fflush(fpout);
//	return 1;

	if (use_pt_diameters) {
		fprintf(fpout,"getDiameters_topo\n");
		fflush(fpout);
		err = getDiameters_topo();
		if (err != 0) {
			printf("Error in getDiameters_topo\n");
			fprintf(fperr,"Error in getDiameters_topo\n");
			return 9;
		}

		fprintf(fpout,"did getDiameters_topo, now fixDiameters_topo\n");
//		fprintf(fpout,"Early exit\n");
		fflush(fpout);
//		return 0;

		err = FixDiameters_topo();
		if (err != 0) {
			printf("Error: FixDiameters_topo\n");
			fprintf(fperr,"Error: FixDiameters_topo\n");
			fclose(fperr);
			return 12;
		}
		fprintf(fpout,"did FixDiameters_topo\n");
		fflush(fpout);
		for (int k=1; k<=nlit; k++) {
			Vlist[k].diameter = avediameter[k];
		}
		fprintf(fpout,"Set Vlist diameters\n");
		fflush(fpout);
	} else {
		err = getDiameters();	
		if (err != 0) {
			printf("Error in getDiameters\n");
			fprintf(fperr,"Error in getDiameters\n");
			return 9;
		}
	}
	printf("Estimated segment lengths and average diameters\n");
	fprintf(fpout,"Estimated segment lengths and average diameters\n");
	fflush(fpout);

	printf("Create vertex list\n");
	fprintf(fpout,"Create vertex list\n");
	fflush(fpout);
	createVertexList();
	printf("Created vertex list\n");
	fprintf(fpout,"Created vertex list\n");
	fflush(fpout);

	int nused = 0;
	for (int k=1; k<=nlit; k++) {
		if (Vlist[k].nbrs > 0) nused++;
	}
	printf("Number of skeleton voxels used: %d\n",nused);
	fprintf(fpout,"Number of skeleton voxels used: %d\n",nused);
	fflush(fpout);

	if (use_pt_diameters) {
		fprintf(fpout,"CreateDistributions_topo\n");
		fflush(fpout);
		err = CreateDistributions_topo();		// scaling for voxelsize now done in the distance calculations
		fprintf(fpout,"did CreateDistributions_topo: err: %d\n",err);
		fflush(fpout);
	} else
		err = CreateDistributions();		// scaling for voxelsize now done in the distance calculations
	if (err != 0) {
		printf("Error: CreateDistributions\n");
		fprintf(fperr,"Error: CreateDistributions\n");
		fclose(fperr);
		return 12;
	} else {
		printf("Created distributions of length and diameter\n");
		fprintf(fpout,"Created distributions of length and diameter\n");
		fflush(fpout);
	}
	fflush(fpout);


	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,net,origin_shift);
		if (err != 0) {
			printf("Error in WriteCmguiData\n");
			fprintf(fperr,"Error in WriteCmguiData\n");
			return 10;
		} else {
			printf("Wrote CMGUI files\n");
		}
	}
	fflush(fpout);
	err = WriteAmiraFile(amfilename, vessFile, skelFile, origin_shift);
	if (err != 0) {
		printf("Error in WriteAmiraFile\n");
		fprintf(fpout,"Error in WriteAmiraFile\n");
		return 11;
	} else {
		printf("Wrote AMIRA file\n");
	}
	fflush(fpout);

/* Moved above
	if (use_pt_diameters) 
		err = CreateDistributions_topo();		// scaling for voxelsize now done in the distance calculations
	else
		err = CreateDistributions();		// scaling for voxelsize now done in the distance calculations
	if (err != 0) {
		printf("Error: CreateDistributions\n");
		fprintf(fperr,"Error: CreateDistributions\n");
		fclose(fperr);
		return 12;
	} else {
		printf("Created distributions of length and diameter\n");
	}
*/

//	showEdges();

	printf("Terminated normally\n");
	fprintf(fperr,"Terminated normally\n");
	fclose(fperr);
	fclose(fpout);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Create Amira SpatialGraph file
// This assumes that all edges are used, i.e. the network has been squeezed
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile(char *amOutFile, char *vessFile, char *skelFile, float *origin_shift)
{
	int k, j, ie, iv, npts, npts_used;
	int *pos;
	VOXEL *pv;

	printf("\nWriteAmiraFile: %s\n",amOutFile);
	fprintf(fpout,"\nWriteAmiraFile: %s\n",amOutFile);


	npts = 0;
//	npts_used = 0;
	for (ie=0;ie<net->ne;ie++) {
//		if (!edgeList[i].used) continue;
		npts += edgeList[ie].npts;
//		npts_used += edgeList[i].npts_used;
	}
	npts_used = npts;

	FILE *fpam = fopen(amOutFile,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
	fprintf(fpam,"# Created by topo2.exe from: %s and %s\n",vessFile,skelFile);
	fprintf(fpam,"\n");
	fprintf(fpam,"define VERTEX %d\n",nv);
	fprintf(fpam,"define EDGE %d\n",ne_used);		// ne -> ne_used
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
	for (iv=0;iv<nv;iv++) {
		pos = Vlist[vertexList[iv].voxel_index].pos;
		fprintf(fpam,"%6.1f %6.1f %6.1f\n",
			vsize[0]*pos[0] - origin_shift[0],
			vsize[1]*pos[1] - origin_shift[1],
			vsize[2]*pos[2] - origin_shift[2]);
	}
	fprintf(fpam,"\n@2\n");
	for (ie=0;ie<ne;ie++) {
		if (!edgeList[ie].used) continue;
		fprintf(fpam,"%d %d\n",edgeList[ie].vert[0],edgeList[ie].vert[1]);
		if (edgeList[ie].vert[0]== edgeList[ie].vert[1]) {
			fprintf(fperr,"Error: WriteAmiraFile: repeated vertices: i: %d  %d\n",ie,edgeList[ie].vert[0]);
		}
	}
	fprintf(fpam,"\n@3\n");
	for (ie=0;ie<ne;ie++) {
		if (!edgeList[ie].used) continue;
		fprintf(fpam,"%d\n",edgeList[ie].npts);
//		fprintf(fpam,"%d\n",edgeList[i].npts_used);
	}
	fprintf(fpam,"\n@4\n");
	for (ie=0;ie<ne;ie++) {
		if (!edgeList[ie].used) continue;
		for (k=0;k<edgeList[ie].npts;k++) {
			j = edgeList[ie].pt[k];
			pv = &Vlist[j];
			fprintf(fpam,"%6.1f %6.1f %6.1f\n",
				vsize[0]*pv->pos[0] - origin_shift[0],
				vsize[1]*pv->pos[1] - origin_shift[1],
				vsize[2]*pv->pos[2] - origin_shift[2]);
		}
	}
	fprintf(fpam,"\n@5\n");
	double diam;
	for (ie=0;ie<ne;ie++) {
		if (!edgeList[ie].used) continue;
		if (!use_object || fixed_diam_flag) {
			diam = FIXED_DIAMETER;
		} else {
			diam = edgeList[ie].segavediam;
			if (diam < 1.0) {
				fprintf(fperr,"d < 1.0: edge: %d npts: %d pt: %d %d  %f\n",ie,edgeList[ie].npts,k,j,diam);
				diam = 1.0;
			}
		}
		for (k=0;k<edgeList[ie].npts;k++) {
			fprintf(fpam,"%6.2f\n",diam);
		}
	}
	fclose(fpam);
	return 0;
}

void freeEdgeList()
{
	EDGE *pe;
	for (int ie=0; ie<ne; ie++) {
		pe = &edgeList[ie];
		free(pe->pt);
	}
	free(edgeList);
}

//--------------------------------------------------------------------------------------------------------
// For METHOD2
// This uses the traced segments to determine point diameters in a sequential procedure.
// First, diameters are estimated (in the usual way with getDiameter()) for interior points on each edge.
// Then other points are filled in by an iterative procedure.
//--------------------------------------------------------------------------------------------------------
int getPointDiameters()
{
	EDGE *pe;
	VOXEL *pv;
	int npin, npts, ie, k, kp0, kp1, kp2, i, kdiam, kdiam2, n2, nsum;
	double p1[3], p2[3], p0[3];
	double r2_ave, r2_min, dpave, dpsum, deave, desum, dave, dsum, len, dlen, lensum, area, dvol, volsum, d;
	bool zero;

	printf("getPointDiameters: nlit: %d\n",nlit);
	npin = 0;
	for (k=1; k<=nlit; k++) {
		if (Vlist[k].nbrs > 0) {
			npin++;
			Vlist[k].diameter = 0;
		}
	}
	printf("excluding short ends: npin: %d\n",npin);
	kdiam = 0;
	n2 = 0;
	dpsum = 0;
	lensum = 0;
	volsum = 0;
	for (ie=0; ie<ne; ie++) {
		pe = &edgeList[ie];
		npts = pe->npts;
		desum = 0;
		len = 0;
		if (npts < 3) {
			n2++;
			kp0 = pe->pt[0];
			kp1 = pe->pt[1];
			len = dist_um(Vlist[kp0].initial_pos,Vlist[kp1].initial_pos);
			continue;
		}
		for (k=1;k<npts-1;k++) {
			kp0 = pe->pt[k-1];
			kp1 = pe->pt[k];
			kp2 = pe->pt[k+1];
			// Basic check for insideness
			pv = &Vlist[kp1];
			if (V3D(pv->pos[0],pv->pos[1],pv->pos[2]) == 0) {
				printf("In getDiameters: edge: %d voxel: %d is not inside the object\n",ie,kp1);
				fprintf(fpout,"In getDiameters: edge: %d voxel: %d is not inside the object\n",ie,kp1);
				exit(1);
			}
			for (i=0; i<3; i++) {
				p0[i] = vsize[i]*(Vlist[kp0].initial_pos[i] + 0.5);		// Note: these pos are guaranteed to fall within the object image
				p1[i] = vsize[i]*(Vlist[kp1].initial_pos[i] + 0.5);
				p2[i] = vsize[i]*(Vlist[kp2].initial_pos[i] + 0.5);
			}
			// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
			EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min,&zero);
			if (zero) {
				printf("getPointDiameters: got a zero for ie: %d  %d %d %d\n",ie,kp0,kp1,kp2);
				fprintf(fpout,"getPointDiameters: got a zero for ie: %d  %d %d %d\n",ie,kp0,kp1,kp2);
				exit(1);
			}
			pv->diameter = 2*sqrt(r2_ave);
			desum += pv->diameter;
			dpsum += pv->diameter;
			kdiam++;
			dlen = dist_um(Vlist[kp0].initial_pos,Vlist[kp1].initial_pos);
			len += dlen;
			if (k == npts-2) {
				len += dist_um(Vlist[kp1].initial_pos,Vlist[kp2].initial_pos);
			}
		}
		deave = desum/(npts-2);
		area = PI*pow(deave/2,2);
		dvol = area*len;
		volsum += dvol;
		lensum += len;
	}
	dpave = dpsum/kdiam;
	printf("Stage 1: pt diam for %d out of %d points, average: %f\n",kdiam,npin,dpave);
	fprintf(fpout,"Stage 1: pt diam for %d out of %d points, average: %f\n",kdiam,npin,dpave);
	printf("approx total vol (from ave edge diam and len): %f\n",volsum);
	fprintf(fpout,"BAD approx total vol (from ave edge diam and len): %f\n",volsum);


// Now need to fill in missing pts
	for (int iter=0; iter<20; iter++) {
		kdiam2 = 0;
		for (k=1; k<=nlit; k++) {
			pv = &Vlist[k];
			if (pv->nbrs == 0 || pv->diameter > 0) continue;
			dsum = 0;
			int nsum = 0;
			for (i=0; i<pv->nbrs; i++) {
				int kk = pv->nbr[i];
				if (Vlist[kk].diameter > 0) {
					nsum++;
					dsum += Vlist[kk].diameter;
				}
			}
			if (nsum > 0) {
				kdiam2++;
				pv->diameter = dsum/nsum;
				dpsum += pv->diameter;
			}
		}
		kdiam += kdiam2;
		printf("Stage 2: iteration: %d filled in %d pt diams\n",iter,kdiam2);
		printf("missing diams at %d pts\n",npin-kdiam);
		if (kdiam2 == 0) break;
	}
// Stage 3: check edges
	kdiam2 = 0;
	for (ie=0; ie<ne; ie++) {
		pe = &edgeList[ie];
		npts = pe->npts;
		desum = 0;
		if (npts == 2) {
			kp0 = pe->pt[0];
			kp1 = pe->pt[1];
			if (Vlist[kp0].diameter == 0 && Vlist[kp1].diameter == 0) {
				for (i=0; i<3; i++) {
					p0[i] = vsize[i]*(Vlist[kp0].initial_pos[i] + 0.5);		// Note: these pos are not guaranteed to fall within the object image
					p1[i] = vsize[i]*(Vlist[kp1].initial_pos[i] + 0.5);
					p2[i] = 2*p1[i] - p0[i];
				}
				EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min,&zero);
				if (zero) {
					printf("getPointDiameters: got a zero for ie: %d  npts=2: %d %d\n",ie,kp0,kp2);
					exit(1);
				}
				d = 2*sqrt(r2_ave);
				Vlist[kp0].diameter = d;
				Vlist[kp1].diameter = d;
				dpsum += 2*pv->diameter;
				kdiam2 += 2;
			} else if (Vlist[kp0].diameter > 0 && Vlist[kp1].diameter == 0) {
				Vlist[kp1].diameter = Vlist[kp0].diameter;
				kdiam2++;
			} else if (Vlist[kp1].diameter > 0 && Vlist[kp0].diameter == 0) {
				Vlist[kp0].diameter = Vlist[kp1].diameter;
				kdiam2++;
			}
		} else {
			if (Vlist[pe->pt[0]].diameter == 0) {
				kdiam2++;
				Vlist[pe->pt[0]].diameter = Vlist[pe->pt[1]].diameter;
				dpsum += Vlist[pe->pt[0]].diameter;
			}
			if (Vlist[pe->pt[npts-1]].diameter == 0) {
				kdiam2++;
				Vlist[pe->pt[npts-1]].diameter = Vlist[pe->pt[npts-2]].diameter;
				dpsum += Vlist[pe->pt[npts-1]].diameter;
			}
		}
	}
	kdiam += kdiam2;
	printf("Stage 3: filled in %d pt diams\n",kdiam2);
	printf("missing diams at %d pts\n",npin-kdiam);

// Are there pts not on an edge?
	int nzero = 0;
	for (ie=0; ie<ne; ie++) {
		pe = &edgeList[ie];
		npts = pe->npts;
		for (i=0; i<npts; i++) {
			k = pe->pt[i];
			if (Vlist[k].diameter == 0) nzero++;
		}
	}
	printf("zero diam pts on edges: %d\n",nzero);
// Are there zero diam pts somewhere?
	nzero = 0;
	for (k=1; k<=nlit; k++) {
		if (Vlist[k].nbrs == 0) continue;
		if (Vlist[k].diameter == 0) nzero++;
	}
	printf("zero diam pts: %d\n",nzero);
	printf("Average pt diameter: %f\n",dpsum/npin);
	return 0;
}

