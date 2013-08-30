// To quantify vessel density per area, and compute tissue volume

//#include <cstdio>
//#include <vector>

//#include <algorithm>
#include <math.h>
//#include <string.h>
//#include <string>
//#include <sstream>

#include "MainWindow.h"
#include "quantify.h"

/*
#include <bitset>
#include <stdio.h>

#include "network.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

struct segment_str {
	POINT end1, end2;
};
typedef segment_str SEGMENT_TYPE;

#define STR_LEN 128
#define BIG 1.0e6

FILE *fperr=NULL, *fpout=NULL;

int nxc, nyc, nzc, nx8, ny8, nz8, nbytes;
float voxelsize[3];
unsigned char *closedata = NULL;
int nsegments;
SEGMENT_TYPE *segment = NULL;
NETWORK *NP0 = NULL;

bool is_setup = false;
*/

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int MainWindow::ReadAmiraFile(char *amFile, NETWORK *net)
{
	int i, j, k, kp, npts;
	int np_used, ne_used;
	EDGE edge;
	char line[STR_LEN];

	printf("ReadAmiraFile: %s\n",amFile);
	fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
	FILE *fpam = fopen(amFile,"r");

	nsegments = 0;
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
					net->edgeList[i].pt[0] = net->edgeList[i].vert[0];							// This was not really necessary -
					net->edgeList[i].pt[net->edgeList[i].npts-1] = net->edgeList[i].vert[1];	// the point coordinates are repeated in @4
					nsegments += net->edgeList[i].npts - 1;
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
						if (k > 0 && k<edge.npts-1) {											// See note above
							sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[k] = kp;
							net->edgeList[i].pt_used[k] = kp;
							kp++;
						}
//						if (k > 0) {
//							len = len + dist(net,net->edgeList[i].pt[k-1],net->edgeList[i].pt[k]);
//						}
					}
//					net->edgeList[i].length_vox = len;
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
						if (j == 13) {
							printf("j=13: d: %f\n",net->point[j].d);
						}
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
	
	printf("nsegments: %d\n",nsegments);

	segment = (SEGMENT_TYPE *)malloc(nsegments*sizeof(SEGMENT_TYPE));

	int ip;
	int iseg = 0;
	for  (int ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
		for (i=0; i<edge.npts-1; i++) {
			ip = edge.pt[i];
			segment[iseg].end1 = net->point[ip];
			segment[iseg].end2 = net->point[ip+1];
            segment[iseg].diam = (net->point[ip].d + net->point[ip+1].d)/2;
			iseg++;
		}
	}
	nsegments = iseg;
	printf("nsegments: %d\n",nsegments);

	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int MainWindow::ReadCloseFile(char *filename)
{
	FILE *fpdata;

    fprintf(fpout,"Reading close data file: %s\n",filename);
	fpdata = fopen(filename,"rb");
	if (fpdata == NULL) {
		printf("fopen failed\n");
		return 1;
	} else {
		printf("fopen OK\n");
	}
	fread(&nxc,4,1,fpdata);
	printf("nxc: %d\n",nxc);
	fread(&nyc,4,1,fpdata);
	printf("nyc: %d\n",nyc);
	fread(&nzc,4,1,fpdata);
	printf("nzc: %d\n",nzc);
	fread(&nx8,4,1,fpdata);
	fread(&ny8,4,1,fpdata);
	fread(&nz8,4,1,fpdata);
	fread(&nbytes,4,1,fpdata);
	printf("nxc,nyc,nzc  nx8,ny8,nz8: %d %d %d  %d %d %d\n",nxc,nyc,nzc,nx8,ny8,nz8);
	if (nbytes != nx8*ny8*nz8/8) {
		printf("Error: this is not a compressed close data file: %s\n",filename);
		return 3;
	}
	closedata = (unsigned char *)malloc(nbytes*sizeof(unsigned char));
	fread(closedata,1,nbytes,fpdata);
	printf("read closedata: nbytes: %d\n",nbytes);
	return 0;
}

//-----------------------------------------------------------------------------------------
// Test if the point p[] corresponds to a lit voxel in closedata(:,:,:)
// Now using compressed close data, with each bit corresponding to a voxel
// Need to check indexing.  Note that p[] uses 1-based indexing.
// Consider a test case with nxc = nyc = nzc = 64
// In this case nx8 = ny8 = 64
// p[] = (1,1,1) should give kbyte=0, kbit=0
//   nb = 0
//   kbyte = 0
//   kbit = 0
//-----------------------------------------------------------------------------------------
bool in_close(int p[3]) 
{
	int kbyte, kbit, nb;
	unsigned char mbyte;
	int in;
	typedef std::bitset<sizeof(unsigned char)> ByteBits;

	// Need to convert to the byte and bit
    nb = p[0]-1 + (p[1]-1)*nx8 + (p[2]-1)*nx8*ny8;    // nb = 0,1,2,3,4,5,6,7,   8,9,10,11,12,13,14,15,
                                                      // kbit 0,1,2,3,4,5,6,7    0,1, 2, 3, 4, 5, 6, 7
                                                      //      kbyte = 0          kbyte = 1
    kbyte = nb/8;
    kbit = nb  - 8*kbyte;
    if (kbyte > nbytes-1) {
        printf("Error: kbyte > nbytes-1: %d %d %d  %d %d %d\n",p[0],p[1],p[2],nb,kbyte,nbytes-1);
		exit(1);
	}
	mbyte = closedata[kbyte];
//	printf("p, kbyte, mbyte, kbit: %d %d %d  %d  %d  %d\n",p[0], p[1], p[2], kbyte, mbyte, kbit);
//	in = ByteBits(mbyte).test(kbit);
	in = ((1 << kbit) & mbyte);
	return (in != 0);
}
    

//--------------------------------------------------------------------
// This uses 1-based indexing!
//--------------------------------------------------------------------
int MainWindow::getArea(int axis, int islice, float *area)
{
	int ix, iy, iz, p[3];
	double darea, total;

	total = 0;
	if (axis == 0) {
		darea = voxelsize[1]*voxelsize[2];
		for (iy=1; iy<=nyc; iy++) {	// this is taken care of in in_close() anyway
			for (iz=1; iz<=nzc; iz++) {
                p[0] = islice;
				p[1] = iy;
				p[2] = iz;
				if (in_close(p)) {
					total += darea;
				}
			}
		}
	}
	if (axis == 1) {
		darea = voxelsize[0]*voxelsize[2];
		for (ix=1; ix<=nxc; ix++) {	// this is taken care of in in_close() anyway
			for (iz=1; iz<=nzc; iz++) {
				p[0] = ix;
                p[1] = islice;
				p[2] = iz;
				if (in_close(p)) {
					total += darea;
				}
			}
		}
	}
	if (axis == 2) {
		darea = voxelsize[0]*voxelsize[1];
		for (ix=1; ix<=nxc; ix++) {	// this is taken care of in in_close() anyway
			for (iy=1; iy<=nyc; iy++) {
				p[0] = ix;
				p[1] = iy;
                p[2] = islice;
				if (in_close(p)) {
					total += darea;
				}
			}
		}
	}
	*area = (float)total;
	return 0;
}

//--------------------------------------------------------------------
// This uses 1-based indexing!
//--------------------------------------------------------------------
int MainWindow::getVolume(float *volume, int *ntvoxels)
{
    int ix, iy, iz, p[3], nt;
	double total, dvol;

	dvol = voxelsize[0]*voxelsize[1]*voxelsize[2];
	total = 0;
    nt = 0;
	for (ix=1; ix<=nxc; ix++) {
		for (iy=1; iy<=nyc; iy++) {
			for (iz=1; iz<=nzc; iz++) {
				p[0] = ix;
				p[1] = iy;
				p[2] = iz;
				if (in_close(p)) {
                    nt++;
					total += dvol;
				}
			}
		}
	}
    *ntvoxels = nt;
	*volume = (float)total;
	return 0;
}

//--------------------------------------------------------------------
// For image creation, the axes are permuted such that the slize plane is always a z plane
// If the case is really a x-slice, x -> z, y -> x, z -> y
// If the case is really a y-slice, y -> z, x -> y, z -> x
//--------------------------------------------------------------------
int MainWindow::histology(int axis, int islice, int *np, float *area)
{
	int iseg, cnt;
	float d;
	float zmin, zmax;
	POINT pos1, pos2;
    float z0, diam, S1[3], S2[3], vsize[2];

	zmin = 1.0e10;
	zmax = -zmin;
	d = islice*voxelsize[axis];
	cnt = 0;
	for (iseg=0; iseg<nsegments;iseg++) {
		pos1 = segment[iseg].end1;
		pos2 = segment[iseg].end2;
        diam = segment[iseg].diam;
		zmin = MIN(pos1.z,zmin);
		zmin = MIN(pos2.z,zmin);
		zmax = MAX(pos1.z,zmax);
		zmax = MAX(pos2.z,zmax);
		if (axis == 0 && ((pos1.x <= d  && d <= pos2.x) || (pos2.x <= d && d <= pos1.x))) {
			cnt++;
            S1[0] = pos1.y;
            S1[1] = pos1.z;
            S1[2] = pos1.x;
            S2[0] = pos2.y;
            S2[1] = pos2.z;
            S2[2] = pos2.x;
            z0 = d;
            vsize[0] = voxelsize[1];
            vsize[1] = voxelsize[2];
        }
		if (axis == 1 && ((pos1.y <= d  && d <= pos2.y) || (pos2.y <= d && d <= pos1.y))) {
			cnt++;
            S1[0] = pos1.z;
            S1[1] = pos1.x;
            S1[2] = pos1.y;
            S2[0] = pos2.z;
            S2[1] = pos2.x;
            S2[2] = pos2.y;
            z0 = d;
            vsize[0] = voxelsize[2];
            vsize[1] = voxelsize[0];
        }
		if (axis == 2 && ((pos1.z <= d  && d <= pos2.z) || (pos2.z <= d && d <= pos1.z))) {
			cnt++;
            S1[0] = pos1.x;
            S1[1] = pos1.y;
            S1[2] = pos1.z;
            S2[0] = pos2.x;
            S2[1] = pos2.y;
            S2[2] = pos2.z;
            z0 = d;
            vsize[0] = voxelsize[0];
            vsize[1] = voxelsize[1];
        }
        fillEllipse(z0,S1,S2,diam,vsize);
	}
	printf("nsegments, cnt: %d %d\n",nsegments,cnt);
    printf("axis, d, zmin, zmax: %d %f %f %f\n",axis,d,zmin,zmax);
	*np = cnt;
	getArea(axis,islice,area);
	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int MainWindow::average_histology(int *np, float *area)
{
    int islice, iseg, cnt, axis;
    float d, darea, totarea;
    POINT pos1, pos2;

    cnt = 0;
    totarea = 0;
    for (axis=0; axis<3; axis++) {
        fprintf(fpout,"axis: %d\n",axis);
        fflush(fpout);
        for (islice=range[axis][0]; islice <=range[axis][1]; islice++) {
            d = islice*voxelsize[axis];
//            fprintf(fpout,"islice: %d  d: %f\n",islice,d);
//            fflush(fpout);
            for (iseg=0; iseg<nsegments;iseg++) {
                pos1 = segment[iseg].end1;
                pos2 = segment[iseg].end2;
                if (axis == 0 && ((pos1.x <= d  && d <= pos2.x) || (pos2.x <= d && d <= pos1.x))) {
                    cnt++;
                }
                if (axis == 1 && ((pos1.y <= d  && d <= pos2.y) || (pos2.y <= d && d <= pos1.y))) {
                    cnt++;
                }
                if (axis == 2 && ((pos1.z <= d  && d <= pos2.z) || (pos2.z <= d && d <= pos1.z))) {
                    cnt++;
                }
            }
            getArea(axis,islice,&darea);
            totarea += darea;
//            fprintf(fpout,"%d %d  %6.1f\n",axis,islice,darea);
//            fflush(fpout);
        }
    }
    fprintf(fpout,"cnt: %d  totarea: %f\n",cnt,totarea);
    fflush(fpout);
    *np = cnt;
    *area = totarea;
    return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int MainWindow::getCloseSize(int nvoxels[])
{
    nvoxels[0] = nxc;
    nvoxels[1] = nyc;
    nvoxels[2] = nzc;
    return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
void free_all()
{
	int ie;
	NETWORK *net = NP0;

	if (net != NULL) {
		for (ie=0; ie < net->ne; ie++) {
			free(net->edgeList[ie].pt);
			free(net->edgeList[ie].pt_used);
		}
		free(net->vertex);
		free(net->edgeList);
		free(net->point);
	}
    if (closedata != NULL) free(closedata);
	if (segment != NULL) free(segment);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int MainWindow::setup(char *input_amfile, char *close_file, char *result_file, float vsize[])
{
	int err;

    fprintf(fpout,"setup\n");
	free_all();
	printf("did free_all\n");
	if (fperr == NULL) fperr = fopen("quantify_error.out","w");
	if (fpout == NULL) fpout = fopen(result_file,"w");
	printf("opened output files\n");
	for (int i=0; i<3; i++) {
		voxelsize[i] = vsize[i];
	}
	printf("set voxelsize\n");
	NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	printf("allocated network\n");
	err = ReadAmiraFile(input_amfile,NP0);
	if (err != 0) return 2;
	printf("read Amira file\n");
	err = ReadCloseFile(close_file);
	if (err != 0) return 3;
	is_setup = true;
	printf("did setup\n");
	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
bool MainWindow::isSetup()
{
	return is_setup;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
void MainWindow::reset()
{
	is_setup = false;
}

// Consider a cylinder intersecting a plane.
// The cylinder has radius R, the centreline C intersects the plane at P0 (x0,y0)
// the line C makes an angle alpha with the plane, the projection of C onto
// the plane is the line L
// First define the rectangle that bounds the ellipse E of intersection.
// Then look at all grid points (pixels) within the rectangle, and select
// those that are within E
// Consider a cylinder cross-section S that touches the plane.
// The point where the perimeter of the disk touches the plane is P1.
// The furthest point of intersection of the cylinder and the plane from P1 is P2.
// The contact points midway between P1 and P2 are P3 and P4.
// P1, P2, P3 and P4 define a rectangle (they are midpoints of the sides).
// The distance P0-P1 = a = R/sin(alpha) = P0-P2
// The distance P0-P3 = b = R = P0-P4
// The focal points F1 and F2 are a distance f from P0, f^2 = a^2 - b^2
// The perimeter of E is points P such that P-F1 + P-F2 = 2a
// Therefore a point P is within E if P-F1 + P-F2 <= 2a
//
// Case of z-slice
// The segment that intersects the plane z=z0 has ends S1(x1,y1,z1) and S2(x2,y2,z2)
// The line C is (x,y,z)' = (x1,y1,z1)' + s*(x2-x1,y2-y1,z2-z1)'
// z = z0 when z0 = z1 + s0*(z2-z1), i.e. when s = s0 = (z0-z1)/(z2-z1)
// and x0 = x1 + s0*(x2-x1), y0 = y1 + s0*(y2-y1)
// The line L is the projection of line C onto the plane z=z0
// (x,y)' = (x1,y1)' + s*(x2-x1,y2-y1)'
// (x,y)' = (x0,y0)' + (s-s0)*(x2-x1,y2-y1)'
// Let t = s-s0, then (x-x0,y-y0)' = t*(x2-x1,y2-y1)'
// and t=0 at P0(x0,y0,z0)
// The unit vector along the line L is u = N(x2-x1,y2-y1)
// The unit normal to the plane is the vector (0,0,1)
// Therefore the angle between the line C and the normal is given by:
// cos(gamma) = (0,0,1).N(x2-x1,y2-y1,z2-z1)'  where Nv is the normalised v
// alpha = PI/2 - gamma.  If gamma > PI/2, gamma = PI - gamma, alpha = gamma - PI/2
// a = R/sin(alpha), b = R, f = sqrt(a*a - b*b)
// P1 is (x0,y0) + a*u
// P2 is (x0,y0) - a*u
// The unit vector normal to u (in the plane) is:
// v = |0 -1|u
//     |1  0|
// The ellipse focal points are at F1 = (x0,y0)' + f*u, F2 = (x0,y0)' - f*u
// The four corners of the rectangle bounding the intersection ellipse E are:
// (x0,y0)' + a*u + b*v
// (x0,y0)' + a*u - b*v
// (x0,y0)' - a*u + b*v
// (x0,y0)' - a*u - b*v
// Therefore the range of x we need to evaluate points over is:
// min(x values) - max(x values)
// and similarly for y:
// min(y values) - max(y values)
// For a candidate grid point P(x,y), evaluate the distance from F1,
// d1 = P-F1, and from F2, d2 = P-F2
// then P is inside the ellipse if (d1+d2) < 2a
// This pixel is lit.
//----------------------------------------------------------------------------
// The axes have been permuted such that the slize plane is always a z plane
// If the case is really a x-slice, x -> z, y -> x, z -> y
// If the case is really a y-slice, y -> z, x -> y, z -> x
//----------------------------------------------------------------------------
void MainWindow::fillEllipse(float z0, float S1[], float S2[], float diam, float vsize[])
{
    float s0, P0[2], C[3], u[2], v[2], d, cosgamma, gamma, alpha;
    float a, b, f, F1[2], F2[2], xc[4], yc[4], xmin, xmax, ymin, ymax;
    float sum, x, y, d1, d2;
    int i, ix, iy, ixmin, ixmax, iymin, iymax;

    s0 = (z0 - S1[2])/(S2[2] - S1[2]);
    P0[0] = S1[0] + s0*(S2[0] - S1[0]);
    P0[1] = S1[1] + s0*(S2[1] - S1[1]);
    u[0] = S2[0] - S1[0];
    u[1] = S2[1] - S1[1];
    d = sqrt(u[0]*u[0] + u[1]*u[1]);
    u[0] /= d;
    u[1] /= d;
    v[0] = -u[1];
    v[1] = u[0];
    sum = 0;
    for (i=0; i<3; i++) {
        C[i] = S2[i] - S1[i];
        sum += C[i]*C[i];
    }
    d = sqrt(sum);
    cosgamma = C[2]/sqrt(d);
    gamma = acos(cosgamma);
    if (gamma < PI/2)
        alpha = PI/2 - gamma;
    else
        alpha = gamma - PI/2;
    a = (diam/2)/sin(alpha);
    b = diam/2;
    f = sqrt(a*a - b*b);
    F1[0] = P0[0] + f*u[0];
    F1[1] = P0[1] + f*u[1];
    F2[0] = P0[0] - f*u[0];
    F2[1] = P0[1] - f*u[1];
    xmin = 1.0e10;
    xmax = -xmin;
    ymin = 1.0e10;
    ymax = -ymin;
    xmin = MIN(xmin,P0[0] + a*u[0] + b*v[0]);
    xmin = MIN(xmin,P0[0] + a*u[0] - b*v[0]);
    xmin = MIN(xmin,P0[0] - a*u[0] + b*v[0]);
    xmin = MIN(xmin,P0[0] - a*u[0] - b*v[0]);
    xmax = MAX(xmax,P0[0] + a*u[0] + b*v[0]);
    xmax = MAX(xmax,P0[0] + a*u[0] - b*v[0]);
    xmax = MAX(xmax,P0[0] - a*u[0] + b*v[0]);
    xmax = MAX(xmax,P0[0] - a*u[0] - b*v[0]);
    ymin = MIN(ymin,P0[1] + a*u[1] + b*v[1]);
    ymin = MIN(ymin,P0[1] + a*u[1] - b*v[1]);
    ymin = MIN(ymin,P0[1] - a*u[1] + b*v[1]);
    ymin = MIN(ymin,P0[1] - a*u[1] - b*v[1]);
    ymax = MAX(ymax,P0[1] + a*u[1] + b*v[1]);
    ymax = MAX(ymax,P0[1] + a*u[1] - b*v[1]);
    ymax = MAX(ymax,P0[1] - a*u[1] + b*v[1]);
    ymax = MAX(ymax,P0[1] - a*u[1] - b*v[1]);
    // These values (xmin, xmax, ymin, ymax) bound the region within which the ellipse E lies
    ixmin = xmin/vsize[0];
    ixmax = xmax/vsize[0];
    iymin = ymin/vsize[1];
    iymax = ymax/vsize[1];
    for (ix=ixmin; ix<ixmax; ix++) {
        for (iy=iymin; iy<iymax; iy++) {
            x = ix*vsize[0];
            y = iy*vsize[1];
            d1 = sqrt((x-F1[0])*(x-F1[0]) + (y-F1[1])*(y-F1[1]));
            d2 = sqrt((x-F2[0])*(x-F2[0]) + (y-F2[1])*(y-F2[1]));
            if (d1+d2 <= 2*a) {
                // inside ellipse, lit voxel
            }
        }
    }
}

/*

//--------------------------------------------------------------------
// For testing
//--------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	int islice, axis, np;
	float area, volume;
	float vsize[3];
	char *input_amfile, *close_file, *result_file;

	if (argc != 7) {
		printf("Usage: quantify input_amfile close_file result_file vx vy vz\n");	// axislabel distance range\n");
		fperr = fopen("quantify_error.log","w");
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	close_file = argv[2];
	result_file = argv[3];
	sscanf(argv[4],"%f",&vsize[0]);
	sscanf(argv[5],"%f",&vsize[1]);
	sscanf(argv[6],"%f",&vsize[2]);

//	sscanf(argv[7],"%f",&z);
	//_splitpath(input_amfile,drive,dir,filename,ext);
	//strcpy(output_basename,drive);
	//strcat(output_basename,dir);
	//strcat(output_basename,basename);
	//sprintf(errfilename,"%s_am_cut.log",output_basename);
	//sprintf(result_file,"%s_am_cut.out",output_basename);
//	fperr = fopen("quantify_error.out","w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	//fpout = fopen(result_file,"w");	
	//NP0 = (NETWORK *)malloc(sizeof(NETWORK));
	//err = ReadAmiraFile(input_amfile,NP0);
	//if (err != 0) return 2;
	//err = ReadCloseFile(close_file);
	//if (err != 0) return 3;

	err = setup(input_amfile,close_file,result_file,vsize);
	if (err != 0) return err;

	getVolume(&volume);
	printf("Tissue volume (um3): %f\n",volume);

	axis = 2;
	islice = nzc/2;

	getArea(axis, islice,&area);
	area *= (float)1.0e-6;	// convert um2 -> mm2
	printf("Area: %8.4f\n",area);

//	iz = z/voxelsize[2] + 1;	// sticking with Fortran indexing for now, for consistency with vdistance.f90
	histology(axis,islice,&np,&area);
	area *= (float)1.0e-6;	// convert um2 -> mm2
	printf("Histology: axis, islice, np, area, np/area: %d %d %d %8.4f %8.1f\n",axis, islice, np, area, np/area);
	return 0;
}

*/
