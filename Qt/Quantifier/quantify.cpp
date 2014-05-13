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
double voxelsize[3];
unsigned char *closedata = NULL;
int nsegments;
SEGMENT_TYPE *segment = NULL;
NETWORK *NP0 = NULL;

bool is_setup = false;
*/

double pointDist(POINT p1, POINT p2)
{
    double dx, dy, dz;

    dx = p1.x - p2.x;
    dy = p1.y - p2.y;
    dz = p1.z - p2.z;
    return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
//-----------------------------------------------------------------------------------------------------
int MainWindow::ReadAmiraFile(char *amFile, NETWORK *net)
{
	int i, j, k, kp, npts;
	int np_used, ne_used;
	EDGE edge;
	char line[STR_LEN];

    printf("Amira file: %s\n",amFile);
    fprintf(fpout,"Amira file: %s\n",amFile);
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
                    sscanf(line,"%lf %lf %lf\n",&(net->vertex[i].point.x),&(net->vertex[i].point.y),&(net->vertex[i].point.z));
//					kp = i;
					net->vertex[i].point.d = 0;
//					net->point[kp] = net->vertex[i].point;
				}
//				kp++;
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
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
//						if (k > 0 && k<edge.npts-1) {											// See note above
                            sscanf(line,"%lf %lf %lf",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[k] = kp;
							net->edgeList[i].pt_used[k] = kp;
							kp++;
//						}
//						if (k > 0) {
//							len = len + dist(net,net->edgeList[i].pt[k-1],net->edgeList[i].pt[k]);
//						}
					}
//					net->edgeList[i].length_vox = len;
				}
			} else if (k == 5) {
				for (i=0;i<net->ne;i++) {
					edge = net->edgeList[i];
                    double dave = 0;
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @5\n");
							return 1;
						}
						j = edge.pt[k];
                        sscanf(line,"%lf",&net->point[j].d);
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
	
    printf("Total vessel segments: %d\n",nsegments);
    fprintf(fpout,"Total vessel segments: %d\n",nsegments);

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
            segment[iseg].len = pointDist(net->point[ip],net->point[ip+1]);
            if (DEBUG) {
                POINT p1 = segment[iseg].end1;
                POINT p2 = segment[iseg].end2;
//                fprintf(fpout,"seg: %6d  ip: %6d %6.1f %6.1f %6.1f --> %6.1f %6.1f %6.1f\n",iseg,ip,p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
            }
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

    fprintf(fpout,"Close data file: %s\n",filename);
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
int MainWindow::getArea(int axis, int islice, int *npixels, double *area)
{
    int ix, iy, iz, p[3], count;
    int ixmin, ixmax, iymin, iymax, izmin, izmax;
    double darea;

    ixmin = 1;
    ixmax = nxc;
    iymin = 1;
    iymax = nyc;
    izmin = 1;
    izmax = nzc;
    if (is_block) {
        ixmin = MAX(ixmin,range[0][0]);
        ixmax = MIN(ixmax,range[0][1]);
        iymin = MAX(iymin,range[1][0]);
        iymax = MIN(iymax,range[1][1]);
        izmin = MAX(izmin,range[2][0]);
        izmax = MIN(izmax,range[2][1]);
    }
    count = 0;
	if (axis == 0) {
		darea = voxelsize[1]*voxelsize[2];
        for (iy=iymin; iy<=iymax; iy++) {	// this is taken care of in in_close() anyway
            for (iz=izmin; iz<=izmax; iz++) {
                p[0] = islice+1;
				p[1] = iy;
				p[2] = iz;
				if (in_close(p)) {
                    count++;
				}
			}
		}
	}
	if (axis == 1) {
		darea = voxelsize[0]*voxelsize[2];
        for (ix=ixmin; ix<=ixmax; ix++) {	// this is taken care of in in_close() anyway
            for (iz=izmin; iz<=izmax; iz++) {
				p[0] = ix;
                p[1] = islice+1;
				p[2] = iz;
				if (in_close(p)) {
                    count++;
                }
			}
		}
	}
	if (axis == 2) {
		darea = voxelsize[0]*voxelsize[1];
        for (ix=ixmin; ix<=ixmax; ix++) {	// this is taken care of in in_close() anyway
            for (iy=iymin; iy<=iymax; iy++) {
				p[0] = ix;
				p[1] = iy;
                p[2] = islice+1;
				if (in_close(p)) {
                    count++;
                }
			}
		}
	}
    *npixels = count;
    *area = count*darea;
	return 0;
}

//--------------------------------------------------------------------
// This uses 1-based indexing!
//--------------------------------------------------------------------
int MainWindow::getVolume(double *volume, int *ntvoxels)
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
    *volume = total;
	return 0;
}

//--------------------------------------------------------------------
// For image creation, the axes are permuted such that the slice plane is always a z plane
// If the case is really a x-slice, x -> z, y -> x, z -> y
// If the case is really a y-slice, y -> z, x -> y, z -> x
//--------------------------------------------------------------------
int MainWindow::SliceHistology(int axis, int islice, int *nvessels, int *nvesselpixels, int *ntissuepixels, double *tissuearea)
{
    int iseg, cnt, nv[2], npixels, ntpixels, kx, ky, rng_x[2], rng_y[2];
    double d;
    double zmin, zmax;
	POINT pos1, pos2;
    double z0, diam, S1[3], S2[3], vsize[2];
    bool hit;

	zmin = 1.0e10;
	zmax = -zmin;
	d = islice*voxelsize[axis];
    ntpixels = 0;
	cnt = 0;
	for (iseg=0; iseg<nsegments;iseg++) {
		pos1 = segment[iseg].end1;
		pos2 = segment[iseg].end2;
        diam = segment[iseg].diam;
		zmin = MIN(pos1.z,zmin);
		zmin = MIN(pos2.z,zmin);
		zmax = MAX(pos1.z,zmax);
		zmax = MAX(pos2.z,zmax);
        hit = false;
        if (axis == 0 && ((pos1.x <= d  && d < pos2.x) || (pos2.x < d && d <= pos1.x))) {   // (x,y,z) -> (z,x,y)
            S1[0] = pos1.y;
            S1[1] = pos1.z;
            S1[2] = pos1.x;
            S2[0] = pos2.y;
            S2[1] = pos2.z;
            S2[2] = pos2.x;
            z0 = d;
            kx = 1;
            ky = 2;
            hit = true;
        }
        if (axis == 1 && ((pos1.y <= d  && d < pos2.y) || (pos2.y < d && d <= pos1.y))) {   // (x,y,z) -> (y,z,x)
            S1[0] = pos1.z;
            S1[1] = pos1.x;
            S1[2] = pos1.y;
            S2[0] = pos2.z;
            S2[1] = pos2.x;
            S2[2] = pos2.y;
            z0 = d;
            kx = 2;
            ky = 0;
//            if (S1[0] < 0.5*nv[0]*vsize[0]) //Note: odd image for y axis, 36
            hit = true;
        }
        if (axis == 2 && ((pos1.z <= d  && d < pos2.z) || (pos2.z < d && d <= pos1.z))) {   // (x,y,z) -> (x,y,z)
            S1[0] = pos1.x;
            S1[1] = pos1.y;
            S1[2] = pos1.z;
            S2[0] = pos2.x;
            S2[1] = pos2.y;
            S2[2] = pos2.z;
            z0 = d;
            kx = 0;
            ky = 1;
            hit = true;
        }
        if (hit) {
            vsize[0] = voxelsize[kx];
            vsize[1] = voxelsize[ky];
            nv[0] = nvoxels[kx];
            nv[1] = nvoxels[ky];
            if (is_block) {
                rng_x[0] = range[kx][0];
                rng_x[1] = range[kx][1];
                rng_y[0] = range[ky][0];
                rng_y[1] = range[ky][1];
            }
            if (DEBUG) {
                fprintf(fpout,"rng_x: %d %d rng_y: %d %d\n",rng_x[0],rng_x[1],rng_y[0],rng_y[1]);
                fprintf(fpout,"\nsegment: %6d S1: %6.1f %6.1f %6.1f S2: %6.1f %6.1f %6.1f diam: %4.1f\n",
                        iseg,S1[0],S1[1],S1[2],S2[0],S2[1],S2[2],diam);
                fflush(fpout);
            }
            fillEllipse(z0,S1,S2,diam,vsize,nv,rng_x,rng_y,&npixels);
            cnt++;
            ntpixels += npixels;
        }
	}
    if (DEBUG) {
        fprintf(fpout,"Vessel count: %6d\nTotal pixels: %8d\n",cnt,ntpixels);
        fflush(fpout);
    }
//    printf("Vessel count: %6d\nTotal vessel pixels: %8d\n",cnt,ntpixels);
//    printf("axis, d, zmin, zmax: %d %f %f %f\n",axis,d,zmin,zmax);
    *nvessels = cnt;
    *nvesselpixels = ntpixels;
    getArea(axis,islice,ntissuepixels,tissuearea);
	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int MainWindow::VolumeHistology(bool *use_axis, int *nvessels, int *nvesselpixels, int *ntissuepixels, double *tissuearea)
{
    int islice, axis, nslices, err;
    int nvess, nvesselpix, ntissuepix;
    int nvess_tot, nvesselpix_tot, ntissuepix_tot;
    double tissarea, tissarea_tot;

//    fprintf(fpout,"VolumeHistology\n");
    fflush(fpout);
    for (axis=0; axis<3; axis++) {
        if (!use_axis[axis]) continue;
        fprintf(fpout,"axis: %d\n",axis);
        fflush(fpout);
        nvess_tot = 0;
        nvesselpix_tot = 0;
        ntissuepix_tot = 0;
        tissarea_tot = 0;
        nslices = range[axis][1] - range[axis][0] + 1;
        for (islice=range[axis][0]; islice <=range[axis][1]; islice++) {
//            fprintf(fpout,"axis,slice: %d %d\n",axis,islice);
            fflush(fpout);
            err = SliceHistology(axis, islice, &nvess, &nvesselpix, &ntissuepix, &tissarea);
            if (err != 0) {
                fprintf(fpout,"error: %d",err);
                fflush(fpout);
            }
            nvess_tot += nvess;
            nvesselpix_tot += nvesselpix;
            ntissuepix_tot += ntissuepix;
            tissarea_tot += tissarea;
        }
        nvessels[axis] = (double)nvess_tot/nslices + 0.5;
        nvesselpixels[axis] = (double)nvesselpix_tot/nslices + 0.5;
        ntissuepixels[axis] = (double)ntissuepix_tot/nslices + 0.5;
        tissuearea[axis] = tissarea_tot/nslices;
//        fprintf(fpout,"nslices: %d  %d %d %d %f\n",nslices,nvess_tot,nvesselpix_tot,ntissuepix_tot,tissarea);
        fflush(fpout);
    }
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
int MainWindow::setup(char *input_amfile, char *close_file, char *result_file, double vsize[])
{
	int err;

    fprintf(fpout,"Quantifier\n----------\n");
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

//--------------------------------------------------------------------
// Network topology is in NETWORK *NP0
// Vertices (branch points) are NP0->vertex[]
//--------------------------------------------------------------------
int MainWindow::branching(int *nbranchpts, double *totlen, double *totvol)
{
    int k, nb;
    double xmin, xmax, ymin, ymax, zmin, zmax, tlen;
    POINT p;

    xmin = range[0][0]*voxelsize[0];
    xmax = range[0][1]*voxelsize[0];
    ymin = range[1][0]*voxelsize[1];
    ymax = range[1][1]*voxelsize[1];
    zmin = range[2][0]*voxelsize[2];
    zmax = range[2][1]*voxelsize[2];
    nb = 0;
    for (k=0; k<NP0->nv; k++) {
        p = NP0->vertex[k].point;
        if (p.x >= xmin && p.x < xmax && p.y >= ymin && p.y < ymax && p.z >= zmin && p.z < zmax) nb++;
    }
    *nbranchpts = nb;
    tlen = 0;
    for (k=0; k<nsegments; k++) {
        p = segment[k].end1;
        if (p.x >= xmin && p.x < xmax && p.y >= ymin && p.y < ymax && p.z >= zmin && p.z < zmax) {
            p = segment[k].end2;
            if (p.x >= xmin && p.x < xmax && p.y >= ymin && p.y < ymax && p.z >= zmin && p.z < zmax) {
                tlen += segment[k].len;
            }
        }
    }
    *totlen = tlen;
    *totvol = rangeVolume();
    return 0;
}

//--------------------------------------------------------------------
// Volume of tissue region defined by range of slices
//--------------------------------------------------------------------
double MainWindow::rangeVolume()
{
    int ix, iy, iz, p[3], nt;

    nt = 0;
    for (ix=range[0][0]; ix<=range[0][1]; ix++){
        for (iy=range[1][0]; iy<=range[1][1]; iy++){
            for (iz=range[2][0]; iz<=range[2][1]; iz++){
                p[0] = ix;
                p[1] = iy;
                p[2] = iz;
                if (in_close(p)) nt++;
            }
        }
    }
    return nt*voxelsize[0]*voxelsize[1]*voxelsize[2];
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
void MainWindow::fillEllipse(double z0, double S1[], double S2[], double diam, double vsize[], int nv[], int rng_x[], int rng_y[], int *npixels)
{
    double s0, P0[2], C[3], u[2], v[2], d, gamma, alpha;
    double a, b, f, F1[2], F2[2], xmin, xmax, ymin, ymax;
    double sum, x, y, d1, d2;
    int i, ix, iy, ixmin, ixmax, iymin, iymax, npix;
    bool circle;

    npix = 0;
    if (S1[2] == S2[2]) S2[2] += 1;
    s0 = (z0 - S1[2])/(S2[2] - S1[2]);
    P0[0] = S1[0] + s0*(S2[0] - S1[0]);
    P0[1] = S1[1] + s0*(S2[1] - S1[1]);
    u[0] = S2[0] - S1[0];
    u[1] = S2[1] - S1[1];
    d = sqrt(u[0]*u[0] + u[1]*u[1]);
    if (d == 0) {
        circle = true;
        xmin = P0[0] - diam/2 - 1;
        xmax = P0[0] + diam/2 + 1;
        ymin = P0[1] - diam/2 - 1;
        ymax = P0[1] + diam/2 + 1;
        if (DEBUG) {
            fprintf(fpout,"circle: centre: %6.1f %6.1f radius: %6.1f\n",P0[0],P0[1],diam/2);
        }
    } else {
        circle = false;
        u[0] /= d;
        u[1] /= d;
        v[0] = -u[1];
        v[1] = u[0];
        if (DEBUG) {
            fprintf(fpout,"s0: %8.3f u: %6.3f %6.3f v: %6.3f %6.3f\n",s0,u[0],u[1],v[0],v[1]);
        }
        sum = 0;
        for (i=0; i<3; i++) {
            C[i] = S2[i] - S1[i];
            sum += C[i]*C[i];
        }
        d = sqrt(sum);
        if (DEBUG) {
            fprintf(fpout,"C: %6.1f %6.1f %6.1f  d: %6.1f\n",C[0],C[1],C[2],d);
        }
        for (i=0; i<3; i++) {
            C[i] /= d;
        }
        gamma = acos(C[2]);
        if (gamma < PI/2)
            alpha = PI/2 - gamma;
        else
            alpha = gamma - PI/2;

        alpha = MAX(alpha,ALPHAMAX);

        a = (diam/2)/sin(alpha);
        b = diam/2;
        f = sqrt(a*a - b*b);
        if (DEBUG) {
            fprintf(fpout,"gamma: %6.1f alpha: %6.1f a,b,f: %6.1f %6.1f %6.1f\n",gamma*180/PI,alpha*180/PI,a,b,f);
        }
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
    }
    // These values (xmin, xmax, ymin, ymax) bound the region within which the ellipse E lies
    // The um values must be converted to voxels
    ixmin = xmin/vsize[0];
    ixmax = xmax/vsize[0];
    iymin = ymin/vsize[1];
    iymax = ymax/vsize[1];
    ixmin = MAX(ixmin,0);
    ixmax = MIN(ixmax,nv[0]);
    iymin = MAX(iymin,0);
    iymax = MIN(iymax,nv[1]);
    if (DEBUG) {
        fprintf(fpout,"ix range: %6d %6d iy range: %6d %6d\n",ixmin,ixmax,iymin,iymax);
    }
    if (is_block) {
        ixmin = MAX(ixmin,rng_x[0]);
        iymin = MAX(iymin,rng_y[0]);
        ixmax = MIN(ixmax,rng_x[1]);
        iymax = MIN(iymax,rng_y[1]);
        if (DEBUG) {
            fprintf(fpout,"is_block -> ix range: %6d %6d iy range: %6d %6d\n",ixmin,ixmax,iymin,iymax);
        }
    }
    for (ix=ixmin; ix<ixmax; ix++) {
        for (iy=iymin; iy<iymax; iy++) {
            x = ix*vsize[0];
            y = iy*vsize[1];
            if (circle) {
                d1 = sqrt((x-P0[0])*(x-P0[0]) + (y-P0[1])*(y-P0[1]));
                if (d1 <= diam/2) {
                    // inside circle, lit voxel
                    if (DEBUG) fprintf(fpout,"lit: %d %d\n",ix,iy);
                    if (is_slice) {
                        imageViewer->myQtImage->setPixel(ix,iy,255);
                    }
                    npix++;
                }
            } else {
                d1 = sqrt((x-F1[0])*(x-F1[0]) + (y-F1[1])*(y-F1[1]));
                d2 = sqrt((x-F2[0])*(x-F2[0]) + (y-F2[1])*(y-F2[1]));
                if (d1+d2 <= 2*a) {
                    // inside ellipse, lit voxel
                    if (is_slice) {
                        imageViewer->myQtImage->setPixel(ix,iy,255);
                    }
                    npix++;
                }
            }
        }
    }
    if (DEBUG) {
        fprintf(fpout,"npixels: %6d\n",npix);
    }
    *npixels = npix;
}

/*

//--------------------------------------------------------------------
// For testing
//--------------------------------------------------------------------
int main(int argc, char **argv)
{
	int err;
	int islice, axis, np;
    double area, volume;
    double vsize[3];
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
    area *= (double)1.0e-6;	// convert um2 -> mm2
	printf("Area: %8.4f\n",area);

//	iz = z/voxelsize[2] + 1;	// sticking with Fortran indexing for now, for consistency with vdistance.f90
	histology(axis,islice,&np,&area);
    area *= (double)1.0e-6;	// convert um2 -> mm2
	printf("Histology: axis, islice, np, area, np/area: %d %d %d %8.4f %8.1f\n",axis, islice, np, area, np/area);
	return 0;
}

*/
