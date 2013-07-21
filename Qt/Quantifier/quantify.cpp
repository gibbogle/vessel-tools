// To quantify vessel density per area, and compute tissue volume

//#include <cstdio>
//#include <vector>

//#include <algorithm>
//#include <math.h>
//#include <string.h>
//#include <string>
//#include <sstream>

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
int ReadAmiraFile(char *amFile, NETWORK *net)
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
			iseg++;
		}
	}
	nsegments = iseg;
	printf("nsegments: %d\n",nsegments);

	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int ReadCloseFile(char *filename)
{
	FILE *fpdata;

	printf("Reading close data file: %s\n",filename);
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
// Test if the point xyz(:) corresponds to a lit voxel in closedata(:,:,:)
// Now using compressed close data, with each bit corresponding to a voxel
// Need to check indexing (now 0-based)
// This need to be tested by comparison with vdistance.f90
//-----------------------------------------------------------------------------------------
bool in_close(int p[3]) 
{
	int kbyte, kbit, nb;
	unsigned char mbyte;
	int in;
	typedef std::bitset<sizeof(unsigned char)> ByteBits;

	// Need to convert to the byte and bit
	nb = p[0] + (p[1]-1)*nx8 + (p[2]-1)*nx8*ny8;
//	nb--;	// for 0-based indexing
	if (nb%8 == 0) {
		kbyte = nb/8;
		kbit = 0;
	} else {
		kbyte = nb/8 + 1;
		kbit = nb  - 8*(kbyte-1);
	}
	if (kbyte > nbytes) {
		printf("Error: kbyte > nbytes: %d %d %d  %d %d %d\n",p[0],p[1],p[2],nb,kbyte,nbytes);
		exit(1);
	}
	kbyte--;	// for 0-based indexing
	if (kbyte >= nbytes) {
		printf("kbyte >= nbytes: %d %d\n",kbyte, nbytes);
		exit(1);
	}
	mbyte = closedata[kbyte];
	//if (btest(mbyte,kbit))	// need equiv of btest()
	//	return true;
	//else
	//	return false;
//	printf("p, kbyte, mbyte, kbit: %d %d %d  %d  %d  %d\n",p[0], p[1], p[2], kbyte, mbyte, kbit);
//	in = ByteBits(mbyte).test(kbit);
	in = ((1 << kbit) & mbyte);
//	printf("in: %d\n",in);
	return (in != 0);
}
    

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int getArea(int axis, int islice, float *area)
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
//--------------------------------------------------------------------
int getVolume(float *volume)
{
	int ix, iy, iz, p[3];
	double total, dvol;

	printf("getVolume\n");
	dvol = voxelsize[0]*voxelsize[1]*voxelsize[2];
	total = 0;
	for (ix=1; ix<=nxc; ix++) {
		for (iy=1; iy<=nyc; iy++) {
			for (iz=1; iz<=nzc; iz++) {
				p[0] = ix;
				p[1] = iy;
				p[2] = iz;
				if (in_close(p)) {
					total += dvol;
				}
			}
		}
	}
	*volume = (float)total;
	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int histology(int axis, int islice, int *np, float *area)
{
	int iseg, cnt;
	float d;
	float zmin, zmax;
	POINT pos1, pos2;

	zmin = 1.0e10;
	zmax = -zmin;
	d = islice*voxelsize[axis];
	cnt = 0;
	for (iseg=0; iseg<nsegments;iseg++) {
		pos1 = segment[iseg].end1;
		pos2 = segment[iseg].end2;
		zmin = MIN(pos1.z,zmin);
		zmin = MIN(pos2.z,zmin);
		zmax = MAX(pos1.z,zmax);
		zmax = MAX(pos2.z,zmax);
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
	printf("nsegments, cnt: %d %d\n",nsegments,cnt);
	printf("axis, d, zmin, zmax: %d %f %f %f\n",axis,d,zmin,zmax);
	*np = cnt;
	getArea(axis,islice,area);
	return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
int getCloseSize(int nvoxels[])
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
int setup(char *input_amfile, char *close_file, char *result_file, float vsize[])
{
	int err;

	printf("setup\n");
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
bool isSetup()
{
	return is_setup;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
void reset()
{
	is_setup = false;
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
