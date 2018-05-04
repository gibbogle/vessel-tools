// To generate b&w image from network description (CMGUI or Amira SpatialGraph file)
// Create a b&w image from the network, where a voxel is white if it falls within
// any conical vessel segment or any nodal sphere.
// This can then be transformed into an STL file in Amira: image -> isosurface -> STL

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSize.h"

#include "network.h"

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
unsigned char *p;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define STR_LEN 128

NETWORK *NP;
FILE *fperr, *fpout;
int iesel;
float delta, dmin, dbox;
double netmin[3], netmax[3];
int N[3];
int Nbox[3];
int ****vlist;
int ***nvlist;
int maxvlist;

bool dbug;

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double netdist(NETWORK *net, int k1, int k2)
{
	double dx = net->point[k2].x - net->point[k1].x;
	double dy = net->point[k2].y - net->point[k1].y;
	double dz = net->point[k2].z - net->point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double dx = x2 - x1;
	double dy = y2 - y1;
	double dz = z2 - z1;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
// Length of a vector
//-----------------------------------------------------------------------------------------------------
double vnorm(double v[3])
{
	double sum=0;
	for (int i=0; i<3; i++) {
		sum += v[i]*v[i];
	}
	return sqrt(sum);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double cosangle(double v1[3], double v2[3])
{
	double d1 = 0;
	double d2 = 0;
	double cosum = 0;
	for (int i=0; i<3; i++) {
		cosum += v1[i]*v2[i];
		d1 += v1[i]*v1[i];
		d2 += v2[i]*v2[i];
	}
	return cosum/sqrt(d1*d2);
}

//-----------------------------------------------------------------------------------------------------
// Crossproduct of u and v, u x v, is the determinant:
//
//  | i  j  k  |
//  | ux uy uz |
//  | vx vy vz |
//
// = i(uy.vz - vy.uz) - j(ux.vz - vx.uz) + k(ux.vy - vx.uy)
//-----------------------------------------------------------------------------------------------------
void crossproduct(double u[3], double v[3], double cross[3])
{
	cross[0] =  u[1]*v[2] - u[2]*v[1];
	cross[1] = -u[0]*v[2] + u[2]*v[0];
	cross[2] =  u[0]*v[1] - u[1]*v[0];
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void SetEdgeBnds(double netmin[], double netmax[])
{
	int i, iseg;
	float r, dave;
	NDPOINT p1, p2;
	EDGE *edge;

	for (i=0; i<3; i++) {
		netmin[i] = 1.0e10;
		netmax[i] = -1.0e10;
	}
	for (int ie=0; ie<NP->ne; ie++) {
		edge = NP->edgeList + ie;
		edge->smin = new POS[edge->npts];
		edge->smax = new POS[edge->npts];
		for (i=0; i<3; i++) {
			edge->emin.pos[i] = 1.0e10;
			edge->emax.pos[i] = -1.0e10;
		}
		dave = 0;
		for (iseg=0; iseg<edge->npts-1; iseg++) {
			p1 = NP->point[edge->pt[iseg]];
			p2 = NP->point[edge->pt[iseg+1]];
			r = MAX(p1.d/2,p2.d/2);
			edge->smin[iseg].pos[0] = MIN(p1.x,p2.x) - r;
			edge->smin[iseg].pos[1] = MIN(p1.y,p2.y) - r;
			edge->smin[iseg].pos[2] = MIN(p1.z,p2.z) - r;
			edge->smax[iseg].pos[0] = MAX(p1.x,p2.x) + r;
			edge->smax[iseg].pos[1] = MAX(p1.y,p2.y) + r;
			edge->smax[iseg].pos[2] = MAX(p1.z,p2.z) + r;
			for (i=0; i<3; i++) {
				edge->emin.pos[i] = MIN(edge->emin.pos[i],edge->smin[iseg].pos[i]);
				edge->emax.pos[i] = MAX(edge->emax.pos[i],edge->smax[iseg].pos[i]);
			}
			dave += p1.d;
		}
		for (i=0; i<3; i++) {
			netmin[i] = MIN(edge->emin.pos[i],netmin[i]);
			netmax[i] = MAX(edge->emax.pos[i],netmax[i]);
		}
//		printf("edge: %d npts: %d diam: %6.2f\n",ie,edge->npts,dave/edge->npts);
	}
}

//-----------------------------------------------------------------------------------------------------
// vx, vy, vz define the tube axes (vx along the centreline), d is length, r is the radius
// check if pos is inside by finding tpos (coordinates in tube axes)
//-----------------------------------------------------------------------------------------------------
bool InsideTube(double vx[3], double vy[3], double vz[3], double d, double r, double pos[3])
{
	int i, j;
	double a[3][3], sum, tpos[3];

	for (i=0; i<3; i++) {
		a[0][i] = vx[i];
		a[1][i] = vy[i];
		a[2][i] = vz[i];
	}
	for (i=0; i<3; i++) {
		sum = 0;
		for (j=0; j<3; j++) {
			sum += a[i][j]*pos[j];
		}
		tpos[i] = sum;
	}
	if (tpos[0] < 0 || tpos[0] > d)
		return false;
	if (tpos[1]*tpos[1]+tpos[2]*tpos[2] > r*r)
		return false;
	if (dbug) printf("tpos: %6.1f %6.1f %6.1f   %6.1f\n",tpos[0],tpos[1],tpos[2],sqrt(tpos[1]*tpos[1]+tpos[2]*tpos[2]));
	return true;
}

//-----------------------------------------------------------------------------------------------------
// Testing
//-----------------------------------------------------------------------------------------------------
void getapex(void)
{
	double x0, y0, z0, x1, y1, z1, d1, x2, y2, z2, d2;
	double d01, d12, h, a0, cosa, cone_cosa, pos[3], v[3], vcone[3];
	double vx[3], vy[3], vz[3];

	x1 = 0;
	y1 = 0;
	z1 = 0;
	x2 = 1.0;
	y2 = 0.0;
	z2 = 0.0;
	d1 = 1.0;
	d2 = 2.0;
	pos[0] = -0.9;
	pos[1] = 1.3;
	pos[2] = 0.5;

	// Find the apex of the cone
	d12 = dist(x1,y1,z1,x2,y2,z2);
	d01 = d12*d1/(d2 - d1);
	a0 = -d01/d12;
	x0 = x1 + a0*(x2-x1);
	y0 = y1 + a0*(y2-y1);
	z0 = z1 + a0*(z2-z1);
	printf("d1,d2: %f %f  d01,d12: %f %f\n",d1,d2,d01,d12);
	printf("apex: %f %f %f\n",x0,y0,z0);
	h = sqrt(d2*d2+(d01+d12)*(d01+d12));
	cone_cosa = (d01+d12)/h;
	vcone[0] = x1 - x0;
	vcone[1] = y1 - y0;
	vcone[2] = z1 - z0;
	v[0] = pos[0] - x0;
	v[1] = pos[1] - y0;
	v[2] = pos[2] - z0;
	cosa = cosangle(v,vcone);
	printf("cone_cosa, cosa: %f %f\n",cone_cosa,cosa);
	d1 = 2;
	vcone[0] = x2 - x1;
	vcone[1] = y2 - y1;
	vcone[2] = z2 - z1;
	for (int i=0; i<3; i++)
		vx[i] = vcone[i]/d12;
	v[0] = 1; v[1] = 0; v[2] = 0;
	if (cosangle(vx,v) > 0.9) {
		v[0] = 0; v[1] = 1; v[2] = 0;
	}
	crossproduct(vx,v,vy);
	crossproduct(vx,vy,vz);
	if (InsideTube(vx,vy,vz,d12,d1/2,pos))
		printf("inside\n");
	else
		printf("outside\n");
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool InsideSphere(NDPOINT *p, double pos[3])
{
	double xv, yv, zv, r, d;

	xv = p->x;
	yv = p->y;
	zv = p->z;
	r = p->d/2;
	d = dist(xv,yv,zv,pos[0],pos[1],pos[2]);
	if (d <= r)
		return true;
	return false;
}

//-----------------------------------------------------------------------------------------------------
// Check if pos[] is inside a truncated cone with end points (x1,y1,z1), (x2,y2,z2) and diameters d1, d2.
//-----------------------------------------------------------------------------------------------------
int InsideSegment(int ie, int iseg, double pos[3])
{
	int ip1, ip2, i;
	double x0, y0, z0, x1, y1, z1, r1, x2, y2, z2, r2;
	double d12, d01, h, a0, cone_cosa, cosa, v[3], vcone[3], d, dcone;
	double vx[3], vy[3], vz[3], rpos[3];
	NDPOINT *p1, *p2;

	ip1 = NP->edgeList[ie].pt[iseg];
	ip2 = NP->edgeList[ie].pt[iseg+1];
	p1 = &NP->point[ip1];
	p2 = &NP->point[ip2];
	// Check spheres
	if (InsideSphere(p1,pos)) 
		return 1;
	if (InsideSphere(p2,pos)) 
		return 1;
	if (p1->d == p2->d) {
		// Special case, cone is a tube
		// Define local axes for the tube:
		// vx along P1-P2
		// vy normal to vx and some reference direction
		// vz normal to vx and vy, vz = vx cross vy
		x1 = p1->x;
		y1 = p1->y;
		z1 = p1->z;
		x2 = p2->x;
		y2 = p2->y;
		z2 = p2->z;
		d12 = dist(x1,y1,z1,x2,y2,z2);
		if (dbug) printf("P1,P2,d12: %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f  %6.2f\n",x1,y1,z1,x2,y2,z2,d12);
		vcone[0] = x2 - x1;
		vcone[1] = y2 - y1;
		vcone[2] = z2 - z1;
		for (i=0; i<3; i++)
			vx[i] = vcone[i]/d12;
		v[0] = 1; v[1] = 0; v[2] = 0;
		if (cosangle(vx,v) > 0.9) {
			v[0] = 0; v[1] = 1; v[2] = 0;
		}
		if (dbug) printf("vx,v: %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f\n",vx[0],vx[1],vx[2],v[0],v[1],v[2]);
		crossproduct(vx,v,vy);
		d = vnorm(vy);
		if (dbug) printf("vy: %6.1f %6.1f %6.1f  length of vy: %6.1f\n",vy[0],vy[1],vy[2],d);
		for (i=0; i<3; i++)
			vy[i] = vy[i]/d;
		crossproduct(vx,vy,vz);
		if (dbug) printf("vy,vz:  %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f\n",vy[0],vy[1],vy[2],vz[0],vz[1],vz[2]);
		// Determine point position relative to P1
		rpos[0] = pos[0] - x1;
		rpos[1] = pos[1] - y1;
		rpos[2] = pos[2] - z1;
		if (dbug) printf("rpos: %6.1f %6.1f %6.1f\n",rpos[0],rpos[1],rpos[2]);
		if (InsideTube(vx,vy,vz,d12,p1->d/2,rpos)) {
			if (dbug) printf("inside\n");
			return 2;
		}
	} else {
		if (p2->d < p1->d) {
			p1 = &NP->point[ip2];
			p2 = &NP->point[ip1];
		}
		x1 = p1->x;
		y1 = p1->y;
		z1 = p1->z;
		r1 = p1->d/2;
		x2 = p2->x;
		y2 = p2->y;
		z2 = p2->z;
		r2 = p2->d/2;
		// Find the apex of the cone
		d12 = dist(x1,y1,z1,x2,y2,z2);
		d01 = d12*r1/(r2 - r1);
		a0 = -d01/d12;
		x0 = x1 + a0*(x2-x1);
		y0 = y1 + a0*(y2-y1);
		z0 = z1 + a0*(z2-z1);
		// First requirement is for pos[] to fall within the cone
		// i.e. angle between (P - P0) and (P1 - P0) less than cone angle.
		// or cos(angle) >= cone_cosa
		h = sqrt(r2*r2+(d01+d12)*(d01+d12));
		cone_cosa = (d01+d12)/h;
		vcone[0] = x1 - x0;
		vcone[1] = y1 - y0;
		vcone[2] = z1 - z0;
		v[0] = pos[0] - x0;
		v[1] = pos[1] - y0;
		v[2] = pos[2] - z0;
		cosa = cosangle(v,vcone);
		if (cosa < cone_cosa)
			return 0;
		// Now need to check that the point falls within the segment section of the cone.
		// d = distance from the apex
		d = dist(x0,y0,z0,pos[0],pos[1],pos[2]);
		// dcone = distance along the centreline
		dcone = cosa*d;
//		printf("d12: %6.1f  d01: %6.1f a0: %6.3f apex: %6.1f %6.1f %6.1f\n",d12,d01,a0,x0,y0,z0);
//		printf("h: %6.1f  cone_cosa: %8.5f cosa: %8.5f dcone: %6.1f\n",h,cone_cosa,cosa,dcone);
		if (dcone >= d01 && dcone <= d01+d12) {
			return 2;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool InsideVertex(int ie, int iv, double pos[3])
{
	double xv, yv, zv, r, d;
	NDPOINT *p;

	p = &NP->point[NP->edgeList[ie].vert[iv]];
	xv = p->x;
	yv = p->y;
	zv = p->z;
	r = p->d/2;
	d = dist(xv,yv,zv,pos[0],pos[1],pos[2]);
	if (d <= r)
		return true;
	return false;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int InsideVessel(int ie, double pos[3])
{
	int iseg, i;
	bool inseg;
	EDGE *edge;

	edge = &(NP->edgeList[ie]);
	/*
	for (i=0; i<3; i++) {
		if (pos[i] < edge->emin.pos[i] || pos[i] > edge->emax.pos[i]) {
			printf("outside edge bnds: \n");
			printf("x: %6.1f  bnds: %6.1f  %6.1f\n",pos[0],edge->emin.pos[0],edge->emax.pos[0]);
			printf("y: %6.1f  bnds: %6.1f  %6.1f\n",pos[1],edge->emin.pos[1],edge->emax.pos[1]);
			printf("z: %6.1f  bnds: %6.1f  %6.1f\n",pos[2],edge->emin.pos[2],edge->emax.pos[2]);
			exit(1);
			return false;
		}
	}
	*/

	for (iseg=0; iseg<edge->npts-1; iseg++) {
		inseg = true;
		for (i=0; i<3; i++) {
			if (pos[i] < edge->smin[iseg].pos[i] || pos[i] > edge->smax[iseg].pos[i]) {
				inseg = false;
				break;
			}
		}
		if (!inseg) continue;
		// pos[] may fall within segment iseg
		int res = InsideSegment(ie, iseg, pos);
		if (res > 0) {
			return res;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool InsideNetwork_old(double pos[3])
{
	int ie;

	for (ie=0; ie<NP->ne; ie++) {
		if (InsideVessel(ie, pos)) {
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int InsideNetwork(double pos[3])
{
	int i, ie, xb, yb, zb, res;

	xb = int((pos[0] - netmin[0])/dbox);
	yb = int((pos[1] - netmin[1])/dbox);
	zb = int((pos[2] - netmin[2])/dbox);
	for (i=0; i<nvlist[xb][yb][zb]; i++) {
		ie = vlist[xb][yb][zb][i];
		res = InsideVessel(ie, pos);
		if (res > 0) {
			return res;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
// A minimum diameter is imposed.
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile(char *amFile, NETWORK *net, float dmin)
{
	int i, j, k, kp, npts;
	int np_used, ne_used;
	EDGE edge;
	char line[STR_LEN];

	printf("ReadAmiraFile: %s\n",amFile);
	fprintf(fpout,"ReadAmiraFile: %s\n",amFile);
	FILE *fpam = fopen(amFile,"r");
	if (fpam == NULL) {
		printf("Open failure\n");
		return 1;
	}
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
	net->point = (NDPOINT *)malloc(net->np*sizeof(NDPOINT));
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
					net->edgeList[i].pt[0] = net->edgeList[i].vert[0];
					net->edgeList[i].pt[net->edgeList[i].npts-1] = net->edgeList[i].vert[1];
				}
				printf("Got edge npts, total: %d\n",npts);
			} else if (k == 4) {
				for (i=0;i<net->ne;i++) {
					edge = net->edgeList[i];
					double len = 0;
					for (k=0;k<edge.npts;k++) {
						if (fgets(line, STR_LEN, fpam) == NULL) {
							printf("ERROR reading section @4\n");
							return 1;
						}
						if (k > 0 && k<edge.npts-1) {
							sscanf(line,"%f %f %f",&net->point[kp].x,&net->point[kp].y,&net->point[kp].z);
							net->edgeList[i].pt[k] = kp;
							net->edgeList[i].pt_used[k] = kp;
							kp++;
						}
						if (k > 0) {
							len = len + netdist(net,net->edgeList[i].pt[k-1],net->edgeList[i].pt[k]);
						}
					}
					net->edgeList[i].length_vox = len;
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
						sscanf(line,"%f",&net->point[j].d);
						if (net->point[j].d == 0) {
							printf("Error: ReadAmiraFile: zero diameter: i: %d npts: %d k: %d j: %d\n",i,edge.npts,k,j);
							return 1;
						}
						net->point[j].d = MAX(net->point[j].d,dmin);
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
//-----------------------------------------------------------------------------------------------------
void ShowEdge(int ie)
{
	int npts;
	EDGE *edge;
	NDPOINT p;

	edge = &(NP->edgeList[ie]);
	npts = edge->npts;
	printf("edge: %d  npts: %d\n",ie,npts);
	for (int iseg=0; iseg<npts; iseg++) {
		p = NP->point[edge->pt[iseg]];
		printf("iseg: %d  %6.1f %6.1f %6.1f  diam: %6.2f\n",iseg,p.x,p.y,p.z,p.d);
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool inlist(int ie, int *list, int nlist)
{
	if (nlist == 0) 
		return false;
	for (int i=0; i<nlist; i++) {
		if (list[i] == ie) return true;
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void CreateBoxList()
{
	int i, xb, yb, zb, ie, iseg, kbox, k1[3], k2[3];
	EDGE *edge;
	float range[3], minrange, maxrange;

	dbox = 100;
	maxvlist = 100;
	SetEdgeBnds(netmin,netmax);
	/*
	iesel = 42;
	int ie = iesel;
	for (int i=0; i<3; i++) {
		netmin[i] = NP->edgeList[ie].emin.pos[i];
		netmax[i] = NP->edgeList[ie].emax.pos[i];
	}
	ShowEdge(ie);
	*/

	for (i=0; i<3; i++) {
		netmin[i] -= 5;	// int(netmin[i] - 5);
		netmax[i] += 5;	// int(netmax[i] + 5);
	}

	printf("Network bounds:\n");
	printf("x:  %6.1f  %6.1f\n",netmin[0],netmax[0]);
	printf("y:  %6.1f  %6.1f\n",netmin[1],netmax[1]);
	printf("z:  %6.1f  %6.1f\n",netmin[2],netmax[2]);
	minrange = 1.0e10;
	maxrange = -1.0e10;
	for (i=0; i<3; i++) {
		range[i] = netmax[i] - netmin[i];
		minrange = MIN(minrange,range[i]);
		maxrange = MAX(maxrange,range[i]);
	}
	for (i=0; i<3; i++) {
		N[i] = int(range[i]/delta + 1);
		Nbox[i] = range[i]/dbox + 1;
		printf("N, Nbox: %d  %6d  %4d delta: %f\n",i,N[i],Nbox[i],delta);
	}

	// Allocate arrays:
	nvlist = new int**[Nbox[0]];
	vlist = new int***[Nbox[0]];
	for (xb=0; xb<Nbox[0]; xb++) {
		nvlist[xb] = new int*[Nbox[1]];
		vlist[xb] = new int**[Nbox[1]];
		for (yb=0; yb<Nbox[1]; yb++) {
			nvlist[xb][yb] = new int[Nbox[2]];
			vlist[xb][yb] = new int*[Nbox[2]];
			for (zb=0; zb<Nbox[2]; zb++) {
				vlist[xb][yb][zb] = new int[maxvlist];
			}
		}
	}
	for (xb=0; xb<Nbox[0]; xb++) {
		for (yb=0; yb<Nbox[1]; yb++) {
			for (zb=0; zb<Nbox[2]; zb++) {
				nvlist[xb][yb][zb] = 0;
			}
		}
	}

	// Set up vessel lists
	// Note that segment iseg of edge spans edge->smin[iseg].pos[i] -> edge->smax[iseg].pos[i]
	// We need to create lists of edges that intersect boxes
	int nsum = 0;
	for (ie=0; ie<NP->ne; ie++) {
		edge = &(NP->edgeList[ie]);
		for (i=0; i<3; i++) {
			k1[i] = 999;
			k2[i] = -999;
		}
		for (iseg=0; iseg<edge->npts-1; iseg++) {
			for (i=0; i<3; i++) {
				kbox = int((edge->smin[iseg].pos[i] - netmin[i])/dbox);
				k1[i] = MIN(k1[i],kbox);
				kbox = int((edge->smax[iseg].pos[i] - netmin[i])/dbox);
				k2[i] = MAX(k2[i],kbox);
			}
		}
		for (xb=k1[0]; xb<=k2[0]; xb++) {
			for (yb=k1[1]; yb<=k2[1]; yb++) {
				for (zb=k1[2]; zb<=k2[2]; zb++) {
					if (!inlist(ie,vlist[xb][yb][zb],nvlist[xb][yb][zb])) {
						vlist[xb][yb][zb][nvlist[xb][yb][zb]] = ie;
						nvlist[xb][yb][zb]++;
						nsum++;
					}
				}
			}
		}
	}
	xb = Nbox[0]/2;
	yb = Nbox[1]/2;
//	zb = Nbox[2]/2;
	for (zb=0; zb<Nbox[2]; zb++) {
		printf("%6d %4d\n",zb,nvlist[xb][yb][zb]);
	}
	printf("nsum: %d\n",nsum);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int err;
	int x, y, z, val, res;
	long long width, height, depth, xysize;
	double pos[3];
	char *input_SGfile;
	char drive[32], dir[128],filename[256], ext[32];
	char errfilename[256], output_tiff[256], result_file[256];
	char output_basename[256];
	bool use_compression = true;

	if (argc != 5) {
		printf("Usage: SG2tiff input_SGfile output_tiff delta dmin\n");
		fperr = fopen("SG2tiff_error.log","w");
		fprintf(fperr,"Usage: SG2tiff input_SGfile output_tiff delta dmin\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_SGfile = argv[1];
	strcpy(output_tiff,argv[2]);
	sscanf(argv[3],"%f",&delta);
	sscanf(argv[4],"%f",&dmin);
	_splitpath(output_tiff,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_SG2tiff.log",output_basename);
	sprintf(result_file,"%s_SG2tiff.out",output_basename);
	fperr = fopen(errfilename,"w");

	printf("delta, dmin: %f %f\n",delta,dmin);

	fpout = fopen(result_file,"w");	
	NP = (NETWORK *)malloc(sizeof(NETWORK));
	err = ReadAmiraFile(input_SGfile,NP,dmin);
	if (err != 0) {
		printf("Error: reading SpatialGraph file: %s\n",input_SGfile);
		return 2;
	}

	CreateBoxList();
	width = N[0];
	height = N[1];
	depth = N[2];

	//width = 578;
	//height = 563;
	//depth = 571;

	printf("width: %d height: %d depth: %d\n",width,height,depth);
	ImageType::Pointer im = ImageType::New();
	ImageType::SizeType imsize; 
	imsize[0] = width;
	imsize[1] = height;
	imsize[2] = depth;
	ImageType::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	ImageType::RegionType imregion; 
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	im->SetRegions(imregion);
	im->Allocate();
	p = (unsigned char *)(im->GetBufferPointer());
	xysize = width*height;

	printf("Created buffer\n");

	dbug = false;
	int n = 0;
	for (x=0; x<N[0]; x++) {
		pos[0] = netmin[0] + x*delta;
		printf("x: %d  %6.1f\n",x,pos[0]);
		for (y=0; y<N[1]; y++) {
			pos[1] = netmin[1] + y*delta;
//			printf(".");
			for (z=0; z<N[2]; z++) {
				pos[2] = netmin[2] + z*delta;
//	dbug = (x == 40 && y == 8 && z == 10);
	if (dbug) printf("x,y,z: %d %d %d  %6.1f %6.1f %6.1f\n",x,y,z,pos[0],pos[1],pos[2]);
				res = InsideNetwork(pos);
				if (res > 0) {
					n++;
					val = 255;
					//if (res == 1)
					//	val = 255;		// sphere
					//else if (res == 2)
					//	val = 128;		// tube
				} else {
					val = 0;
				}
				V(x,y,z) = val;
			}
		}
//		printf("\n");
	}
	printf("Inside voxels: %d\n",n);

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(output_tiff);
	writer->SetInput(im);
	if (use_compression) {
		writer->UseCompressionOn();
	}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 3;
	}
	return 0;
}
