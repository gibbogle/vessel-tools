/*
 * To cut an image into two pieces using a triangulated (non-planar) surface at the cut surface
 */

//using namespace std;

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
#include <algorithm>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkAndImageFilter.h"

/*
#define MAXPAIRS 100
#define MAXLABELS 1000000	// 400000
#define MAXLABLIST 1000		// 100
#define MAXOBJECTS 200000
#define NBIGOBJECTS 5

#define USE_SHORT false
*/

#define MAXEXP 10000000
#define NBRS 26
typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im, im_mask;
int width, height, depth, xysize;
unsigned char *p, *p_mask;
//unsigned short *label;
unsigned int *label;

int **lablist;
int *nlablist;
int maxnlablist;
int nbr[26][3];
int nbr2D[8][2];

struct triangle_str 
{
	double v1[3], v2[3], v3[3];		// v[0] = x, v[1] = y, v[2] = z
};
typedef triangle_str TRIANGLE;

struct cutline_str
{
	int npts;
	int *x, *y;
	int z;
	double *f;
};
typedef cutline_str CUTLINE;

TRIANGLE *triangle;
int ntriangles;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]			// V is used to access the voxel array p[]
#define M(a,b,c)  p_mask[(c)*xysize+(b)*width+(a)]		// M is used to access the voxel array p_mask[]

int max16 = pow(2.,16)-1;

FILE *fp, *fpout;

//------------------------------------------------------------------------------------------------
// Returns true if j is in list[], otherwise false.
//------------------------------------------------------------------------------------------------
bool inlist(int list[], int nlist, int j)
{
	int i;
	bool in = false;
	for (i=0; i<nlist; i++)
	{
		if (list[i] == j)
		{
			in = true;
			break;
		}
	}
	return in;
}


//------------------------------------------------------------------------------------------------
// The value of NBRS determines how two voxels are to be considered neighbours.
// NBRS = 6  The von Neumann neighbourhood
//       18  The Moore18 neighbourhood (excludes longest diagonals)
//       26  The Moore26 neighbourhood - all voxels in the cube are neighbours
// nbr[k][:] is the offset of the kth neighbour
//------------------------------------------------------------------------------------------------
void setupNeighbours()
{
	int k=0;

	if (NBRS == 6)
	{
		nbr[0][0] = -1;	nbr[0][1] =  0;	nbr[0][2] = 0;
		nbr[1][0] =  1;	nbr[1][1] =  0;	nbr[1][2] = 0;
		nbr[2][0] =  0;	nbr[2][1] = -1;	nbr[2][2] = 0;
		nbr[3][0] =  0;	nbr[3][1] =  1;	nbr[3][2] = 0;
		nbr[4][0] =  0;	nbr[4][1] =  0;	nbr[4][2] = -1;
		nbr[5][0] =  0;	nbr[5][1] =  0;	nbr[5][2] =  1;
	}
	else
	{
		for (int x=-1; x<=1; x++)
		{
			for (int y=-1; y<=1; y++)
			{
				for (int z=-1; z<=1; z++)
				{	
					if (x == 0 && y == 0 && z == 0) continue;
					if (NBRS == 18 && abs(x) == 1 && abs(y) == 1 && abs(z) == 1) continue;
					nbr[k][0] = x;
					nbr[k][1] = y;
					nbr[k][2] = z;
					k++;
				}
			}
		}
	}
//	printf("Number of 3D neighbours: %d\n",k);
	k = 0;
	for (int x=-1; x<=1; x++)
	{
		for (int y=-1; y<=1; y++)
		{
			if (x == 0 && y == 0) continue;
			nbr2D[k][0] = x;
			nbr2D[k][1] = y;
			k++;
		}
	}
//	printf("Number of 2D neighbours: %d\n",k);
}

//------------------------------------------------------------------------------------------------
// We know that V(x,y,z) > 0, i.e. the voxel (x,y,z) is lit
// Look for a neighbour with a label to use.
// klist[] is the list of lit neighbour voxels not previously labelled.
// If (x,y,z) has a label, it is given to voxels in klist[].
// If a neighbour voxel (xx,yy,zz) already has a label, then
//    if (x,y,z) has no label, it gets given the same label as (xx,yy,zz)
//    else, if the label of (x,y,z) and the neighbour (xx,yy,zz) are different, the two labels 
//    are added to the pair lists. 
//------------------------------------------------------------------------------------------------
void explore(int x, int y, int z, int *n, int klist[])
{
	int c1, c2, k, dx, dy, dz, xx, yy, zz, v2, nn;

//	c1 = M(x,y,z);		// the label of voxel (x,y,z) which must be 255
	nn = 0;
	for (k=0; k<NBRS; k++) {
		dx = nbr[k][0];
		dy = nbr[k][1];
		dz = nbr[k][2];
		xx = x+dx;
		if (xx < 0 || xx >= width) continue;
		yy = y+dy;
		if (yy < 0 || yy >= height) continue;
		zz = z+dz;
		if (zz < 0 || zz >= depth) continue;
	
		c2 = M(xx,yy,zz);
		if (c2 == 128)	{ // lit neighbour voxel at (xx,yy,zz)
			M(xx,yy,zz) = 255;
			klist[nn] = k;
			nn++;
		}
	}
	*n = nn;
}
/*
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void connecter(void)
{
	int x, y, z, xx, yy, zz, x0, y0, z0, c1, k, kk, n, iexp, nexp;
	int klist[26];
	int **explist;
	unsigned char v1;

	setupNeighbours();
	printf("connecter:\n");
	fprintf(fpout,"connecter:\n");
	fflush(fpout);
	explist = (int **)malloc(MAXEXP*sizeof(int *));
	for (k=0; k<MAXEXP; k++)
		explist[k] = (int *)malloc(3*sizeof(int));

}
*/
/*
//------------------------------------------------------------------------------------------------
// First create a lookup table for all used label values, such that those connected to lab map
// to 255, and the rest map to 0.  This requires calling combine, then using nlist and list[].
// Then the construction of the desired image is simple: each non-zero C(x,y,z) maps to 255 or 0.
//------------------------------------------------------------------------------------------------
int createObjectImage(int lab, int cmap[], bool taken[])
{
	int x, y, z, n, k, total;
	int nlist;
	int list[MAXLABELS];
	unsigned char *labelmap;

	fprintf(fpout,"createObjectImage: %d\n",lab);
	fflush(fpout);
	combine(lab, list, &nlist, cmap, taken);
	total = 0;
	for (k=0; k<nlist; k++)
		total += cnt[list[k]];
	printf("lab: %d nlist: %d total: %d\n",lab,nlist,total);
	fprintf(fpout,"lab: %d nlist: %d total: %d\n",lab,nlist,total);
	fflush(fpout);

	labelmap = (unsigned char *)malloc((cmax+1)*sizeof(unsigned char));
	for (k=0; k<=cmax; k++)
		labelmap[k] = 0;
	for (k=0; k<nlist; k++)
		labelmap[list[k]] = 255;

	n = 0;
	for (x=0; x<width; x++)
	{
		for (y=0; y<height; y++)
		{
			for (z=0; z<depth; z++)
			{
				if (C(x,y,z) != 0)
				{
					V(x,y,z) = labelmap[C(x,y,z)];
					if (V(x,y,z) != 0)
						n++;
				} 
				else
				{
					V(x,y,z) = 0;
				}
			}
		}
	}
	free(labelmap);
	return n;
}
*/

//------------------------------------------------------------------------------------------------
// Flood fill mask V starting from (x0,y0,z0), setting inside to 255
// The start point must be a corner of the image! 
// x0 must be either 0 or width-1,
// y0 must be either 0 or height-1, 
// z0 must be either 0 or depth-1
//------------------------------------------------------------------------------------------------
void flood3D(int x0, int y0, int z0)
{
	int x, y, z, dx, dy, dz, xx, yy, zz, k,kk,  n;
	int delx, dely, delz;
	int klist[26];
	short **explist;
	int kexp, kadd, nexp;

	setupNeighbours();
	printf("flood3D:\n");
	fprintf(fpout,"flood3D:\n");
	fflush(fpout);
	explist = (short **)malloc(MAXEXP*sizeof(short *));
	for (k=0; k<MAXEXP; k++)
		explist[k] = (short *)malloc(3*sizeof(short));

	nexp = 0;
	M(x0,y0,z0) = 255;		// This is the seed
	explist[0][0] = x0;
	explist[0][1] = y0;
	explist[0][2] = z0;
	nexp++;
	kexp = 0;
	kadd = 1;
	for (;;) {
		if (kexp%1000000 == 0) printf("kexp: %d kadd: %d nexp: %d\n",kexp,kadd,nexp);
		x = explist[kexp][0];
		y = explist[kexp][1];
		z = explist[kexp][2];
		kexp++;
		if (kexp == MAXEXP) kexp = 0;
		explore(x,y,z,&n,klist);
		nexp--;
		if (n > 0) {
			for (k=0; k<n; k++) {
				kk = klist[k];
				xx = x + nbr[kk][0];
				yy = y + nbr[kk][1];
				zz = z + nbr[kk][2];
				explist[kadd][0] = xx;
				explist[kadd][1] = yy;
				explist[kadd][2] = zz;
				kadd++;
				if (kadd == MAXEXP) kadd = 0;
				if (kexp == kadd) {
					printf("Error: flood3D: need to increase MAXEXP\n");
					exit(1);
				}
				nexp++;
			}
		}
		if (nexp == 0) break;
	}
	printf("Completed flood3D\n");
}

void oldflood3D(int x0, int y0, int z0)
{
	int x, y, z, dx, dy, dz, xx, yy, zz, k;
	int delx, dely, delz;

	if (x0 == 0) {
		delx = 1;
	} else if (x0 == width-1) {
		delx = -1;
	} else {
		printf("Error: x0 = %d.  Must be 0 or width-1\n",x0);
		exit(1);
	}
	if (y0 == 0) {
		dely = 1;
	} else if (y0 == height-1) {
		dely = -1;
	} else {
		printf("Error: y0 = %d.  Must be 0 or height-1\n",y0);
		exit(1);
	}
	if (z0 == 0) {
		delz = 1;
	} else if (z0 == depth-1) {
		delz = -1;
	} else {
		printf("Error: z0 = %d.  Must be 0 or depth-1\n",z0);
		exit(1);
	}

	M(x0,y0,z0) = 255;
//	for (x=0; x<width; x++) {
	x = x0;
	for (;;) {
//		for (y=0; y<height; y++) {
		y = y0;
		for (;;) {
//			for (z=0; z<depth; z++) {
			z = z0;
			for (;;) {
				if (M(x,y,z) == 255) {
					for (k=0; k<NBRS; k++) {
						dx = nbr[k][0];
						dy = nbr[k][1];
						dz = nbr[k][2];
						xx = x+dx;
						if (xx < 0 || xx >= width) continue;
						yy = y+dy;
						if (yy < 0 || yy >= height) continue;
						zz = z+dz;
						if (zz < 0 || zz >= depth) continue;

						if (M(xx,yy,zz) == 128) {
							// light neighbour voxel at (xx,yy,zz)
							M(xx,yy,zz) = 255;
						}
					}
				}
				z += delz;
				if (z < 0 || z >= depth) break;
			}
			y += dely;
			if (y < 0 || y >= height) break;
		}
		x += delx;
		if (x < 0 || x >= width) break;
	}
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void binarise()
{
	int x, y, z;

	for (x=0; x<width; x++) {
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				if (M(x,y,z) != 255) M(x,y,z) = 0;
			}
		}
	}
}

void normalise(double *u)
{
	double sum;

	sum = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
	sum = sqrt(sum);
	u[0] /= sum;
	u[1] /= sum;
	u[2] /= sum;
}

//------------------------------------------------------------------------------------------------
// Generate the rotation matrix R for rotation by alpha about the vector v[]
//------------------------------------------------------------------------------------------------
void rotationMatrix(double alpha, double v[], double R[][3])
{
	double c, s, c1;

	c = cos(alpha);
	s = sin(alpha);
	c1 = 1 - c;
	R[0][0] = c + v[0]*v[0]*c1;
	R[1][0] = v[1]*v[0]*c1 + v[2]*s;
	R[2][0] = v[2]*v[0]*c1 - v[1]*s;
	R[0][1] = v[0]*v[1]*c1 - v[2]*s;
	R[1][1] = c + v[1]*v[1]*c1;
	R[2][1] = v[2]*v[1]*c1 + v[0]*s;
	R[0][2] = v[0]*v[2]*c1 + v[1]*s;
	R[1][2] = v[1]*v[2]*c1 - v[0]*s;
	R[2][2] = c + v[2]*v[2]*c1;
}

//------------------------------------------------------------------------------------------------
// u2[] and u3[3] are unit vectors
// We find the rotation matrix A that takes u2 to the z axis, and takes u3 to the yz plane
// This fails if u2 is already aligned with the z axis
//------------------------------------------------------------------------------------------------
void getTransform(double *u2, double *u3, double A[][3])
{
	int i, j, k;
	double R1[3][3], R2[3][3];
	double zaxis[3], v[3], u3r[3], vamp, sum;
	double alpha, cosa, sina;
	double PI;
	double tol = 0.001;

	PI = 4.0*atan(1.0);
//	printf("\n");
//	printf("u2:  %6.3f %6.3f %6.3f\n",u2[0],u2[1],u2[2]);
//	printf("u3:  %6.3f %6.3f %6.3f\n",u3[0],u3[1],u3[2]);
	// First need to find the transformation matrix R1 to rotate u2 to lie on the z axis
	zaxis[0] = 0;
	zaxis[1] = 0;
	zaxis[2] = 1;
	if (fabs(u2[2] - 1.0) < tol) {	// u2 is aligned with the +z axis
		v[0] = 1;
		v[1] = 0;
		v[2] = 0;
		alpha = 0;
	} else if (fabs(u2[2] + 1.0) < tol) {	// u2 is aligned with -z axis
		v[0] = 1;
		v[1] = 0;
		v[2] = 0;
		alpha = PI;
	} else {
	//	crossproduct(u2,zaxis,v);	// this gives the vector v that is the axis of rotation.
		v[0] = u2[1];
		v[1] = -u2[0];
		v[2] = 0;
		normalise(v);
	//	printf("v:  %6.3f %6.3f %6.3f\n",v[0],v[1],v[2]);
	//	cosa = dotproduct(u2,zaxis);
		cosa = u2[2];
		alpha = acos(cosa);
	//	printf("cosa, alpha: %f %f %f\n",cosa,alpha,alpha*180/PI);
	}
	// R1 = rotation matrix to rotate by alpha about v
	rotationMatrix(alpha,v,R1);
	// check!
	v[0] = R1[0][0]*u2[0] + R1[0][1]*u2[1] + R1[0][2]*u2[2];
	v[1] = R1[1][0]*u2[0] + R1[1][1]*u2[1] + R1[1][2]*u2[2];
	v[2] = R1[2][0]*u2[0] + R1[2][1]*u2[1] + R1[2][2]*u2[2];
//	printf("Rotated u2: %f %f %f\n",v[0],v[1],v[2]);
	// u3 -> u3r
	u3r[0] = R1[0][0]*u3[0] + R1[0][1]*u3[1] + R1[0][2]*u3[2];
	u3r[1] = R1[1][0]*u3[0] + R1[1][1]*u3[1] + R1[1][2]*u3[2];
	u3r[2] = R1[2][0]*u3[0] + R1[2][1]*u3[1] + R1[2][2]*u3[2];
//	printf("Rotated u3 -> u3r: %6.3f %6.3f %6.3f\n",u3r[0],u3r[1],u3r[2]);
	// Now we need to find the transformation matrix R2 to rotate u3r about the z axis to lie in the yz plane
	cosa = u3r[1]/sqrt(u3r[0]*u3r[0] + u3r[1]*u3r[1]);
	alpha = acos(cosa);
	if (u3r[0] < 0) alpha *= -1;
//	printf("cosa, alpha: %f %f %f\n",cosa,alpha,alpha*180/PI);
	// R2 = rotation matrix to rotate by alpha about zaxis
	rotationMatrix(alpha,zaxis,R2);
	// check!
//	v[0] = R2[0][0]*u3r[0] + R2[0][1]*u3r[1] + R2[0][2]*u3r[2];
//	v[1] = R2[1][0]*u3r[0] + R2[1][1]*u3r[1] + R2[1][2]*u3r[2];
//	v[2] = R2[2][0]*u3r[0] + R2[2][1]*u3r[1] + R2[2][2]*u3r[2];
//	printf("Rotated u3r: %f %f %f\n",v[0],v[1],v[2]);
	// The total rotation A = R2.R1
	for (i=0; i<3; i++) {
		for (j=0; j<3; j++) {
			sum = 0;
			for (k=0; k<3; k++) {
				sum += R2[i][k]*R1[k][j];
			}
			A[i][j] = sum;
		}
	}
}

//------------------------------------------------------------------------------------------------
// Check if p lies within the triangle with vertices P1, P2, P3.
// We need to determine if the vertex ordering is anticlockwise or clockwise.
// Since P1 and P2 both lie on the z axis, the order is clockwise if P3[1] > 0, anticlockwise otherwise.
// Note that instead of (x,y) we now have (y,z)
// If A is (y1,z1) and B is (y2,z2) then a point (y,z) lies on AB if
// E = (y-y1)(z2-z1) - (z-z1)(y2-y1) = 0
// In the anticlockwise case: (y,z) is to the left (inside) if E < 0
// In the clockwise case: (y,z) is to the right (inside) if E > 0
//------------------------------------------------------------------------------------------------
bool inside(double p[], double P1[], double P2[], double P3[])
{
	bool clockwise;
	double E;

	clockwise = (P3[1] > 0);
	E = (p[1] - P1[1])*(P2[2] - P1[2]) - (p[2] - P1[2])*(P2[1] - P1[1]);
	if ((clockwise && E < 0) || (!clockwise && E > 0)) return false;
	E = (p[1] - P2[1])*(P3[2] - P2[2]) - (p[2] - P2[2])*(P3[1] - P2[1]);
	if ((clockwise && E < 0) || (!clockwise && E > 0)) return false;
	E = (p[1] - P3[1])*(P1[2] - P3[2]) - (p[2] - P3[2])*(P1[1] - P3[1]);
	if ((clockwise && E < 0) || (!clockwise && E > 0)) return false;
	return true;
}

//------------------------------------------------------------------------------------------------
// All mask voxels "near" the surface of the triangle are set to 0.
// Problem: triangle ordering affects the result!
//------------------------------------------------------------------------------------------------
void processTriangle(TRIANGLE *t) 
{
	int x, y, z, i, k, dx, dy, dz, xx, yy, zz;
	int xmin, xmax, ymin, ymax, zmin, zmax;
	double u2[3], u3[3], v[3], P1[1], P2[3], P3[3], p[3], d;
	double shift[3];	// this is the translation vector to take P1 -> (0,0,0)
	double A[3][3];		// this is the rotation matrix
	double dlim = 0.75;
	bool dbug = false;;

//	printf("v1: %6.1f %6.1f %6.1f\n",t->v1[0],t->v1[1],t->v1[2]);
//	printf("v2: %6.1f %6.1f %6.1f\n",t->v2[0],t->v2[1],t->v2[2]);
//	printf("v3: %6.1f %6.1f %6.1f\n",t->v3[0],t->v3[1],t->v3[2]);
	for (i=0; i<3; i++) {
		shift[i] = -t->v1[i];
		u2[i] = t->v2[i] + shift[i];
		u3[i] = t->v3[i] + shift[i];
	}
//	printf("shift: %6.1f %6.1f %6.1f\n",shift[0],shift[1],shift[2]);
//	printf("u2: %6.1f %6.1f %6.1f\n",u2[0],u2[1],u2[2]);
//	printf("u3: %6.1f %6.1f %6.1f\n",u3[0],u3[1],u3[2]);
	// we need to add shift[] to a point before rotating it
	normalise(u2);
	normalise(u3);
//	printf("normalised\n");
//	printf("u2: %6.3f %6.3f %6.3f\n",u2[0],u2[1],u2[2]);
//	printf("u3: %6.3f %6.3f %6.3f\n",u3[0],u3[1],u3[2]);
	// Now u2 and u3 are relative to P1 as origin
	// We now find the rotation matrix A that takes u2 to the z axis and u3 to the yz plane
	getTransform(u2,u3,A);
	// Find the transformed triangle vertices
	for (i=0; i<3; i++) {
		u2[i] = t->v2[i] + shift[i];
		u3[i] = t->v3[i] + shift[i];
	}
	P1[0] = 0;
	P1[1] = 0;
	P1[2] = 0;
	P2[0] = A[0][0]*u2[0] + A[0][1]*u2[1] + A[0][2]*u2[2];
	P2[1] = A[1][0]*u2[0] + A[1][1]*u2[1] + A[1][2]*u2[2];
	P2[2] = A[2][0]*u2[0] + A[2][1]*u2[1] + A[2][2]*u2[2];
//	printf("Rotated u2 -> P2: %f %f %f\n",P2[0],P2[1],P2[2]);
	P3[0] = A[0][0]*u3[0] + A[0][1]*u3[1] + A[0][2]*u3[2];
	P3[1] = A[1][0]*u3[0] + A[1][1]*u3[1] + A[1][2]*u3[2];
	P3[2] = A[2][0]*u3[0] + A[2][1]*u3[1] + A[2][2]*u3[2];
//	printf("Rotated u3 -> P3: %f %f %f\n",P3[0],P3[1],P3[2]);

	xmin = ymin = zmin = 9999999;
	xmax = ymax = zmax = -1;
	xmin = MIN(xmin,int(t->v1[0]));
	xmin = MIN(xmin,int(t->v2[0]));
	xmin = MIN(xmin,int(t->v3[0]));
	xmax = MAX(xmax,int(t->v1[0]));
	xmax = MAX(xmax,int(t->v2[0]));
	xmax = MAX(xmax,int(t->v3[0]));
	xmax = xmax+1;
	xmax = MIN(xmax,width-1);
	ymin = MIN(ymin,int(t->v1[1]));
	ymin = MIN(ymin,int(t->v2[1]));
	ymin = MIN(ymin,int(t->v3[1]));
	ymax = MAX(ymax,int(t->v1[1]));
	ymax = MAX(ymax,int(t->v2[1]));
	ymax = MAX(ymax,int(t->v3[1]));
	ymax = ymax+1;
	ymax = MIN(ymax,height-1);
	zmin = MIN(zmin,int(t->v1[2]));
	zmin = MIN(zmin,int(t->v2[2]));
	zmin = MIN(zmin,int(t->v3[2]));
	zmax = MAX(zmax,int(t->v1[2]));
	zmax = MAX(zmax,int(t->v2[2]));
	zmax = MAX(zmax,int(t->v3[2]));
	zmax = zmax+1;
	zmax = MIN(zmax,depth-1);
//	printf("x: %6d %6d\n",xmin,xmax);
//	printf("y: %6d %6d\n",ymin,ymax);
//	printf("z: %6d %6d\n",zmin,zmax);
	// Now search within these bounds
	for (x=xmin; x<=xmax; x++) {
		for (y=ymin; y<=ymax; y++) {
			for (z=zmin; z<=zmax; z++) {
		//		dbug = (z == 80 && x == 400);
				// see if (x,y,z) is close to the triangle
				// First translate the point
				v[0] = x + shift[0];
				v[1] = y + shift[1];
				v[2] = z + shift[2];
				// Next carry out the rotation with A
				p[0] = A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
				p[1] = A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
				p[2] = A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2];
				if (dbug) {
					printf("x,y,z: %6d %6d %6d  p: %6.1f %6.1f %6.1f\n",x,y,z,p[0],p[1],p[2]);
				}
				// p is the transformed point 
				if (inside(p,P1,P2,P3)) {
					d = fabs(p[0]);
					if (dbug) {
						printf("inside: %6.3f\n",d);
					}
					if (d <= dlim) {		// this is a voxel close to the surface
						M(x,y,z) = 0;
						// Zero neighbours
						for (k=0; k<NBRS; k++) {
							dx = nbr[k][0];
							dy = nbr[k][1];
							dz = nbr[k][2];
							xx = x+dx;
							if (xx < 0 || xx >= width) continue;
							yy = y+dy;
							if (yy < 0 || yy >= height) continue;
							zz = z+dz;
							if (zz < 0 || zz >= depth) continue;
							M(xx,yy,zz) = 0;
						}
					}
				}
			}
		}
	}
}

//------------------------------------------------------------------------------------------------
// Cut lines, each a list of 1-based voxel locations, are read from an input file.
// Triangles are created from these locations.
//------------------------------------------------------------------------------------------------
void createTriangles(char *cutfile, int prnflag)
{
	int i, j, nlines, npts, x, y, z, n;
	int i1, i2, ip1, ip2, np1, np2, ntri, itri, ip, last;
	int nxval, nzval, *xval, *yval;
	double dx, dy, d[100], sum, dsum;
	double *f1, *f2;
	double *v;
	CUTLINE *cutline;
	FILE *fplines = NULL;
	char line[1024];
	TRIANGLE tri;
	int vert[100][3];
	bool prn, done;

	prn = (prnflag == 1);
	fplines = fopen(cutfile,"r");
	if (prn) {
		n = fscanf(fplines,"%d",&nxval);
		npts = nxval;
		xval = (int *)malloc(nxval*sizeof(int));
		yval = (int *)malloc(nxval*sizeof(int));
		// read xval[];
		for (i=0; i<nxval; i++) {
			n = fscanf(fplines,"%d",&xval[i]);
		}
	}
	// read nlines
	n = fscanf(fplines,"%d",&nlines);
	printf("nlines: %d  %d\n",n,nlines);
	cutline = (CUTLINE *)malloc(nlines*sizeof(CUTLINE));

	for (i=0; i<nlines; i++) {
		// read npts, z (1-based)
		if (!prn) {
			n = fscanf(fplines,"%d %d",&npts,&z);
			printf("npts,z: %d  %d %d\n",n,npts,z);
			cutline[i].npts = npts;
			cutline[i].z = z-1;
			cutline[i].x = (int *)malloc(npts*sizeof(int));
			cutline[i].y = (int *)malloc(npts*sizeof(int));
			cutline[i].f = (double *)malloc(npts*sizeof(double));
			for (j=0; j<npts; j++) {
				n = fscanf(fplines, "%d %d",&x, &y);
//				printf("x,y: %d  %d %d %d\n",n,j,x,y);
				cutline[i].x[j] = x-1;
				cutline[i].y[j] = y-1;
			}
		} else {
			n = fscanf(fplines,"%d",&z);
			cutline[i].npts = npts;
			cutline[i].z = z-1;
			cutline[i].x = (int *)malloc(npts*sizeof(int));
			cutline[i].y = (int *)malloc(npts*sizeof(int));
			cutline[i].f = (double *)malloc(npts*sizeof(double));
			for (j=0; j<nxval; j++) {
				n = fscanf(fplines,"%d",&yval[j]);
				cutline[i].x[j] = xval[j]-1;
				cutline[i].y[j] = yval[j]-1;
			}
		}
		sum = 0;
		d[0] = 0;
		for (j=1; j<npts; j++) {
			dx = cutline[i].x[j] - cutline[i].x[j-1];
			dy = cutline[i].y[j] - cutline[i].y[j-1];
			d[j] = sqrt(dx*dx+dy*dy);
			sum += d[j];
		}		
		dsum = 0;
		for (j=0; j<npts-1; j++) {
			dsum += d[j];
			cutline[i].f[j] = dsum/sum;
		}
		cutline[i].f[npts-1] = 1;
	}
	printf("read cutlines\n");
	fclose(fplines);

	triangle = (TRIANGLE *)malloc(100000*sizeof(TRIANGLE));
	printf("\nTriangles:\n");
	ntriangles = 0;
//	i1 = 0;
//	i2 = i1+1;
	for (i1=0; i1<nlines-1; i1++) {
		i2 = i1+1;
		//printf("cutline: %d npts: %d\n",i1,cutline[i1].npts);
		//for (j=0; j<cutline[i1].npts; j++) {
		//	printf("%d %6d %6d %6d  %6.3f\n",j,cutline[i1].x[j],cutline[i1].y[j],cutline[i1].z,cutline[i1].f[j]);
		//}
		//printf("\n");
		//printf("cutline: %d npts: %d\n",i2,cutline[i1].npts);
		//for (j=0; j<cutline[i2].npts; j++) {
		//	printf("%d %6d %6d %6d  %6.3f\n",j,cutline[i2].x[j],cutline[i2].y[j],cutline[i2].z,cutline[i2].f[j]);
		//}
		//printf("\n");
		np1 = cutline[i1].npts;
		np2 = cutline[i2].npts;
		f1 = cutline[i1].f;
		f2 = cutline[i2].f;
		ip1 = 0;
		ip2 = 0;
		ntri = 0;
		last = 1;
		done = false;
		for (;;) {
			if (ip1 == np1-1) {
				vert[ntri][0] = 1000 + ip2;		// ip>=1000 => line i2
				vert[ntri][1] = ip1;
				vert[ntri][2] = 1000 + (ip2+1);
//				printf("vert (b) done: %d  %d %d\n",ntri,ip1,ip2);
				ip2++;
				ntri++;
				done = true;
			} else if (ip2 == np2-1) {
				vert[ntri][0] = ip1;
				vert[ntri][1] = 1000 + ip2;
				vert[ntri][2] = ip1+1;
//				printf("vert (a) done: %d  %d %d\n",ntri,ip1,ip2);
				ip1++;
				ntri++;
				done = true;
			} else if (f1[ip1+1] < f2[ip2+1]) {
				vert[ntri][0] = ip1;
				vert[ntri][1] = 1000 + ip2;
				vert[ntri][2] = ip1+1;
//				printf("vert (a): %d  %d %d  %f %f\n",ntri,ip1,ip2,f1[ip1+1],f2[ip2+1]);
				ip1++;
				ntri++;
				last = 1;
			} else if (f1[ip1+1] > f2[ip2+1]) {
				vert[ntri][0] = 1000 + ip2;		// -sign => line i2
				vert[ntri][1] = ip1;
				vert[ntri][2] = 1000 + (ip2+1);
//				printf("vert (b): %d  %d %d  %f %f\n",ntri,ip1,ip2,f1[ip1+1],f2[ip2+1]);
				ip2++;
				ntri++;
				last = 2;
			} else {	// a tie (maybe reached the end)
				if (last == 1) {
					vert[ntri][0] = 1000 + ip2;		// -sign => line i2
					vert[ntri][1] = ip1;
					vert[ntri][2] = 1000 + (ip2+1);
//					printf("vert (b): %d  %d %d  %f %f\n",ntri,ip1,ip2,f1[ip1+1],f2[ip2+1]);
					ip2++;
					ntri++;
					last = 2;
				} else {
					vert[ntri][0] = ip1;
					vert[ntri][1] = 1000 + ip2;
					vert[ntri][2] = ip1+1;
//					printf("vert (a): %d  %d %d  %f %f\n",ntri,ip1,ip2,f1[ip1+1],f2[ip2+1]);
					ip1++;
					ntri++;
					last = 1;
				}
			}
			if (done) break;
		}
		for (itri=0; itri<ntri; itri++) {
			printf("triangle: %4d\n",ntriangles);
			for (j=0; j<3; j++) {
				if (j == 0) 
					v = triangle[ntriangles].v1;
				else if (j == 1)
					v = triangle[ntriangles].v2;
				else if (j == 2)
					v = triangle[ntriangles].v3;
				ip = vert[itri][j];
				if (ip >= 1000) {
					i = i2;
					ip -= 1000;
				} else {
					i = i1;
				}
				v[0] = cutline[i].x[ip];
				v[1] = cutline[i].y[ip];
				v[2] = cutline[i].z;
//				printf("  vertex: %d  %8.1f %8.1f %8.1f\n",j,v[0],v[1],v[2]);
			}
			ntriangles++;
		}
	}
	free(cutline);
}

//------------------------------------------------------------------------------------------------
// Note:  The first and last voxels on each cut line drawn on a slice need to lie on the image 
// boundary, i.e. with voxel coord = 0 or width-1, height-1, depth-1.
//------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	time_t t1;
	t1 = time(NULL);
	char maskfile[]= "mask.tif";
	char surfacefile[]= "surface.tif";
	char *infile, *cutfile, *outfile;
	char *baseName;
	char numstr[2];
	char drive[32], dir[128],filename[64], ext[32];
	char errorFile[128], outputFile[128], outputPath[128];
	int prnflag;
	FILE *fperr;
	bool use_triangles = true;
	int itri;
	TRIANGLE *t;

//	createTriangles();
//	return 1;

	if (argc != 5) {
		fperr = fopen("makecut_error.log","w");
		printf("Usage: makecut input_tiff cutline_file prnflag output_tiff\n");
		fprintf(fperr,"Usage: makecut input_tiff cutline_file prnflag output_tiff\n");
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line 
	}

	infile = argv[1];
	cutfile = argv[2];
	sscanf(argv[3],"%d",&prnflag);
	outfile = argv[4];
	fpout = fopen("maskcut.out","w");

	/*
	baseName = argv[2]; 
	_splitpath(infile,drive,dir,filename,ext);
	strcpy(outputPath,drive);
	strcat(outputPath,dir);
	sprintf(errorFile,"%s%s%s",outputPath,filename,"_connect.log");
	printf("Error file: %s\n",errorFile);
	sprintf(outputFile,"%s%s%s",outputPath,filename,"_connect.out");
	printf("Output file: %s\n",outputFile);
	fperr = fopen(errorFile,"w");
	fpout = fopen(outputFile,"w");
	printf("Input image file: %s\n",infile);
	fprintf(fpout,"Input image file: %s\n",infile);
	printf("Base for output image files: %s\n",baseName);
	fprintf(fpout,"Base for output image files: %s\n",baseName);
	sscanf(argv[3],"%d",&minvoxels);
	printf("Applying minvoxels: %d\n",minvoxels);
	fprintf(fpout,"Applying minvoxels: %d\n",minvoxels);
	*/
	/*
	lablist = (int **)malloc(MAXLABELS * sizeof(int *));
	if (lablist == NULL)
	{
		printf("out of memory: lablist\n");
		fprintf(fp,"out of memory: lablist\n");
		fclose(fp);
		return 4;
	}
	*/

	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

	reader->SetFileName(infile);
	try
	{
		printf("Loading input TIFF: %s\n",infile);
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Read error on input file\n");
		fclose(fp);
		return 2;	// Read error on input file
	}

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	p = (unsigned char *)(im->GetBufferPointer());

	im_mask = ImageType::New();
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
	im_mask->SetRegions(imregion);
	im_mask->Allocate();
	p_mask = (unsigned char *)(im_mask->GetBufferPointer());

	setupNeighbours();
// Set all voxels to 128
	memset(p_mask,128,xysize*depth);

	createTriangles(cutfile, prnflag);

	/*
// Create 2 test triangles
	ntriangles = 2;
	triangle = (TRIANGLE *)malloc(ntriangles*sizeof(TRIANGLE));
	t = &triangle[0];
	t->v2[0] = 0;	// swapped 1 <-> 2
	t->v2[1] = 0;
	t->v2[2] = depth/4;
	t->v1[0] = width-1;
	t->v1[1] = 0;
	t->v1[2] = depth/4;
	t->v3[0] = width-1;
	t->v3[1] = height;
	t->v3[2] = depth/2;	// /2
	t = &triangle[1];
	t->v1[0] = 0;
	t->v1[1] = 0;
	t->v1[2] = depth/4;
	t->v2[0] = width-1;
	t->v2[1] = height;
	t->v2[2] = depth/2;	// /2
	t->v3[0] = 0;
	t->v3[1] = height;
	t->v3[2] = depth/2;	// /2
	
	ntriangles = 4;
	triangle = (TRIANGLE *)malloc(ntriangles*sizeof(TRIANGLE));
	t = &triangle[0];
	t->v2[0] = 0;	// swapped vertices 1 <-> 2
	t->v2[1] = 0;
	t->v2[2] = 20;
	t->v1[0] = width-1;
	t->v1[1] = 0;
	t->v1[2] = 50;
	t->v3[0] = width-1;
	t->v3[1] = height/2;
	t->v3[2] = 50;

	t = &triangle[1];
	t->v1[0] = 0;
	t->v1[1] = 0;
	t->v1[2] = 20;
	t->v2[0] = width-1;
	t->v2[1] = height/2;
	t->v2[2] = 50;
	t->v3[0] = 0;
	t->v3[1] = height/2;
	t->v3[2] = 20;

	t = &triangle[2];
	t->v1[0] = 0;
	t->v1[1] = height/2;
	t->v1[2] = 20;
	t->v2[0] = width-1;
	t->v2[1] = height/2;
	t->v2[2] = 50;
	t->v3[0] = width-1;
	t->v3[1] = height-1;
	t->v3[2] = 50;

	t = &triangle[3];
	t->v1[0] = 0;
	t->v1[1] = height/2;
	t->v1[2] = 20;
	t->v2[0] = width-1;
	t->v2[1] = height-1;
	t->v2[2] = 50;
	t->v3[0] = 0;
	t->v3[1] = height-1;
	t->v3[2] = 20;
	*/
	for (itri=0; itri<ntriangles; itri++) {
		printf("triangle: %4d\n",itri);
		printf("  vertex: %8.1f %8.1f %8.1f\n",triangle[itri].v1[0],triangle[itri].v1[1],triangle[itri].v1[2]);
		printf("  vertex: %8.1f %8.1f %8.1f\n",triangle[itri].v2[0],triangle[itri].v2[1],triangle[itri].v2[2]);
		printf("  vertex: %8.1f %8.1f %8.1f\n",triangle[itri].v3[0],triangle[itri].v3[1],triangle[itri].v3[2]);
	}

	if (use_triangles) {
		for (int itri=0; itri<ntriangles; itri++) {
//			itri = 4;
			printf("processTriangle: %d\n",itri);
			processTriangle(&triangle[itri]);
		}
	} else {
	// Create a simple planar cut
		int zcut = depth/2;
		for (int x=0; x<width; x++) {
			for (int y=0; y<height; y++) {
				M(x,y,zcut) = 0;
				M(x,y,zcut+1) = 0;
			}
		}
	}

	// Write out surface file
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetInput(im_mask);
	writer->UseCompressionOn();

	writer->SetFileName(surfacefile);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Write error on surface file\n");
		fclose(fp);
		return 3;	// Write error on surface file
	}
	printf("Created surface image file: %s\n",surfacefile);

	// We need to specify a corner that is within the "to remove" region
	printf("Flooding\n");
	flood3D(0,0,0);	// the UI might accept 1-based numbers (1, width, height, depth) which would then be decremented.
	binarise();

	// Write out mask file
//	typedef itk::ImageFileWriter<ImageType> FileWriterType;
//	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetInput(im_mask);
	writer->UseCompressionOn();

	writer->SetFileName(maskfile);
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Write error on mask file\n");
		fclose(fp);
		return 3;	// Write error on mask file
	}
	printf("Created mask image file: %s\n",maskfile);

	typedef itk::AndImageFilter <ImageType>	AndImageFilterType;
 
	AndImageFilterType::Pointer andFilter = AndImageFilterType::New();
	andFilter->SetInput(0, im);
	andFilter->SetInput(1, im_mask);
	andFilter->Update();
 
	// Write out masked output file
	writer->SetFileName(outfile);
	writer->SetInput(andFilter->GetOutput());
	//typedef itk::ImageFileWriter<ImageType> FileWriterType;
	//FileWriterType::Pointer writer = FileWriterType::New();
	writer->UseCompressionOn();

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Write error on output file\n");
		fclose(fp);
		return 3;	// Write error on output file
	}
	printf("Created masked image file: %s\n",outfile);

	return 0;
}

