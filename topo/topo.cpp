/*
 * To extract topology from a skeletonized .tif  
 */

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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//-------------------------- 
#define REVISED_VERSION true
//--------------------------

#define USE_NEW false
#define USE_HEALING false
#define NEW_PRUNE true

#define MAXVERTICES 200
#define MAXBRANCHES 200000
#define MAXPOINTS 1
#define MAXEDGEPTS 1000
#define MAXEDGES 100000
#define MAXNBRS 20		// was 10
//#define MINSTUB 3		// remove an edge that has an unconnected end if npts <= MINSTUB 
#define MINSEGLEN 3		// segments less than MINSEGLEN are ignored
#define NBOX 800

// If the number of neighbours is not 2, the voxel is a vertex.
// nbrs = 1 => an unconnected vessel end
// nbrs > 2 => a vertex with nbrs links
struct voxel_str {
	int pos[3];
	int nbrs;
	int nbr[MAXNBRS][3];
	int nid[MAXNBRS];
	int ivert;
};
typedef voxel_str VOXEL;

struct vec_str
{
	int dx, dy, dz;
};
typedef vec_str VEC;

struct vecset_str
{
	VEC vector[8];
};
typedef vecset_str VECSET;

struct xyz_str
{
	double x, y, z;
};
typedef xyz_str XYZ;

struct vertex_str
{
	int ivox;	// index into the voxel list voxel[]
	int ivisit;	// index into the visit list visited[] (-1 if not visited yet)
	int pos[3];
	int nlinks;
	int nlinks_used;
	int nfollowed;
	bool *followed;
	int *edge;
	int pt[10];
	bool used;
};
typedef vertex_str VERTEX;

struct edge_str
{
	int vert[2];
	int npts;
	int npts_used;
	int *pt;
	int *pt_used;
//	float *avediameter;
//m	double mindiameter[MAXEDGEPTS];	// try removing mindiameter to save memory
	double segavediam;
	double segmindiam;
	double length_um;	// voxel-voxel length
	bool used;
	bool ok;
};
typedef edge_str EDGE;

struct pair_str
{
	int i1, i2;
};
typedef pair_str PAIR;

struct loosend_str
{
	int iedge;			// edge
	int iend;			// free end (0, 1)
	double dir[3];
	double len;
	double dave;
	int joined;
};
typedef loosend_str LOOSEND;

#define V3Dskel(a,b,c)  pskel[(c)*imsize_xy+(b)*width+(a)]
#define V3D(a,b,c)  p[(c)*imsize_xy+(b)*width+(a)]
typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im, imskel;
long long width, height, depth, imsize_xy;
unsigned char *pskel, *p;
int nv, ne, negmin, np, np_used;
int ne_max, np_max;
int edge, node;
int maxnbrs;
EDGE *edgeList;
VERTEX *vertex;
VOXEL *voxel;
int vertEdge[MAXVERTICES][10], nve[MAXVERTICES];
int point[MAXPOINTS][3];
bool *point_used;
float *avediameter, *mindiameter, *radius;
int nb, branchlist[MAXBRANCHES];
//double voxelsize_xy, voxelsize_z;	// in um
double calib_param;
double vsize[3];	// in um
float FIXED_DIAMETER;
//int volume;
VECSET plane[9];
FILE *fperr, *fpout;
FILE *exelem, *exnode, *dotcom;
char output_basename[512];
bool uniform_diameter;
bool junction_max;
double max_ratio;
bool dbug = false;

#define PI 3.14159

int adjoinEdges();

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float zdist(int k1, int k2)
{
	float dx = vsize[0]*(voxel[k2].pos[0] - voxel[k1].pos[0]);
	float dy = vsize[0]*(voxel[k2].pos[1] - voxel[k1].pos[1]);
	float dz = vsize[0]*(voxel[k2].pos[2] - voxel[k1].pos[2]);
	return sqrt(dx*dx+dy*dy+dz*dz);
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
//-----------------------------------------------------------------------------------------------------
void showvertex(int iv)
{
	int k, ie, npts;
	EDGE *edge;

	printf("vertex: %d nlinks: %d nlinks_used: %d\n",iv,vertex[iv].nlinks, vertex[iv].nlinks_used);
	fprintf(fpout,"vertex: %d nlinks: %d nlinks_used: %d\n",iv,vertex[iv].nlinks, vertex[iv].nlinks_used);
	if (!vertex[iv].used) {
		printf("NOT USED\n");
		fprintf(fpout,"NOT USED\n");
		return;
	}
	for (int k=0; k<vertex[iv].nlinks_used; k++) {
		ie = vertex[iv].edge[k];
		edge = &edgeList[ie];
		npts = edge->npts;
		printf("link: %d edge: %d npts: %d vert: %d %d endpts: %d %d\n",k,ie,npts,
			edge->vert[0], edge->vert[1], edge->pt[0], edge->pt[npts-1]);
		fprintf(fpout,"link: %d edge: %d npts: %d vert: %d %d endpts: %d %d\n",k,ie,npts,
			edge->vert[0], edge->vert[1], edge->pt[0], edge->pt[npts-1]);
	}
}

//-----------------------------------------------------------------------------------------------------
// mode = 'D' to print diameters, 'N' to print node numbers
//-----------------------------------------------------------------------------------------------------
void showedge(int ie, char mode) {
	EDGE *edge;
	int npts, k, kp, kv0, kv1;

	edge = &edgeList[ie];
	npts = edge->npts;
	kv0 = edge->vert[0];
	kv1 = edge->vert[1];
//	printf("\nedge: %5d  vert: %d %d -> %d %d npts: %2d: ",ie,kv0,kv1,vertex[kv0].ivox, vertex[kv1].ivox,npts);
	fprintf(fpout,"\nedge: %5d  vert: %d %d -> %d %d npts: %2d: ",ie,kv0,kv1,vertex[kv0].ivox, vertex[kv1].ivox,npts);
	if (!edge->used) {
//		printf("NOT USED\n");
		fprintf(fpout,"NOT USED\n");
		return;
	}
	for (k=0; k<npts; k++) {
		kp = edge->pt[k];
		if (mode == 'N') {
//			printf("%6d ",kp);
			fprintf(fpout,"%6d ",kp);
		} else {
//			printf("%6.1f",avediameter[kp]);
			fprintf(fpout,"%6.1f",avediameter[kp]);
		}
	}
//	printf("\n");
	fprintf(fpout,"\n");
}

//-----------------------------------------------------------------------------------------------------
// Check consistency between edge vert[] and pt[]
//-----------------------------------------------------------------------------------------------------
int checkEdgeEndPts()
{
	EDGE *edge;
	int ie, npts, k, kp0, kp1, kv0, kv1, kvp0, kvp1, err;

	err = 0;
	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		npts = edge->npts;
		kp0 = edge->pt[0];
		kp1 = edge->pt[npts-1];
		kv0 = edge->vert[0];
		kv1 = edge->vert[1];
		kvp0 = vertex[kv0].ivox;
		kvp1 = vertex[kv1].ivox;
		if ((kp0==kvp0 && kp1==kvp1) || (kp0==kvp1 && kp1==kvp0)) continue;
		err = 1;
		fprintf(fpout,"checkEdgeEndPts: edge: %d vert: %d %d -> %d %d pt: %d %d\n",ie,kv0,kv1,kvp0,kvp1,kp0,kp1);
		fflush(fpout);
	}
	return err;
}

//-----------------------------------------------------------------------------------------------------
// If the value k is found in list (list[0] ... list[n-1]) then the index is returned,
// otherwise -1 is returned.
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

//-----------------------------------------------------------------------------------------------------
// Reorder element and node identifiers
// writeCMGUIfile() does not use any vertex data
// writeAMIRAfile() uses edge.vert[] and vertex[].pos[]
// Need to check for verticies with only two links, join two edges into one.
//-----------------------------------------------------------------------------------------------------
int squeezer(void)
{
	int i, j, k, kv, ivox, kp, kv0, kv1, err;
	int ne_x, nv_x, np_x, knew, i_x;
	EDGE *edge;
	EDGE *edgeList_x;
	VERTEX *vertex_x;
	VOXEL *voxel_x;
	float *avediameter_x;
	int *oldpt;
	bool dbug;

	printf("squeezer: np: %d\n",np);
	vertex_x = (VERTEX *)malloc(nv*sizeof(VERTEX));
	edgeList_x = (EDGE *)malloc(2*ne*sizeof(EDGE));		// 2* for added joining edges
	voxel_x = (VOXEL *)malloc(2*np*sizeof(VOXEL));
	oldpt = (int *)malloc(2*np*sizeof(int));
	avediameter_x = (float *)malloc(2*np*sizeof(float));
	for (kv=0; kv<nv; kv++) {
		vertex_x[kv].edge = (int *)malloc(MAXNBRS*sizeof(int));
	}

	err = 0;
	ne_x = 0;
	nv_x = 0;
	np_x = 0;
	// First process the vertices
	for (i=0; i<ne; i++) {
		if (!edgeList[i].used) continue;
		if (edgeList[i].vert[0] == edgeList[i].vert[1]) {
			edgeList[i].used = false;
			continue;	// repeated pt
		}
		dbug = false;
		edge = &edgeList_x[ne_x];
		edge->used = true;
		for (j=0; j<2; j++) {
			kv = edgeList[i].vert[j];
			ivox = vertex[kv].ivox;
			knew = inlist(oldpt,np_x,ivox);
			if (knew == -1) {
				oldpt[np_x] = ivox;		//edgeList[i].vert[j];						// index in old voxel list
				for (k=0; k<3; k++) {
					vertex_x[np_x].pos[k] = voxel[oldpt[np_x]].pos[k];
				}
				knew = np_x;
				np_x++;
			}
			edge->vert[j] = knew;
			vertex_x[knew].ivox = -1;
		}
		ne_x++;
	}
	nv_x = np_x;
	printf("nv_x: %d\n",nv_x);
	fprintf(fpout,"nv_x: %d\n",nv_x);

// Now add the edge points to the list.

// TRY THIS
	np_x = 0;
//
	i_x = 0;
	for (i=0; i<ne; i++) {
		if (!edgeList[i].used) continue;
//		fprintf(fpout,"old edge: %d npts: %d  new edge: %d npts: %d\n",i,edgeList[i].npts,i_x,edgeList[i].npts_used);
		int npts = edgeList[i].npts_used;
		if (npts < 1) {
			printf("squeezer: i: %d npts: %d\n",i,npts);
			fprintf(fpout,"squeezer: i: %d npts: %d\n",i,npts);
			return 1;
		}
		edgeList_x[i_x].pt = (int *)malloc(npts*sizeof(int));
		edgeList_x[i_x].pt_used = (int *)malloc(npts*sizeof(int));
		edge = &edgeList_x[i_x];
		edge->npts = npts;
		kv0 = edge->vert[0];
		kv1 = edge->vert[1];
		for (k=0; k<npts; k++) {
			j = edgeList[i].pt_used[k];
			if (j > np || j < 0) {
				printf("i: %d k: %d  j: %d\n",i,k,j);
				return 1;
			}
			voxel_x[np_x] = voxel[j];
//			avediameter_x[np_x] = avediameter[j];
//			edge->pt[k] = np_x;

			if (k == 0) {
				if (vertex_x[kv0].ivox < 0) {
					vertex_x[kv0].ivox = np_x;
					edge->pt[k] = np_x;
					avediameter_x[np_x] = avediameter[j];
					np_x++;
				} else {
					edge->pt[k] = vertex_x[kv0].ivox;
				}
			} else if (k == npts-1) {
				if (vertex_x[kv1].ivox < 0) {
					vertex_x[kv1].ivox = np_x;
					edge->pt[k] = np_x;
					avediameter_x[np_x] = avediameter[j];
					np_x++;
				} else {
					edge->pt[k] = vertex_x[kv1].ivox;
				}
			} else {
				edge->pt[k] = np_x;
				avediameter_x[np_x] = avediameter[j];
				np_x++;
			}

			// Check for repeated point
			if (edge->pt[k] == edge->pt[k-1]) {
				printf("Error: squeezer: repeated point: old edge: %d new edge: %d old point: %d new point: %d\n",i,i_x,j,np_x);
				fprintf(fperr,"Error: squeezer: repeated point: old edge: %d new edge: %d old point: %d new point: %d\n",i,i_x,j,np_x);
				err = 1;
			}
		}
		//kv0 = edge->vert[0];
		//kv1 = edge->vert[1];
		//fprintf(fpout,"edge: %d vert: %d %d -> %d %d\n",i_x,kv0,kv1,vertex_x[kv0].ivox,vertex_x[kv1].ivox);
		i_x++;
	}
	printf("Added edge points\n");
	printf("ne, ne_x: %d %d  nv, nv_x: %d %d  np, np_x: %d %d\n",ne,ne_x,nv,nv_x,np,np_x);
	fprintf(fpout,"Added edge points\n");
	fprintf(fpout,"ne, ne_x: %d %d  nv, nv_x: %d %d  np, np_x: %d %d\n",ne,ne_x,nv,nv_x,np,np_x);
	fflush(fpout);

	// Now copy the revised data back into the original arrays
	nv = nv_x;
	for (i=0; i<nv; i++) {
		vertex[i] = vertex_x[i];
	}
	ne = ne_x;
	for (i=0; i<ne; i++) {
		edgeList[i] = edgeList_x[i];
		edge = &edgeList[i];
		edge->npts_used = edge->npts;			// The output functions need to work with or without the squeezer step,
		for (k=0; k<edge->npts; k++)			// therefore need data in the form of "_used"
			edge->pt_used[k] = edge->pt[k];
		// Check for repeated point
		if (edge->pt[1] == edge->pt[0]) {
			printf("Error: squeezer: repeated point: new edge: %d new point: %d\n",i,edge->pt[0]);
			fprintf(fperr,"Error: squeezer: repeated point: new edge: %d new point: %d\n",i,edge->pt[0]);
			err = 2;
		}
		// Set nlinks for verticies
		for (k=0; k<2; k++) {
			kv = edge->vert[k];
			vertex[kv].nlinks++;
			vertex[kv].nlinks_used++;
		}
	}
	np = np_x;

	for (i=0; i<np; i++) {
		voxel[i] = voxel_x[i];
		avediameter[i] = avediameter_x[i];
	}
	free(vertex_x);
	free(edgeList_x); 
	free(voxel_x);
	free(oldpt);
	free(avediameter_x);
	return err;
}

//-----------------------------------------------------------------------------------------------------
// Create CMGUI .com file
//-----------------------------------------------------------------------------------------------------
void write_com(char *fileName)
{
    fprintf(dotcom, "#Lymph node structure for section %s\n\n", fileName);
    fprintf(dotcom, "# Create a material in addition to the default.\n");
    fprintf(dotcom, "gfx cre mat gold ambient 1 0.7 0 diffuse 1 0.7 0 specular 0.5 0.5 0.5 shininess 0.8\n");
    fprintf(dotcom, "# Read in the reticular mesh (group vessels) and hide the axes.\n");
    fprintf(dotcom, "gfx read nodes %s.exnode\n", fileName);
    fprintf(dotcom, "gfx read elements %s.exelem\n", fileName);
    fprintf(dotcom, "# Open the graphics window and turn on perspective.\n");
//    fprintf(dotcom, "gfx mod win 1 view perspective\n");
    fprintf(dotcom, "# Destroy the default lines.\n");
    fprintf(dotcom, "gfx modify g_element vessels lines delete\n");
    fprintf(dotcom, "# The radius of the vessel is stored in component 1 of field\n");
    fprintf(dotcom, "# 'vessel_radius', defined over the elements in the vessels group.\n");
    fprintf(dotcom, "# Now draw spheres using these radii with the following command.\n");
    fprintf(dotcom, "gfx destroy node all\n");
    fprintf(dotcom, "gfx modify g_element vessels general clear\n");
    fprintf(dotcom, "gfx modify g_element vessels cylinders coordinate coordinates tessellation default local circle_discretization 12 radius_scalar vessel_radius scale_factor 1.0 native_discretization NONE select_on material gold selected_material default_selected render_shaded\n");
    fprintf(dotcom, "gfx modify g_element vessels node_points coordinate coordinates local glyph sphere general size \"0*0*0\" centre 0,0,0 font default orientation vessel_radius scale_factors \"2*2*2\" select_on material gold selected_material default_selected\n");
    fprintf(dotcom, "gfx cre win 1\n");
    fprintf(dotcom, "gfx mod win 1 view perspective\n");
}

//-----------------------------------------------------------------------------------------------------
// Write initial section of .exnode file
//-----------------------------------------------------------------------------------------------------
void write_exnode(void)
{
    fprintf(exnode, "Group name: vessels\n");
    fprintf(exnode, " #Fields=2\n");
    fprintf(exnode, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exnode, "  x.  Value index=1, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  y.  Value index=2, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  z.  Value index=3, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exnode, "  1.  Value index=4, #Derivatives=0, #Versions=1\n");
}

//-----------------------------------------------------------------------------------------------------
// Write initial section of .exelem file
//-----------------------------------------------------------------------------------------------------
void write_exelem(void)
{
    fprintf(exelem, "Group name: vessels\n");
    fprintf(exelem, " Shape.  Dimension=1\n");
    fprintf(exelem, " #Scale factor sets= 1\n");
    fprintf(exelem, "  l.Lagrange, #Scale factors= 2\n");
    fprintf(exelem, " #Nodes= 2\n #Fields=2\n");
    fprintf(exelem, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exelem, "   x.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   y.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   z.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
}

//-----------------------------------------------------------------------------------------------------
// Create CMGUI files.
// Required squeezed network description
//-----------------------------------------------------------------------------------------------------
int WriteCmguiData(void)
{
	int k, ie, ip, npts, err;
//	int kv0, kv1;
	EDGE edge;
	char dotcomname[256], exelemname[256], exnodename[256];

	printf("WriteCmguiData: np: %d\n",np);
	fprintf(fperr,"WriteCmguiData: np: %d\n",np);
	fflush(fperr);
	err = 0;
	sprintf(dotcomname,"%s.com.txt",output_basename);
	printf("dotcomname: %s\n",dotcomname);
	fprintf(fperr,"dotcomname: %s\n",dotcomname);
	fflush(fperr);
	sprintf(exelemname,"%s.exelem",output_basename);
	printf("exelemname: %s\n",exelemname);
	fprintf(fperr,"exelemname: %s\n",exelemname);
	fflush(fperr);
	sprintf(exnodename,"%s.exnode",output_basename);
	fprintf(fperr,"Com file: %s exelem file: %s exnode file: %s\n",dotcomname,exelemname,exnodename);
	fflush(fperr);
	dotcom = fopen(dotcomname,"w");
	exelem = fopen(exelemname,"w");
	exnode = fopen(exnodename,"w");
	write_com(output_basename);
	write_exelem();
	printf("created exelem header\n");
	fprintf(fperr,"created exelem header\n");
	write_exnode();
	printf("created exnode header\n");
	fprintf(fperr,"created exnode header\n");
	// First set up radius values for each point.  Note that junction points occur on multiple segments,
	// with different radius values.  Select the maximum in this case.
	radius = (float *)malloc(np*sizeof(float));
	point_used = (bool *)malloc(np*sizeof(bool));
	for (k=0; k<np; k++) {
		point_used[k] = false;
	}
	int kelem = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		npts = edge.npts;
		/*
		int kfrom = edge.pt_used[0];
		int kto = edge.pt_used[npts-1];
		if (kfrom == kto) {
			printf("Error: writecmguidata: repeated node: element: %d npts: %d kv0,kv1: %d %d kfrom: %d\n",ie,npts,edge.vert[0],edge.vert[1],kfrom);
			fprintf(fperr,"repeated node: element: %d npts: %d kv0.kv1: %d %d kfrom: %d\n",ie,npts,edge.vert[0],edge.vert[1],kfrom);
			for (k=0; k<npts; k++)
				fprintf(fperr,"  k: %3d  pt: %6d\n",k,edge.pt_used[k]);
			edgeList[ie].used = false;
			continue;
		}
		radius[kfrom] = MAX(radius[kfrom],avediameter[kfrom]/2);
		radius[kto] = MAX(radius[kto],avediameter[kto]/2);
		point_used[kfrom] = true;
		*/
		for (ip=1; ip<npts; ip++) {
			int k1 = edge.pt_used[ip-1];
			int k2 = edge.pt_used[ip];
			point_used[k1] = true;
			point_used[k2] = true;
			kelem++;
	        fprintf(exelem, "Element: %d 0 0\n", kelem);
	        fprintf(exelem, "  Nodes: %d %d\n", k1+1, k2+1);
	        fprintf(exelem, "  Scale factors: 1 1\n");
		}
	}
	for (k=0; k<np; k++) {
		if (point_used[k]) {	// Note: after squeezing, all points are used.
			fprintf(exnode, "Node: %d\n", k+1);
//			if (REVISED_VERSION) {
				fprintf(exnode, "%6.1f %6.1f %6.1f\n", vsize[0]*voxel[k].pos[0],vsize[1]*voxel[k].pos[1],vsize[2]*voxel[k].pos[2]);
//			} else {
//				fprintf(exnode, "%6.1f %6.1f %6.1f\n", vsize[0]*point[k][0],vsize[1]*point[k][1],vsize[2]*point[k][2]);
//			}
			fprintf(exnode, "%6.2f\n", avediameter[k]/2);
		}
	}
	fclose(dotcom);
	fclose(exelem);
	fclose(exnode);
	free(radius);
	free(point_used);
	return err;
}

//-----------------------------------------------------------------------------------------------------
// Create Amira SpatialGraph file
// This assumes that all edges are used, i.e. the network has been squeezed
//-----------------------------------------------------------------------------------------------------
int WriteAmiraFile(char *outFile, char *vessFile, char *skelFile)
{
	int i, k, j, npts, npts_used;

	printf("\nWriteAmiraFile: %s\n",outFile);
	fprintf(fpout,"\nWriteAmiraFile: %s\n",outFile);
	npts = 0;
	npts_used = 0;
	for (i=0;i<ne;i++) {
		if (!edgeList[i].used) continue;
		npts += edgeList[i].npts;
		npts_used += edgeList[i].npts_used;
	}

	FILE *fpam = fopen(outFile,"w");
	fprintf(fpam,"# AmiraMesh 3D ASCII 2.0\n");
	fprintf(fpam,"# Created by topo.exe from: %s and %s\n",vessFile,skelFile);
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
		fprintf(fpam,"%6.1f %6.1f %6.1f\n",vsize[0]*vertex[i].pos[0],vsize[1]*vertex[i].pos[1],vsize[2]*vertex[i].pos[2]);
	}
	fprintf(fpam,"\n@2\n");
	for (i=0;i<ne;i++) {
		fprintf(fpam,"%d %d\n",edgeList[i].vert[0],edgeList[i].vert[1]);
		if (edgeList[i].vert[0]== edgeList[i].vert[1]) {
			fprintf(fperr,"Error: WriteAmiraFile: repeated vertices: i: %d  %d\n",i,edgeList[i].vert[0]);
		}
	}
	fprintf(fpam,"\n@3\n");
	for (i=0;i<ne;i++) {
		fprintf(fpam,"%d\n",edgeList[i].npts);
//		fprintf(fpam,"%d\n",edgeList[i].npts_used);
	}
	fprintf(fpam,"\n@4\n");
	for (i=0;i<ne;i++) {
		for (k=0;k<edgeList[i].npts;k++) {
			j = edgeList[i].pt[k];
			fprintf(fpam,"%6.1f %6.1f %6.1f\n",vsize[0]*voxel[j].pos[0],vsize[1]*voxel[j].pos[1],vsize[2]*voxel[j].pos[2]);
		}
	}
	fprintf(fpam,"\n@5\n");
	double diam;
	for (i=0;i<ne;i++) {
		for (k=0;k<edgeList[i].npts;k++) {
			j = edgeList[i].pt[k];
			diam = avediameter[j];
			if (diam < 1.0) {
//				printf("d < 1.0: edge: %d npts: %d pt: %d %d  %f\n",i,edgeList[i].npts,k,j,diam);
				fprintf(fperr,"d < 1.0: edge: %d npts: %d pt: %d %d  %f\n",i,edgeList[i].npts,k,j,diam);
				diam = 1.0;
			}
			fprintf(fpam,"%6.2f\n",avediameter[j]);
		}
	}
	fclose(fpam);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double dist1(int k1, int k2)
{
	double p1[3], p2[3], d[3], d2;

	d2 = 0;
	for (int k=0; k<3; k++) {
		p1[k] = point[k1][k];
		p2[k] = point[k2][k];
		d[k] = p1[k] - p2[k];
		d2 += d[k]*d[k];
	}
	return sqrt(d2);
}

//-----------------------------------------------------------------------------------------------------
// Distance between voxels k1 and k2.  Previously the distance was returned in units of voxels,
// but to account for dx=dy!=dz the distance is now in um.
//-----------------------------------------------------------------------------------------------------
double dist_um(int k1, int k2)
{
	double del, d2;

	d2 = 0;
	for (int k=0; k<3; k++) {
		del = vsize[k]*(voxel[k1].pos[k] - voxel[k2].pos[k]);
		d2 += del*del;
	}
	return sqrt(d2);
}
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
double dotProduct(double *v1, double *v2)
{
	int i;
	double sum;

	sum = 0;
	for (i=0;i<3;i++) {
		sum += v1[i]*v2[i];
	}
	return sum;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];

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

/*
//-----------------------------------------------------------------------------------------------------
// The point of interest, p0, lies on the line connecting p1 and p2.
// We want to estimate the mean square radius of the vessel at p0.
// We do this by shooting lines equispaced through 360 degrees from p0,
// in the plane normal to p1-p2, to find distance to the first black voxel.
//-----------------------------------------------------------------------------------------------------
float MeanSquareRadius(int *p0, int *p1, int *p2)
{
	int i;
	double N[3], sum, dp, r=1;
	double v0[3], v[3], p0r[3];
	double angle;
	int nrays = 16;

	sum = 0;
	for (i=0; i<3; i++) {
		N[i] = p2[i] - p1[i];	// axis of vessel, normal to the plane
		sum += N[i]*N[i];
		p0r[i] = p0[i];
	}
	sum = sqrt(sum);
	for (i=0; i<3; i++) {
		N[i] = N[i]/sum;		// unit normal
	}
	// N is normal to the plane that approximates the slice through the vessel at p0
	// Need to determine a a base vector through p0 in this plane
	v[0] = 1; v[1] = 0; v[2] = 0;	// try this
	dp = dotProduct(v,N);
	dp = abs(dp);
	if (dp > 0.9) {					// too close to N, use this
		v[0] = 0; v[1] = 1; v[2] = 0;
	}
	// Define v0 as the vector normal to both v and N, i.e. lies in the plane
	crossProduct(v0,v,N);
	dp = dotProduct(v0,v0);
	dp = sqrt(dp);
	for (i=0; i<3; i++) {
		v0[i] = v0[i]/dp;		// unit vector in the plane
	}
//	printf("N:   %f %f %f\n",N[0],N[1],N[2]);
//	printf("v0:  %f %f %f\n",v0[0],v0[1],v0[2]);
	sum = 0;
	for (i=0;i<nrays;i++) {
		angle = i*(2.*PI/nrays);
		Rotate(v0,N,angle,v);
//		r = GetRadius(p0r,v);
		sum += r*r;
//		printf("%d  %f %f %f  %f\n",i,v[0],v[1],v[2],r);
	}
	return sum/nrays;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int VertexDiameters(void)
{
	int i, j, k, kp, kv;
	EDGE edge;

	// Need to handle vertices.  Use max of connected points
	for (i=0;i<nv;i++)
		nve[i] = 0;
	for (int ie=0;ie<ne;ie++) {
		edge = edgeList[ie];
		for (j=0;j<2;j++) {
			kv = edge.vert[j];
			if (kv == 0) printf("edge: %d zero vertex: %d\n",ie,j);
			bool inlist = false;
			if (nve[kv] > 0) {
				for (k=0;k<nve[kv];k++) {
					if (vertEdge[kv][k] == ie) {
						inlist = true;
						break;
					}
				}
			}
			if (!inlist) {
				vertEdge[kv][nve[kv]] = ie;
				nve[kv]++;
			}
		}
	}
	// These are the stubs that have been dropped
//	for (kv=0;kv<nv;kv++) {
//		if (nve[kv] == 0) {
//			printf("vertex with no edges: %d\n",kv);
//			exit(1);
//		}
//	}

	// We now know which edges are connected to each vertex
	k = 0;
	for (;;) {
		k++;
		int missing = 0;
		for (kv=0;kv<nv;kv++) {
			if (nve[kv] == 0) continue;
			double d, dmax = 0;
			for (i=0;i<nve[kv];i++) {
				bool dbug = (kv == -1);
				edge = edgeList[vertEdge[kv][i]];
				if (dbug) {
					printf("vertex kv: %d i: %d edge: %d  %d %d  %d %d\n",kv,i,
						vertEdge[kv][i],
						edge.vert[0],edge.vert[1],
						edge.pt[0],edge.pt[edge.npts-1]);
				}
				if (edge.vert[0] == kv) {
					kp = edge.pt[0];
					d = avediameter[edge.pt[1]];
					if (dbug) printf("vert 0: kp: %d d: %f\n",kp,d);
				} else {
					kp = edge.pt[edge.npts-1];
					d = avediameter[edge.pt[edge.npts-2]];
					if (dbug) printf("vert 1: kp: %d d: %f\n",kp,d);
				}
				if (d > dmax) {
					dmax = d;
					if (dbug) printf("kp: %d dmax: %f\n",kp,dmax);
				}
			}
			if (dmax > 0) {
				avediameter[kp] = dmax;
			} else {
				missing++;
			}
		}
		if (missing == 0) break;
		printf("Vertex diameter iteration: %d  missing: %d\n",k,missing);
//		exit(1);
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Get distance from p1[] in direction of unit vector v[] to the first black voxel.
// This is not necessarily the best method.  It leads to odd oscillation in the min diameter dist. plot.
// Add delta = 0.5 because a voxel at (i,j,k) spans i-0.5 <= x < i+0.5 etc.
//-----------------------------------------------------------------------------------------------------
double GetRadius(double p1[3], double v[3])
{
	int i, x, y, z;
	double r;
	double delta = 0.5;

	i = 0;
	for(;;) {
		r = i*0.1;
		x = int(p1[0] + r*v[0] + delta);
		y = int(p1[1] + r*v[1] + delta);
		z = int(p1[2] + r*v[2] + delta);
		if (x <= 0 || x >= width-1 || y <= 0 || y >= height-1 || z <= 0 || z >= depth-1) {
			r = 0.5;
			break;
		}
		if (V3D(x,y,z) == 0) {
			if (r == 0) {
				printf("GetRadius (a): r=0: i,x,y,z,V: %d  %d %d %d  %d\n",i,x,y,z,V3D(x,y,z));
				fprintf(fperr,"GetRadius (a): r=0: i,x,y,z,V: %d  %d %d %d  %d\n",i,x,y,z,V3D(x,y,z));
			}
			break;
		}
		i++;
	}
	if (r == 0) {
		printf("GetRadius (b): r=0: i,x,y,z: %d  %d %d %d\n",i,x,y,z);
		fprintf(fperr,"GetRadius (b): r=0: i,x,y,z: %d  %d %d %d\n",i,x,y,z);
	}
	return r;
}
*/

//-----------------------------------------------------------------------------------------------------
// Get distance from p1[] in direction of unit vector v[] to the first black voxel.
//-----------------------------------------------------------------------------------------------------
double GetRadius2(double p1[3], double v[3])
{
	int i, k, ixyz[3], ixyzmax[3];
	double r, dr, xyz;
	bool out;

//	printf("GetRadius2\n");
	dr = vsize[2]/10;
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
			if (ixyz[k] < 0) out = true;;
			if (ixyz[k] > ixyzmax[k]) out = true;
		}
//		printf("ixyz: %d %d %d\n",ixyz[0],ixyz[1],ixyz[2]);
		if (out) break;
		if (V3D(ixyz[0],ixyz[1],ixyz[2]) == 0) {
			if (r == 0) {
				printf("GetRadius2 (a): r2=0: i,x,y,z,V: %d  %d %d %d  %d\n",i,ixyz[0],ixyz[1],ixyz[2],V3D(ixyz[0],ixyz[1],ixyz[2]));
				fprintf(fperr,"GetRadius2 (a): r2=0: i,x,y,z,V: %d  %d %d %d  %d\n",i,ixyz[0],ixyz[1],ixyz[2],V3D(ixyz[0],ixyz[1],ixyz[2]));
				return 1.0;
			}
			break;
		}
		i++;
	}
	return r*r;
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
int EstimateDiameter(double p0[3], double p1[3], double p2[3], double *r2ave, double *r2min)
{
	int i;
	double N[3], sum, dp, r2, r2sum;
	double v0[3], v[3];
	double angle;
	int nrays = 16;

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
	dp = dotProduct(v,N);
	dp = abs(dp);
	if (dp > 0.9) {					// too close to N, use this
		v[0] = 0; v[1] = 1; v[2] = 0;
	}
	// Define v0 as the vector normal to both v and N, i.e. lies in the plane
	crossProduct(v0,v,N);
	dp = dotProduct(v0,v0);
	dp = sqrt(dp);
	for (i=0; i<3; i++) {
		v0[i] = v0[i]/dp;		// unit vector in the plane
	}
	r2sum = 0;
	*r2min = 1.0e10;
	for (i=0;i<nrays;i++) {
		angle = i*(PI/nrays);
//		printf("angle: %2d %6.2f %6.3f %6.3f %6.3f\n",i,angle,v0[0],v0[1],v0[2]);
		Rotate(v0,N,angle,v);
//		printf("angle: %2d %6.2f %6.3f %6.3f %6.3f\n",i,angle,v[0],v[1],v[2]);
		r2 = GetRadius2(p1,v);
		if (r2 <= 1.0) {
//			printf("EstimateDiameter: r2<=1.0: p1: %f %f %f\n",p1[0],p1[1],p1[2]);
			fprintf(fperr,"EstimateDiameter: r2<=1.0: p1: %f %f %f\n",p1[0],p1[1],p1[2]);
			break;;
		}
		r2sum += r2;
		if (r2 < *r2min) *r2min = r2;
		angle += PI;
		Rotate(v0,N,angle,v);
		r2 = GetRadius2(p1,v);
		if (r2 <= 1.0) {
//			printf("EstimateDiameter: r2<=1.0: p1: %f %f %f\n",p1[0],p1[1],p1[2]);
			fprintf(fperr,"EstimateDiameter: r2<=1.0: p1: %f %f %f\n",p1[0],p1[1],p1[2]);
			break;;
		}
		r2sum += r2;
		if (r2 < *r2min) *r2min = r2;
	}
	*r2ave = r2sum/(2*nrays);
	*r2ave = MAX(1.0,*r2ave);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Last ditch estimate of diameter at a point.
// Try a number of planes (9).  The average square radius is computed for each cutting plane,
// and the minimum is chosen to estimate the average diameter.
// Modified to account for dx=dy!=dz
// The square incremental distance dr2 now accounts for voxelsize_xy and voxelsize_z,
// and the value returned by the function is now the diameter in um (needs no scaling later).
//-----------------------------------------------------------------------------------------------------
double pointdiameter(int kp)
{
	int xyz[3];
	int i, iplane, ivec, k, dx, dy, dz, x, y, z;
	double r2min, r2sum, dr2, diam;
	float delta = 1.3;	// increment to subtract from distance to zero voxel

	for (i=0; i<3; i++)
		xyz[i] = voxel[kp].pos[i];
//	printf("pointdiameter: kp: %d  position: %d %d %d\n",kp,xyz[0],xyz[1],xyz[2]);
	if (kp == 447832)
		fprintf(fperr,"pointdiameter: kp: %d  position: %d %d %d\n",kp,xyz[0],xyz[1],xyz[2]);
	fflush(fperr);
	r2min = 1.0e10;
	for (iplane=0; iplane<9; iplane++) {
		r2sum = 0;
		for (ivec=0; ivec<8; ivec++) {
			dx = plane[iplane].vector[ivec].dx;
			dy = plane[iplane].vector[ivec].dy;
			dz = plane[iplane].vector[ivec].dz;

			dr2 = dx*dx*vsize[0]*vsize[0] + dy*dy*vsize[1]*vsize[1] + dz*dz*vsize[2]*vsize[2];

//			printf("plane: %d vector: %d  %d %d %d  %f\n",iplane,ivec,dx,dy,dz,dr2);
//			fprintf(fperr,"plane: %d vector: %d  %d %d %d  %f\n",iplane,ivec,dx,dy,dz,dr2);
//			fflush(fperr);
			k = 0;
			for (;;) {
				k++;
				x = xyz[0] + k*dx;
				y = xyz[1] + k*dy;
				z = xyz[2] + k*dz;
				if (x<0 || x>width-1 || y<0 || y>height-1 || z<0 || z>depth-1) {
					r2sum += k*k*dr2;
					break;
				}
				if (V3D(x,y,z) == 0) {
					r2sum += (k-delta)*(k-delta)*dr2;
					break;
				}
			}
		}
		if (r2sum < r2min)
			r2min = r2sum;
	}
	if (r2min <= 0) {
		return 1.0;
	}
	diam = 2*sqrt(r2min/8);
	return diam;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int CrudeGetDiameters(void)
{
	int k, kmax;
	float dmax;

	printf("GetDiameters: np: %d\n",np);
	fprintf(fperr,"GetDiameters: np: %d\n",np);
	dmax = 0;
	for (k=0; k<np; k++) {
		avediameter[k] = pointdiameter(k);
		if (avediameter[k] > dmax) {
			dmax = avediameter[k];
			kmax = k;
		}
	}
	fprintf(fperr,"Max. diameter at: kp = %d  x,y,z: %d %d %d d: %6.2f\n",kmax,
		voxel[kmax].pos[0],voxel[kmax].pos[2],voxel[kmax].pos[2],dmax);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// This method uses the preceding and following pts on an edge to estimate the centreline direction
// vector.  To improve the estimate of the cutting plane perpendicular to the vessel, the pts need 
// to be well separated, i.e. this must be applied to the simplified network.
// It can be used only for the interior pts, therefore not when npts = 2.
//-----------------------------------------------------------------------------------------------------
double getDiameter(int kp0, int kp1, int kp2)
{
	int i;
	double p1[3], p2[3], p0[3];
	double r2_ave, r2_min, diam, factor;
//	double alpha = 0.0;
	double dlim = 50.0;

	for (i=0; i<3; i++) {
		p0[i] = vsize[i]*(voxel[kp0].pos[i] + 0.5);		// Note: centres of voxel cubes
		p1[i] = vsize[i]*(voxel[kp1].pos[i] + 0.5);
		p2[i] = vsize[i]*(voxel[kp2].pos[i] + 0.5);
	}
	// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
	EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min);
	diam = 2*sqrt(r2_ave);
	if (calib_param != 0) {
		//factor = 1.0 + calib_param*diam/dlim;
		//diam = factor*diam;
		diam *= calib_param;
	}
	return diam;
}

//-----------------------------------------------------------------------------------------------------
// Estimates diameters for all interior points.
// Note that an edge with npts=2 (no intermediate points) is not processed.
//-----------------------------------------------------------------------------------------------------
int GetDiameters(void)
{
	int ie, npts, k, kp, kp0, kp1, kp2;
	float ad;
	EDGE *edge;

	printf("GetDiameters\n");
	fprintf(fpout,"GetDiameters\n");
	for (k=0; k<np; k++) {
		avediameter[k] = 0;
	}
	for (ie=0;ie<ne;ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
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
	return 0;
}

/*
//-----------------------------------------------------------------------------------------------------
// This version estimates the average diameter, and minimum diameter, at edge midpoints.
// Note that the diameters are edge properties.
// Get edge lengths too. There are two methods:
// (a) add voxel-voxel distances (this will overestimate length)
// (b) compute distance from start to end of the edge (segment).  (this will underestimate length)
//-----------------------------------------------------------------------------------------------------
int OldGetDiameters(void)
{
	int i, j, k, k0, k1, k2;
	double p1[3], p2[3], p0[3];
	double diam_ave, diam_min, d2, dsum;
	EDGE edge;
	int nout = 0;

	printf("GetDiameters\n");
	for (i=0;i<ne;i++) {
		edge = edgeList[i];
		for (j=0;j<edge.npts;j++) {
			k = edge.pt[j];
			avediameter[k] = 0;
		}
		if (!edge.ok) continue;
		k0 = edge.pt[0];
		k2 = edge.pt[edge.npts-1];
		edgeList[i].length_um = dist_um(k0,k2);
		dsum = 0;
		for (j=0;j<edge.npts-1;j++) {
			k0 = edge.pt[j];
			k2 = edge.pt[j+1];
			dsum += dist_um(k0,k2);
			// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
//			EstimateDiameter(p0,p1,p2,&diam_ave,&diam_min);
//			edgeList[i].avediameter[j] = diam_ave;
//			edgeList[i].mindiameter[j] = diam_min;
		}
		edgeList[i].length_um = dsum;
		edge = edgeList[i];
//		printf("edge: %d  %d\n",i,edge.npts);
//		fprintf(fperr,"edge: %d  %d\n",i,edge.npts);
		// Try to improve the estimation of diameters
		// (a) Use skeleton points spaced 5 voxels apart (if possible) to provide estimate of centreline p0,p2
		// (b) Compute average values of ave and min diameter for the segment
		if (!edge.ok) continue;
		if (edge.npts <= 3) continue;	// do not compute anything for very short segments
		for (j=1;j<edge.npts-1;j++) {
			int j0 = MAX(j-2,0);
			int j2 = MIN(j+2,edge.npts-1);
			k0 = edge.pt[j0];
			k1 = edge.pt[j];
			k2 = edge.pt[j2];
			d2 = 0;
			for (k=0; k<3; k++) {
				p0[k] = point[k0][k];
				p1[k] = point[k1][k];
				p2[k] = point[k2][k];
//				p1[k] = (p0[k] + p2[k])/2.;
			}
//			printf("i,j0,j,j2,k0,k1,k2,p1: %d %d %d %d %d %d %d  %d %d %d\n",i,j0,j,j2,k0,k1,k2, point[k1][0], point[k1][1], point[k1][2]);
			if (V3D(point[k1][0], point[k1][1], point[k1][2]) == 0) {
				printf("Skeleton point outside vessel: %d %d  %d %d %d\n",i,j, point[k1][0], point[k1][1], point[k1][2]);
				fprintf(fperr,"Skeleton point outside vessel: %d %d  %d %d %d\n",i,j, point[k1][0], point[k1][1], point[k1][2]);
				return 1;
			}
			EstimateDiameter(p0,p1,p2,&diam_ave,&diam_min);
			edgeList[i].avediameter[j] = diam_ave;
//m			edgeList[i].mindiameter[j] = diam_min;
			mindiameter[j] = diam_min;
		}
		edgeList[i].avediameter[0] = edgeList[i].avediameter[1];
		edgeList[i].avediameter[edge.npts-1] = edgeList[i].avediameter[edge.npts-2];
//m		edgeList[i].mindiameter[0] = edgeList[i].mindiameter[1];
//m		edgeList[i].mindiameter[edge.npts-1] = edgeList[i].mindiameter[edge.npts-2];
		mindiameter[0] = mindiameter[1];
		mindiameter[edge.npts-1] = mindiameter[edge.npts-2];
		// Now compute segment averages
		double asum = 0;
		double msum = 0;
		for (j=0;j<edge.npts;j++) {
			k = edge.pt[j];
			avediameter[k] = edgeList[i].avediameter[j];
			asum += edgeList[i].avediameter[j];
//m			msum += edgeList[i].mindiameter[j];
			msum += mindiameter[j];
		}
		edgeList[i].segavediam = asum/edge.npts;
		edgeList[i].segmindiam = msum/edge.npts;
		if (edgeList[i].segmindiam < 0.5) {
			printf("segmindiam < 0.5: %d %d %f\n",i,edge.npts,edgeList[i].segmindiam);
			fprintf(fperr,"segmindiam < 1: %d %d %f %f\n",i,edge.npts,edgeList[i].segmindiam,edgeList[i].segavediam);
		}
		if (edgeList[i].segavediam < 0.5) {
			printf("segavediam < 0.5: %d %d %f\n",i,edge.npts,edgeList[i].segavediam);
			fprintf(fperr,"segavediam < 0.5: %d %d  min, ave: %f %f\n",i,edge.npts,edgeList[i].segmindiam,edgeList[i].segavediam);
		}
	}
	// Now need to treat edges with <=3 points
	for (int it=0; it<10; it++) {
		int nmiss = 0;
		float d[3];
		for (i=0;i<ne;i++) {
			edge = edgeList[i];
			if (!edge.ok) continue;
			if (edge.npts > 3) continue;
			float dsum = 0;
			int nsum = 0;
			for (j=0;j<edge.npts;j++) {
				k = edge.pt[j];
				d[j] = avediameter[k];
				if (d[j] > 0) {
					nsum++;
					dsum += d[j];
				}
			}
			if (nsum == 0) {
				nmiss++;
				for (j=0;j<edge.npts;j++) {
					k = edge.pt[j];
					d[j] = pointdiameter(k);
					nsum++;
					dsum += d[j];
				}
				for (j=0;j<edge.npts;j++) {
					k = edge.pt[j];
					avediameter[k] = dsum/nsum;
				}
			}  else {
				for (j=0;j<edge.npts;j++) {
					k = edge.pt[j];
					if (d[j] == 0)
						avediameter[k] = dsum/nsum;
				}
			}
		}
		fprintf(fperr,"No diameter data for: %d\n",nmiss);
		fflush(fperr);
		if (nmiss == 0) break;
	}
	printf("nout: %d\n",nout);
	return 0;
}
*/

//-----------------------------------------------------------------------------------------------------
// Initialize vectors for the 9 planes.
// For each plane there is a set of 8 vectors, corresponding to the 8 neighbour voxels.
// Since in general dx=dy but dx != dz, distance calculations need to take account of this asymmetry.
// How to generate more than 9 planes?
//-----------------------------------------------------------------------------------------------------
int InitVector(void)
{
	int iplane, k, dx, dy, dz;

	iplane = 0;		// XY plane
	dz = 0;
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dy=-1; dy<=1; dy++) {
			if (dx==0 && dy==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = dz;
			k++;
		}
	}
	iplane = 1;		// YZ plane
	dx = 0;
	k = 0;
	for (dz=-1; dz<=1; dz++) {
		for (dy=-1; dy<=1; dy++) {
			if (dz==0 && dy==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = dz;
			k++;
		}
	}
	iplane = 2;		// ZX plane
	dy = 0;
	k = 0;
	for (dz=-1; dz<=1; dz++) {
		for (dx=-1; dx<=1; dx++) {
			if (dz==0 && dx==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = dz;
			k++;
		}
	}
	iplane = 3;		// parallel to x axis (dz = dy)
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dy=-1; dy<=1; dy++) {
			if (dx==0 && dy==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = dy;
			k++;
		}
	}
	iplane = 4;		// parallel to x axis (dz = -dy)
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dy=-1; dy<=1; dy++) {
			if (dx==0 && dy==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = -dy;
			k++;
		}
	}
	iplane = 5;		// parallel to y axis (dz = dx)
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dy=-1; dy<=1; dy++) {
			if (dx==0 && dy==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = dx;
			k++;
		}
	}
	iplane = 6;		// parallel to y axis (dz = -dx)
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dy=-1; dy<=1; dy++) {
			if (dx==0 && dy==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dy;
			plane[iplane].vector[k].dz = -dx;
			k++;
		}
	}
	iplane = 7;		// parallel to z axis (dy = dx)
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dz=-1; dz<=1; dz++) {
			if (dx==0 && dz==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = dx;
			plane[iplane].vector[k].dz = dz;
			k++;
		}
	}
	iplane = 8;		// parallel to z axis (dy = -dx)
	k = 0;
	for (dx=-1; dx<=1; dx++) {
		for (dz=-1; dz<=1; dz++) {
			if (dx==0 && dz==0) continue;
			plane[iplane].vector[k].dx = dx;
			plane[iplane].vector[k].dy = -dx;
			plane[iplane].vector[k].dz = dz;
			k++;
		}
	}
	//for (iplane=0; iplane<9; iplane++) {
	//	fprintf(fperr,"plane: %d\n",iplane);
	//	for (k=0; k<8; k++)
	//		fprintf(fperr,"%d  %d %d %d\n",k,plane[iplane].vector[k].dx,plane[iplane].vector[k].dy,plane[iplane].vector[k].dz);
	//}
	return 0;
}


/*
//-----------------------------------------------------------------------------------------------------
// This version estimates the average diameter at node points. 
//-----------------------------------------------------------------------------------------------------
int GetDiameters1(void)
{
	int i, j, k, k0, k1, k2, kv, kp;
	int p1[3], p2[3], p0[3];
	double d, r2;
	EDGE edge;

	printf("GetDiameters\n");
	for (i=0;i<ne;i++) {
		edge = edgeList[i];
		if (edge.npts == 2) continue;
		for (j=1;j<edge.npts-1;j++) {
//		for (j=0;j<edge.npts;j++) {
			k0 = edge.pt[j];
			k1 = edge.pt[j-1];
			k2 = edge.pt[j+1];
			for (k=0; k<3; k++) {
				p0[k] = point[k0][k];
				p1[k] = point[k1][k];
				p2[k] = point[k2][k];
			}
			if (V3D(p0[0],p0[1],p0[2]) == 0) {
				printf("Error: skeleton point is outside: %d %d %d  %d %d %d\n",i,j,k0,p0[0],p0[1],p0[2]);
				exit(1);
			}
			r2 = MeanSquareRadius(p0,p1,p2);
			avediameter[k0] = 2*sqrt(r2);
		}
	}
	// Need to handle vertices.  Use max of connected points
	for (i=0;i<nv;i++)
		nve[i] = 0;
	for (int ie=0;ie<ne;ie++) {
		edge = edgeList[ie];
		for (j=0;j<2;j++) {
			kv = edge.vert[j];
			if (kv == 0) printf("edge: %d zero vertex: %d\n",ie,j);
			bool inlist = false;
			if (nve[kv] > 0) {
				for (k=0;k<nve[kv];k++) {
					if (vertEdge[kv][k] == ie) {
						inlist = true;
						break;
					}
				}
			}
			if (!inlist) {
				vertEdge[kv][nve[kv]] = ie;
				nve[kv]++;
			}
		}
	}
	// These are the stubs that have been dropped
//	for (kv=0;kv<nv;kv++) {
//		if (nve[kv] == 0) {
//			printf("vertex with no edges: %d\n",kv);
//			exit(1);
//		}
//	}

	// We now know which edges are connected to each vertex
	k = 0;
	for (;;) {
		k++;
		int missing = 0;
		for (kv=0;kv<nv;kv++) {
			if (nve[kv] == 0) continue;
			double dmax = 0;
			for (i=0;i<nve[kv];i++) {
				bool dbug = (kv == -1);
				edge = edgeList[vertEdge[kv][i]];
				if (dbug) {
					printf("vertex kv: %d i: %d edge: %d  %d %d  %d %d\n",kv,i,
						vertEdge[kv][i],
						edge.vert[0],edge.vert[1],
						edge.pt[0],edge.pt[edge.npts-1]);
				}
				if (edge.vert[0] == kv) {
					kp = edge.pt[0];
					d = avediameter[edge.pt[1]];
					if (dbug) printf("vert 0: kp: %d d: %f\n",kp,d);
				} else {
					kp = edge.pt[edge.npts-1];
					d = avediameter[edge.pt[edge.npts-2]];
					if (dbug) printf("vert 1: kp: %d d: %f\n",kp,d);
				}
				if (d > dmax) {
					dmax = d;
					if (dbug) printf("kp: %d dmax: %f\n",kp,dmax);
				}
			}
			if (dmax > 0) {
				avediameter[kp] = dmax;
			} else {
				missing++;
			}
		}
		if (missing == 0) break;
		printf("Vertex diameter iteration: %d  missing: %d\n",k,missing);
//		exit(1);
	}

	return 0;
}
*/

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int GetNumberOfNeighbors(int x0, int y0, int z0)
{
	int dx, dy, dz, x, y, z;

	int n = -1;	// to exclude (x0,y0,z0)
	for (dx=-1; dx<=1; dx++) {
		x = x0 + dx;
		if (x < 0 || x > width-1) continue;
		for (dy=-1; dy<=1; dy++) {
			y = y0 + dy;
			if (y < 0 || y > height-1) continue;
			for (dz=-1; dz<=1; dz++) {
				z = z0 + dz;
				if (z < 0 || z > depth-1) continue;
				if (V3Dskel(x,y,z) > 0) n++;
			}
		}
	}
	return n;
}

// Look for another member of tmplist[] that is a neighbour of k
// If such a member has smaller d2, k is dominated
bool dominated1(int k, int ntmp, int tmplist[][3], bool drop[], int d2[]) {
	int j;

	for (j=0; j<ntmp; j++) {
		if (j == k) continue;
		if (abs(tmplist[k][0]-tmplist[j][0]) <= 1 &&
		    abs(tmplist[k][1]-tmplist[j][1]) <= 1 &&
		    abs(tmplist[k][2]-tmplist[j][2]) <= 1) {
				{
					if (d2[j] < d2[k]) return true;
				}
		}
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------
// A neighbour V1 of V0 is dominated if:
// (1) V1 is not N6 of V0
// (2) V1 is N6 of V2, which is N26 of V0
//-----------------------------------------------------------------------------------------------------
bool dominated(int k, int pos[], int ntmp, int tmplist[][3], bool drop[]) {
	int i, j, neq, nd, offset;

	// Is V1 = tmplist[k][] N6 of pos?
	nd = 0;
	neq = 0;
	for (i=0; i<3; i++) {
		offset = tmplist[k][i] - pos[i];
		if (offset == 1 || offset == -1) nd++;
		if (offset == 0) neq++;
	}
	if (neq == 2 && nd == 1) return false;	// V0 is N6 of pos, therefore not dominated

	// Is V1 N6 of another neighbour V2 tmplist[j][]?
	for (j=0; j<ntmp; j++) {
		if (j==k) continue;
		if (drop[j]) continue;
		nd = 0;
		neq = 0;
		for (i=0; i<3; i++) {
			offset = tmplist[k][i] - tmplist[j][i];
			if (offset == 1 || offset == -1) nd++;
			if (offset == 0) neq++;
		}
		if (neq == 2 && nd == 1) 
			return true;		// V0 is N6 of V2, therefore dominated
	}
	return false;
}

//-----------------------------------------------------------------------------------------------------
// Old code
//-----------------------------------------------------------------------------------------------------
int GetNeighborsOld(int *pos, int nbrlist[][3])
{
	int dx, dy, dz, x, y, z, nbrs;

	nbrs = 0;
	for (dz=-1; dz<=1; dz++) {
		z = pos[2] + dz;
		if (z < 0 || z > depth-1)
			continue;
		for (dy=-1; dy<=1; dy++) {
			y = pos[1] + dy;
			if (y < 0 || y > height-1)
				continue;
			for (dx=-1; dx<=1; dx++) {
				if (dx == 0 && dy == 0 && dz == 0) continue;
				x = pos[0] + dx;
				if (x < 0 || x > width-1)
					continue;
				if (V3Dskel(x,y,z) != 0) {
					nbrlist[nbrs][0] = x;
					nbrlist[nbrs][1] = y;
					nbrlist[nbrs][2] = z;
					nbrs++;
				}
			}
		}
	}
	return nbrs;
}

//-----------------------------------------------------------------------------------------------------
// New version
//-----------------------------------------------------------------------------------------------------
int GetNeighborsNew(int *pos, int nbrlist[][3])
{
	int dx, dy, dz, x, y, z, nbrs, i, k, ntmp;
	int tmplist[20][3], d2[20];
	bool drop[20];
	bool dbug=false;;

//	if (dbug) printf("pos: %3d %3d %3d\n",pos[0],pos[1],pos[2]);
	nbrs = 0;
	for (dz=-1; dz<=1; dz++) {
		z = pos[2] + dz;
		if (z < 0 || z > depth-1)
			continue;
		for (dy=-1; dy<=1; dy++) {
			y = pos[1] + dy;
			if (y < 0 || y > height-1)
				continue;
			for (dx=-1; dx<=1; dx++) {
				x = pos[0] + dx;
				if (dx == 0 && dy == 0 && dz == 0) continue;
				if (x < 0 || x > width-1)
					continue;
				if (V3Dskel(x,y,z) != 0) {
					nbrlist[nbrs][0] = x;
					nbrlist[nbrs][1] = y;
					nbrlist[nbrs][2] = z;
//					if (dbug) printf("nbrs: %d  %3d %3d %3d\n",nbrs,x,y,z);
					nbrs++;
				}
			}
		}
	}
	if (nbrs <= 2) return nbrs;
	// We need to check that there are no redundant neighbour links
	if (dbug) printf("pos: %3d %3d %3d  nbrs: %d\n",pos[0],pos[1],pos[2],nbrs);
	for (k=0; k<nbrs; k++) {
		for (i=0; i<3; i++) {
			tmplist[k][i] = nbrlist[k][i];
			drop[k] = false;
			//dx = tmplist[k][0] - pos[0];
			//dy = tmplist[k][1] - pos[1];
			//dz = tmplist[k][2] - pos[2];
			//d2[k] = dx*dx + dy*dy + dz*dz;
		}
		if (dbug) printf("%d %d %d %d\n",k,tmplist[k][0],tmplist[k][1],tmplist[k][2]);
	}
	ntmp = nbrs;
	for (k=0; k<ntmp; k++) {
//		if (dominated1(k,ntmp,tmplist,drop,d2)) {
		if (dominated(k,pos,ntmp,tmplist,drop)) {
			drop[k] = true;
			if (dbug) printf("drop: %d\n",k);
		}
	}
	nbrs = 0;
	for (k=0; k<ntmp; k++) {
		if (!drop[k]) {
			for (i=0; i<3; i++) {
				nbrlist[nbrs][i] = tmplist[k][i];
			}
//			printf("nbrs: %d  %3d %3d %3d\n",nbrs,nbrlist[nbrs][0],nbrlist[nbrs][1],nbrlist[nbrs][2]);
			nbrs++;
		}
	}
	return nbrs;
}

/*
//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void EndVertex(int *pos)
{
	VERTEX *vp = &vertex[nv];
	for (int i=0;i<3;i++)
		vp->pos[i] = pos[i];
	nv++;
	vp->nlinks = 1;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void AddVertex(int *pos, int *prev, int nbrlist[][3], int nbrs, int *j)
{
	VERTEX *vp = &vertex[nv];
	int i;

	for (i=0;i<3;i++)
		vp->pos[i] = pos[i];
	nv++;
	if (nv == MAXVERTICES) {
		printf("ERROR: nv == MAXVERTICES\n");
		fprintf(fperr,"ERROR: nv == MAXVERTICES\n");
		fclose(fperr);
		exit(1);
	}
	vp->nlinks = nbrs;
	vp->nfollowed = 1;
	*j = -1;
	for (i=0;i<nbrs;i++) {
		if (nbrlist[i][0] != prev[0] || nbrlist[i][1] != prev[1] || nbrlist[i][2] != prev[2]) {
			vp->followed[i] = false;
			if (*j < 0) {
				*j = i;
			}
		} else {
			vp->followed[i] = true;
		}
	}
	if (nbrs > maxnbrs)
		maxnbrs = nbrs;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void RevisitVertex(int kv, int *prev,  int nbrlist[][3], int nbrs, int *j)
{
	int i;

	VERTEX *vp = &vertex[kv];
	*j = -1;
	for (i=0;i<nbrs;i++) {
		if (nbrlist[i][0] == prev[0] && nbrlist[i][1] == prev[1] && nbrlist[i][2] == prev[2]) {
			vp->followed[i] = true;
		} else if (*j == -1 && !vp->followed[i]) {
			*j = i;
		}
	}
	vp->nfollowed++;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void AddPoint(int *pos)
{
	for(int i=0;i<3;i++)
		point[np][i] = pos[i];
	np++;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int PointIndex(int *pos)
{
	int kp;

	for (kp=0;kp<np;kp++) {
		if (point[kp][0] == pos[0] && point[kp][1] == pos[1] && point[kp][2] == pos[2])
			return kp;
	}
	return -1;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int TraceSkeleton1(void)
{
	int x, y, z, i, kp, kv, nbrs, nbrlist[10][3];
	EDGE edge;
//	int  xprev, yprev, zprev, xnext, ynext, znext;
	int prev[3], next[3];
	bool start, end;
	
	dbug = false;

	printf("TraceSkeleton\n");

	maxnbrs = 0;
	nv = 0;
	np = 0;
	nb = 0;
	ne = 0;
	negmin = 0;
	node = 0;
	// Find a starting point - an end of a line
	for (z=0; z<depth; z++) {
		printf("z: %d\n",z);
		for (y=0; y<height; y++) {
			for (x=0; x<width; x++) {
				if (V3Dskel(x,y,z) == 0) continue;
				next[0] = x; next[1] = y; next[2] = z;
				nbrs = GetNeighborsOld(next,nbrlist);
				if (nbrs != 1) continue;	// an end point has just one neighbor
				EndVertex(next);
				AddPoint(next);
				// Start a new edge
				edge.npts = 0;
				edge.pt[0] = np-1;
				edge.npts++;
				edge.vert[0] = nv-1;

				node++;
				break;
			}
			if (nv != 0) break;
		}
		if (nv != 0) break;
	}
	start = true;
	end = false;
	printf("Start point: %d %d %d\n",next[0],next[1],next[2]);
	int it = 0;
	for (;;) {
		nbrs = GetNeighborsOld(next,nbrlist);
		it++;
		if (dbug) fprintf(fperr,"it: %d  next: %d %d %d  nbrs: %d\n",it,next[0],next[1],next[2],nbrs);
		if (nbrs == 1) {
			if (start) {
				start = false;
				for (i=0;i<3;i++) {
					prev[i] = next[i];
					next[i] = nbrlist[0][i];
				}
			} else {
				AddPoint(next);
				EndVertex(next);
				// Complete an edge
				edge.pt[edge.npts] = np-1;
				edge.npts++;
				edge.vert[1] = nv-1;
				if (edge.vert[1] == edge.vert[0]) {
					fprintf(fperr,"Error: TraceSkeleton (a): repeated vertex: edge: %d\n",ne);
					return 1;
				}
//				if (edge.npts > MINSTUB) {
					edgeList[ne] = edge;
					if (dbug) fprintf(fperr,"Completed edge (a): %d npts: %d vertices: %d %d\n",ne,edge.npts,edge.vert[0],edge.vert[1]);
					ne++;
					negmin++;
					if (ne == MAXEDGES) {
						printf("ERROR: ne == MAXEDGES\n");
						fprintf(fperr,"ERROR: ne == MAXEDGES\n");
						fclose(fperr);
						return 2;
					}
//				}

				node++;
				if (dbug) fprintf(fperr,"            Reached the end of a vessel, need to start again at branch: %d\n",nb);
				end = true;
			}
		} else if (nbrs == 2) {		// normal non-vertex point
			AddPoint(next);
			if (nbrlist[0][0] != prev[0] || nbrlist[0][1] != prev[1] || nbrlist[0][2] != prev[2]) {
				for (i=0;i<3;i++) {
					prev[i] = next[i];
					next[i] = nbrlist[0][i];
				}
			} else {
				for (i=0;i<3;i++) {
					prev[i] = next[i];
					next[i] = nbrlist[1][i];
				}
			}
			// Add an edge point
			edge.pt[edge.npts] = np-1;
			edge.npts++;
			if (edge.npts == MAXEDGEPTS) {
				printf("ERROR: MAXEDGEPTS exceeded\n");
				return 3;
			}

			node++;
		} else if (nbrs >= 3) {		// vertex
			fprintf(fperr,"Vertex nbrlist:\n");
			if (dbug) {
				for (int k=0;k<nbrs;k++) {
					fprintf(fperr,"%d  %d %d %d\n",k,nbrlist[k][0],nbrlist[k][1],nbrlist[k][2]);
				}
			}
			// Need to check if it has already been visited by scanning the list
			kv = -1;
			for (int k=0;k<nv;k++) {
				VERTEX *vp = &vertex[k];
				if (vp->pos[0] == next[0] && vp->pos[1] == next[1] && vp->pos[2] == next[2]) {
					// The vertex is in the list
					kv = k;
					break;
				}
			}
			if (kv >= 0) {
				// This vertex has been visited before
				// Need to set link as followed
				// May need to go back to the last vertex (branchlist[nb]) if no link out is unfollowed
				// Find point in the list as kp
				kp = PointIndex(next);
				if (dbug) fprintf(fperr,"      Revisiting vertex #: %d  %d %d %d kp: %d\n",kv,next[0],next[1],next[2],kp);
				int j;	// this is the link to take away from the vertex
				RevisitVertex(kv,prev,nbrlist,nbrs,&j);
				node++;
				if (j < 0) {
					if (dbug) fprintf(fperr,"                  All links have been followed, go to branch: %d\n",nb);
					end = true;
				} else {
					node = 0;
					for (i=0;i<3;i++) {
						prev[i] = next[i];
						next[i] = nbrlist[j][i];
					}
					vertex[kv].followed[j] = true;
					vertex[kv].nfollowed++;
					if (vertex[kv].nfollowed < vertex[kv].nlinks) {
						branchlist[nb] = kv;
						if (dbug) fprintf(fperr,"Add visited vertex: %d to branchlist as nb: %d\n",kv,nb);
						nb++;
					}
				}
				// Find point in the list as kp
//				kp = PointIndex(next);
			} else {
				// This is a new vertex
				AddPoint(next);
				kp = np-1;
				node++;
				int j;	// this is the link to take away from the vertex
				AddVertex(next,prev,nbrlist,nbrs,&j);
				kv = nv-1;
				if (dbug) {
					fprintf(fperr,"New point kp: %d  %d %d %d\n",kp,point[kp][0],point[kp][1],point[kp][2]);
					fprintf(fperr,"New vertex kv: %d kp: %d next: %d %d %d\n",kv,kp,next[0],next[1],next[2]);
				}
				node = 0;
				for (i=0;i<3;i++) {
					prev[i] = next[i];
					next[i] = nbrlist[j][i];
				}
				if (dbug) fprintf(fperr,"From j: %d set next: %d %d %d\n",j,next[0],next[1],next[2]);
				vertex[kv].followed[j] = true;
				vertex[kv].nfollowed++;
				branchlist[nb] = kv;
				if (dbug) fprintf(fperr,"Add new vertex: %d to branchlist as nb: %d\n",kv,nb);
				nb++;
			}
			// Complete an edge
			edge.pt[edge.npts] = kp;
			edge.npts++;
			if (edge.npts == MAXEDGEPTS) {
				printf("ERROR: MAXEDGEPTS exceeded\n");
				fprintf(fperr,"ERROR: MAXEDGEPTS exceeded\n");
				return 4;
			}
			if (edge.npts < 2) {
				printf("Error: edge.npts: %d\n",edge.npts);
				fprintf(fperr,"Error: edge.npts: %d\n",edge.npts);
				return 5;
			}
			edge.vert[1] = kv;
			if (edge.vert[1] == edge.vert[0]) {
				fprintf(fperr,"Error: TraceSkeleton (b): loop: repeated vertex: edge: %d\n",ne);
//				return 1;
			} else {
				edgeList[ne] = edge;
				if (dbug) {
					fprintf(fperr,"Completed edge (b): %d npts: %d vertices: %d %d  %d %d\n",ne,edge.npts,edge.vert[0],edge.vert[1],
						edge.pt[0],edge.pt[edge.npts-1]);
				}
				ne++;
					if (ne == MAXEDGES) {
						printf("ERROR: ne == MAXEDGES\n");
						fprintf(fperr,"ERROR: ne == MAXEDGES\n");
						return 6;
					}
				if (edge.npts > MINSEGLEN) {
					negmin++;
				}
			}
			if (!end) {
				// Start an edge
				edge.npts = 0;
				edge.pt[0] = kp;
				edge.npts++;
				edge.vert[0] = kv;
				if (dbug) fprintf(fperr,"Start edge: %d at kp: %d  kv: %d  next: %d %d %d\n",ne,kp,kv,next[0],next[1],next[2]);
			}
		}
		if (end) {
			if (dbug) fprintf(fperr,"End\n\n");

//			//Testing
//			if (nv >= 10) break;

			for(;;) {
				if (nb == 0) {
					printf("No more branches\n");
					break;
				}
				kv = branchlist[nb-1];
				if (vertex[kv].nfollowed < vertex[kv].nlinks) break;
				nb--;
			}
			if (nb == 0) break;

			node = 0;
			if (dbug) fprintf(fperr,"                  Redirected to branch: %d  vertex: %d\n",nb,kv);
			nbrs = GetNeighborsOld(vertex[kv].pos,nbrlist);
			for (i=0;i<nbrs;i++) {
				if (!vertex[kv].followed[i]) {
					for (int j=0;j<3;j++) {
						next[j] = nbrlist[i][j];
						prev[j] = vertex[kv].pos[j];
					}
					break;
				}
			}
			// Need to start afresh from vertex kv
			fprintf(fperr,"                  Starting on link #: %d of vertex: %d\n",i,kv);
			vertex[kv].followed[i] = true;
			vertex[kv].nfollowed++;
			kp = PointIndex(vertex[kv].pos);
			// Start a new edge at a previously saved point
			edge.npts = 0;
			edge.pt[0] = kp;
			edge.npts++;
			edge.vert[0] = kv;
			if (dbug) fprintf(fperr,"Restart from vertex: %d  %d\n",kv,kp);
			end = false;
		}
	}
	for (int ie=0;ie<ne;ie++) {
		edgeList[ie].ok = true;
//		int n = edgeList[ie].npts;
//		int k1 = edgeList[ie].pt[0];
//		int k2 = edgeList[ie].pt[n-1];
//		if (n <= 3 || k1 == k2)
//			edgeList[ie].ok = false;
//		else
//			edgeList[ie].ok = true;
	}

	printf("maxnbrs: %d\n",maxnbrs);
	return 0;
}
*/

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
	double ad, len, ddiam, dlen, ltot, lsum, dsum, dvol, r2, r2prev, lsegdtot;
	double ave_len, volume, d95;
	double ave_pt_diam, ave_seg_diam;
	int ie, ip, k, ka, kp, kpprev, ndpts, nlpts, ndtot, nsegdtot, nptstot, nptsusedtot;
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
//		printf("ie: %d npts: %d\n",ie,edge.npts);
		fprintf(fperr,"ie: %d npts: %d npts_used: %d\n",ie,edge.npts,edge.npts_used);
		fflush(fperr);
		nptstot += edge.npts;
		nptsusedtot += edge.npts_used;
		dbug = false;
		kpprev = 0;
		r2prev = 0;
		dsum = 0;
		lsum = 0;
		dvol = 0;
		for (ip=0; ip<edge.npts; ip++) {
			kp = edge.pt[ip];
//			ad = point[kp].d;
			ad = avediameter[kp];
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
				dlen = dist_um(kp,kpprev);
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
		k = int(len/dlen);
		if (len <= lenlimit) continue;
		if (k >= NBOX) {
			printf("Edge too long: ie: %d  k: %d len: %f  k: %d\n",ie,k,len);
			fprintf(fperr,"Edge too long: ie: %d  k: %d len: %f  k: %d\n",ie,k,len);
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
	printf("Total vessel volume: %10.0f\n\n",volume);
	fprintf(fpout,"Total vessel volume: %10.0f\n\n",volume);

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
	printf("nptstot: %d\n",nptstot);
	printf("nptsusedtot: %d\n",nptsusedtot);
	return 0;
}


//-----------------------------------------------------------------------------------------------------
// Create CMGUI files.
// To adjust for the reduced number of nodes in the simplified network, we need a way to handle the 
// node numbers.  The nodes are currently numbered 1...np, but not all these node numbers are used.
// Does the node file need all nodes from 1 to np?
// A quick test seems to show that neither element nor node numbers need to be consecutive.
// This should be made to work with squeezed and unsqueezed data.
// SOME PROBLEM HERE
//-----------------------------------------------------------------------------------------------------
int oldWriteCmguiData(void)
{
	int k, ie, ip, npts, err;
//	int kv0, kv1;
	EDGE edge;
	char dotcomname[256], exelemname[256], exnodename[256];

	printf("WriteCmguiData: np: %d\n",np);
	fprintf(fperr,"WriteCmguiData: np: %d\n",np);
	fflush(fperr);
	err = 0;
	sprintf(dotcomname,"%s.com.txt",output_basename);
	printf("dotcomname: %s\n",dotcomname);
	fprintf(fperr,"dotcomname: %s\n",dotcomname);
	fflush(fperr);
	sprintf(exelemname,"%s.exelem",output_basename);
	printf("exelemname: %s\n",exelemname);
	fprintf(fperr,"exelemname: %s\n",exelemname);
	fflush(fperr);
	sprintf(exnodename,"%s.exnode",output_basename);
	fprintf(fperr,"Com file: %s exelem file: %s exnode file: %s\n",dotcomname,exelemname,exnodename);
	fflush(fperr);
	dotcom = fopen(dotcomname,"w");
	exelem = fopen(exelemname,"w");
	exnode = fopen(exnodename,"w");
	write_com(output_basename);
	write_exelem();
	printf("created exelem header\n");
	fprintf(fperr,"created exelem header\n");
	write_exnode();
	printf("created exnode header\n");
	fprintf(fperr,"created exnode header\n");
	// First set up radius values for each point.  Note that junction points occur on multiple segments,
	// with different radius values.  Select the maximum in this case.
	radius = (float *)malloc(np*sizeof(float));
	point_used = (bool *)malloc(np*sizeof(bool));
	for (k=0; k<np; k++) {
		radius[k] = vsize[0];
		point_used[k] = false;
	}
	int kelem = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;	// Note: after squeezing, all edges are used
//		npts = edge.npts;
		npts = edge.npts_used;

		int kfrom = edge.pt_used[0];
		int kto = edge.pt_used[npts-1];
		if (kfrom == kto) {
			printf("Error: writecmguidata: repeated node: element: %d npts: %d kv0,kv1: %d %d kfrom: %d\n",ie,npts,edge.vert[0],edge.vert[1],kfrom);
			fprintf(fperr,"repeated node: element: %d npts: %d kv0.kv1: %d %d kfrom: %d\n",ie,npts,edge.vert[0],edge.vert[1],kfrom);
			for (k=0; k<npts; k++)
				fprintf(fperr,"  k: %3d  pt: %6d\n",k,edge.pt_used[k]);
			edgeList[ie].used = false;
			continue;
		}
		radius[kfrom] = MAX(radius[kfrom],avediameter[kfrom]/2);
		radius[kto] = MAX(radius[kto],avediameter[kto]/2);
		point_used[kfrom] = true;
//		fprintf(fperr,"edge: %6d npts: %3d kfrom,kto: %6d %6d  %6d\n",ie,npts,kfrom,kto,kelem);
		for (ip=1; ip<npts; ip++) {
			int k2 = edge.pt_used[ip];
			point_used[k2] = true;
			radius[k2] = MAX(radius[k2],avediameter[k2]/2);
			double d = dist_um(kfrom,k2);
// Since simplify() is executed, no need to drop points
//			if (FIXED_DIAMETER > 0 || d > 0.5*(radius[kfrom]+radius[k2]) || ip == npts-1) {
				kelem++;
	            fprintf(exelem, "Element: %d 0 0\n", kelem);
	            fprintf(exelem, "  Nodes: %d %d\n", kfrom+1, k2+1);
	            fprintf(exelem, "  Scale factors: 1 1\n");
				kfrom = k2;
//			}
		}
	}
	for (k=0; k<np; k++) {
		if (point_used[k]) {	// Note: after squeezing, all points are used.
			fprintf(exnode, "Node: %d\n", k+1);
			if (REVISED_VERSION) {
				fprintf(exnode, "%6.1f %6.1f %6.1f\n", vsize[0]*voxel[k].pos[0],vsize[1]*voxel[k].pos[1],vsize[2]*voxel[k].pos[2]);
			} else {
				fprintf(exnode, "%6.1f %6.1f %6.1f\n", vsize[0]*point[k][0],vsize[1]*point[k][1],vsize[2]*point[k][2]);
			}
//			fprintf(exnode, "%6.2f\n", vsize*avediameter[k]/2);
			fprintf(exnode, "%6.2f\n", avediameter[k]/2);
		}
	}
	fclose(dotcom);
	fclose(exelem);
	fclose(exnode);
	free(radius);
	free(point_used);
	return err;
}


//-----------------------------------------------------------------------------------------------------
// A simplified network is produced by dropping close pts
// Edges and vertices are unchanged, but the number of points on an edge is reduced.
// The simplification is applied only if npts >= 5
//-----------------------------------------------------------------------------------------------------
int simplify()
{
	int ie, iv, kp1, kp2, kp3, ip0, ip1, nptot, npts, k;
	int *pt;
	double d12, d23, dmin, diam1, diam2;
	EDGE edge;
	bool dbug;
	double fraction = 0.5;

	printf("simplify\n");

	pt = (int *)malloc(1000*sizeof(int));
	nptot = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		if (edge.npts < 5) continue;
		dbug = false;
		kp1 = edge.pt[0];
		kp3 = edge.pt[edge.npts-1];
//		edgeList[ie].pt[0] = kp1;
		pt[0] = kp1;
		diam1 = getDiameter(edge.pt[0],edge.pt[1],edge.pt[2]);
		if (dbug) printf("diam1: %f  %d %d %d\n",diam1,edge.pt[0],edge.pt[1],edge.pt[2]);
		ip1 = 0;
		ip1++;
		for (ip0=1; ip0<edge.npts; ip0++) {
			kp2 = edge.pt[ip0];
			if (ip0 == edge.npts-1)
				diam2 = getDiameter(edge.pt[ip0-2],edge.pt[ip0-1],edge.pt[ip0]);	//only approx. because pt spacing is small
			else
				diam2 = getDiameter(edge.pt[ip0-1],edge.pt[ip0],edge.pt[ip0+1]);	//only approx. because pt spacing is small
			d12 = zdist(kp1,kp2);
			d23 = zdist(kp2,kp3);
			dmin = fraction*(diam1 + diam2)/2;
			if (dbug) printf("ip0: %d kp2: %d diam2: %f d12,d23: %f %f dmin: %f\n",ip0,kp2,diam2,d12,d23,dmin);
			if ((d12 > dmin &&  d23 > dmin) || ip0 == edge.npts-1) {
//				edgeList[ie].pt[ip1] = kp2;
//				edgeList[ie].pt_used[ip1] = kp2;
				pt[ip1] = kp2;
				if (dbug) printf("pt: %d %d\n",ip1,kp2);
				ip1++;
				kp1 = kp2;
				diam1 = diam2;
				nptot++;
			}
		}
//		edgeList[ie].npts = ip1;
//		edgeList[ie].npts_used = edgeList[ie].npts;
		npts = ip1;
		
		if (npts == 2) {
			kp1 = edge.pt[0];
			kp2 = edge.pt[edge.npts/2];
			kp3 = edge.pt[edge.npts-1];
			pt[0] = kp1;
			pt[1] = kp2;
			pt[2] = kp3;
			npts = 3;
		}
		
		edgeList[ie].npts = npts;
		edgeList[ie].npts_used = npts;
		for (k=0; k<npts; k++) {
			edgeList[ie].pt[k] = pt[k];
			edgeList[ie].pt_used[k] = pt[k];
		}
		if (npts == 2) {
			printf("simplify: ie: %d npts: %d -> %d\n",ie,edge.npts,edgeList[ie].npts);
		}
	}
	printf("done: nptot: %d\n",nptot);
	return 0;
}

/*
//-----------------------------------------------------------------------------------------------------
// Reduce the number of points on an edge by imposing a lower bound DMIN on the distance between points.
//-----------------------------------------------------------------------------------------------------
int simplify1(void)
{
	int ie, j1, j2, k, kk, n, nused;
	double pdist[MAXEDGEPTS], dave, distsum, del, DFAC;

	printf("simplify\n");
	fprintf(fperr,"simplify\n");
	fflush(fperr);
//	DFAC = 1.0;
	DFAC = 3.0;
	np_used = 0;
	for (ie=0;ie<ne;ie++) {
		if (!edgeList[ie].used) continue;
//		printf("  edge: %6d  npts: %4d\n",ie,edgeList[ie].npts);
//		fprintf(fperr,"  edge: %6d  npts: %4d\n",ie,edgeList[ie].npts);
		if (edgeList[ie].npts > MAXEDGEPTS) exit(1);
		distsum = 0;
		dave = 0;
		for (k=1;k<edgeList[ie].npts;k++) {
			j1 = edgeList[ie].pt[k-1];
			j2 = edgeList[ie].pt[k];
			if (j1 == j2) {
				printf("Error: simplify: ie: %d j1=j2: %d\n",ie,j1);
				fprintf(fperr,"Error: simplify: ie: %d j1=j2: %d\n",ie,j1);
				return 1;
			}
			pdist[k] = dist_um(j1,j2);		// distance from pt[k-1] to pt[k]
			distsum += pdist[k];
			if (k == 1) {
				dave += avediameter[j1];
			}
			dave += avediameter[j2];
//			printf("    k: %6d  pdist: %f  %f\n",k,pdist[k],distsum);
		}
//		printf("    distsum: %f dave: %f\n",distsum,dave);
		if (FIXED_DIAMETER > 0) {
			n = edgeList[ie].npts/3. + 0.5;
		} else {
			dave /= edgeList[ie].npts;
			n = distsum/(DFAC*dave) + 0.5;
		}
		n = MAX(n,1);
		del = distsum/n;
//		printf("ie: %d npts: %d n: %d del: %f\n",ie,edgeList[ie].npts,n,del);
		nused = 1;
		edgeList[ie].pt_used[0] = edgeList[ie].pt[0];
		if (n > 2) {
			// Now step along the edge and choose the points closest to multiples of del
			distsum = 0;
			kk = 1;
			for (k=1;k<edgeList[ie].npts-1; k++) {
				distsum += pdist[k];
				if (distsum > kk*del) {
					nused++;
					edgeList[ie].pt_used[kk] = edgeList[ie].pt[k];
					if (ie == 20) printf("edge: %d pt_used: %d  %d\n",kk,k,edgeList[ie].pt[k]);
					kk++;
//					if (kk == n) break;
				}
			}
//			edgeList[ie].pt_used[n] = edgeList[ie].pt[edgeList[ie].npts-1];
//			edgeList[ie].npts_used = n+1;
			edgeList[ie].pt_used[kk] = edgeList[ie].pt[edgeList[ie].npts-1];
			edgeList[ie].npts_used = kk+1;
		} else {
			for (k=1; k<edgeList[ie].npts; k++) {
				edgeList[ie].pt_used[k] = edgeList[ie].pt[k];
			}
			edgeList[ie].npts_used = edgeList[ie].npts;
		}
		np_used += edgeList[ie].npts_used;
		if (edgeList[ie].npts_used != edgeList[ie].npts) {
			printf("simplify: edge: %d dropped: %d\n",edgeList[ie].npts - edgeList[ie].npts_used);
			fprintf(fperr,"simplify: edge: %d dropped: %d\n",edgeList[ie].npts - edgeList[ie].npts_used);
		}
	}
	return 0;
}
*/

//-----------------------------------------------------------------------------------------------------
// When an edge (edrop) is removed and a vertex (kv) is reduced to two used edges, the two edges are 
// joined to create a single edge.
//-----------------------------------------------------------------------------------------------------
int joiner(int kv)
{
	int i, ie, ive[2], n, nep[2], npts, newvert[2], kv1;
	int *ptemp;
	bool reverse[2];
	EDGE *newedge;
	bool dbug = false;

	if (dbug) printf("\njoiner: kv: %d\n",kv);
	n = 0;
	for (i=0; i<vertex[kv].nlinks; i++) {
		ie = vertex[kv].edge[i];
		if (edgeList[ie].used) {
			ive[n] = ie;
			nep[n] = edgeList[ie].npts;
			if (dbug) printf("link: %d  %d npts: %d\n",n,ie,nep[n]);
			n++;
			if (n > 2) {
				printf("Error: joiner: n: %d\n",n);
				exit(1);
			}
		}
	}
	// ive[2] are the two connected used edges.  Keep ive[0]
	newedge = &edgeList[ive[0]];
	// Save the points of ive[0] in ptemp
	ptemp = (int *)malloc(nep[0]*sizeof(int));
	for (i=0; i<nep[0]; i++)
		ptemp[i] = newedge->pt[i];
	free(newedge->pt);
	free(newedge->pt_used);
//	free(newedge->avediameter);
	npts = nep[0] + nep[1] - 1;
	newedge->pt = (int *)malloc(npts*sizeof(int));
	newedge->pt_used = (int *)malloc(npts*sizeof(int));
//	newedge->avediameter = (float *)malloc(npts*sizeof(float));
	newedge->npts = npts;
	newedge->npts_used = npts;
	if (kv == newedge->vert[0]) {
		newvert[0] = newedge->vert[1];
		reverse[0] = true;
	} else {
		newvert[0] = newedge->vert[0];
		reverse[0] = false;
	}
	if (kv == edgeList[ive[1]].vert[0]) {
		newvert[1] = edgeList[ive[1]].vert[1];
		reverse[1] = false;
	} else {
		newvert[1] = edgeList[ive[1]].vert[0];
		reverse[1] = true;
	}
	for (i=0; i<2; i++)
		newedge->vert[i] = newvert[i];
	// Copy pts from edge ive[0]
	if (reverse[0]) {
		for (i=0; i<nep[0]; i++)
			newedge->pt[i] = ptemp[nep[0]-i-1];
	} else {
		for (i=0; i<nep[0]; i++)
			newedge->pt[i] = ptemp[i];
	}
	// Copy pts from edge ive[1]
	if (reverse[1]) {
		for (i=0; i<nep[1]-1; i++)
			newedge->pt[nep[0]+i] = edgeList[ive[1]].pt[nep[1]-i-2];
	} else {
		for (i=0; i<nep[1]-1; i++)
			newedge->pt[nep[0]+i] = edgeList[ive[1]].pt[i+1];
	}
	for (i=0; i<newedge->npts; i++) {
		newedge->pt_used[i] = newedge->pt[i];
//		newedge->avediameter[i] = FIXED_DIAMETER;
	}
	if (dbug) {
		printf("newedge: %d npts: %d vert: %d %d\n",ive[0],newedge->npts,newedge->vert[0],newedge->vert[1]);
		printf("pts: ");
		for (int j=0; j<newedge->npts; j++)
			printf("%d ",newedge->pt[j]);
		printf("Remove edge: %d  vert: %d %d\n",ive[1],edgeList[ive[1]].vert[0],edgeList[ive[1]].vert[1]);
	}
	edgeList[ive[1]].used = false;
	vertex[kv].used = false;
	kv1 = newedge->vert[1];
	// Need to locate edge ive[1] in the list of edges for vertex newedge->vert[1], and replace it with ive[0].
	for (i=0; i<vertex[kv1].nlinks; i++) {
		if (vertex[kv1].edge[i] == ive[1]) {
			vertex[kv1].edge[i] = ive[0];
			break;
		}
	}
	free(ptemp);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Look for vertices with two links, and join the two connected edges into a single edge.
// This creates spurious vessels!!!!
//-----------------------------------------------------------------------------------------------------
void checkVerticies(bool join)
{
	int kv, ie, nlinks, nlinks_used, n2, j, k;
	EDGE *edge;

	printf("checkVerticies\n");
	fprintf(fpout,"checkVerticies\n");
	for (kv=0; kv < nv; kv++) {
		vertex[kv].used = false;
		vertex[kv].nlinks = 0;			//---
		vertex[kv].nlinks_used = 0;		//---
	}
//	printf("zeroed nlinks\n");
	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
//		printf("ie: %d vert: %d %d\n",ie,edge->vert[0],edge->vert[1]);
		kv = edge->vert[0];
		vertex[kv].edge[vertex[kv].nlinks_used] = ie;	//---
		vertex[kv].nlinks++;			//---
		vertex[kv].nlinks_used++;		//---
		vertex[kv].used = true;	
//		if (kv == 0) fprintf(fpout,"vertex 0 on edge: %d\n",ie);
		kv = edge->vert[1];
		vertex[kv].edge[vertex[kv].nlinks_used] = ie;	//---
		vertex[kv].used = true;
		vertex[kv].nlinks++;			//---
		vertex[kv].nlinks_used++;		//---
//		if (kv == 0) fprintf(fpout,"vertex 0 on edge: %d\n",ie);
	}
	if (!join) return;
	n2 = 0;
	for (kv=0; kv < nv; kv++) {
		if (!vertex[kv].used) continue;
		nlinks_used = vertex[kv].nlinks;
		if (nlinks_used == 2) {
			n2++;
//			fprintf(fpout,"vertex: %d nlinks: %d\n",kv,nlinks_used);
			for (j=0; j<2; j++) {
				ie = vertex[kv].edge[j];
//				fprintf(fpout,"edge: %d %d npts: %d vert: %d %d ",j,ie,edgeList[ie].npts,edgeList[ie].vert[0],edgeList[ie].vert[1]);
//				for (k=0; k<edgeList[ie].npts; k++) {
//					fprintf(fpout,"%d ",edgeList[ie].pt[k]);
//				}
//				fprintf(fpout,"\n");
			}
			joiner(kv);					//---
		}
	}
	printf("Number of 2-verticies: %d\n",n2);
}

//-----------------------------------------------------------------------------------------------------
// Check consistency of vertex nlinks_used and edges
//-----------------------------------------------------------------------------------------------------
int checker(void)
{
	int ie, kv, k, ivox;
	int *nlused;
	bool err;

	printf("checker: nv: %d\n",nv);
	fprintf(fperr,"checker: nv: %d\n",nv);
	err = false;
	nlused = (int *)malloc(nv*sizeof(int));
	for (kv=0; kv<nv; kv++) {
		nlused[kv] = 0;
	}
	for (ie=0; ie<ne; ie++) {
		if (!edgeList[ie].used) continue;
		for (k=0; k<2; k++) {
			kv = edgeList[ie].vert[k];
			ivox = vertex[kv].ivox;
			nlused[kv]++;
		}
	}
	for (kv=0; kv<nv; kv++) {
		if (!vertex[kv].used) continue;
		if (nlused[kv] > 0 && nlused[kv] != vertex[kv].nlinks_used) {
			printf("checker: vertex: %d nlinks: %d nlinks_used: %d used: %d\n",kv,vertex[kv].nlinks,vertex[kv].nlinks_used,nlused[kv]);
			fprintf(fperr,"checker: vertex: %d nlinks: %d nlinks_used: %d used: %d\n",kv,vertex[kv].nlinks,vertex[kv].nlinks_used,nlused[kv]);
			err = true;
		}
		//if (nlused[kv] == 2) {
		//	printf("checker: vertex: %d nlinks: %d\n",kv,2);
		//	err = true;
		//}
		ivox = vertex[kv].ivox;
	}
	free(nlused);
	if (err) return 1;
	printf("OK\n");
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int checkEdge(int ie)
{
	int i, kv0, kv1;
	bool connected = false;;

	kv0 = edgeList[ie].vert[0];
	kv1 = edgeList[ie].vert[1];
	for (i=0; i<ne; i++) {
		if (!edgeList[i].used) continue;
		if (i == ie) continue;
		if (edgeList[i].vert[0] == kv0 || edgeList[i].vert[0] == kv1) {
			connected = true;
			return 0;
		}
		if (edgeList[i].vert[1] == kv0 || edgeList[i].vert[1] == kv1) {
			connected = true;
			return 0;
		}
	}
	printf("Edge: %d not connected\n",ie);
	exit(1);
}

//-----------------------------------------------------------------------------------------------------
// A loop is broken if its edges are still used.
// Edges j1, j2 and j3 form a loop.
// Each vertex must be checked to see if it has external connections.  There are three cases:
// (1) One vertex is connected externally (this is a dead-end).
// (2) Two vertices are connected externally (side loop)
// (3) Three vertices are connected externally (3-way junction loop)
//-----------------------------------------------------------------------------------------------------
int fixloop(int j1, int j2, int j3)
{
	int e[3];
	int vert[3];
	bool used[3];
	int nconn[3];
	int i, j, ie,  kv, kv0, kv1, ncmin, nc2, imax, n, edrop;
	double d, dmax;
	bool dbug;

	dbug = false;
	if (dbug) {
		printf("fixloop: edges: %d %d %d\n",j1,j2,j3);
		fprintf(fperr,"fixloop: edges: %d %d %d\n",j1,j2,j3);
	}
	e[0] = j1;
	e[1] = j2;
	e[2] = j3;
	for (i=0; i<3; i++) {
		used[i] = edgeList[e[i]].used;
		nconn[i] = 0;
		kv0 = edgeList[e[i]].vert[0];
		kv1 = edgeList[e[i]].vert[1];
		if (dbug) {
			printf("fixloop: edge: %d vertices: %d %d  nlinks: %d %d\n",e[i],kv0,kv1,vertex[kv0].nlinks,vertex[kv1].nlinks);
			fprintf(fperr,"fixloop: edge: %d vertices: %d %d  nlinks: %d %d\n",e[i],kv0,kv1,vertex[kv0].nlinks,vertex[kv1].nlinks);
		}
	}
	vert[0] = edgeList[e[0]].vert[0];
	vert[1] = edgeList[e[0]].vert[1];
	if (edgeList[e[1]].vert[0] == vert[0] || edgeList[e[1]].vert[0] == vert[1]) {
		vert[2] = edgeList[e[1]].vert[1];
	} else {
		vert[2] = edgeList[e[1]].vert[0];
	}
	// nconn[j] is the number of used edges connected to vertex vert[j]
	for (i=0; i<3; i++) {
		n = 0;
		for (j=0; j<vertex[vert[i]].nlinks; j++) {
			ie = vertex[vert[i]].edge[j];
			if (edgeList[ie].used) nconn[i]++;
		}
	}
	if (dbug) {
		printf("nconn: %d %d %d\n",nconn[0],nconn[1],nconn[2]);
		fprintf(fperr,"nconn: %d %d %d\n",nconn[0],nconn[1],nconn[2]);
	}
	nc2 = 0;
	ncmin = 999;
	for (i=0; i<3; i++) {
		if (nconn[i] < ncmin) ncmin = nconn[i];
		if (nconn[i] == 2) nc2++;
	}
	if (nc2 > 0) {
		edrop = -1;
		for (i=0; i<3; i++) {
			kv0 = edgeList[e[i]].vert[0];
			kv1 = edgeList[e[i]].vert[1];
			for (j=0; j<3; j++) {
				if (nconn[j] == 2 && (vert[j] == kv0 || vert[j] == kv1)) {
					edrop = e[i];
					edgeList[edrop].used = false;
					if (dbug) {
						printf("Drop edge (a): %d  vert: %d %d  nlinks: %d %d\n",edrop,kv0,kv1,vertex[kv0].nlinks,vertex[kv1].nlinks);
						fprintf(fperr,"Drop edge (a): %d  vert: %d %d  nlinks: %d %d\n",edrop,kv0,kv1,vertex[kv0].nlinks,vertex[kv1].nlinks);
					}
					break;
				}
			}
			if (edrop >= 0) break;
		}
	} else if (nc2 == 0) {	// DO NOT arbitrarily choose edge j1, choose longest
		dmax = 0;
		for (i=0; i<3; i++) {
			kv0 = edgeList[e[i]].vert[0];
			kv1 = edgeList[e[i]].vert[1];
			d = dist_um(kv0,kv1);
			if (d > dmax) {
				imax = i;
				dmax = d;
			}
		}
		edrop = e[imax];
		edgeList[edrop].used = false;
		if (dbug) {
			printf("Drop edge (b): %d\n",edrop);
			fprintf(fperr,"Drop edge (b): %d\n",edrop);
		}
	}
//	return edrop;
	kv0 = edgeList[edrop].vert[0];
	kv1 = edgeList[edrop].vert[1];
	if (dbug) {
		printf("Dropped edge: %d  vert: %d %d  nlinks: %d %d\n",edrop,kv0,kv1,vertex[kv0].nlinks_used,vertex[kv1].nlinks_used);
		fprintf(fperr,"Dropped edge: %d  vert: %d %d  nlinks: %d %d\n",edrop,kv0,kv1,vertex[kv0].nlinks_used,vertex[kv1].nlinks_used);
	}

	// Now need to adjust edges and vertices to account for the removal of edge edrop.
	// Possibly two linked edges at a vertex need to be joined, in which case the vertex goes unused.
	// One edge is extended, the other goes unused.
	
	for (i=0; i<2; i++) {
		kv = edgeList[edrop].vert[i];
		n = 0;
		for (j=0; j<vertex[kv].nlinks; j++) {
			ie = vertex[kv].edge[j];
			if (edgeList[ie].used) {
				n++;
			}
		}
		vertex[kv].nlinks_used = n;
	}
	
	return edrop;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int rejoin(int ndropped, int dropped[])
{
	int i, k, edrop, kv0, kv1, kv, n, j, ie;
	bool dbug = false;

	printf("rejoin\n");
	for (k=0; k<ndropped; k++) {
		edrop = dropped[k];
		kv0 = edgeList[edrop].vert[0];
		kv1 = edgeList[edrop].vert[1];
		if (dbug) {
			printf("Dropped edge: %d  vert: %d %d  nlinks: %d %d\n",edrop,kv0,kv1,vertex[kv0].nlinks_used,vertex[kv1].nlinks_used);
			fprintf(fperr,"Dropped edge: %d  vert: %d %d  nlinks: %d %d\n",edrop,kv0,kv1,vertex[kv0].nlinks_used,vertex[kv1].nlinks_used);
		}

		// Now need to adjust edges and vertices to account for the removal of edge edrop.
		// Possibly two linked edges at a vertex need to be joined, in which case the vertex goes unused.
		// One edge is extended, the other goes unused.
		for (i=0; i<2; i++) {
			kv = edgeList[edrop].vert[i];
			n = 0;
			for (j=0; j<vertex[kv].nlinks; j++) {
				ie = vertex[kv].edge[j];
				if (edgeList[ie].used) {
					n++;
				}
			}
			vertex[kv].nlinks_used = n;
			if (n == 2) {
				if (dbug) {
					printf("joiner: kv: %d\n",kv);
					fprintf(fperr,"joiner: kv: %d\n",kv);
				}
				joiner(kv);
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int checkUnconnected(void)
{
	int i, kv0, kv1, n0, n1, ii, err;
	EDGE edge;

	printf("Check for unconnected edges\n");
	// Step through edges, looking for unconnected edges.
	err = 0;
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
		if (n0+n1 == 0 && ne > 1) {
			printf("Error: edge: %d is unconnected: npts: %d vert: %d %d\n",i,edge.npts,kv0,kv1);
			fprintf(fperr,"Error: edge: %d is unconnected: npts: %d vert: %d %d\n",i,edge.npts,kv0,kv1);
			err = 1;
		}
	}
	return err;
}


//-----------------------------------------------------------------------------------------------------
// Search for loops
// deloop works fine.  Should remove longest side when all three have nconn = 3
// Not perfect: sometimes creates a vertex with only two links
//-----------------------------------------------------------------------------------------------------
int deloop(int iter)
{
	int i, ii, kv0, kv1, kkv0, kkv1, kv2, kv3, edrop;
	int npairs, ne2, j1, j2, j3, nloops, ndropped, NE2MAX;
//	int dropped[1000000];
	int *dropped;
	EDGE edge, eedge;
	PAIR *pair;
	int *e2;
	bool dup;
	bool NEW_dup_delete = true;		// try this!

	printf("deloop: iter: %d\n",iter);
	printf("deloop: iter: %d\n",iter);
	printf("deloop: iter: %d\n",iter);
	fprintf(fperr,"deloop: iter: %d\n",iter);

	ne2 = 0;
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		if (edge.used && edge.npts == 2) {
			ne2++;
		}
	}
	NE2MAX = 10*ne2;
	e2 = (int *)malloc(NE2MAX*sizeof(int));
	pair = (PAIR *)malloc(NE2MAX*sizeof(PAIR));
	ne2 = 0;
	for (i=0; i<ne; i++) {
		edge = edgeList[i];
		if (!edge.used) continue;
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		if (edge.npts == 2) {
			if (ne2 == NE2MAX) {
				printf("Error: deloop: array dimension e2 exceeded\n");
				fprintf(fperr,"Error: deloop: array dimension e2 exceeded\n");
				return -1;
			}
//			printf("2-point edge: %6d vert: %6d %6d\n",i,kv0,kv1);
			e2[ne2] = i;
			ne2++;
			// Weed out duplicates
			dup = false;
			for (ii=i+1; ii<ne; ii++) {
				eedge = edgeList[ii];
				if (!eedge.used) continue;
				kkv0 = eedge.vert[0];
				kkv1 = eedge.vert[1];
				// This finds a number of edges
				if ((kkv0 == kv0 && kkv1 == kv1) || (kkv0 == kv1 && kkv1 == kv0)) {
					dup = true;
					if (eedge.npts == 2) {
						printf("Error: deloop: not using point\n");
						fprintf(fperr,"Error: deloop: not using point\n");
						return -2;
					}
					if (NEW_dup_delete) {
						edgeList[ii].used = false;
						vertex[kkv0].nlinks_used--;
						vertex[kkv1].nlinks_used--;
						printf("dup: dropped edge: %d  vert: %d %d\n",ii,kkv0,kkv1);
					}
					break;
				}
			}
			if (dup && !NEW_dup_delete) {
				ne2--;
				edgeList[i].used = false;
				vertex[kv0].nlinks_used--;
				vertex[kv1].nlinks_used--;
				printf("dup: dropped edge: %d  vert: %d %d\n",i,kv0,kv1);
				if (iter >= 0) {
					fprintf(fperr,"dup: dropped edge: %d  kv0,kv1: %d %d\n",i,kv0,kv1);
				}
			}
		}
	}
	printf("Number of 2-edges (edges with two points): %d\n",ne2);
	fprintf(fperr,"Number of 2-edges (edges with two points): %d\n",ne2);

	// Find all pairs of connected 2-edges
	/*
	npairs = 0;
	for (j1=0; j1<ne2; j1++) {
		edge = edgeList[e2[j1]];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		for (j2=j1+1; j2<ne2; j2++) {
			eedge = edgeList[e2[j2]];
			if (eedge.vert[0] == kv1) {			// -->  -->
				pair[npairs].i1 = e2[j1];
				pair[npairs].i2 = e2[j2];
				npairs++;
			} else if (eedge.vert[1] == kv1) {	// -->  <--
				pair[npairs].i1 = e2[j1];
				pair[npairs].i2 = -e2[j2];
				npairs++;
			} else if (eedge.vert[0] == kv0) {	// <--  -->
				pair[npairs].i1 = -e2[j1];
				pair[npairs].i2 = e2[j2];
				npairs++;
			} else if (eedge.vert[1] == kv0) {	// <--  <--
				pair[npairs].i1 = -e2[j1];
				pair[npairs].i2 = -e2[j2];
				npairs++;
			}
			if (npairs == NE2MAX) {
				printf("Error: npairs == NE2MAX: %d  %d\n",npairs,NE2MAX);
				fprintf(fperr,"Error: npairs == NE2MAX: %d  %d\n",npairs,NE2MAX);
				return -3;
			}
		}
	}
	
	printf("Number of 2-edge pairs (connections): %d\n",npairs);
	fprintf(fperr,"Number of 2-edge pairs (connections): %d\n",npairs);
	
	bool hit;
	nloops = 0;
	for (i=0; i<npairs; i++) {
		for (ii=i+1; ii<npairs; ii++) {
			hit = false;
			if (abs(pair[i].i2) == abs(pair[ii].i1)) {
				j1 = pair[i].i1;
				j2 = abs(pair[i].i2);
				j3 = pair[ii].i2;
				hit = true;
			} else if (abs(pair[i].i1) == abs(pair[ii].i2)) {
				j1 = pair[i].i2;
				j2 = abs(pair[i].i1);
				j3 = pair[ii].i1;
				hit = true;
			} else if (abs(pair[i].i2) == abs(pair[ii].i2)) {
				j1 = pair[i].i1;
				j2 = abs(pair[i].i2);
				j3 = pair[ii].i1;
				hit = true;
			} else if (abs(pair[i].i1) == abs(pair[ii].i1)) {
				j1 = pair[i].i2;
				j2 = abs(pair[i].i1);
				j3 = pair[ii].i2;
				hit = true;
			}
			if (hit) {
				if (j1 > 0)
					k1 = edgeList[j1].vert[0];
				else
					k1 = edgeList[-j1].vert[1];
				if (j3 > 0)
					k3 = edgeList[j3].vert[1];
				else
					k3 = edgeList[-j3].vert[0];
				if (k1 == k3) {
					nloops++;
					if (iter > 0) {
						printf("fixloop: nloops: %d j1,j2,j3: %d %d %d  k1,k3: %d %d\n",nloops,abs(j1),abs(j2),abs(j3),k1,k3);
						fprintf(fperr,"fixloop: nloops: %d j1,j2,j3: %d %d %d  k1,k3: %d %d\n",nloops,abs(j1),abs(j2),abs(j3),k1,k3);
					}
					fixloop(abs(j1),abs(j2),abs(j3));
				}
			}
		}
	}
	*/

	// Change the method.  The edges in pair, pair[].i1 and pair[].i2, are given signs:
	//		+ if the edge is directed towards the common vertex
	//      - if the edge is directed away from the common vertex
	npairs = 0;
	for (j1=0; j1<ne2; j1++) {
		edge = edgeList[e2[j1]];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		for (j2=j1+1; j2<ne2; j2++) {
			eedge = edgeList[e2[j2]];
			if (eedge.vert[0] == kv1) {			// -->  -->
				pair[npairs].i1 = e2[j1];
				pair[npairs].i2 = -e2[j2];
				npairs++;
			} else if (eedge.vert[1] == kv1) {	// -->  <--
				pair[npairs].i1 = e2[j1];
				pair[npairs].i2 = e2[j2];
				npairs++;
			} else if (eedge.vert[0] == kv0) {	// <--  -->
				pair[npairs].i1 = -e2[j1];
				pair[npairs].i2 = -e2[j2];
				npairs++;
			} else if (eedge.vert[1] == kv0) {	// <--  <--
				pair[npairs].i1 = -e2[j1];
				pair[npairs].i2 = e2[j2];
				npairs++;
			}
			if (npairs == NE2MAX) {
				printf("Error: npairs == NE2MAX: %d  %d\n",npairs,NE2MAX);
				fprintf(fperr,"Error: npairs == NE2MAX: %d  %d\n",npairs,NE2MAX);
				return -3;
			}
		}
	}
	
	printf("Number of 2-edge pairs (connections): %d\n",npairs);
	fprintf(fperr,"Number of 2-edge pairs (connections): %d\n",npairs);

	// For a triangle (loop) there must be corresponding edges in each pair with opposite sign

	dropped = (int *)malloc(1000000*sizeof(int));

	bool hit;
	nloops = 0;
	ndropped = 0;
	for (i=0; i<npairs; i++) {
		int ie1 = pair[i].i1;
		int ie2 = pair[i].i2;
		if (!edgeList[abs(ie1)].used || !edgeList[abs(ie2)].used) continue;
		for (ii=i+1; ii<npairs; ii++) {
			int iie1 = pair[ii].i1;
			int iie2 = pair[ii].i2;
			if (!edgeList[abs(iie1)].used || !edgeList[abs(iie2)].used) continue;
			hit = false;
			if (ie1 == -iie1) {
				j1 = ie1;
				j2 = ie2;
				j3 = iie2;
				hit = true;
			} else if (ie1 == -iie2) {
				j1 = ie1;
				j2 = ie2;
				j3 = iie1;
				hit = true;
			} else if (ie2 == -iie1) {
				j1 = ie2;
				j2 = ie1;
				j3 = iie2;
				hit = true;
			} else if (ie2 == -iie2) {
				j1 = ie2;
				j2 = ie1;
				j3 = iie1;
				hit = true;
			}
			if (hit) {		// For a loop the other verticies of edges j2 and j3 must be the same
				if (!edgeList[abs(j1)].used || !edgeList[abs(j2)].used  || !edgeList[abs(j3)].used) continue;
				if (j2 > 0) 
					kv2 = edgeList[j2].vert[0];
				else
					kv2 = edgeList[-j2].vert[1];
				if (j3 > 0) 
					kv3 = edgeList[j3].vert[0];
				else
					kv3 = edgeList[-j3].vert[1];
				if (kv2 == kv3) {
					nloops++;
					if (iter > 0) {
						printf("fixloop: nloops: %d j1,j2,j3: %d %d %d\n",nloops,abs(j1),abs(j2),abs(j3));
						fprintf(fperr,"fixloop: nloops: %d j1,j2,j3: %d %d %d\n",nloops,abs(j1),abs(j2),abs(j3));
					}
					edrop = fixloop(abs(j1),abs(j2),abs(j3));
					dropped[ndropped] = edrop;
					ndropped++;
				}
			}
		}
	}
	free(e2);
	free(pair);
	printf("Number of loops removed: %d\n",nloops);
	fprintf(fperr,"Number of loops removed: %d\n",nloops);
	rejoin(ndropped,dropped);
	printf("did rejoin: ndropped: %d\n",ndropped);
	// Now check for any verticies with only two links, rejoin edges
	//for (int kv=0; kv<nv; kv++) {
	//	if (vertex[kv].used) continue;
	//	if (vertex[kv].nlinks_used == 2) {
	//		joiner(kv);
	//	}
	//}
	free(dropped);
	return nloops;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int edgePtDist(void)
{
	int i, ie, n, ne_used;
	int ecount[MAXEDGEPTS];

	printf("edgePtDist: ne: %d\n",ne);
	fprintf(fperr,"Vertices: %d\n",nv);
	//for (i=0;i<nv;i++)
	//	fprintf(fperr,"vertex: %4d  ivox: %6d  nlinks: %d\n",i,vertex[i].ivox,vertex[i].nlinks);
	for (i=0; i<MAXEDGEPTS; i++)
		ecount[i] = 0;
	ne_used = 0;
	for (ie=0; ie<ne; ie++) {
		if (!edgeList[ie].used) continue;
		ne_used++;
		n = edgeList[ie].npts_used;
		ecount[n]++;
	}
	printf("Number of edges used: %d\n",ne_used);
	//for (i=0; i<20; i++) {
	//	printf("%4d %8d\n",i,ecount[i]);
	//	if (i > 10 && ecount[i] == 0) break;
	//}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Follow the branch with ID = link to the next vertex, either branching or dead-end.
//-----------------------------------------------------------------------------------------------------
int follower(int kvert, int link, int *link_next, int *nevox, int evox[])
{
	int ivox1, ivox2, link1, nbrs, i;

//	printf("follower: kvert: %d  link: %d\n",kvert,link);
	*nevox = 0;
	link1 = link;
	ivox1 = vertex[kvert].ivox;
	evox[*nevox] = ivox1;
	(*nevox)++;
	for (;;) {
		ivox2 = voxel[ivox1].nid[link1];		// index to first voxel on this branch
		nbrs = voxel[ivox2].nbrs;				// number of voxel neighbours
//		printf("follower: link1, ivox1: %d %d ivox2: %d nbrs: %d\n",link1,ivox1,ivox2,nbrs);
		evox[*nevox] = ivox2;
		(*nevox)++;
		if (*nevox == MAXEDGEPTS) {
			printf("Error: follower: MAXEDGEPTS exceeded\n");
			exit(1);
		}
		if (nbrs == 1) {						// this is a dead-end vertex
			*link_next = 0;
			return voxel[ivox2].ivert;			// return the vertex ID
		} else if (nbrs > 2) {					// this is a branching vertex, need to identify the link number
			for (i=0; i<nbrs; i++) {
				if (voxel[ivox2].nid[i] == ivox1) {
					*link_next = i;
					return voxel[ivox2].ivert;	// return the vertex ID
				}
			}
			printf("Error: follower: we should never get here!\n");
			exit(1);
		} else {								// this is a simple link voxel
			if (voxel[ivox2].nid[0] == ivox1)	// choose the link number for the outgoing link
				link1 = 1;
			else
				link1 = 0;
			ivox1 = ivox2;
		}
	}
}


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int tracer(void)
{
	int i, kvert0, kvert, link, kvisit, ivox, nshort, nloop, isweep;
	int kvert_next, link_next;
	int *visited;
	int nvisited, ndeadend;
	int nevox, evox[MAXEDGEPTS];
	int ecount[MAXEDGEPTS];
	bool todo;

	for (i=0; i<MAXEDGEPTS; i++)
		ecount[i] = 0;
	visited = (int *)malloc(nv*sizeof(int));
	// First choose a starting vertex, one with nlinks = 1
	for (i=0; i<nv; i++) {
		if (vertex[i].nlinks == 1) {
			kvert0 = i;
			break;
		}
	}

	for (isweep=0; isweep<2; isweep++) {
		ivox = vertex[kvert0].ivox;
		printf("tracer: starting vertex: %d  %d %d %d\n",kvert0,voxel[ivox].pos[0],voxel[ivox].pos[1],voxel[ivox].pos[2]);
		printf("sweep: %d\n",isweep);
		nshort = 0;
		nloop = 0;
		ne = 0;
		ndeadend = 0;
		nvisited = 0;
		for (kvert=0; kvert<nv; kvert++) {
			visited[kvert] = 0;
			vertex[kvert].nfollowed = 0;
			vertex[kvert].ivisit = -1;
			for (i=0; i<vertex[kvert].nlinks; i++) {
				vertex[kvert].followed[i] = false;
			}
		}
		visited[nvisited] = kvert0;			// the starting vertex kvert0 is the only entry in the visited list
		vertex[kvert0].ivisit = nvisited;	// the visit list position of kvert0 is nvisited = 0
		nvisited++;

		for (;;) {
			todo = false;
			for (kvisit=nvisited-1; kvisit>=0; kvisit--) {	// traverse the visited list looking for the last with nfollowed < nlinks
				kvert = visited[kvisit];
				if (vertex[kvert].nfollowed < vertex[kvert].nlinks) {
					todo = true;
					break;
				}
			}
			if (!todo) break;
			// Choose the first unfollowed link to follow
			for (link=0; link<vertex[kvert].nlinks; link++) {
				if (!vertex[kvert].followed[link]) break;
			}
			if (ne == 0) printf("unfollowed link of vertex: %d %d nlinks: %d\n",kvert,link,vertex[kvert].nlinks);
			// Now follow this link, to another vertex or to a dead-end
			// The follower returns the next vertex ID, and the link number on this vertex
			kvert_next = follower(kvert,link,&link_next,&nevox,evox);
			if (ne == 0) printf("kvert, kvert_next: %d %d link_next: %d nevox: %d\n",kvert,kvert_next,link_next,nevox);

			if (kvert_next == kvert) {
				nloop++;
			}
			if (vertex[kvert_next].nlinks == 1 && nevox < 0) {
				nshort++;
				if (ne == 0) printf("nshort: %d\n",nshort);
				if (isweep == 1)
					printf("short end: vertices: %d  %d nevox: %d\n",kvert,kvert_next,nevox);
			} else {
				if (isweep == 1) {		// store this edge
					edgeList[ne].used = true;
					edgeList[ne].vert[0] = kvert;
					edgeList[ne].vert[1] = kvert_next;
					edgeList[ne].npts = nevox;
					edgeList[ne].npts_used = nevox;
					edgeList[ne].pt = (int *)malloc(nevox*sizeof(int));
					edgeList[ne].pt_used = (int *)malloc(nevox*sizeof(int));
//					edgeList[ne].avediameter = (float *)malloc(nevox*sizeof(float));
					for (i=0; i<nevox; i++) {
						edgeList[ne].pt[i] = evox[i];
						edgeList[ne].pt_used[i] = evox[i];
//						edgeList[ne].avediameter[i] = FIXED_DIAMETER;
					}
				}
				ne++;
				ecount[nevox]++;
				if (ne <= 1) printf("nevox: %d\n",nevox);
			}

			vertex[kvert].edge[link] = ne-1;
			vertex[kvert].followed[link] = true;
			vertex[kvert].nfollowed++;
			vertex[kvert_next].edge[link_next] = ne-1;
			vertex[kvert_next].followed[link_next] = true;
			vertex[kvert_next].nfollowed++;
			if (ne <= 1) printf("vertex: %d link: %d edge: %d\n",kvert,link,ne-1);
			if (vertex[kvert_next].ivisit < 0) {	// If this is the first visit to kvert_next
				visited[nvisited] = kvert_next;
				vertex[kvert_next].ivisit = nvisited;
				if (ne <= 1) printf("first visit to: %d ivisit: %d\n",kvert_next,nvisited);
				nvisited++;
			}
		}
		if (isweep == 0) {
			printf("Number of edges: %d  nshort: %d  nloop: %d\n",ne,nshort,nloop);
			//for (i=0; i<20; i++) {
			//	printf("%4d %8d\n",i,ecount[i]);
			//	if (i > 10 && ecount[i] == 0) break;
			//}
			ne_max = 1.5*ne;
			edgeList = (EDGE *)malloc(ne_max*sizeof(EDGE));
			printf("Allocated edgeList: %d\n",ne);
		}
	}
	free(visited);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Point kp0 is the end of an edge that has a final segment with orientation of the unit vector e0,
// point kp1 is the end of an edge that has a final segment with orientation of the unit vector e1.
// We seek a point P equidistant from P0 and P1 such that the sum of squares of angles along the
// connection P0-P-P1 is minimised.  The three angles are:
//   the angle e0 makes with P0-P = alpha1
//   the angle P0-P makes with P-P1 = alpha2
//   the angle P-P1 makes with -e1 = alpha3
//   Q = alpha1^2 + alpha2^2 + alpha3^2 is to be minimised.
// 
// At the end P0 we start by defining three orthogonal unit vectors that become the axes of
// a local coordinate system, represented conveniently in spherical coordinates.
// Define Vx as the unit vector corresponding to P0-P1.
// Define Vy as any unit vector perpendicular to Vx.
// Define Vz as the unit vector perpendicular to Vx and Vy, Vz = Vx x Vy (cross-product).
//
// Then a general unit vector is given by:
//   Vu = x.Vx + y.Vy + z.Vz
// subject to:
//   x^2 + y^2 + z^2 = 1 (which implies that -1<=x<=1, -1<=y<=1, -1<=z<=1)
// A general unit vector in this coordinate system can be parametrised by theta and phi,
// where theta is the angle with the Z axis (Vx) 0 <= theta <= PI
// and phi is the angle about the Z axis, starting from 0 at the X axis, 0 <= phi <= 2PI
// where:
//   x = sin(theta).cos(phi)
//   y = sin(theta).sin(phi)
//   z = cos(theta)
// For phi = 0, a line from P0 in the direction Vu, with L = d/(2.sin(theta)), 
// ends at a point P that is equidistant from P0 and P1.
// If the line segment length is L, for any phi, x = d/2 ==> L = d/(2.sin(theta).cos(phi))
// Then the parameters (theta, phi) give the point P, offset by the corresponding vector from P0, as:
//   P = P0 + L.(x.Vx + y.Vy + z.Vz)
// From e0, P-P0, P1-P, and e1 the angles alpha1, alpha2, and alpha3 can be computed, hence Q.
// We need to determine the (theta, phi) pair that minimises Q.
// Let t = theta, p = phi.
// Q = Q(t,p)
// The minimum occurs where dQ/dt = 0 and dQ/dp = 0.
// Let F1(t,p) = dQ/dt,
//     F2(t,p) = dQ/dp.
// We now have a two-dimensional root-finding problem, where the evaluation of F1 and F2 is numerical.
// Using Newton-Raphson, the derivatives dFi/dt and dFi/dp (i = 1,2) are numerical.
// dF1/dt = d/dt(dQ/dt) = (Q(t+dt,p) - 2.Q(t,p) + Q(t-dt,p))/(dt^2)
// dF1/dp = d/dp(dQ/dt) = ((Q(t+dt,p+dp) - Q(t,p+dp))/dt - (Q(t+dt,p) - Q(t,p))/dt)/dp
// dF2/dp = d/dp(dQ/dp) = (Q(t,p+dp) - 2.Q(t,p) + Q(t,p-dp))/(dp^2)
// dF2/dt = d/dt(dQ/dp) = ((Q(t+dt,p+dp) - Q(t+dt,p))/dp - (Q(t,p+dp) - Q(t,p))/dp)/dt
//
// We solve:
//   |J11 J12| |dt| = |-F1|
//   |J21 J22| |dp|   |-F2|
//-----------------------------------------------------------------------------------------------------
int joining_edge(double *e0, double *e1, int kv0, int kv1)		// was kp0, kp1
{
	double p0[3], p1[3], p[3];		// locations of P0 and P1
	double Vx[3], Vy[3], Vz[3];		// orthogonal unit vectors
	double d;	// distance from P0 -> P1
	double theta;	// theta
	double phi;	// phi
	double L;	// distance from P0 -> P, and from P -> P1
	double J11, J12, J21, J22;		// dF1/dt, dF1/dp, dF2/dt, dF2/dp
	double Qzz, Qzm, Qmz, Qzp, Qpz, Qpp;	// z = 0, m = -1, p = +1
	double F1, F2;
	double del, dt, dp, len;
	int i, k, kp, kp0, kp1, nlinks;
	int *enew;
	double epsilon = 0.00001;

	kp0 = vertex[kv0].ivox;
	kp1 = vertex[kv1].ivox;
//	p0[0] = point[kp0].x; p0[1] = point[kp0].y; p0[2] = point[kp0].z; 
//	p1[0] = point[kp1].x; p1[1] = point[kp1].y; p1[2] = point[kp1].z; 
	for (k=0; k<3; k++) {
		p0[k] = voxel[kp0].pos[k];
		p1[k] = voxel[kp1].pos[k];
	}
	makeAxes(e0,p0,p1,Vx,Vy,Vz);

	d = 0;
	for (i=0; i<3; i++) {
//		del = p0[i]-p1[i];
		del = vsize[i]*(p0[i]-p1[i]);
		d += del*del;
	}
	d = sqrt(d);
	theta = PI/2;	// initial guess of theta 
	phi = 0;		// and phi: point midway between P0 and P1
//	printf("P0: %8.4f %8.4f %8.4f P1: %8.4f %8.4f %8.4f\n",p0[0],p0[1],p0[2],p1[0],p1[1],p1[2]);
//	printf("d = %8.4f\n",d);
	del = 0.01;

	k = 0;
	for (;;) {
		k++;
		Qzz = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta,phi);
		Qmz = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta-del,phi);
		Qzm = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta,phi-del);
		Qpp = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta+del,phi+del);
		Qpz = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta+del,phi);
		Qzp = Q(e0,e1,p0,p1,Vx,Vy,Vz,d,theta,phi+del);
		F1 = (Qpz - Qzz)/del;
		F2 = (Qzp - Qzz)/del;
		/*
		printf("Qzz: %10.6f\n",Qzz);
		printf("Qmz: %10.6f\n",Qmz);
		printf("Qzm: %10.6f\n",Qzm);
		printf("Qpp: %10.6f\n",Qpp);
		printf("Qpz: %10.6f\n",Qpz);
		printf("Qzp: %10.6f\n",Qzp);
		
		printf("Qzz: %10.6f F1: %10.6f  F2: %10.6f  fabs: %10.6f %10.6f\n",Qzz,F1,F2,fabs(F1),fabs(F2));
		*/
		if (fabs(F1) < epsilon && fabs(F2) < epsilon) break;
		J11 = (Qpz - 2*Qzz + Qmz)/(del*del);
		J22 = (Qzp - 2*Qzz + Qzm)/(del*del);
		J12 = (Qpp - Qzp - Qpz + Qzz)/(del*del);
		J21 = (Qpp - Qpz - Qzp + Qzz)/(del*del);
		/*
		printf("J11: %10.6f\n",J11);
		printf("J12: %10.6f\n",J12);
		printf("J21: %10.6f\n",J21);
		printf("J22: %10.6f\n",J22);
		*/
		dt = (J12*F2 - J22*F1)/(J11*J22 - J12*J21);
		dp = (-F1 - J11*dt)/J12;
		theta += dt;
		phi += dp;
		L = d/(2*sin(theta)*cos(phi));
		if (k == 20) {
			printf("Too many iterations\n");
			exit(1);
		}
		if (L > 0.75*d || (fabs(dt) < epsilon && fabs(dp) < epsilon)) break;
	}

	getP(d,theta,phi,p0,Vx,Vy,Vz,p);

	kp = np;
	//point[kp].x = p[0];
	//point[kp].y = p[1];
	//point[kp].z = p[2];
	//point[kp].d = (point[kp0].d + point[kp1].d)/2;
	//point[kp].used = true;
	voxel[kp].pos[0] = p[0];		// size of voxel array ???????????????????????????
	voxel[kp].pos[1] = p[1];
	voxel[kp].pos[2] = p[2];
	avediameter[kp] = FIXED_DIAMETER;
//	voxel[kp].d = (voxel[kp0].d + voxel[kp1].d)/2;
//	voxel[kp].used = true;
	np++;
	if (np == np_max) {		// voxel array size exceeded
		printf("voxel array size exceeded\n");
		fprintf(fperr,"voxel array size exceeded\n");
		return 1;
	}

	edgeList[ne].npts = 3;			// size of edgeList array ????????????????????????
	edgeList[ne].pt = (int *)malloc(edgeList[ne].npts*sizeof(int));
	edgeList[ne].pt_used = (int *)malloc(edgeList[ne].npts*sizeof(int));
//	edgeList[ne].avediameter = (float *)malloc(edgeList[ne].npts*sizeof(float));
	edgeList[ne].npts_used = 3;
	edgeList[ne].pt[0] = kp0;
	edgeList[ne].pt[1] = kp;
	edgeList[ne].pt[2] = kp1;
	len = dist_um(kp0,kp) + dist_um(kp,kp1);
	if (len > 300) {
		printf("Too long! %d %d %d  %f\n",kp0,kp,kp1,len);
		printf("P0, P, P1\n");
		for (i=0; i<3; i++) {
			printf("%6.3f  %6.3f  %6.3f\n",p0[i],p[i],p1[i]);
		}
		exit(1);
	}
	edgeList[ne].pt_used[0] = kp0;
	edgeList[ne].pt_used[1] = kp;
	edgeList[ne].pt_used[2] = kp1;
//	edgeList[ne].vert[0] = kp0;
//	edgeList[ne].vert[1] = kp1;
	edgeList[ne].vert[0] = kv0;
	edgeList[ne].vert[1] = kv1;
//	for (i=0; i<3; i++)
//		edgeList[ne].avediameter[i] = FIXED_DIAMETER;
	edgeList[ne].used = true;
// Need to modify the two vertexes, kv0 and kv1	
// increment nlinks, nlinks_used
// resize edge array, add edge ne to the list

	nlinks = vertex[kv0].nlinks;
	nlinks++;
	vertex[kv0].nlinks_used = vertex[kv0].nlinks_used + 1;
	enew = (int *)malloc(nlinks*sizeof(int));
	for (i=0;i<nlinks-1;i++)
		enew[i] = vertex[kv0].edge[i];
	free(vertex[kv0].edge);
	vertex[kv0].edge = (int *)malloc(nlinks*sizeof(int));
	for (i=0;i<nlinks-1;i++)
		vertex[kv0].edge[i] = enew[i];
	free(enew);

	vertex[kv0].edge[nlinks-1] = ne;
	nlinks = vertex[kv1].nlinks;
	nlinks++;
	vertex[kv1].nlinks_used = vertex[kv1].nlinks_used + 1;
	enew = (int *)malloc(nlinks*sizeof(int));
	for (i=0;i<nlinks-1;i++)
		enew[i] = vertex[kv1].edge[i];
	free(vertex[kv1].edge);
	vertex[kv1].edge = (int *)malloc(nlinks*sizeof(int));
	for (i=0;i<nlinks-1;i++)
		vertex[kv1].edge[i] = enew[i];
	vertex[kv1].edge[nlinks-1] = ne;
	free(enew);

	ne++;
	if (ne == ne_max) {		// edgeList array size exceeded
		printf("edgeList array size exceeded\n");
		fprintf(fperr,"edgeList array size exceeded\n");
		return 1;
	}
	return 0;
}
//-----------------------------------------------------------------------------------------------------
// Try to join loose ends.  There are several criteria making up the joining score for a loose end pair:
// (1) closeness - close ends get a high score
// (2) angle - small angle made a by straight joining line at each end gets a high score
// (3) diameter match - closely matched diameters get a high score (NOT USED)
//
// Currently a single segment is used to join the two vertices.  Need to use intermediate points
// to create a curve.
//
// Needs modifying to account for the variable voxel dimensions vsize[:]
//-----------------------------------------------------------------------------------------------------
int healer(LOOSEND end[], int nloose)
{
	int i0, i1, k, kp0, kp1, njoined, imax, kv0, kv1, kv0c, kv1c, err;
	double *dir0, *dir1;
	double dave0, dave1, d, vec[3], ang0, ang1, s, smax;
	EDGE edge0, edge1;
#define MAX_END_SEPARATION 50.
#define THRESHOLD_SCORE 2./MAX_END_SEPARATION
#define MAX_DIAMETER 14

	printf("healer: %d\n",nloose);
	njoined = 0;
	for (i0=0; i0<nloose; i0++) {
		if (end[i0].joined >= 0) continue;
		edge0 = edgeList[end[i0].iedge];
//		kp0 = edge0.vert[end[i0].iend];
		kv0 = edge0.vert[end[i0].iend];
		if (end[i0].iend == 0)
			kv0c = edge0.vert[1];
		else
			kv0c = edge0.vert[0];
		kp0 = vertex[kv0].ivox;
		dir0 = end[i0].dir;
		dave0 = end[i0].dave;
		if (dave0 > MAX_DIAMETER) continue;
		imax = -1;
		smax = 0;
		for (i1=i0; i1<nloose; i1++) {
			if (i0 == i1) continue;
			if (end[i1].joined >= 0) continue;

			edge1 = edgeList[end[i1].iedge];
//			kp1 = edge1.vert[end[i1].iend];
			kv1 = edge1.vert[end[i1].iend];
			if (end[i1].iend == 0)
				kv1c = edge1.vert[1];
			else
				kv1c = edge1.vert[0];
			kp1 = vertex[kv1].ivox;
			if (kv1c == kv0c) continue;
			dir1 = end[i1].dir;
			dave1 = end[i1].dave;
			if (dave1 > MAX_DIAMETER) continue;
			d = dist_um(kp0,kp1);
			if (d > MAX_END_SEPARATION) continue;
			for (k=0; k<3; k++) 
				vec[k] = voxel[kp1].pos[k] - voxel[kp0].pos[k];
			ang0 = angle(dir0,vec);
			if (ang0 < -1.001*PI || ang0 > 1.001*PI) {
				printf("Bad ang0: %f\n",ang0);
				return 1;
			}
			for (k=0; k<3; k++)
				vec[k] = -vec[k];
			ang1 = angle(dir1,vec);
			if (ang1 < -1.001*PI || ang1 > 1.001*PI) {
				printf("Bad ang1: %f\n",ang1);
				return 1;
			}
			s = fang(ang0)*fang(ang1)*fdist(d);
			if (s < THRESHOLD_SCORE) continue;
			if (s > smax) {
				imax = i1;
				smax = s;
			}
		}
		if (imax >= 0) {
			njoined += 2;
			end[i0].joined = imax;
			end[imax].joined = i0;
			edge1 = edgeList[end[imax].iedge];
//			kp1 = edge1.vert[end[imax].iend];
			kv1 = edge1.vert[end[imax].iend];
			kp1 = vertex[kv1].ivox;
			dir1 = end[imax].dir;

			printf("Joining loose ends: new edge: %6d  %6d %6d  %6d %6d  %8.4f\n",ne,end[i0].iedge,end[imax].iedge,kv0,kv1,smax);
			fprintf(fperr,"Joining loose ends: new edge: %6d  %6d %6d  %6d %6d  %8.4f\n",ne,end[i0].iedge,end[imax].iedge,kv0,kv1,smax);
//			joining_edge(dir0,dir1,kp0,kp1);
			err = joining_edge(dir0,dir1,kv0,kv1);
			if (err != 0) 
				return err;
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// List the loose ends
//-----------------------------------------------------------------------------------------------------
void showEnds(LOOSEND end[], int nloose)
{
	int i, ie, iend,kv0, kv1, kv;

	for (i=0; i<nloose; i++) {
		ie = end[i].iedge;
		iend = end[i].iend;
		kv0 = edgeList[ie].vert[0];
		kv1 = edgeList[ie].vert[1];
		kv = edgeList[ie].vert[iend];
//		printf("end: %d  edge: %d  vert: %d %d  iend: %d  kv: %d\n",i,ie,kv0,kv1,iend,kv);
		fprintf(fperr,"end: %d  edge: %d  vert: %d %d  iend: %d  kv: %d\n",i,ie,kv0,kv1,iend,kv);
	}
}

//-----------------------------------------------------------------------------------------------------
// Replaces oldprune()
//-----------------------------------------------------------------------------------------------------
void prune() 
{
	int iv, ie, nlinks, iecon, ivcon, k, kp0, kp1, count, n, etmp[20];
	float len, dave;
	bool deadend;
	EDGE *edge;
	int *vtmp;

	printf("NO prune\n");
	return;

	printf("prune\n");
	// Check nlinks_used
	//vtmp = (int *)malloc(nv*sizeof(int));
	//for (iv=0; iv<nv; iv++) {
	//	vtmp[iv] = vertex[iv].nlinks_used;
	//}
	checkVerticies(false);	// To ensure that .nlinks = .nlinks_used are consistent with edgeList.
	//for (iv=0; iv<nv; iv++) {
	//	if (vtmp[iv] != vertex[iv].nlinks_used) {
	//		fprintf(fpout,"iv: %d vtmp: %d nlinks_used: %d\n",iv,vtmp[iv],vertex[iv].nlinks_used);
	//	}
	//}
	//exit(1);

	count = 0;
	for (iv=0; iv<nv; iv++) {
		if (!vertex[iv].used) continue;
		if (vertex[iv].nlinks_used == 1) {
			count++;
			deadend = true;
		} else {
			deadend = false;
		}
		nlinks = 0;
		for (ie=0; ie<ne; ie++) {
			if (!edgeList[ie].used) continue;
			if (edgeList[ie].vert[0] == iv) {
				nlinks++;
				if (deadend) {
					iecon = ie;
					ivcon = edgeList[ie].vert[1];
				}
			}
			if (edgeList[ie].vert[1] == iv) {
				nlinks++;
				if (deadend) {
					iecon = ie;
					ivcon = edgeList[ie].vert[0];
				}
			}
		}
		//if (nlinks != vertex[iv].nlinks_used) {
		//	printf("vertex: %d nlinks_used: %d actual: %d\n",iv,vertex[iv].nlinks_used,nlinks);
		//}
		// Now find the length if nlinks == 1
		if (deadend) {
			edge = &edgeList[iecon];
			len = 0;
			dave = 0;
			for (k=0; k<edge->npts; k++) {
				if (k > 0) {
					kp0 = edge->pt[k-1];
					kp1 = edge->pt[k];
					len += zdist(kp0,kp1);
				}
				dave += pointdiameter(edge->pt[k]);
			}
			dave/= edge->npts;
//			printf("loose end: %d %d %6.1f %6.1f %6.3f\n",iecon,iv,len,dave,len/dave);
			if (len/dave < 5) {
				//if (ivcon == 10016) {
				//	fprintf(fpout,"Remove the loose end: %d edge: %d on vertex: %d\n",iv,iecon,ivcon);
				//	showvertex(10016);
				//}
				edge->used = false;
				// remove iecon from the edge list for vertex ivcon
				//n = 0;
				//for (k=0; k<vertex[ivcon].nlinks; k++) {
				//	ie = vertex[ivcon].edge[k];
				//	if (ie != iecon) {
				//		vertex[ivcon].edge[n] = ie;
				//		n++;
				//	}
				//}
				//vertex[ivcon].nlinks--;
				vertex[ivcon].nlinks_used--;
				//if (ivcon == 10016) {
				//	fprintf(fpout,"Removed the edge: %d\n",iecon);
				//	showvertex(10016);
				//}
				// Now the non-deadend vertex may have nlinks=2, and the edges may need to be concatenated
				if (vertex[ivcon].nlinks == 2) {
					//if (ivcon == 10016) {
					//	fprintf(fpout,"Now nlinks = 2\n");
					//	showvertex(ivcon);
					//}
					joiner(ivcon);
					//if (ivcon == 10016) {
					//	fprintf(fpout,"did joiner: %d\n",ivcon);
					//	showedge(18818,'N');
					//	showedge(18824,'N');
					//}
				}
			}
		}
	}
	printf("ne: %d loose ends: %d\n",ne,count);
}

//-----------------------------------------------------------------------------------------------------
// Prune short dead-end edges
// If not USE_HEALING, iter=1 makes no changes to the network.
//-----------------------------------------------------------------------------------------------------
int oldprune(int iter)
{
	int ie, i, k, kv0, kv1, k1, k2, j1, j2, npruned, nlused[2], nloose, last, err;
	float diam, len;
	EDGE edge, edge0;
	LOOSEND *end;
	int MINTWIGFACTOR = 3;

	printf("prune: iter: %d min twig length factor: %d\n", iter, MINTWIGFACTOR);
	fprintf(fperr,"prune: iter: %d min twig length factor: %d\n", iter, MINTWIGFACTOR);
	if (iter == 0) {	// Need to count loose ends, to allocate array end[]
		nloose = 0;
		for (ie=0; ie<ne; ie++) {
			edge = edgeList[ie];
			if (!edge.used) continue;
			kv0 = edge.vert[0];
			kv1 = edge.vert[1];
			nlused[0] = vertex[kv0].nlinks_used;
			nlused[1] = vertex[kv1].nlinks_used;
			if (nlused[0] != 1 && nlused[1] != 1) continue;
			nloose++;
		}
		printf("nloose: %d\n",nloose);
		fprintf(fperr,"nloose: %d\n",nloose);
		end = (LOOSEND *)malloc(nloose*sizeof(LOOSEND));
	}
	
	npruned = 0;
	nloose = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
		if (!edge.used) continue;
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		nlused[0] = vertex[kv0].nlinks_used;
		nlused[1] = vertex[kv1].nlinks_used;
		if (nlused[0] != 1 && nlused[1] != 1) continue;
//		printf("edge: %d npts: %d kv0,kv1: %d %d nlused: %d %d\n",ie,edge.npts,kv0,kv1,nlused[0],nlused[1]);
		// One end or the other is loose
		diam = 0;
		len = 0;
		for (i=0; i<edge.npts; i++) {
			diam = diam + pointdiameter(edge.pt[i]);
			if (i > 0)
				len += dist_um(last,edge.pt[i]);	// now length in um
			last = edge.pt[i];
		}
		diam = diam/edge.npts;					// now diam in um
		if (iter == 0) {
			end[nloose].iedge = ie;
			if (nlused[0] == 1) {
				end[nloose].iend = 0;
				j2 = 0;
				j1 = 1;
			} else {
				end[nloose].iend = 1;
				j2 = edge.npts-1;
				j1 = j2 - 1;
			}
			end[nloose].dave = diam;	//dave;
			end[nloose].len = len;
			k1 = edge.pt[j1];
			k2 = edge.pt[j2];
			len = dist_um(k1,k2);
			for (k=0; k<3; k++) {
//				end[nloose].dir[k] = (voxel[k2].pos[k] - voxel[k1].pos[k])/len;
				end[nloose].dir[k] = vsize[k]*(voxel[k2].pos[k] - voxel[k1].pos[k])/len;
			}
			end[nloose].joined = -1;
			edge0 = edgeList[end[nloose].iedge];
			printf("loose edge: %d npts: %d\n",end[nloose].iedge,edge0.npts);
			fprintf(fperr,"loose edge: %d npts: %d\n",end[nloose].iedge,edge0.npts);
			nloose++;
		} else {		// Is there a way to process the list of loose ends, as in prune.cpp?
			if (len <= MINTWIGFACTOR*diam) {
				if (nlused[0] == 1) {
					vertex[kv0].used = false;
					edgeList[ie].used = false;
					vertex[kv0].nlinks_used--;
					vertex[kv1].nlinks_used--;
					npruned++;
					printf("nlused[0]=1: ie: %d kv0,kv1: %d %d\n",ie,kv0,kv1);
					fprintf(fperr,"nlused[0]=1: ie: %d kv0,kv1: %d %d\n",ie,kv0,kv1);
				}
				if (nlused[1] == 1) {
					vertex[kv1].used = false;
					edgeList[ie].used = false;
					vertex[kv0].nlinks_used--;
					vertex[kv1].nlinks_used--;
					npruned++;
					printf("nlused[1]=1: ie: %d kv0,kv1: %d %d\n",ie,kv0,kv1);
					fprintf(fperr,"nlused[1]=1: ie: %d kv0,kv1: %d %d\n",ie,kv0,kv1);
				}
			}	
		}
	}
	if (iter == 0 && USE_HEALING) {
		err = 0;
		showEnds(end,nloose);
		err = healer(end,nloose);
		free(end);
		return err;
	}
	free(end);
	printf("Number of short twigs pruned: %d\n",npruned);
	return npruned;
}

//-----------------------------------------------------------------------------------------------------
// Check that all edges are completed, i.e. no vertex appears within an edge.
//-----------------------------------------------------------------------------------------------------
int adjoinEdges(void)
{
	int *nvrefs;
	PAIR *pair;
	int ie, iv, kv0, kv1, ivmax, nloose, i1, i2, n1, n2, k1, k2, i, kp0, kp1, ntwo, nvunused, neunused;
	int temp[1000];	// should be big enough
	EDGE edge1, edge2;

	printf("adjoinEdges: ne: %d nv: %d np: %d\n",ne,nv,np);
	fprintf(fpout,"adjoinEdges\n");

	nvrefs = (int *)malloc(np*sizeof(int));
	pair = (PAIR *)malloc(np*sizeof(PAIR));
	for (;;) {
		for (iv=0; iv<nv; iv++) {
			nvrefs[iv] = 0;
			pair[iv].i1 = 0;
			pair[iv].i2 = 0;
		}
		neunused = 0;
		ivmax = 0;
		for (ie=0; ie<ne; ie++) {
			edge1 = edgeList[ie];
			if (!edge1.used) {
				neunused++;
//				printf("Unused edge: %d  %d\n",neunused,ie);
//				fprintf(fpout,"Unused edge: %d  %d\n",neunused,ie);
				continue;
			}
			kv0 = edge1.vert[0];
			kp0 = vertex[kv0].ivox;
			kv1 = edge1.vert[1];
			kp1 = vertex[kv1].ivox;
			if (kp0 >= np || kp1 >= np) {
				printf("kp >= np: %d %d %d\n",kp0,kp1,np);
				exit(1);
			}
			if (nvrefs[kv0] == 0)
				pair[kv0].i1 = ie;
			else if (nvrefs[kv0] == 1)
				pair[kv0].i2 = ie;
			nvrefs[kv0]++;
			if (nvrefs[kv1] == 0)
				pair[kv1].i1 = ie;
			else if (nvrefs[kv1] == 1)
				pair[kv1].i2 = ie;
			nvrefs[kv1]++;
			if (kv0 > ivmax) ivmax = kv0;
			if (kv1 > ivmax) ivmax = kv1;
		}
//		printf("Number of unused edges: %d\n",neunused);
		fprintf(fperr,"Number of unused edges: %d\n",neunused);
		ntwo = 0;
		nloose = 0;
		nvunused = 0;
		for (iv=0; iv<=ivmax; iv++) {
			if (nvrefs[iv] == 0) {
	//			printf("Unused: %d\n",iv);
				nvunused++;
			} else if (nvrefs[iv] == 1) {
				nloose++;
			} else if (nvrefs[iv] == 2) {
				i1 = pair[iv].i1;
				i2 = pair[iv].i2;
				if (i1 == i2) {
					printf("adjoinEdges: error: i1 = i2: %d\n",i1);
					fprintf(fperr,"adjoinEdges: error: i1 = i2: %d\n",i1);
					exit(1);
				}
				edge1 = edgeList[i1];
				edge2 = edgeList[i2];
				if (!edge1.used || !edge2.used) continue;	// already processed, adjoined
				ntwo++;
				// The two edges are edge1 and edge2
				// They will be combined into one, and replace edge1, edge2 will be deleted
				n1 = edge1.npts;
				kv0 = edge1.vert[0];
				kv1 = edge1.vert[1];
				if (kv0 == iv)
					k1 = 0;
				else
					k1 = n1-1;
//				printf("edge1: %d npts: %d vert: %d %d\n",i1,n1,kv0,kv1);
				n2 = edge2.npts;
				kv0 = edge2.vert[0];
				kv1 = edge2.vert[1];
				if (kv0 == iv)
					k2 = 0;
				else
					k2 = n2-1;
//				printf("edge2: %d npts: %d vert: %d %d\n",i2,n2,kv0,kv1);
				if (k1 == 0 && k2 == 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[n1-i-1];
					kv0 = edge1.vert[1];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[i];
					kv1 = edge2.vert[1];
				} else if (k1 == 0 && k2 != 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[n1-i-1];
					kv0 = edge1.vert[1];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[n2-i-1];
					kv1 = edge2.vert[0];
				} else if (k1 != 0 && k2 == 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[i];
					kv0 = edge1.vert[0];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[i];
					kv1 = edge2.vert[1];
				} else if (k1 != 0 && k2 != 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[i];
					kv0 = edge1.vert[0];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[n2-i-1];
					kv1 = edge2.vert[0];
				}
				edge1.npts = n1 + n2 - 1;
				edge1.npts_used = edge1.npts;
				free(edge1.pt);
				free(edge1.pt_used);
				edge1.pt = (int *)malloc(edge1.npts*sizeof(int));
				edge1.pt_used = (int *)malloc(edge1.npts*sizeof(int));
				for (i=0; i<edge1.npts; i++) {
					edge1.pt[i] = temp[i];
					edge1.pt_used[i] = temp[i];
				}
				edge1.vert[0] = kv0;
				edge1.vert[1] = kv1;
				edge1.used = true;
				edgeList[i1] = edge1;
				edgeList[i2].used = false;
//				printf("Adjoined edge: %6d to edge: %6d  npts: %3d  vert: %6d %6d\n",i2,i1,edge1.npts,kv0,kv1);
//				fprintf(fperr,"Adjoined edge: %6d to edge: %6d  npts: %3d  vert: %6d %6d\n",i2,i1,edge1.npts,kv0,kv1);
				if (kv0 == kv1) {
					printf("Loop: kv0 = kv1: %d  nlinks: %d nlinks_used: %d\n",kv0,vertex[kv0].nlinks,vertex[kv0].nlinks_used);
					// drop this edge
					edgeList[i1].used = false;
					vertex[kv0].nlinks -= 2;
					vertex[kv0].nlinks_used -= 2;
				}
			}
		}
		printf("Number of two-healings: %d\n",ntwo);
//		printf("Number of unused edges: %d\n",neunused);
//		printf("Number of unused vertices: %d\n",nvunused);
		if (ntwo == 0) break;
	}
	free(nvrefs);
	free(pair);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Note: n_prune_cycles not used now
//-----------------------------------------------------------------------------------------------------
int TraceSkeleton(int n_prune_cycles)
{
	int x, y, z, ie, i, j, k, kv, nbrs=0, nbrlist[MAXNBRS][3], maxnbrs, npts;
	int next[3];
	int nbcount[MAXNBRS];
	int nhit, k1, xnb, ynb, znb, nloops, iter, npruned, err;
	bool hit;

	printf("CountNbrs\n");
	for (i=0;i<MAXNBRS;i++)
		nbcount[i] = 0;
	npts = 0;
	maxnbrs = 0;
	for (z=0; z<depth; z++) {
//		printf("z: %d\n",z);
		for (y=0; y<height; y++) {
			for (x=0; x<width; x++) {
				if (V3Dskel(x,y,z) == 0) continue;
				npts++;
				next[0] = x; next[1] = y; next[2] = z;
				if (USE_NEW)
					nbrs = GetNeighborsNew(next,nbrlist);
				else
					nbrs = GetNeighborsOld(next,nbrlist);
//				printf("nbrs: %d\n",nbrs);
//				if (nbrs != 1) continue;	// an end point has just one neighbor
				nbcount[nbrs] = nbcount[nbrs] + 1;
				if (nbrs > maxnbrs) {
//					printf(" nbrs: %d\n",z,nbrs);
					maxnbrs = nbrs;
				}
			}
		}
	}
	printf("npts: %d  maxnbrs: %d\n",npts,maxnbrs);
//	exit(0);

	for (i=0; i<=maxnbrs; i++) {
		printf("%4d %8d\n",i,nbcount[i]);
	}
	np = npts;
	np_max = 1.1*npts;
	voxel = (VOXEL *)malloc(np_max*sizeof(VOXEL));
	avediameter = (float *)malloc(np_max*sizeof(float));
	for (i=0; i<npts; i++)
		avediameter[i] = FIXED_DIAMETER;

	// The number of vertices is the number of voxels with nbrs != 2.
	nv = npts - nbcount[2];
	printf("Number of vertices: %d\n",nv);
	vertex = (VERTEX *)malloc(nv*sizeof(VERTEX));
	k = 0;
	kv = 0;
	for (z=0; z<depth; z++) {
		for (y=0; y<height; y++) {
			for (x=0; x<width; x++) {
				if (V3Dskel(x,y,z) == 0) continue;
				next[0] = x; next[1] = y; next[2] = z;
				if (USE_NEW)
					nbrs = GetNeighborsNew(next,nbrlist);
				else
					nbrs = GetNeighborsOld(next,nbrlist);
				if (nbrs > MAXNBRS) {
					printf("Error: TraceSkeleton: nbrs > MAXNBRS: %d  %d\n",nbrs,MAXNBRS);
					fprintf(fperr,"Error: TraceSkeleton: nbrs > MAXNBRS: %d  %d\n",nbrs,MAXNBRS);
					return 1;
				}
				voxel[k].nbrs = nbrs;
				voxel[k].pos[0] = x;
				voxel[k].pos[1] = y;
				voxel[k].pos[2] = z;
				for (i=0;i<nbrs;i++) {
					for (j=0; j<3; j++) {
						voxel[k].nbr[i][j] = nbrlist[i][j];
					}
				}
				if (nbrs != 2) {	// initialise the vertex
					voxel[k].ivert = kv;
					vertex[kv].ivox = k;
					vertex[kv].pos[0] = x;
					vertex[kv].pos[1] = y;
					vertex[kv].pos[2] = z;
					vertex[kv].ivisit = -1;
					vertex[kv].nlinks = nbrs;
					vertex[kv].nlinks_used = nbrs;
					vertex[kv].followed = (bool *)malloc(nbrs*sizeof(bool));
					vertex[kv].edge = (int *)malloc(nbrs*sizeof(int));
					for (i=0; i<nbrs; i++)
						vertex[kv].followed[i] = false;
					vertex[kv].nfollowed = 0;
					vertex[kv].used = true;
					kv++;
				}
				k++;
			}
		}
	}
	printf("Made voxel list\n");

	// Now find the neighbour IDs
	for (k=0; k<npts; k++) {
		nbrs = voxel[k].nbrs;
		nhit = 0;
		for (j=0; j<nbrs; j++) {	// look at the neighbour points one by one
			xnb = voxel[k].nbr[j][0];
			ynb = voxel[k].nbr[j][1];
			znb = voxel[k].nbr[j][2];
			hit = false;
			if (k > 0) {	// look at preceding voxels in the list
				k1 = k;
				for (;;) {
					k1--;
					if (k1 < 0) break;
					if (voxel[k1].pos[2] < voxel[k].pos[2]-1) break;
					if (voxel[k1].pos[0]==xnb && voxel[k1].pos[1]==ynb && voxel[k1].pos[2]==znb) {
						voxel[k].nid[j] = k1;
						hit = true;
						nhit++;
						break;
					}
				}
			}
			if (hit) continue;
			if (k < npts-1) {	// look at succeeding voxels in the list
				k1 = k;
				for (;;) {
					k1++;
					if (k1 > npts-1) break;
					if (voxel[k1].pos[2] > voxel[k].pos[2]+1) break;
					if (voxel[k1].pos[0]==xnb && voxel[k1].pos[1]==ynb && voxel[k1].pos[2]==znb) {
						voxel[k].nid[j] = k1;
						hit = true;
						nhit++;
						break;
					}
				}
			}
		}
		if (nhit != nbrs) {
			printf("nhit != nbrs: %d  %d %d\n",k,nhit,nbrs);
			return 1;
		}
	}
	printf("Determined all neighbour IDs\n");

	// Check for disconnected fragments
	for (i=0; i<npts; i++) {
		nbrs = voxel[i].nbrs;
		if (nbrs == 1) {		// look at dead-end vertices
			j = voxel[i].nid[0];
			if (voxel[j].nbrs == 1) {
				printf("Solitary pair: %d %d  %d %d\n",i,j,voxel[i].nid[0],voxel[j].nid[0]);
			}
		}
	}

	// I would like to eliminate tight loops here
	// It should be possible to find loops using a recursive function

	tracer();
	printf("did tracer\n");
	fprintf(fpout,"did tracer\n");
	err = checker();
	if (err != 0) return 1;

	printf("Do deloop\n");
	for (iter=0; iter<3; iter++) {
		nloops = deloop(iter);
		if (nloops < 0) return 1;
		printf("did deloop: nloops: %d\n",nloops);
		fprintf(fpout,"did deloop: nloops: %d\n",nloops);
		err = checker();
		if (err != 0) return 1;
	}

	if (NEW_PRUNE) {
		prune();
		fprintf(fpout,"did prune\n");
	} else {
		for (iter=0; iter<n_prune_cycles; iter++) {
			printf("call prune: %d\n",iter);
			npruned = oldprune(iter);
			err = checker();
			if (err != 0) return 1;
			printf("iter: %d npruned: %d\n",iter,npruned);
			if (iter > 0 && npruned < 10) break;
		}
	
		adjoinEdges();
		printf("did adjoinEdges\n");
	}
	// Drop  loops created by pruning
	for (ie=0; ie<ne; ie++) {
		npts = edgeList[ie].npts;
		int kfrom = edgeList[ie].pt_used[0];
		int kto = edgeList[ie].pt_used[npts-1];
		if (kfrom == kto) edgeList[ie].used = false;
	}
	err = checkUnconnected();
	if (err != 0) return 1;
	edgePtDist();
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
int FixDiameters()
{
	EDGE *edge;
	int ie, k, kp0, kp1, kp, npts, n0, n1, ipass, nzero, iv, err;
	double d0, d1, diam, d2ave, dave, vsum, dmax, dmin;
	double len, len0, len1, lsum;
	bool done;
	double alpha = 0.3;
	bool FIX_JUNCTIONS = true;

	printf("FixDiameters\n");
	//for (ie=0; ie<100; ie++) {
	//	edge = &edgeList[ie];
	//	if (!edge->used) continue;
	//	showedge(ie,'D');
	//}

	for (ie=0; ie<ne; ie++) {
		edge = &edgeList[ie];
		if (!edge->used) continue;
		npts = edge->npts;
		kp0 = edge->pt[0];
		kp1 = edge->pt[npts-1];
//		edge->avediameter[0] = 0;
//		edge->avediameter[npts-1] = 0;
		for (k=1; k<npts-1; k++) {
			kp = edge->pt[k];
			diam = avediameter[kp];
			d0 = dist_um(kp,kp0);
			d1 = dist_um(kp,kp1);
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
//		if (ie == 4) printf("ie: %d npts: %d n0: %d d0: %f\n",ie,npts,n0,d0);
		d1 = 0;
		n1 = npts-1;
		for (k=npts-2; k>0; k--) {
			kp = edge->pt[k];
			d1 = avediameter[kp];
			if (d1 > 0) break;
			n1 = k;
		}
//		if (ie == 4) printf("ie: %d npts: %d n1: %d d1: %f\n",ie,npts,n1,d1);
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
			//for (k=0; k<npts; k++) {
			//	kp = edge->pt[k];
			//	printf("%6.1f",avediameter[kp]);
			//}
			//printf("\n");
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
				len = zdist(kp0,kp1);
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
	// Need to identify the point corresponding to the vertex!
	// This is edge->pt[0] if iv == edge->vert[0], edge->pt[npts-1] if iv = edge->vert[1]
	for (iv=0; iv<nv; iv++) {
		if (!vertex[iv].used || vertex[iv].nlinks_used == 0) continue;
		dmin = 1.0e10;
		dmax = 0;
		dave = 0;
		n1 = 0;
		for (k=0; k<vertex[iv].nlinks; k++) {
			ie = vertex[iv].edge[k];
			edge = &edgeList[ie];
			if (!edge->used) continue;
			n1++;
			npts = edge->npts;
//			if (npts > 2) {
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
//			}
		}
		dave /= n1;
//		if (dave == 0) continue;
		for (k=0; k<vertex[iv].nlinks; k++) {
			ie = vertex[iv].edge[k];
			edge = &edgeList[ie];
			if (!edge->used) continue;
			npts = edge->npts;
			if (iv == edge->vert[0]) {
				kp = edgeList[ie].pt[0];
			} else {
				kp = edgeList[ie].pt[npts-1];
			}
//			printf("iv: %6d  %3d  %6d %6d %6.1f %6.1f\n",iv,n1,k,kp,avediameter[kp],dmax);
			if (junction_max || dmax/dmin > max_ratio) {
				avediameter[kp] = dmax;
			} else {
				avediameter[kp] = dave;
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// To track down the problem with r2=0
//-----------------------------------------------------------------------------------------------------
int checkDiameter()
{
	int i, x, y, z, dx, dy, dz;
	double p1[3], p2[3], p0[3];
	double r2_ave, r2_min, diam, factor;
//	double alpha = 0.0;
	double dlim = 50.0;
	unsigned char val;
	int v0[]={1332,744,57};
	int v1[]={1330,744,57};
	int v2[]={1334,744,57};

	for (dx=-2; dx<=2; dx++) {
		for (dy=-2; dy<=2; dy++) {
			for (dz=-2; dz<=2; dz++) {
				x = v0[0]+dx;
				y = v0[1]+dy;
				z = v0[2]+dz;
				val = V3D(x,y,z);
				fprintf(fpout,"val: %d %d %d  %d\n",x,y,z,val);
			}
		}
	}
	for (i=0; i<3; i++) {
		p0[i] = vsize[i]*(v0[i] + 0.5);		// Note: centres of voxel cubes
		p1[i] = vsize[i]*(v1[i] + 0.5);
		p2[i] = vsize[i]*(v2[i] + 0.5);
	}
	// This estimates the average and minimum diameter at the point p1, centreline p0 -> p2
	EstimateDiameter(p0,p1,p2,&r2_ave,&r2_min);
	diam = 2*sqrt(r2_ave);
	if (calib_param != 0) {
		//factor = 1.0 + calib_param*diam/dlim;
		//diam = factor*diam;
		diam *= calib_param;
	}
	fprintf(fpout,"diam: %f\n",diam);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int uniform_flag, junction_max_flag;
	int n_prune_cycles;
	double voxelsize_x, voxelsize_y, voxelsize_z;
	double volume;
	char *vessFile, *skelFile;
	int count, err;
	char *outfilename;
	char drive[32], dir[256],filename[256], ext[32];
	char errfilename[256], amfilename[256];
	bool squeeze = true;

	if (argc != 13) {
		printf("Usage: topo skel_tiff object_tiff output_file voxelsize_x voxelsize_y voxelsize_z calib_param fixed_diam uniform_flag junction_max_flag max_ratio n_prune_cycles\n");
		fperr = fopen("topo_error.log","w");
		fprintf(fperr,"Usage: topo skel_tiff object_tiff output_file voxelsize_x voxelsize_y voxelsize_z calib_param fixed_diam uniform_flag junction_max_flag max_ratio n_prune_cycles\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	outfilename = argv[3];
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s%s",output_basename,"_topo.log");
	fperr = fopen(errfilename,"w");
	printf("errfilename: %s\n",errfilename);

	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
	fprintf(fperr,"Basename: %s\n",output_basename);
	sprintf(amfilename,"%s.am",output_basename);

	skelFile = argv[1];
	vessFile = argv[2];
	sscanf(argv[4],"%lf",&voxelsize_x);
	sscanf(argv[5],"%lf",&voxelsize_y);
	sscanf(argv[6],"%lf",&voxelsize_z);
	sscanf(argv[7],"%lf",&calib_param);
	sscanf(argv[8],"%f",&FIXED_DIAMETER);
	sscanf(argv[9],"%d",&uniform_flag);
	sscanf(argv[10],"%d",&junction_max_flag);
	sscanf(argv[11],"%lf",&max_ratio);
	sscanf(argv[12],"%d",&n_prune_cycles);

	vsize[0] = voxelsize_x;
	vsize[1] = voxelsize_y;
	vsize[2] = voxelsize_z;

	fpout = fopen(outfilename,"w");
	printf("Input skeleton file: %s\n",skelFile);
	fprintf(fpout,"Input skeleton file: %s\n",skelFile);
	printf("Input object file: %s\n",vessFile);
	fprintf(fpout,"Input object file: %s\n",vessFile);
	printf("Output basename: %s\n",output_basename);
	fprintf(fpout,"Output basename: %s\n",output_basename);
	printf("Voxel size: x,y,z: %f %f %f\n",voxelsize_x, voxelsize_y,voxelsize_z);
	fprintf(fpout,"Voxel size: x,y,z: %f %f %f\n",voxelsize_x, voxelsize_y,voxelsize_z);
	printf("Uniform diameter flag: %d\n",uniform_flag);
	fprintf(fperr,"Uniform diameter flag: %d\n",uniform_flag);
	printf("Junction max diameter flag: %d\n",junction_max_flag);
	fprintf(fperr,"Junction max diameter flag: %d\n",junction_max_flag);
	if (FIXED_DIAMETER > 0) {
		printf("Using a fixed vessel diameter: %f\n",FIXED_DIAMETER);
		fprintf(fpout,"Using a fixed vessel diameter: %f\n",FIXED_DIAMETER);
	}
	uniform_diameter = (uniform_flag==1);
	junction_max = (junction_max_flag==1);
	fprintf(fpout,"Amira file: %s\n",amfilename);

	if (USE_HEALING) {
		fprintf(fpout, "Using healing\n");
	} else {
		fprintf(fpout, "NOT using healing\n");
	}

	typedef itk::ImageFileReader<ImageType> FileReaderType;

	// Read skeleton file
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

	width = imskel->GetLargestPossibleRegion().GetSize()[0];
	height = imskel->GetLargestPossibleRegion().GetSize()[1];
	depth = imskel->GetLargestPossibleRegion().GetSize()[2];
	imsize_xy = width*height;

	printf("Image dimensions: width, height, depth: %lld %lld %lld\n",width,height,depth);
	pskel = (unsigned char *)(imskel->GetBufferPointer());
	int nt = 0;
	for (long long i=0; i<width*height*depth; i++) {
		if (pskel[i] > 0) nt++;
	}

//	if (FIXED_DIAMETER == 0) {
	
	// Read original vessel file (binary)
	FileReaderType::Pointer reader = FileReaderType::New();
	reader->SetFileName(vessFile);
	try
	{
		printf("Reading input vessel file: %s\n",vessFile);
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fperr,"Read error on vessel file\n");
		fclose(fperr);
		return 3;	// Read error on input file
	}

	im = reader->GetOutput();

	int w = im->GetLargestPossibleRegion().GetSize()[0];
	int h = im->GetLargestPossibleRegion().GetSize()[1];
	int d = im->GetLargestPossibleRegion().GetSize()[2];

	if (w != width || h != height || d != depth) {
		printf("Error: skeleton and vessel files differ in size\n");
		fprintf(fperr,"Error: skeleton and vessel files differ in size\n");
		fclose(fperr);
		return 4;
	}

	p = (unsigned char *)(im->GetBufferPointer());

	count = 0;
	for (long long i=0; i<width*height*depth; i++) {
		if (p[i] > 0) count++;
	}
	//for (int x=0; x<width; x++) {
	//	for (int y=0; y<height; y++) {
	//		for (int z=0; z<depth; z++) {
	//			if (V3D(x,y,z) > 0) {
	//				count++;
	//				if (z == 57) fprintf(fpout,"%d %d %d\n",x,y,z);
	//				printf("%d %d %d\n",x,y,z);
	//			}
	//		}
	//	}
	//}
	volume = count*vsize[0]*vsize[1]*vsize[2];

	fprintf(fpout,"lit voxel count: %d\n",count);

	/*
	checkDiameter();

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName("zzz.tif");
	writer->SetInput(im);
	writer->UseCompressionOn();
	try
	{
		printf("Writing vessel file: %s\n","zzz.tif");
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fperr,"Write error on vessel file\n");
		fclose(fperr);
		return 3;	// Read error on input file
	}
	return 1;
	*/

	InitVector();

	err = TraceSkeleton(1);
	if (err != 0) {
		printf("Error: TraceSkeleton\n");
		fprintf(fperr,"Error: TraceSkeleton\n");
		fclose(fperr);
		return 5;
	}
	printf("did TraceSkeleton\n");
	fprintf(fpout,"did TraceSkeleton\n");

	checkVerticies(true);

	err = simplify();
	if (err != 0) {
		printf("Error: simplify\n");
		fprintf(fperr,"Error: simplify\n");
		fclose(fperr);
		return 7;
	}
	printf("did simplify\n");
	fprintf(fpout,"did simplify\n");

	// Prune again
	for (int i=0; i<n_prune_cycles; i++) {
		prune();
		checkVerticies(true);
		err = checkEdgeEndPts();
		if (err != 0) {
			printf("Error: checkEdgeEndPts\n");
			fprintf(fperr,"Error: checkEdgeEndPts\n");
			fclose(fperr);
			return 13;
		}
	}

	if (FIXED_DIAMETER == 0) {
		err = GetDiameters();
		if (err != 0) {
			printf("Error: GetDiameters\n");
			fprintf(fperr,"Error: GetDiameters\n");
			fclose(fperr);
			return 6;
		}
		err = FixDiameters();
		if (err != 0) {
			printf("Error: FixDiameters\n");
			fprintf(fperr,"Error: FixDiameters\n");
			fclose(fperr);
			return 12;
		}
	}

	printf("Total voxels: %d edges: %d vertices: %d points: %d\n",count,ne,nv,np);
	printf("Voxel size: %6.2f %6.2f %6.2f\n",vsize[0],vsize[1],vsize[2]);
	printf("Total voxel volume: %10.0f\n",volume);
	fprintf(fpout,"Total voxels: %d edges: %d vertices: %d points: %d\n",count,ne,nv,np);
	fprintf(fpout,"Voxel size: %6.2f %6.2f %6.2f\n",vsize[0],vsize[1],vsize[2]);
	fprintf(fpout,"Total voxel volume: %10.0f\n",volume);

//	checker();
	err = checkEdgeEndPts();
	if (err != 0) {
		printf("Error: checkEdgeEndPts\n");
		fprintf(fperr,"Error: checkEdgeEndPts\n");
		fclose(fperr);
		return 13;
	}
//	checkVerticies(true);

	if (squeeze) {
		err = squeezer();
		if (err != 0) {
			printf("Error: squeezer\n");
			fprintf(fperr,"Error: squeezer\n");
			fclose(fperr);
			return 8;
		} else {
			printf("squeezed\n");
		}
	}
	checker();
	err = checkEdgeEndPts();
	if (err != 0) {
		printf("Error: checkEdgeEndPts\n");
		fprintf(fperr,"Error: checkEdgeEndPts\n");
		fclose(fperr);
		return 13;
	}
	err = CreateDistributions();		// scaling for voxelsize now done in the distance calculations
	if (err != 0) {
		printf("Error: CreateDistributions\n");
		fprintf(fperr,"Error: CreateDistributions\n");
		fclose(fperr);
		return 9;
	}

	if (squeeze) {
		err = WriteCmguiData();
		if (err != 0) {
			printf("Error: WriteCmguiData\n");
			fprintf(fperr,"Error: WriteCmguiData\n");
			fclose(fperr);
			return 10;
		}

		err = WriteAmiraFile(amfilename,vessFile,skelFile);
		if (err != 0) {
			printf("Error: WriteAmiraFile\n");
			fprintf(fperr,"Error: WriteAmiraFile\n");
			fclose(fperr);
			return 11;
		}
	} else {
		printf("Amira and CMGUI files not written - data not squeezed\n");
		fprintf(fperr,"Amira and CMGUI files not written - data not squeezed\n");
	}
	printf("Terminated normally\n");
	fprintf(fperr,"Terminated normally\n");
	fclose(fperr);
	fclose(fpout);
	return 0;
}