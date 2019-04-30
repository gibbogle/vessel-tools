// prune.cpp

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

//#include "prune.h"

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


#include "network.h"

FILE *fperr, *fpout;

float ddiam, dlen, dminimum;
bool use_len_limit, use_len_diam_limit;
float len_limit, len_diam_limit;

double ratio_limit;
bool use_ratio;
bool use_fixed_diam;
float fixed_diam;

#define JOIN_LOOSE_ENDS false
#define DROP_UNCONNECTED false

int CheckNetwork(NETWORK *net, char *str);		// temporary here

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
#define crossProduct(a,b,c) \
	(a)[0] = (b)[1] * (c)[2] - (c)[1] * (b)[2]; \
	(a)[1] = (b)[2] * (c)[0] - (c)[2] * (b)[0]; \
	(a)[2] = (b)[0] * (c)[1] - (c)[0] * (b)[1];


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float dist(NETWORK *net, int k1, int k2)
{
	float dx = net->point[k2].x - net->point[k1].x;
	float dy = net->point[k2].y - net->point[k1].y;
	float dz = net->point[k2].z - net->point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
void showEdges(NETWORK *net)
{
	int i, k, j;
	EDGE edge;

	fprintf(fperr,"Edge points:\n");
	for (i=0; i<100; i++) {
		edge = net->edgeList[i];
		fprintf(fperr,"edge: %6d %6d\n",i,edge.npts);
		for (k=0; k<edge.npts; k++) {
			j = edge.pt[k];
			fprintf(fperr,"%6d %6d  %6.1f %6.1f %6.1f\n",k,j,net->point[j].x,net->point[j].y,net->point[j].z);
		}
	}
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
int showEdge(NETWORK *net, int ie)
{
	int i, kv0, kv1;
	EDGE edge;

	printf("Edge: %6d Used: %d Points: ",ie,net->edgeList[ie].used);
	for (i=0; i<net->edgeList[ie].npts; i++) 
		printf("%6d ",net->edgeList[ie].pt[i]);
	printf("\n");
	kv0 = net->edgeList[ie].vert[0];
	kv1 = net->edgeList[ie].vert[1];
	if (kv0 != net->edgeList[ie].pt[0]) {
		printf("kv0 not pt[0]: %d %d\n",kv0,net->edgeList[ie].pt[0]);
		exit(1);
	}
	if (kv1 != net->edgeList[ie].pt[net->edgeList[ie].npts-1]) {
		printf("kv1 not pt[n-1]: %d %d\n",kv1,net->edgeList[ie].pt[net->edgeList[ie].npts-1]);
		exit(1);
	}
	for (i=0;i<net->ne;i++) {
		if (i == ie) continue;
		edge = net->edgeList[i];
		if (!edge.used) continue;
		if (edge.vert[0] == kv0 || edge.vert[1] == kv0)
			printf("  vert[0]: %d connected to edge: %d\n",kv0,i);
		if (edge.vert[0] == kv1 || edge.vert[1] == kv1)
			printf("  vert[1]: %d connected to edge: %d\n",kv1,i);
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool samevalue(float v1, float v2)
{
	float tol = 0.01;

	return (fabs(v1-v2) < tol);
}

//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
bool samepoint(POINT *p1, POINT *p2)
{
	if (!samevalue(p1->x,p2->x)) return false;
	if (!samevalue(p1->y,p2->y)) return false;
	if (!samevalue(p1->z,p2->z)) return false;
	return true;
}

//-----------------------------------------------------------------------------------------------------
// Check consistency between edge vert[] and pt[]
//-----------------------------------------------------------------------------------------------------
int checkEdgeEndPts(NETWORK *net)
{
	EDGE *edge;
	POINT *pv0, *pv1, *pt0, *pt1;
	int ie, npts, k, kp0, kp1, kv0, kv1, kvp0, kvp1, err;

	err = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = &net->edgeList[ie];
		npts = edge->npts;
		pt0 = &net->point[edge->pt[0]];
		pt1 = &net->point[edge->pt[npts-1]];
		kv0 = edge->vert[0];
		kv1 = edge->vert[1];
		pv0 = &net->vertex[kv0].point;
		pv1 = &net->vertex[kv1].point;
//		if ((kp0==kvp0 && kp1==kvp1) || (kp0==kvp1 && kp1==kvp0)) continue;
		if ((samepoint(pt0,pv0) && samepoint(pt1,pv1)) || (samepoint(pt0,pv1) && samepoint(pt1,pv0))) continue;
		err = 1;
		fprintf(fpout,"checkEdgeEndPts: edge: %d vert: %d %d -> %d %d pt: %d %d\n",ie,kv0,kv1,kvp0,kvp1,kp0,kp1);
		fflush(fpout);
	}
	return err;
}

//-----------------------------------------------------------------------------------------------------
// Read Amira SpatialGraph file
// Old version modified to use *net
//-----------------------------------------------------------------------------------------------------
int ReadAmiraFile_old(char *amFile, NETWORK *net)
{
	int i, j, k, kp, npts, nee, npp, err;
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
			nee = 2*net->ne;
			k++;
		}
		if (strncmp(line,"define POINT",12) == 0) {
			sscanf(line+12,"%d",&net->np);
			npp = 2*net->np;
			k++;
		}
	}

	net->vertex = (VERTEX *)malloc(net->nv*sizeof(VERTEX));
	net->edgeList = (EDGE *)malloc(nee*sizeof(EDGE));	// 2* for added joining edges
	net->point = (POINT *)malloc(npp*sizeof(POINT));		// 2* for added joining edges

	// Initialize
	for (i=0; i<nee; i++) {
		net->edgeList[i].used = false;
	}
	for (i=0; i<npp; i++) {
		net->point[i].used = false;
	}

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
					sscanf(line,"%f %f %f\n",&net->vertex[i].point.x,&net->vertex[i].point.y,&net->vertex[i].point.z);
					kp = i;
					net->vertex[i].point.d = 0;
//					vertex[i].point.used = true;
					net->point[kp] = net->vertex[i].point;
				}
				kp++;
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
					float len = 0;
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
							len = len + dist(net,net->edgeList[i].pt[k-1],net->edgeList[i].pt[k]);
						}
					}
					net->edgeList[i].length_um = len;
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
						if (use_fixed_diam) {
							net->point[j].d = fixed_diam;
						} else {
							sscanf(line,"%f",&net->point[j].d);
							if (net->point[j].d == 0) {
								printf("Error: ReadAmiraFile: zero diameter: i: %d npts: %d k: %d j: %d\n",i,edge.npts,k,j);
								return 1;
							}
							if (j < net->nv) {		// because the first nv points are vertices
								net->vertex[j].point.d = net->point[j].d;
							}
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
	if (use_fixed_diam) return 0;

	err = CheckNetwork(net, "ReadAmiraFile");

	return err;
}


//-----------------------------------------------------------------------------------------------------
// Check that all edges are completed, i.e. no vertex appears within an edge.
// This assumes that edge.vert[0] = edge.pt[0] and edge.vert[1] = edge.pt[npts-1]
// This equality was guaranteed by the previous readAmiraFile(), but not by the revised (simple) version.
//-----------------------------------------------------------------------------------------------------
int adjoinEdges(NETWORK *net)
{
	int *nvrefs;
	PAIR *pair;
	int ie, iv, kv0, kv1, ivmax, nloose, i1, i2, n1, n2, k1, k2, i, ntwo, nvunused, neunused;
	int temp[1000];	// should be big enough
	EDGE edge1, edge2;

	printf("adjoinEdges: nv: %d\n",net->nv);
	nvrefs = (int *)malloc(10*net->nv*sizeof(int));
	pair = (PAIR *)malloc(10*net->nv*sizeof(PAIR));
	int cnt = 0;
	for (;;) {
		cnt++;
		for (iv=0; iv<net->nv; iv++) {
			nvrefs[iv] = 0;
			pair[iv].i1 = 0;
			pair[iv].i2 = 0;
		}
		neunused = 0;
		ivmax = 0;
		for (ie=0; ie<net->ne; ie++) {
			edge1 = net->edgeList[ie];
			if (!edge1.used) {
				neunused++;
//				printf("Unused edge: %d  %d\n",neunused,ie);
				continue;
			}
			n1 = edge1.npts;
			kv0 = edge1.pt[0];
			kv1 = edge1.pt[n1-1];
//			printf("edge: %d %d %d %d\n",ie,n1,kv0,kv1);
//			fflush(stdout);
			if (kv0 >= 10*net->nv) {
				printf("kv0 too big: %d %d\n",kv0,net->nv);
				fflush(stdout);
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
		ntwo = 0;
		nloose = 0;
		nvunused = 0;
		for (iv=0; iv<=ivmax; iv++) {
			if (nvrefs[iv] == 0) {
//				printf("Unused: %d\n",iv);
				nvunused++;
			} else if (nvrefs[iv] == 1) {
				nloose++;
			} else if (nvrefs[iv] == 2) {
				i1 = pair[iv].i1;
				i2 = pair[iv].i2;
				edge1 = net->edgeList[i1];
				edge2 = net->edgeList[i2];
				if (!edge1.used || !edge2.used) continue;	// already processed, adjoined
				ntwo++;
				// The two edges are edge1 and edge2
				// They will be combined into one, and replace edge1, edge2 will be deleted
				n1 = edge1.npts;
				kv0 = edge1.pt[0];
				kv1 = edge1.pt[n1-1];
				if (kv0 == iv)
					k1 = 0;
				else
					k1 = n1-1;
				n2 = edge2.npts;
				kv0 = edge2.pt[0];
				kv1 = edge2.pt[n2-1];
	//			fprintf(fpout,"Two occurrences: %d  %d %d\n",iv,i1,i2);
	//			if ((kv0 == edge1.pt[0] && kv1 == edge1.pt[n1-1]) ||
	//				(kv1 == edge1.pt[0] && kv0 == edge1.pt[n1-1])) {	// this is a loop
	//				printf("Duplicated edge: %d\n",i2);
	//				edgeList[i2].used = false;
	//				continue;
	//			}
				if (kv0 == iv)
					k2 = 0;
				else
					k2 = n2-1;
				if (k1 == 0 && k2 == 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[n1-i-1];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[i];
				} else if (k1 == 0 && k2 != 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[n1-i-1];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[n2-i-1];
				} else if (k1 != 0 && k2 == 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[i];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[i];
				} else if (k1 != 0 && k2 != 0) {
					for (i=0; i<n1; i++)
						temp[i] = edge1.pt[i];
					for (i=1; i<n2; i++)
						temp[n1+i-1] = edge2.pt[n2-i-1];
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
				edge1.vert[0] = temp[0];
				edge1.vert[1] = temp[edge1.npts-1];
				edge1.used = true;
				net->edgeList[i1] = edge1;
				net->edgeList[i2].used = false;
//				fprintf(fpout,"Adjoined edge: %d to edge: %d\n",i2,i1);
			}
		}
		printf("Number of two-healings: %d\n",ntwo);
		printf("Number of unused edges: %d\n",neunused);
		printf("Number of unused vertices: %d\n",nvunused);
		if (ntwo == 0) break;
	}
	free(nvrefs);
	free(pair);
	printf("did adjoinEdges\n");
	fflush(stdout);
	return 0;
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
int joining_edge( NETWORK *net, double *e0, double *e1, int kp0, int kp1)
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
	int i, k, kp;
	double epsilon = 0.00001;

	p0[0] = net->point[kp0].x; p0[1] = net->point[kp0].y; p0[2] = net->point[kp0].z; 
	p1[0] = net->point[kp1].x; p1[1] = net->point[kp1].y; p1[2] = net->point[kp1].z; 
	makeAxes(e0,p0,p1,Vx,Vy,Vz);

	d = 0;
	for (i=0; i<3; i++) {
		del = p0[i]-p1[i];
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

	kp = net->np;
	net->point[kp].x = p[0];
	net->point[kp].y = p[1];
	net->point[kp].z = p[2];
	net->point[kp].d = (net->point[kp0].d + net->point[kp1].d)/2;
	net->point[kp].used = true;
	net->np++;
	net->edgeList[net->ne].npts = 3;
	net->edgeList[net->ne].pt = (int *)malloc(net->edgeList[net->ne].npts*sizeof(int));
	net->edgeList[net->ne].pt_used = (int *)malloc(net->edgeList[net->ne].npts*sizeof(int));
	net->edgeList[net->ne].npts_used = 3;
	net->edgeList[net->ne].pt[0] = kp0;
	net->edgeList[net->ne].pt[1] = kp;
	net->edgeList[net->ne].pt[2] = kp1;
	len = dist(net,kp0,kp) + dist(net,kp,kp1);
	if (len > 300) {
		printf("Too long! %d %d %d  %f\n",kp0,kp,kp1,len);
		printf("P0, P, P1\n");
		for (i=0; i<3; i++) {
			printf("%6.3f  %6.3f  %6.3f\n",p0[i],p[i],p1[i]);
		}
		exit(1);
	}
	net->edgeList[net->ne].pt_used[0] = kp0;
	net->edgeList[net->ne].pt_used[1] = kp;
	net->edgeList[net->ne].pt_used[2] = kp1;
	net->edgeList[net->ne].vert[0] = kp0;
	net->edgeList[net->ne].vert[1] = kp1;
	net->edgeList[net->ne].used = true;
	net->ne++;
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
//-----------------------------------------------------------------------------------------------------
int joinLooseEnds(NETWORK *net, LOOSEND end[], int nloose)
{
	int i0, i1, k, kp0, kp1, njoined, imax;
	double *dir0, *dir1;
	double dave0, dave1, d, vec[3], ang0, ang1, s, smax;
	EDGE edge0, edge1;
#define MAX_END_SEPARATION 100.
#define THRESHOLD_SCORE 1./MAX_END_SEPARATION
#define MAX_DIAMETER 14

	njoined = 0;
	for (i0=0; i0<nloose; i0++) {
		if (end[i0].joined >= 0) continue;
		edge0 = net->edgeList[end[i0].iedge];
		kp0 = edge0.vert[end[i0].iend];
		dir0 = end[i0].dir;
		dave0 = end[i0].dave;
		if (dave0 > MAX_DIAMETER) continue;
		imax = -1;
		smax = 0;
		for (i1=i0; i1<nloose; i1++) {
			if (i0 == i1) continue;
			if (end[i1].joined >= 0) continue;

			edge1 = net->edgeList[end[i1].iedge];
			kp1 = edge1.vert[end[i1].iend];
			dir1 = end[i1].dir;
			dave1 = end[i1].dave;
			if (dave1 > MAX_DIAMETER) continue;
			d = dist(net,kp0,kp1);
			if (d > MAX_END_SEPARATION) continue;
			vec[0] = net->point[kp1].x - net->point[kp0].x;
			vec[1] = net->point[kp1].y - net->point[kp0].y;
			vec[2] = net->point[kp1].z - net->point[kp0].z;
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
//			printf("Joined: %6d  %6d %6d  %6d %6d  %8.4f\n",njoined,i0,imax,kp0,kp1,smax);
			edge1 = net->edgeList[end[imax].iedge];
			kp1 = edge1.vert[end[imax].iend];
			dir1 = end[imax].dir;

			joining_edge(net,dir0,dir1,kp0,kp1);
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// To remove obviously unwanted dead-end twigs.
// A twig is the segment from a terminus to a vertex (junction).
// The first criterion is that short fat twigs should be removed, since being fat such a twig
// is very unlikely to have lost its connection because of insufficient label intensity, and
// is much more likely to be an artifact of the thinning algorithm.
// One problem is that we do not want to remove the end of a main artery or vein, which
// may branch soon after the terminus.
//
// A secondary criterion, for thin vessels, is that there is no apparent matching twig nearby,
// i.e. another thin twig within a feasible distance, and possibly with feasible orientation.
//
// I don't believe this is finding all the loose ends.
//-----------------------------------------------------------------------------------------------------
int pruner(NETWORK *net, int iter)
{
	int i, j, k, ii, kv0, kv1, n0, n1, j1, j2, k1, k2, nloose, npruned, kp0, err;
	float len, dave, dsum;
	EDGE edge, edge0;
	LOOSEND *end;
	int NLOOSEMAX = 100000;

	printf("pruner: ne: %d\n",net->ne);
	fflush(stdout);

	nloose = 0;
	npruned = 0;
	end = (LOOSEND *)malloc(NLOOSEMAX*sizeof(LOOSEND));
	// Step through edges, looking for loose ends.
	for (i=0; i<net->ne; i++) {
		if (!net->edgeList[i].used) continue;
		edge = net->edgeList[i];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		n0 = 0;
		n1 = 0;
		for (ii=0; ii<net->ne; ii++) {
			if (!net->edgeList[ii].used) continue;
			if (i == ii) continue;
			if (net->edgeList[ii].vert[0] == kv0 || net->edgeList[ii].vert[1] == kv0) n0=1;
			if (net->edgeList[ii].vert[0] == kv1 || net->edgeList[ii].vert[1] == kv1) n1=1;
		}
		if (n0+n1 == 0 && DROP_UNCONNECTED) {
			printf("Error: pruner: edge: %d is unconnected\n",i);
			fflush(stdout);
			fprintf(fperr,"Error: pruner: edge: %d is unconnected\n",i);
			net->edgeList[i].used = false;
//			return 1;
		} else if ((n0+n1 == 1) || (n0+n1 == 0 && !DROP_UNCONNECTED)) {	// This is a loose end
			// Need the length and average diameter
			len = 0;
			dsum = 0;
			for (k=0; k<edge.npts; k++) {
				j = edge.pt[k];
				dsum += net->point[j].d;
				if (k > 0) {
					len += dist(net,j,edge.pt[k-1]);
				}
			}
			dave = dsum/edge.npts;
			end[nloose].iedge = i;
			if (n0 == 0) {
				end[nloose].iend = 0;
				j2 = 0;
				j1 = 1;
			} else {
				end[nloose].iend = 1;
				j2 = edge.npts-1;
				j1 = j2 - 1;
			}
			end[nloose].dave = dave;
			end[nloose].len = len;
			k1 = edge.pt[j1];
			k2 = edge.pt[j2];
			len = dist(net,k1,k2);
			end[nloose].dir[0] = (net->point[k2].x - net->point[k1].x)/len;
			end[nloose].dir[1] = (net->point[k2].y - net->point[k1].y)/len;
			end[nloose].dir[2] = (net->point[k2].z - net->point[k1].z)/len;
			end[nloose].joined = -1;
			edge0 = net->edgeList[end[nloose].iedge];
			kp0 = edge0.vert[end[nloose].iend];
//			fprintf(fpout,"%6d  %6d %2d %2d %2d %6d %6d %6.3f %6.3f %6.3f  %6d\n",nloose,i,end[nloose].iend,j1,j2,k1,k2,
//				end[nloose].dir[0],end[nloose].dir[1],end[nloose].dir[2],kp0);
			nloose++;
		}
	}
	if (iter == 0 && JOIN_LOOSE_ENDS) {
		err = joinLooseEnds(net,end,nloose);
		if (err != 0) {
			free(end);
			return err;
		}
	}

	npruned = 0;
	for (i=0; i<nloose; i++) {
		if (end[i].joined < 0) {
			if (use_ratio) {
				if (end[i].len/end[i].dave <= ratio_limit) {	// prune this twig 
					net->edgeList[end[i].iedge].used = false;
					npruned++;
				}
			} else {
				if (end[i].len <= ratio_limit) {	// prune this twig 
					net->edgeList[end[i].iedge].used = false;
					npruned++;
				}
			}
		}
	}
	printf("Number of loose ends: %d, number pruned: %d\n",nloose,npruned);
	fprintf(fperr,"Number of loose ends: %d, number pruned: %d\n",nloose,npruned);
	free(end);
	return 0;
}

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
// A loop is broken if its edges are still used.
// Each vertex must be checked to see if it has external connections.  There are three cases:
// (1) One vertex is connected externally (this is a dead-end).
// (2) Two vertices are connected externally (side loop)
// (3) Three vertices are connected externally (3-way junction loop)
//-----------------------------------------------------------------------------------------------------
int fixloop(NETWORK *net, int j1, int j2, int j3)
{
	int e[3];
	int vert[3];
	bool used[3];
	int nconn[3];
	int i, j, kv0, kv1, ncmin, nc2, imax;
	double d, dmax;

	e[0] = j1;
	e[1] = j2;
	e[2] = j3;
	for (i=0; i<3; i++) {
		used[i] = net->edgeList[e[i]].used;
		nconn[i] = 0;
	}
	vert[0] = net->edgeList[e[0]].vert[0];
	vert[1] = net->edgeList[e[0]].vert[1];
	if (net->edgeList[e[1]].vert[0] == vert[0] || net->edgeList[e[1]].vert[0] == vert[1]) {
		vert[2] = net->edgeList[e[1]].vert[1];
	} else {
		vert[2] = net->edgeList[e[1]].vert[0];
	}
	for (i=0; i<net->ne; i++) {
		kv0 = net->edgeList[i].vert[0];
		kv1 = net->edgeList[i].vert[1];
		for (j=0; j<3; j++) {
			if (vert[j] == kv0 || vert[j] == kv1) {
				nconn[j] += 1;
			}
		}
	}
//	printf("fixloop: %6d %6d %6d  used: %2d %2d %2d  nconn: %2d %2d %2d\n",
//		j1,j2,j3,used[0],used[1],used[2],nconn[0],nconn[1],nconn[2]);

	nc2 = 0;
	ncmin = 999;
	for (i=0; i<3; i++) {
		if (nconn[i] < ncmin) ncmin = nconn[i];
		if (nconn[i] == 2) nc2++;
	}
	if (nc2 > 0) {
		for (i=0; i<3; i++) {
			kv0 = net->edgeList[e[i]].vert[0];
			kv1 = net->edgeList[e[i]].vert[1];
			for (j=0; j<3; j++) {
				if (nconn[j] == 2 && (vert[j] == kv0 || vert[j] == kv1)) {
					net->edgeList[e[i]].used = false;
					break;
				}
			}
		}
	} else if (nc2 == 0) {	// DO NOT arbitrarily choose edge j1, choose longest
		dmax = 0;
		for (i=0; i<3; i++) {
			kv0 = net->edgeList[e[i]].vert[0];
			kv1 = net->edgeList[e[i]].vert[1];
			d = dist(net,kv0,kv1);
			if (d > dmax) {
				imax = i;
				dmax = d;
			}
		}
		net->edgeList[e[imax]].used = false;
	}

	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Search for loops
// deloop works fine.  Should remove longest side when all three have nconn = 3
//-----------------------------------------------------------------------------------------------------
int deloop(NETWORK *net)
{
#define NE2MAX 100000
	int i, ii, kv0, kv1, kkv0, kkv1;
	int npairs, ne2, j1, j2, j3, nloops, k1, k3;
	EDGE edge, eedge;
	PAIR *pair;
	int *e2;
	bool dup;

	printf("deloop\n");
	fflush(stdout);
	pair = (PAIR *)malloc(NE2MAX*sizeof(PAIR));
	e2 = (int *)malloc(NE2MAX*sizeof(int));
// Does any non-vertex point occur in more than one edge?  No.
	
	ne2 = 0;
	for (i=0; i<net->ne; i++) {
		edge = net->edgeList[i];
		if (!edge.used) continue;
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		if (edge.npts == 2) {
			if (ne2 == NE2MAX) {
				printf("Array dimension e2 exceeded\n");
				free(pair);
				free(e2);
				return 1;
			}
//			printf("2-edge: %4d %6d  %6d %6d\n",ne2,i,kv0,kv1);
//			fflush(stdout);
			e2[ne2] = i;
			ne2++;
			dup = false;
			for (ii=i+1; ii<net->ne; ii++) {
				eedge = net->edgeList[ii];
				if (!eedge.used) continue;
				kkv0 = eedge.vert[0];
				kkv1 = eedge.vert[1];
				// This finds a number of edges
				if ((kkv0 == kv0 && kkv1 == kv1) || (kkv0 == kv1 && kkv1 == kv0)) {
//					printf("Duplicate edges: %6d  %6d  %4d %4d\n",i,ii,edge.npts,eedge.npts);
//					fprintf(fpout,"%6d  %6d  %4d %4d\n",i,ii,edge.npts,eedge.npts);
					dup = true;
					if (eedge.npts == 2) {
						net->point[kkv0].d = 1.414*net->point[kkv0].d;
						net->point[kkv1].d = 1.414*net->point[kkv1].d;
					}
					break;
				}
			}
			if (dup) {
				ne2--;
				net->edgeList[i].used = false;
			}
		}
	}
//	printf("Number of 2-edges (edges with two points): %d\n",ne2);
//	fflush(stdout);

	npairs = 0;
	for (j1=0; j1<ne2; j1++) {
		edge = net->edgeList[e2[j1]];
		kv0 = edge.vert[0];
		kv1 = edge.vert[1];
		for (j2=j1+1; j2<ne2; j2++) {
			eedge = net->edgeList[e2[j2]];
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
		}
	}
//	printf("Number of 2-edge pairs (connections): %d\n",npairs);
//	fflush(stdout);
	
	bool hit;
	nloops = 0;
	for (i=0; i<npairs; i++) {
		for (ii=i+1; ii<npairs; ii++) {
			if (ii > npairs-1) continue;
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
					k1 = net->edgeList[j1].vert[0];
				else
					k1 = net->edgeList[-j1].vert[1];
				if (j3 > 0)
					k3 = net->edgeList[j3].vert[1];
				else
					k3 = net->edgeList[-j3].vert[0];
				if (k1 == k3) {
					nloops++;
					fixloop(net,abs(j1),abs(j2),abs(j3));
				}
			}
		}
	}
	free(pair);
	free(e2);
	printf("Number of loops: %d\n",nloops);
	fflush(stdout);
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Reorder element and node identifiers, net0 -> net1
//-----------------------------------------------------------------------------------------------------
int squeezer(NETWORK *net0, NETWORK *net1)
{
	int i, j, k;
	int ne_x, nv_x, np_x, knew, i_x;
	int ne, nv, np;
	EDGE edge;
	EDGE *edgeList_x;
	VERTEX *vertex_x;
	POINT *point_x;
	int *oldpt;

	printf("squeezer\n");
	ne = net0->ne;
	nv = net0->nv;
	np = net0->np;
	vertex_x = (VERTEX *)malloc(nv*sizeof(VERTEX));
	edgeList_x = (EDGE *)malloc(2*ne*sizeof(EDGE));		// 2* for added joining edges
	point_x = (POINT *)malloc(2*np*sizeof(POINT));
	net1->vertex = vertex_x;
	net1->edgeList = edgeList_x;
	net1->point = point_x;
	oldpt = (int *)malloc(2*np*sizeof(int));

	ne_x = 0;
	nv_x = 0;
	np_x = 0;
	// First process the vertices
	for (i=0; i<ne; i++) {
		if (!net0->edgeList[i].used) continue;
		edge = edgeList_x[ne_x];
		edge.used = true;
		for (j=0; j<2; j++) {
			knew = inlist(oldpt,np_x,net0->edgeList[i].vert[j]);
			if (knew == 0) {
				oldpt[np_x] = net0->edgeList[i].vert[j];
				vertex_x[np_x].point = net0->point[oldpt[np_x]];
				point_x[np_x] = net0->point[oldpt[np_x]];
				knew = np_x;
				np_x++;
			}
			edge.vert[j] = knew;
		}
		edgeList_x[ne_x] = edge;
		ne_x++;
	}
	nv_x = np_x;

	// Now add the edge interior points to the list.  These occur once only.
	i_x = 0;
	for (i=0; i<ne; i++) {
//		printf("edge: %d used: %d npts:%d\n",i,edgeList[i].used,edgeList[i].npts);
		if (!net0->edgeList[i].used) continue;
		int npts = net0->edgeList[i].npts;
		if (npts < 1) {
			printf("squeezer: i: %d npts: %d\n",i,npts);
			return 1;
		}
		edgeList_x[i_x].pt = (int *)malloc(npts*sizeof(int));
		edge = edgeList_x[i_x];
		edge.npts = npts;
		edge.pt[0] = edge.vert[0];
		for (k=1; k<npts-1; k++) {
			j = net0->edgeList[i].pt[k];
			oldpt[np_x] = j;
			point_x[np_x] = net0->point[j];
			edge.pt[k] = np_x;
			np_x++;
		}
		edge.pt[npts-1] = edge.vert[1];
		edgeList_x[i_x] = edge;
		i_x++;
	}
	printf("Added interior edge points\n");
	printf("ne, ne_x: %d %d  nv, nv_x: %d %d  np, np_x: %d %d\n",ne,ne_x,nv,nv_x,np,np_x);

	net1->ne = ne_x;
	net1->nv = nv_x;
	net1->np = np_x;
	free(oldpt);
	return 0;
}
//-----------------------------------------------------------------------------------------------------
// This program doubles as pruner-cleaner and as a way to create a network with constant vessel diameter
//-----------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int n_prune_cycles, err;
	char *input_amfile;
	char drive[32], dir[1024],filename[256], ext[32];
	char errfilename[1024], output_amfile[1024], outfilename[1024], result_file[1024];
	char output_basename[128];
	int prune_flag, cmgui_flag, ratio_flag;
	float origin_shift[3];
	NETWORK *net0, *net1;

	if (argc != 11 && argc != 4) {
		printf("To generate a network with constant diameters: \n");
		printf("Usage: prune input_amfile output_amfile fixed_diam\n");
		printf("\nTo clean up and prune network:\n");
		printf("Usage: prune input_amfile output_amfile ratio_flag ratio_limit n_prune_cycles prune_flag cmgui_flag delta_diam delta_len dminimum\n");
		printf("       ratio_flag: = 1 if length/diameter is to be used, 0 if length limit\n");
		printf("       ratio_limit: either length/diameter limit or length limit\n");
		printf("       prune_cycle: number of cycles of pruning\n");
		printf("       prune_flag: = 1 to prune, 0 else\n");
		printf("       cmgui_flag: = 1 to generate CMGUI files\n");
		printf("       delta_diam: box size (um) for diameter histogram\n");
		printf("       delta_len: box size (um) for length histogram\n");
		printf("       dminimum: all edges with diameter < dminimum are to be removed\n");
		fperr = fopen("prune_error.log","w");
		fprintf(fperr,"Usage: prune input_amfile output_amfile ratio_flag ratio_limit n_prune_cycles prune_flag cmgui_flag delta_diam delta_len dminimum\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line
	}

	input_amfile = argv[1];
	strcpy(outfilename,argv[2]);
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_prune.log",output_basename);
	sprintf(output_amfile,"%s.am",output_basename);
	sprintf(result_file,"%s.out",output_basename);
	fperr = fopen(errfilename,"w");
	if (argc == 4) {
		sscanf(argv[3],"%f",&fixed_diam);
		use_fixed_diam = true;
		printf("fixed_diam: %f\n",fixed_diam);
		origin_shift[0] = 0;
		origin_shift[1] = 0;
		origin_shift[2] = 0;
		net0 = (NETWORK *)malloc(sizeof(NETWORK));
		fpout = fopen(result_file,"w");	
		err = ReadAmiraFile_old(input_amfile,net0);
		if (err != 0) return 2;
		printf("did ReadAmiraFile_old\n");
		err = WriteAmiraFile(output_amfile,input_amfile,net0,origin_shift);
		if (err != 0) return 8;
		printf("did WriteAmiraFile\n");
		err = WriteCmguiData(output_basename,net0,origin_shift);
		if (err != 0) return 10;
		printf("did WriteCmguiData\n");
		return 0;
	} else {
		sscanf(argv[3],"%d",&ratio_flag);
		sscanf(argv[4],"%lf",&ratio_limit);
		sscanf(argv[5],"%d",&n_prune_cycles);
		sscanf(argv[6],"%d",&prune_flag);
		sscanf(argv[7],"%d",&cmgui_flag);
		sscanf(argv[8],"%f",&ddiam);
		sscanf(argv[9],"%f",&dlen);
		sscanf(argv[10],"%f",&dminimum);
		use_fixed_diam = false;
	}
	use_ratio = (ratio_flag == 1);
	if (prune_flag == 0) n_prune_cycles = 0;
	_splitpath(outfilename,drive,dir,filename,ext);
	strcpy(output_basename,drive);
	strcat(output_basename,dir);
	strcat(output_basename,filename);
	sprintf(errfilename,"%s_prune.log",output_basename);
	sprintf(output_amfile,"%s.am",output_basename);
	sprintf(result_file,"%s.out",output_basename);
	fperr = fopen(errfilename,"w");

//	fprintf(fperr,"drive: %s dir: %s filename: %s ext: %s\n",drive,dir,filename,ext);
//	fprintf(fperr,"Basename: %s\n",output_basename);

	use_len_diam_limit = false;
	use_len_limit = false;

	net0 = (NETWORK *)malloc(sizeof(NETWORK));

	fpout = fopen(result_file,"w");	
	err = ReadAmiraFile_old(input_amfile,net0);
	if (err != 0) return 2;
	if (n_prune_cycles > 0) {
		if (dminimum > 0) {
			err = EdgeDimensions(net0->edgeList,net0->point,net0->ne);
			if (err != 0) return 8;
		}
		err = adjoinEdges(net0);
		if (err != 0) return 3;
		err = deloop(net0);
		if (err != 0) return 4;
		err = adjoinEdges(net0);
		if (err != 0) return 3;
		for (int k=0; k<n_prune_cycles; k++) {
			err = pruner(net0,k);
			if (err != 0) return 5;
			err = adjoinEdges(net0);
			if (err != 0) return 3;
			err = checkEdgeEndPts(net0);
			if (err != 0) return 6;
		}
		net1 = (NETWORK *)malloc(sizeof(NETWORK));
		err = squeezer(net0,net1);	// must squeeze, or SpatialGraph and CMGUI files are not consistent
		if (err != 0) return 7;
	} else {
		net1 = net0;
	}
	origin_shift[0] = 0;
	origin_shift[1] = 0;
	origin_shift[2] = 0;
	err = WriteAmiraFile(output_amfile,input_amfile,net1,origin_shift);
	if (err != 0) return 8;
	if (cmgui_flag == 1) {
		err = WriteCmguiData(output_basename,net1,origin_shift);
		if (err != 0) return 10;
	}
	err = EdgeDimensions(net1->edgeList,net1->point,net1->ne);
	if (err != 0) return 8;
	err = CreateDistributions(net1);
//	err = CreateDistributions(ddiam, dlen);
	if (err != 0) return 9;
	return 0;
}