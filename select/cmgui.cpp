// test
#include <cstdio>

#include "network.h"

extern int nv, ne, np;
extern int nv_used, ne_used, np_used;
extern EDGE *edgeList;
extern VERTEX *vertex;
extern POINT *point;
extern FILE *fperr, *fpout;
extern float dist(NETWORK *net, int k1, int k2);

int testfunction(void)
{
	printf("testfunction\n");
	return 0;
}

//-----------------------------------------------------------------------------------------------------
// Write initial section of .exnode file
//-----------------------------------------------------------------------------------------------------
void WriteExnodeHeader(FILE *exnode)
{
//    fprintf(exnode, "Region: /vessels\n");
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
void WriteExelemHeader(FILE *exelem)
{
 //   fprintf(exelem, "Region: /vessels\n");
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
// To adjust for the reduced number of nodes in the simplified network, we need a way to handle the 
// node numbers.  The nodes are currently numbered 1...np, but not all these node numbers are used.
// Does the node file need all nodes from 1 to np?
// A quick test seems to show that neither element nor node numbers need to be consecutive.
//-----------------------------------------------------------------------------------------------------
int WriteCmguiData(char *basename, NETWORK *net, float origin_shift[])
{
	int k, ie, ip, npts;
	EDGE edge;
	char exelemname[128], exnodename[128];
//	char dotcomname[64];
	FILE *exelem, *exnode;

	printf("WriteCmguiData: %s %d %d\n",basename,net->ne,net->np);
	fprintf(fperr,"WriteCmguiData: %s\n",basename);
	fflush(fperr);
//	sprintf(dotcomname,"%s.com",output_basename);
	sprintf(exelemname,"%s.exelem",basename);
	sprintf(exnodename,"%s.exnode",basename);
	printf("exelem file: %s\n",exelemname);
	printf("exnode file: %s\n",exnodename);
	fprintf(fperr,"exelem file: %s\n",exelemname);
	fprintf(fperr,"exnode file: %s\n",exnodename);
	fflush(fperr);
//	dotcom = fopen(dotcomname,"w");
	exelem = fopen(exelemname,"w");
	exnode = fopen(exnodename,"w");
//	write_com(output_basename);
	WriteExelemHeader(exelem);
	WriteExnodeHeader(exnode);
	printf("wrote headers: ne: %d np: %d\n",net->ne,net->np);
	int kelem = 0;
	for (ie=0; ie<net->ne; ie++) {
		edge = net->edgeList[ie];
//		printf("ie: %d used: %d npts: %d\n",ie,edge.used,edge.npts);
		if (edge.used) {
			npts = edge.npts;
	//		npts = edge.npts_used;
			int kfrom = edge.pt[0];
	//		int kfrom = edge.pt_used[0];
	//		point_used[kfrom] = true;
			for (ip=1; ip<npts; ip++) {
				int k2 = edge.pt[ip];
				double d = dist(net,kfrom,k2);
				if (d > 0.5*(net->point[kfrom].d/2+net->point[k2].d/2) || ip == npts-1) {
					kelem++;
					fprintf(exelem, "Element: %d 0 0\n", kelem);
					fprintf(exelem, "  Nodes: %d %d\n", kfrom+1, k2+1);
					fprintf(exelem, "  Scale factors: 1 1\n");
					kfrom = k2;
				}
			}
		}
	}
	printf("wrote elem data\n");
	for (k=0; k<net->np; k++) {
		if (net->point[k].used) {
			fprintf(exnode, "Node: %d\n", k+1);
			fprintf(exnode, "%6.1f %6.1f %6.1f\n",
				net->point[k].x - origin_shift[0],
				net->point[k].y - origin_shift[1],
				net->point[k].z - origin_shift[2]);
			fprintf(exnode, "%6.2f\n", net->point[k].d/2);
//			if (net->point[k].d > 75) printf("point: %d diameter: %6.1f\n",k,net->point[k].d);
			if (net->point[k].d == 0) {
				printf("Error: zero diameter: %d\n",k);
				return 1;
			}
		}
	}
	printf("wrote node data\n");
//	fclose(dotcom);
	fclose(exelem);
	fclose(exnode);
	return 0;
}