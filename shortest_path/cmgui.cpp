// test
#include <cstdio>
#include <math.h>

#include "shortest_path.h"

extern int nv, ne, np;
extern int nv_used, ne_used, np_used;
extern int vertex_case;
extern bool use_nb;
extern int nbmax[2];
extern float artery_col[3], vein_col[3];
extern float distmin[2], distmax[2];

extern EDGE *edgeList;
extern VERTEX *vertex;
extern POINT *point;
extern FILE *fperr, *fpout;


//-----------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------
float dist(int k1, int k2)
{
	float dx = point[k2].x - point[k1].x;
	float dy = point[k2].y - point[k1].y;
	float dz = point[k2].z - point[k1].z;
	return sqrt(dx*dx+dy*dy+dz*dz);
}

//-----------------------------------------------------------------------------------------------------
// Create CMGUI .com file
//-----------------------------------------------------------------------------------------------------
void write_com(FILE *dotcom, char *fileName)
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
//    fprintf(dotcom, "gfx mod g_e vessels element_points glyph sphere_hires general size \"0*0*0\" ");
//    fprintf(dotcom, "orientation vessel_radius scale_factors \"2*2*2\" discretization \"3*3*3\" material gold");
    fprintf(dotcom, "gfx destroy node all\n");
    fprintf(dotcom, "gfx modify g_element vessels general clear\n");
    fprintf(dotcom, "gfx modify g_element vessels cylinders coordinate coordinates tessellation default local circle_discretization 12 radius_scalar vessel_radius scale_factor 0.5 native_discretization NONE select_on material gold selected_material default_selected render_shaded\n");
    fprintf(dotcom, "gfx modify g_element vessels node_points coordinate coordinates local glyph sphere general size \"0*0*0\" centre 0,0,0 font default orientation vessel_radius scale_factors \"1*1*1\" select_on material gold selected_material default_selected\n");
    fprintf(dotcom, "gfx cre win 1\n");
    fprintf(dotcom, "gfx mod win 1 view perspective\n");
}


//-----------------------------------------------------------------------------------------------------
// Write initial section of .exnode file
//-----------------------------------------------------------------------------------------------------
void WriteExnodeHeader(FILE *exnode)
{
//    fprintf(exnode, "Region: /vessels\n");
    fprintf(exnode, "Group name: vessels\n");
    fprintf(exnode, " #Fields=4\n");
    fprintf(exnode, " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exnode, "  x.  Value index=1, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  y.  Value index=2, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  z.  Value index=3, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 2) vessel_radius, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exnode, "  1.  Value index=4, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 3) branch_count, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exnode, "  1.  Value index=5, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, " 4) RGB_vector, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exnode, "  1.  Value index=6, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  2.  Value index=7, #Derivatives=0, #Versions=1\n");
    fprintf(exnode, "  3.  Value index=8, #Derivatives=0, #Versions=1\n");
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
    fprintf(exelem, " #Nodes= 2\n #Fields=4\n");
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
    fprintf(exelem, " 3) branch_count, coordinate, rectangular cartesian, #Components=1\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, " 4) RGB_vector, coordinate, rectangular cartesian, #Components=3\n");
    fprintf(exelem, "   1.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   2.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
    fprintf(exelem, "      1.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   1\n");
    fprintf(exelem, "      2.  #Values=1\n");
    fprintf(exelem, "       Value indices:     1\n       Scale factor indices:   2\n");
    fprintf(exelem, "   3.  l.Lagrange, no modify, standard node based.\n     #Nodes= 2\n");
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
int WriteCmguiData(char *basename)
{
	int k, ie, ip, npts, iv, i;
	float fraction, rgb[3];
	EDGE edge;
	char exelemname[256], exnodename[256], dotcomname[256];
	FILE *exelem, *exnode, *dotcom;
	bool is_end;

	printf("WriteCmguiData: %s %d %d\n",basename,ne,np);
	fprintf(fperr,"WriteCmguiData: %s\n",basename);
	fflush(fperr);
	sprintf(dotcomname,"%s.com.txt",basename);
	sprintf(exnodename,"%s.exnode",basename);
	sprintf(exelemname,"%s.exelem",basename);
	printf("exelem file: %s exnode file: %s\n",exelemname,exnodename);
	fprintf(fperr,"exelem file: %s exnode file: %s\n",exelemname,exnodename);
	fflush(fperr);
	dotcom = fopen(dotcomname,"w");
	exelem = fopen(exelemname,"w");
	exnode = fopen(exnodename,"w");
	write_com(dotcom, basename);
	WriteExelemHeader(exelem);
	WriteExnodeHeader(exnode);
	printf("wrote headers: ne: %d np: %d\n",ne,np);

	for (k=0; k<np; k++) {
		point[k].used = false;
	}
	int kelem = 0;
	for (ie=0; ie<ne; ie++) {
		edge = edgeList[ie];
//		printf("ie: %d nb: %d\n",ie,edge.nb);
		if (edge.used) {
			npts = edge.npts;
//			if (edge.nb[0] == 0) {
//				printf("edge with nb = 0: %d vertices: %d %d\n",ie,edge.vert[0],edge.vert[1]);
//				if (ie > 10) return 1;
//			}
	//		npts = edge.npts_used;
			int kfrom = edge.pt[0];
			if (edge.vert[0] != kfrom) {
				printf("vert[0] != kfrom: ie: %d  %d %d\n",ie,edge.vert[0],kfrom);
				return 1;
			}
			if (edge.vert[1] != edge.pt[npts-1]) {
				printf("vert[1] != kto: ie: %d  %d %d\n",ie,edge.vert[1],edge.pt[npts-1]);
				return 1;
			}
			point[kfrom].used = true;
	//		int kfrom = edge.pt_used[0];
	//		point_used[kfrom] = true;
			for (ip=1; ip<npts; ip++) {
				int k2 = edge.pt[ip];
//				if (point[k2].nb == 0) {
//					printf("nb=0: edge: %d ip: %d k2: %d\n",ie,ip,k2);
//				}
	//			int k2 = edge.pt_used[ip];
				double d = dist(kfrom,k2);
				if (d > 0.5*(point[kfrom].d/2+point[k2].d/2) || ip == npts-1) {
					point[k2].used = true;
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
	int nbmaxx = nbmax[0];
	if (nbmax[1] > nbmaxx) nbmaxx = nbmax[1];
	for (k=0; k<np; k++) {
		if (point[k].used) {
			is_end = false;
			if (k == 9164 || k == 7562) is_end = true;
			rgb[0] = rgb[1] = rgb[2] = 0;
			fprintf(exnode, "Node: %d\n", k+1);
			fprintf(exnode, "%6.1f %6.1f %6.1f\n", point[k].x,point[k].y,point[k].z);
			fprintf(exnode, "%6.2f\n", point[k].d/2);
			if (point[k].d == 0) {
				printf("Error: zero diameter: %d\n",k);
				return 1;
			}
			if (vertex_case == 1 || vertex_case == 3)
				iv = 0;
			else
				iv = 1;
			fprintf(exnode, "%d\n", point[k].nb[iv]);		// currently only a single nb value is stored
			if (vertex_case == 1 || vertex_case == 3) {
				iv = 0;
				if (use_nb) {
					fraction = float(nbmax[iv] - point[k].nb[iv])/nbmax[iv];
//					fraction = float(nbmaxx - point[k].nb[iv])/nbmaxx;
					fraction = fraction*fraction;
					if (is_end) {
						fprintf(fpout,"end: iv,k,nbmax,nb: %3d %6d %4d %4d %6.3f\n",iv,k,nbmaxx,point[k].nb[iv],fraction);
						printf("end: iv,k,nbmax,nb: %3d %6d %4d %4d %6.3f\n",iv,k,nbmaxx,point[k].nb[iv],fraction);
					}
				} else {
					fraction = (distmax[iv] - point[k].distance[iv])/distmax[iv];
//					fprintf(fpout,"%6d %6.1f %6.1f %6.3f\n",k,point[k].distance[0],distmax[iv],fraction);
				}
				for (i=0; i<3; i++)
					rgb[i] += fraction*artery_col[i];
			}
			if (vertex_case == 2 || vertex_case == 3) {
				iv = 1;
				if (use_nb) {
					fraction = float(nbmax[iv] - point[k].nb[iv])/nbmax[iv];
//					fraction = float(nbmaxx - point[k].nb[iv])/nbmaxx;
					fraction = fraction*fraction;
					if (is_end) {
						fprintf(fpout,"end: iv,k,nbmax,nb: %3d %6d %4d %4d %6.3f\n",iv,k,nbmaxx,point[k].nb[iv],fraction);
						printf("end: iv,k,nbmax,nb: %3d %6d %4d %4d %6.3f\n",iv,k,nbmaxx,point[k].nb[iv],fraction);
					}
				} else {
					fraction = (distmax[iv] - point[k].distance[iv])/distmax[iv];
				}
				for (i=0; i<3; i++)
					rgb[i] += fraction*vein_col[i];
			}
			fprintf(exnode, "%5.2f %5.2f %5.2f\n",rgb[0],rgb[1],rgb[2]);
		}
	}
	printf("wrote node data\n");
	fclose(dotcom);
	fclose(exelem);
	fclose(exnode);
	return 0;
}
