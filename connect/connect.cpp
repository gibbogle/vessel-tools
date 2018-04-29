/*
 * To find small objects in an image 
 */
using namespace std;

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
#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkSize.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMeanImageFilter.h"
//#include "itkBoxMeanImageFilter.h"

#define MAXPAIRS 100
#define MAXLABELS 1000000	// 400000
#define MAXLABLIST 2000		// 100
#define NBRS 26
#define MAXEXP 300
#define MAXOBJECTS 200000
#define NBIGOBJECTS 5

#define USE_SHORT false

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
long long width, height, depth, imsize;
unsigned char *p;
//unsigned short *label;
unsigned int *label;

//int lablist[MAXLABELS][MAXLABLIST];
int **lablist;
//int nlablist[MAXLABELS];
int *nlablist;
int maxnlablist;
int nbr[26][3];
int nbr2D[8][2];

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]		// V is used to access the voxel array p[]
#define C(a,b,c)  label[(c)*imsize+(b)*width+(a)]	// C is used to access the label[] array

int max16 = pow(2.,16)-1;

struct pair_str 
{
	int c1, c2;
};
typedef pair_str PAIR;

PAIR apair[MAXPAIRS];
int cnt[MAXLABELS];
int npairs;
int cmax;
//int bigobject_label[MAXOBJECTS];
//int bigobject_nvoxels[MAXOBJECTS];

struct bigObject {
  int label;
  int nvoxels;
} bigobj[MAXOBJECTS];

FILE *fp, *fpout;

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

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
bool by_nvoxels ( const bigObject &a, const bigObject &b )
{
    return a.nvoxels > b.nvoxels;
}

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
int countvoxels(int nlist, int list[])
{
	int x, y, z, n, k, total;
	unsigned char *labelmap;

	total = 0;
	for (k=0; k<nlist; k++)
		total += cnt[list[k]];
	printf("countvoxels: nlist: %d total: %d\n",nlist,total);

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
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int combine(int lab, int list[], int *nlistresult, int cmap[], bool taken[])
{

	int k, i2, c2;
	bool added, in2;
//	int list[MAXLABELS];
	int nlist;
//	bool taken[MAXLABELS];

//	printf("combine: %d\n",lab);
	taken[lab] = true;
	list[0] = lab;
	nlist = 1;
	for (k=0; k<nlablist[lab]; k++)
	{
		c2 = lablist[lab][k];
		list[nlist] = c2;
		nlist++;
		cmap[c2] = lab;
	}
	added = true;
	while (added)		// I think one sweep is NOT enough
	{
		added = false;
		// Look in all untaken lablists for a label in list[].  If found, take all the labels.
		for (i2=0; i2<cmax; i2++)
		{
			if (taken[i2])
				continue;
			in2 = false;
			for (k=0; k<nlablist[i2]; k++)
			{
				c2 = lablist[i2][k];
				if (inlist(list,nlist,c2))
				{
					in2 = true;
					break;
				}
			}
			if (in2)
			{
				// copy all the labels into list
				if (nlist + lablist[i2][k] >= MAXLABELS)
				{
					printf("Size of list[] in tracer (MAXLABELS) exceeded\n");
					exit(6);
				}
				for (k=0; k<nlablist[i2]; k++)
				{
					c2 = lablist[i2][k];
					if (!inlist(list,nlist,c2))
					{
						list[nlist] = c2;
						nlist++;
						cmap[c2] = lab;
					}
				}
				taken[i2] = true;
				added = true;
			}
		}
	}
	*nlistresult = nlist;
	return 0;
}

//------------------------------------------------------------------------------------------------
// At this stage the image with buffer C() has values 0,1,2,...,cmax
// Based on the info in lablist[][], some of these values are equivalent (i.e. connected).
// Need to determine the mapping from all C-numbers to the unique values
// corresponding to the connected regions, and the associated counts.
// For a starting label lab, list[] holds all the labels that are connected to lab.
//------------------------------------------------------------------------------------------------
int tracer(int minvoxels, bool drop[])
{
	int lab, k;
	int total, nobjects, nsignificant, ndropped;
	int *list;
	int *cmap;
	int nlist;
	bool *taken;

	printf("tracer: cmax: %d\n",cmax);
	fprintf(fpout,"tracer: cmax: %d\n",cmax);
	fflush(fpout);
	list = (int *)malloc(MAXLABELS * sizeof(int));
	cmap = (int *)malloc(MAXLABELS * sizeof(int));
	taken = (bool *)malloc(MAXLABELS * sizeof(bool));
	fflush(fpout);

	nobjects = 0;
	nsignificant = 0;
	for (lab=1; lab<=cmax; lab++) 
	{
		drop[lab] = false;
		cmap[lab] = lab;
		taken[lab] = false;
	}

	for (lab=1; lab<=cmax; lab++) 
	{
		if (cmap[lab] != lab) continue;

		combine(lab, list, &nlist, cmap, taken);
		/*
		// TESTING!!!!!!!!!!!
		if (lab == 99)
		{
			int n = countvoxels(nlist,list);
			printf("lab: %d has %d voxels\n",lab,n);
			exit(1);
		}
		*/
		nobjects++;
		total = 0;
		for (k=0; k<nlist; k++)
			total += cnt[list[k]];
		if (total < minvoxels)
		{
			printf("Dropping label: %d nvoxels: %d\n",lab,total);
			for (k=0; k<nlist; k++)
				drop[list[k]] = true;
		}
		else
		{
			printf("object: %6d label: %6d nvoxels: %8d\n",nobjects,lab,total);
			if (nsignificant == MAXOBJECTS-1)
			{
				printf("Too many objects\n");
				exit(7);
			}
			bigobj[nsignificant].label = lab;
			bigobj[nsignificant].nvoxels = total;
			nsignificant++;
		}
	}
	sort( bigobj, bigobj+nsignificant, by_nvoxels );

	ndropped = 0;
	for (lab=1; lab<=cmax; lab++)
	{
		if (cmap[lab] == lab)
		{
			if (drop[lab])
				ndropped++;
		}
		else 
		{
//			printf("%d --> %d\n",i,cmap[i]);
		}
	}
	printf("nobjects: %d  ndropped: %d  retained: %d\n\n",nobjects,ndropped,nobjects-ndropped);
	fprintf(fpout,"nobjects: %d  ndropped: %d  retained: %d\n\n",nobjects,ndropped,nobjects-ndropped);
	fflush(fpout);
	for (lab=0; lab<NBIGOBJECTS; lab++)
	{
		printf("label: %d  nvoxels: %d\n",bigobj[lab].label,bigobj[lab].nvoxels);
		fprintf(fpout,"label: %d  nvoxels: %d\n",bigobj[lab].label,bigobj[lab].nvoxels);
		fflush(fpout);
	}
	free(list);
	free(cmap);
	free(taken);
	return nobjects-ndropped;
}


//------------------------------------------------------------------------------------------------
// lablist[i][:] is the list of all labels that are paired with (connected to) label i.
// This means that somewhere a voxel labelled i1 is a neighbour of a voxel labelled i2.
// The current number of entries in the list is nlablist[i].
// A list is maintained for every label value, therefore when a pair of labels, i1 and i2, are
// found to be connected, an addition is made to lablist[i1] and to lablist[i2], provided
// i2 is not already in lablist[i1] (and by symmetry i1 is in lablist[i2]).
// npairs is the total number of label pairings.
//------------------------------------------------------------------------------------------------
bool addPair(int i1, int i2)
{
	int k;

	for (k=0; k<nlablist[i1]; k++) {
		if (lablist[i1][k] == i2)
			return false;		// already in the list(s)
	}
	// Add to lists
	lablist[i1][nlablist[i1]] = i2;
	nlablist[i1]++;
	if (nlablist[i1] == MAXLABLIST)
	{
		printf("MAXLABLIST exceeded\n");
		exit(8);
	}
	maxnlablist = MAX(maxnlablist,nlablist[i1]);

	lablist[i2][nlablist[i2]] = i1;
	nlablist[i2]++;
	if (nlablist[i2] == MAXLABLIST)
	{
		printf("MAXLABLIST exceeded\n");
		exit(8);
	}
	maxnlablist = MAX(maxnlablist,nlablist[i2]);

	npairs++;
	return true;
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

	c1 = C(x,y,z);		// the label of voxel (x,y,z)
	nn = 0;
	for (k=0; k<NBRS; k++)
	{
		dx = nbr[k][0];
		dy = nbr[k][1];
		dz = nbr[k][2];
		xx = x+dx;
		if (xx < 0 || xx >= width) continue;
		yy = y+dy;
		if (yy < 0 || yy >= height) continue;
		zz = z+dz;
		if (zz < 0 || zz >= depth) continue;

		v2 = V(xx,yy,zz);
		if (v2 == 0) continue;

		// lit neighbour voxel at (xx,yy,zz)
		c2 = C(xx,yy,zz);
		if (c2 > 0)
		{
			if (c1 == 0)			// c1=0, c2>0
			{
				C(x,y,z) = c2;
				cnt[c2]++;	
				c1 = c2;
			}
			else					// c1>0, c2>0
			{
				if (c1 != c2)
				{
					addPair(c1,c2);
				}
			}
		}
		else 
		{
			if (c1 > 0) 
			{			// c1>0, c2=0
				C(xx,yy,zz) = c1;
				cnt[c1]++;
			}
			klist[nn] = k;
			nn++;
		}
	}
	*n = nn;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
void connecter(void)
{
	int x, y, z, xx, yy, zz, x0, y0, z0, c1, k, kk, n, iexp, nexp;
	int klist[26];
//	int explist[MAXEXP][3];
	int **explist;
	unsigned char v1;

	explist = (int **)malloc(MAXEXP*sizeof(int *));
	for (k=0; k<MAXEXP; k++)
		explist[k] = (int *)malloc(3*sizeof(int));

	maxnlablist = 0;
	setupNeighbours();
	printf("connecter:\n");
	fprintf(fpout,"connecter:\n");
	fflush(fpout);
	npairs = 0;
	cmax = 0;
	for (k=0; k<MAXLABELS; k++)
	{
		nlablist[k] = 0;
		cnt[k] = 0;
	}
//	printf("Initialized nlablist[]\n");
	memset(label,0,imsize*depth*sizeof(unsigned int));	// a label value is an unsigned integer
//	printf("Initialized label[]\n");

	for (z=0; z<depth; z++)
	{
		printf("z: %d\n",z);
		for (y=0; y<height; y++)
		{
//		printf("y: %d\n",y);
			for (x=0; x<width; x++)
			{
//		printf("x: %d\n",x);
				v1 = V(x,y,z);
				if (v1 == 0) continue;

				explore(x,y,z,&n,klist);

				c1 = C(x,y,z);
				if (c1 == 0)	// This is possibly a new label
				{
					cmax++;
					if (cmax == MAXLABELS)
					{
						printf("MAXLABELS exceeded in connecter\n");
						fprintf(fpout,"MAXLABELS exceeded in connecter\n");
						fflush(fpout);
						exit(9);
					}
					if (USE_SHORT && cmax > max16)
					{
						printf("max short size exceeded\n");
						exit(10);
					}
					C(x,y,z) = cmax;
					c1 = cmax;
					cnt[c1]++;
				}

				if (n > 0) 
				{
					// Now give this label to the other lit voxels not yet given a label
					nexp = 0;
					for (k=0; k<n; k++)
					{
						kk = klist[k];
						xx = x + nbr[kk][0];
						yy = y + nbr[kk][1];
						zz = z + nbr[kk][2];
						if (C(xx,yy,zz) == 0) {
							C(xx,yy,zz) = c1;
							cnt[c1]++;
						}
//						printf("base exp: %d  %d %d %d\n",nexp,xx,yy,zz);
						if (nexp == MAXEXP-1)
						{
							printf("explist dimension exceeded (a)\n");
							exit(11);
						}
						explist[nexp][0] = xx;
						explist[nexp][1] = yy;
						explist[nexp][2] = zz;
						nexp++;
					}
					// This is the basis for a ranging exploration, in which new connected lit voxels
					// are added to the explist, up to a specified maximum.  The exploration is restricted
					// to the x-y plane, i.e. to a slice.
					iexp = 0;
					while (iexp < nexp && nexp < MAXEXP-1)
					{
						x0 = explist[iexp][0];
						y0 = explist[iexp][1];
						z0 = explist[iexp][2];
//						printf("    iexp: %d %d   %d %d %d\n",iexp,nexp,x0,y0,z0);
						// look for lit, unlabeled voxels adjacent to explist[iexp].
						// if one is found, label it c1, increment cnt[c1], 
						// and add it to explist[] as long as nexp < MAXEXP
						for (k=0; k<8; k++)
						{
							xx = x0 + nbr2D[k][0];
							if (xx < 0 || xx > width-1) continue;
							yy = y0 + nbr2D[k][1];
							if (yy < 0 || yy > height-1) continue;
							zz = z0;
//							printf("k: %d  %d %d %d  V: %d C: %d\n",k,xx,yy,zz,V(xx,yy,zz),C(xx,yy,zz));
							if (V(xx,yy,zz) == 0) continue;
//							if (C(xx,yy,zz) != 0) continue;		// try this modification
							if (C(xx,yy,zz) != 0) 
							{
								if ((int)C(xx,yy,zz) != c1)
								{
									addPair(C(xx,yy,zz),c1);
								}
								continue;
							}
							C(xx,yy,zz) = c1;
							cnt[c1]++;
							explist[nexp][0] = xx;
							explist[nexp][1] = yy;
							explist[nexp][2] = zz;
//							printf("        add: %d %d   %d %d %d\n",k,nexp,xx,yy,zz);
							nexp++;
							if (nexp == MAXEXP-1)
							{
//								printf("explist dimension exceeded (b)\n");
								break;
							}
						}
						iexp++;
					}
//					printf("Did: %d %d %d nexp: %d  cmax: %d\n",x,y,z,nexp,cmax);
				}
			}
		}
	}
	printf("\ncmax: %d npairs: %d\n",cmax,npairs);
	printf("max nlablist: %d\n",maxnlablist);
	fprintf(fpout,"\ncmax: %d npairs: %d\n",cmax,npairs);
	fprintf(fpout,"max nlablist: %d\n",maxnlablist);
	fflush(fpout);

	for (k=0; k<MAXEXP; k++)
		free(explist[k]);
	free(explist);
}


//------------------------------------------------------------------------------------------------
// Connected objects with fewer than the threshold voxels are removed (set to 0).
//------------------------------------------------------------------------------------------------
void remover(bool drop[])
{
	int x, y, z, lab;

	for (x=0; x<width; x++)
	{
		for (y=0; y<height; y++)
		{
			for (z=0; z<depth; z++)
			{
				lab = C(x,y,z);
				if (lab != 0)
				{
					if (drop[lab])
						V(x,y,z) = 0;
				}
			}
		}
	}
}

//------------------------------------------------------------------------------------------------
// Connected objects with fewer than the threshold voxels (<=2) are removed (set to 0).
//------------------------------------------------------------------------------------------------
void despeckle(void)
{
	int x, y, z, dx, dy, dz, xx, yy, zz, k;

	for (x=0; x<width; x++)
	{
		for (y=0; y<height; y++)
		{
			for (z=0; z<depth; z++)
			{
				if (V(x,y,z) == 0) continue;
				k = 0;
				for (dx=-1;dx<=1;dx++) {
					xx = x+dx;
					if (xx < 0 || xx > width-1) continue;
					for (dy=-1;dy<=1;dy++) {
						yy = y+dy;
						if (yy < 0 || yy > height-1) continue;
						for (dz=-1;dz<=1;dz++) {
							zz = z+dz;
							if (zz < 0 || zz > depth-1) continue;
							if (V(xx,yy,zz) != 0) k++;
						}
					}
				}
//				printf("%d %d %d   %d\n",x,y,z,k);
				if (k <= 2) {
					V(x,y,z) = 0;
				}

			}
		}
	}
}

//------------------------------------------------------------------------------------------------
// First create a lookup table for all used label values, such that those connected to lab map
// to 255, and the rest map to 0.  This requires calling combine, then using nlist and list[].
// Then the construction of the desired image is simple: each non-zero C(x,y,z) maps to 255 or 0.
//------------------------------------------------------------------------------------------------
int createObjectImage(int lab, int cmap[], bool taken[])
{
	int x, y, z, n, k, total;
	int nlist;
//	int list[MAXLABELS];
	int *list;
	unsigned char *labelmap;

	fprintf(fpout,"createObjectImage: %d\n",lab);
	fflush(fpout);
	list = (int *)malloc(MAXLABELS*sizeof(int));
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
	free(list);
	return n;
}


//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
//	bool drop[MAXLABELS];
//	int cmap[MAXLABELS];
//	bool taken[MAXLABELS];
	bool *drop;
	int *cmap;
	bool *taken;
	int minvoxels=0, i, lab, nobjects;
	time_t t1;
	t1 = time(NULL);
	char objectFile[128];
//	char baseName[] = "object";
	char *infile;
	char *baseName;
	char numstr[2];
	char drive[32], dir[128],filename[64], ext[32];
	char errorFile[128], outputFile[128], outputPath[128];
	char *temp_name;
	FILE *fperr;

	if (argc != 4) {
		fperr = fopen("connect_error.log","w");
		printf("Usage: connect input_tiff base minvoxels\n");
		printf("       where: base is the base for output filenames,\n");
		printf("              e.g. if base = 'A' the output files are 'A00.tif', 'A01.tif', ...\n");
		printf("            : minvoxels is the minimum size of a connected object\n");
		fprintf(fperr,"Usage: connect input_tiff base minvoxels\n");
		fprintf(fperr,"       where: base is the base for output filenames,\n");
		fprintf(fperr,"              e.g. if base = 'A' the output files are 'A00.tif', 'A01.tif', ...\n");
		fprintf(fperr,"            : minvoxels is the minimum size of a connected object\n");
		fprintf(fperr,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fperr,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fperr);
		return 1;	// Wrong command line 
	}

	infile = argv[1];
	baseName = argv[2]; 
//	_splitpath(infile,drive,dir,filename,ext);
//	strcpy(outputPath,drive);
//	strcat(outputPath,dir);
#if (defined (_WIN32) || defined (_WIN64))
    // windows code
	_splitpath(infile,drive,dir,filename,ext);
	strcpy(outputPath,drive);
	strcat(outputPath,dir);
//	strcat(output_basename,filename);
#elif (defined (LINUX) || defined (__linux__))
    // linux code
	strcpy(outputPath, basename(infile));
	temp_name = remove(outputPath);
	strcpy(outputPth,temp_name);
#endif
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

	drop = (bool *)malloc(MAXLABELS * sizeof(bool));
	cmap = (int *)malloc(MAXLABELS * sizeof(int));
	taken = (bool *)malloc(MAXLABELS * sizeof(bool));
	nlablist = (int *)malloc(MAXLABELS * sizeof(int));
//	printf("Allocated nlablist\n");
	lablist = (int **)malloc(MAXLABELS * sizeof(int *));
	if (lablist == NULL)
	{
		printf("out of memory: lablist\n");
		fprintf(fp,"out of memory: lablist\n");
		fclose(fp);
		return 4;
	}
//	printf("Allocated lablist **\n");
	for (i = 0; i < MAXLABELS; i++)
	{
		lablist[i] = (int *)malloc(MAXLABLIST * sizeof(int));
		if (lablist[i] == NULL)
		{
			printf("out of memory: lablist[i]\n");
			fprintf(fp,"out of memory: lablist\n");
			fclose(fp);
			return 4;
		}
	}
//	printf("Allocated lablist[][]\n");


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
	imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	p = (unsigned char *)(im->GetBufferPointer());
	if (USE_SHORT) {
//		label = (unsigned short*)malloc(imsize*depth*sizeof(unsigned short));  
	}
	else
		label = (unsigned int*)malloc(imsize*depth*sizeof(unsigned int));
	if (label == NULL) {
		printf("Allocation of label failed: size: %d\n",imsize*depth);
		return 5;
	}
	despeckle();
	printf("Did despeckle\n");
	fprintf(fpout,"Did despeckle\n");
	fflush(fpout);
	fclose(fpout);

	int nvoxels = 0;
	for(i=0; i<imsize*depth; i++) {
		if (p[i] != 0) nvoxels++;
	}
	fprintf(fpout,"After despeckling, total lit voxels = %d\n",nvoxels);
//	fprintf(fpout,"Not doing connecting\n");
//	return 1;

	connecter();
	printf("Did connecter\n\n");
	fprintf(fpout,"Did connecter\n\n");
	fflush(fpout);

	nobjects = tracer(minvoxels, drop);
	printf("Did tracer: nobjects: %d\n",nobjects);

//	remover(drop);

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetInput(im);
	writer->UseCompressionOn();

	for (lab=1; lab<=cmax; lab++) 
	{
		cmap[lab] = lab;
		taken[lab] = false;
	}

	int nfiles = NBIGOBJECTS;
	if (nobjects < nfiles) {
		nfiles = nobjects;
	}
	double totvoxels = 0;
	for (int i=0; i<nfiles; i++)
	{
		totvoxels += bigobj[i].nvoxels;
	}
	for (int i=0; i<nfiles; i++)
	{
		lab = bigobj[i].label;
		createObjectImage(lab,cmap,taken);
		sprintf(numstr,"%02d",i);
		objectFile[0] = '\0';
		strcat(objectFile,outputPath);
		strcat(objectFile, baseName);
		strcat(objectFile, numstr);
		strcat(objectFile, ".tif");
		printf("%s n: %d  fraction: %f\n",objectFile,bigobj[i].nvoxels,bigobj[i].nvoxels/totvoxels);
		fprintf(fpout,"%s n: %d  fraction: %f\n",objectFile,bigobj[i].nvoxels,bigobj[i].nvoxels/totvoxels);

		writer->SetFileName(objectFile);
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
		printf("Created object image file: %s\n",objectFile);
	}

	printf("Elapsed time: %ld seconds \n",(long int)(time(NULL)-t1));

	return 0;
}

