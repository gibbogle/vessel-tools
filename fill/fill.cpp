/*
 * To fill holes in a binary tiff image
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
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

unsigned char *p;

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]
#define KMAX_THRESHOLD 90
#define epsilonw 1.0

#define NPROBES 26
int probe[NPROBES][3];
double unitprobe[NPROBES][3];
double probeLen, projLimit;
double vsize[3];

struct cpoint
{
	short x, y, z;
	short nhits;
	unsigned char val;
};
typedef cpoint CPOINT;

struct litvoxel{
	int x,y,z;
	int plane;
};
typedef litvoxel LITVOXEL;

CPOINT *candidate;
int maxcand, ncand;
int candHitLimit, insideHitLimit;
int width, height, depth;
int imsize;
typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;

#define WINDOW 5
#define MAX_SLICE (2*WINDOW+1)*(2*WINDOW+1)*(2*WINDOW+1)
#define NPLANES 25		// was 21

double planes[NPLANES][3]= {
	{1, 0, -1},
	{1, 0, 1},
	{1, -1, 0},
	{1, 1, 0},
	{0, 1, 1},
	{0, -1, 1},
	{1, 0, 0},
	{0, 1, 0},
	{0, 0, 1},
	{0.5, 0.5, 1},
	{0.5, -0.5, 1},
	{-0.5, -0.5, 1},
	{-0.5, 0.5, 1},
	{1, 0.5, 0.5},
	{1, 0.5, -0.5},
	{1, -0.5, 0.5},
	{1, -0.5, -0.5},
	{0.5, 1, 0.5},
	{-0.5, 1, 0.5},
	{0.5, 1, -0.5},
	{-0.5, 1, -0.5},
	{1, 1, 1},
	{1, -1, 1},
	{-1, 1, 1},
	{-1, -1, 1}
};

int nfore_pts[NPLANES], nback_pts[NPLANES];
int fore_pt[NPLANES][MAX_SLICE][3];
int back_pt[NPLANES][MAX_SLICE][3];

bool HEAL=true;

//----------------------------------------------------------------------------------------------------------------
// The 25 planes are defined by their unit normals.
// For each plane, we look at a cube of grid points 2*WINDOW+1 on the side.
// Grid points that are within epsilonw (1.5) of the plane are recorded as "foreground" in fore_pt[][][],
// while the rest are "background" and are stored in back_pt[][][].  The total numbers of such points are
// nfore_pts[] and nbak_pts[] respectively.
// In effect, "foreground" points are those that lie in (are close to) a planar slice through (0,0,0).
//----------------------------------------------------------------------------------------------------------------
int SetupSlicePts()
{
	for(int k =0; k < NPLANES; k++)
	{
		double sum  = sqrt(planes[k][0]*planes[k][0] + planes[k][1]*planes[k][1] + planes[k][2]*planes[k][2]);
		planes[k][0]/=sum;
		planes[k][1]/=sum;
		planes[k][2]/=sum;
	}
	for(int ip = 0; ip < NPLANES; ip++)
	{
		double a = planes[ip][0];
		double b = planes[ip][1];
		double c = planes[ip][2];
		int kf = 0;
		int kb = 0;
		for(int wz =-WINDOW;wz<=WINDOW;wz++)
		{
			for(int wy=-WINDOW;wy<=WINDOW;wy++)
			{
				for(int wx=-WINDOW;wx<=WINDOW;wx++)
				{
					if(fabs(wz*a+wy*b+wx*c) <= epsilonw)
					{
						kf++;
						fore_pt[ip][kf][0] = wx;
						fore_pt[ip][kf][1] = wy;
						fore_pt[ip][kf][2] = wz;
					}
					else {
						kb++;
						back_pt[ip][kb][0] = wx;
						back_pt[ip][kb][1] = wy;
						back_pt[ip][kb][2] = wz;
					}
				}
			}
		}
		nfore_pts[ip] = kf;
		nback_pts[ip] = kb;
	}
	return 0;
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
int probeSetup(void)
{
	int dir, x, y, z;
	double d, proj, projmin;

	dir = 0;
	for (x=-1;x<=1; x++)
	{
		for (y=-1;y<=1; y++)
		{
			for (z=-1;z<=1; z++)
			{
				if (x==0 && y==0 && z==0)
					continue;
				probe[dir][0] = x;
				probe[dir][1] = y;
				probe[dir][2] = z;
				d = sqrt(double(x*x+y*y+z*z));
				unitprobe[dir][0] = x/d;
				unitprobe[dir][1] = y/d;
				unitprobe[dir][2] = z/d;
				dir++;
			}
		}
	}
	projmin = 10.0;
	for (x=0; x<10; x++)
	{
		for (y=0; y<10; y++)
		{
			for (z=0; z<10; z++)
			{
				d = x*x + y*y + z*z;
				if (d == 0)
					continue;
				d = sqrt(d);
				double projmax = 0;
				for (dir=0; dir<NPROBES; dir++)
				{
					proj = (x*unitprobe[dir][0] + y*unitprobe[dir][1] + z*unitprobe[dir][2])/d;
					if (proj > projmax)
						projmax = proj;
				}
				if (projmax < projmin)
					projmin = projmax;
			}
		}
	}
	projLimit = 0.99*projmin;
	printf("projLimit: %f\n",projLimit);
	return 0;
}

//----------------------------------------------------------------------------------------------------------------
// To account for varying vsize[:] (non-cubic voxel shape) need to use double instead of int
//----------------------------------------------------------------------------------------------------------------
//int prober(int x, int y, int z, int nhlim, int *nhits, int *dirmin, int *d2min, int d2hit[])
int prober(int x, int y, int z, int nhlim, int *nhits, int *dirmin, double *d2min, double d2hit[])
{
	int dx, dy, dz, dir, k, xx, yy, zz, n;	// dd
	double dd;
	double probeLen2 = probeLen*probeLen;

	n = 0;
	*dirmin = -1;
	*d2min = 99999;
	for (dir=0; dir<NPROBES; dir++)
	{
		dx = probe[dir][0];
		if (x+dx < 0 || x+dx >= width) continue;
		dy = probe[dir][1];
		if (y+dy < 0 || y+dy >= height) continue;
		dz = probe[dir][2];
		if (z+dz < 0 || z+dz >= depth) continue;
//		dd = dx*dx + dy*dy + dz*dz;
		dd = dx*dx*vsize[0]*vsize[0] + dy*dy*vsize[1]*vsize[1] + dz*dz*vsize[2]*vsize[2];
		d2hit[dir] = 0;
		for (k=1; k*k*dd<probeLen2; k++)
		{
			xx = x+k*dx;
			if (xx < 0 || xx >= width) continue;
			yy = y+k*dy;
			if (yy < 0 || yy >= height) continue;
			zz = z+k*dz;
			if (zz < 0 || zz >= depth) continue;
//			printf("k: %d V: %d\n",k,V(xx,yy,zz));
			if (V(xx,yy,zz) != 0)
			{
				n++;
				d2hit[dir] = k*dd;
				if (d2hit[dir] < *d2min)
				{
					*dirmin = dir;
					*d2min = d2hit[dir];
				}
				break;
			}
		}
		if (n == nhlim)
			break;
	}
	*nhits = n;
//	printf("nhits: %d\n",*nhits);
	return 0;
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
int setupCandidateList(void)
{
	int x0, y0, z0, dirmin, nhits, nhlim;	// , d2min
//	int d2hit[NPROBES];
	double d2hit[NPROBES], d2min;

	nhlim = candHitLimit;
	ncand = 0;
	for (z0=0; z0<depth; z0++)
	{
		printf("z0: %d ncand: %d\n",z0,ncand);
		for (y0=0; y0<height; y0++)
		{
			for (x0=0; x0<width; x0++)
			{
				if (V(x0,y0,z0) == 255) continue;
				prober(x0,y0,z0,nhlim,&nhits,&dirmin,&d2min,d2hit);
				if (nhits >= candHitLimit)
				{
					candidate[ncand].x = x0;
					candidate[ncand].y = y0;
					candidate[ncand].z = z0;
					candidate[ncand].nhits = nhits;
					candidate[ncand].val = 0;
					ncand++;
					if (ncand == maxcand)
					{
						printf("Error: candidate list size exceeded\n");
						return 1;
					}
				}
			}
		}
	}
	return 0;
}

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
//int filler(int x0, int y0, int z0, int d2min, int *d2hit, int flist[][3])
int filler(int x0, int y0, int z0, double d2min, double *d2hit, int flist[][3])
{
	int dx, dy, dz, dir, x, y, z, n;	// d, dx2, dy2, dz2, d2,
	double d, dx2, dy2, dz2, d2;
	double sd[3];
	int k, id[3];
	double dlen, proj;
	bool inside;
	unsigned char in = 255;
//	ImageType::IndexType index;


//	d = int(sqrt(double(d2min)) + 0.5);
	d = sqrt(d2min);
	for (k=0; k<3; k++) {
		id[k] = d/vsize[k] + 0.5;
		sd[k] = vsize[k]*vsize[k];
	}
	n = 0;
//	for (dx=-d; dx<=d; dx++)
	for (dx=-id[0]; dx<=id[0]; dx++)
	{
		x = x0 + dx;
		if (x < 0 || x >= width)
			continue;
		dx2 = dx*dx*sd[0];
//		for (dy=-d; dy<=d; dy++)
		for (dy=-id[1]; dy<=id[1]; dy++)
		{
			y = y0 + dy;
			if (y < 0 || y >= height)
				continue;
			dy2 = dy*dy*sd[1];
			if (dx2 + dy2 > d2min)
				continue;
//			for (dz=-d; dz<=d; dz++)
			for (dz=-id[2]; dz<=id[2]; dz++)
			{
				z = z0 + dz;
				if (z < 0 || z >= depth)
					continue;
				dz2 = dz*dz*sd[2];
				d2 = dx2 + dy2 + dz2;

				if (d2 > d2min)
					continue;
				if (V(x,y,z) != 0)
					continue;

				if (d2 == 0)
				{
					inside = true;
				} 
				else
				{
//					dlen = sqrt(double(d2));
					dlen = sqrt(d2);
					inside = false;
					for (dir=0; dir<NPROBES; dir++)
					{
						if (d2hit[dir] == 0)
							continue;
						proj = (unitprobe[dir][0]*dx*vsize[0] + unitprobe[dir][1]*dy*vsize[1] 
								+ unitprobe[dir][2]*dz*vsize[2])/dlen;
						if (proj > projLimit)
						{
							inside = true;
							break;
						}
					}
				}

				if (inside) {
					V(x,y,z) = in;
					flist[n][0] = x;
					flist[n][1] = y;
					flist[n][2] = z;
					n++;
				}
			}
		}
	}
	return n;
}

//------------------------------------------------------------------------------------------------
// Connected objects with fewer than the threshold voxels are removed (set to 0).
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


//----------------------------------------------------------------------------------------------------------------
// Method:
// First a list of all candidate sites for filling is created.
// A voxel is in the list if:
// (a) It is not lit (V(x,y,z) = 0)
// (b) The number of probe hits is >= candHitLimit (currently = 0.5*insideHitLimit)
//     There are 26 probes and a hit occurs if a probe encounters a lit voxel less than probeLen distant
//
// The candidate list is traversed several times (nsweeps).  Each candidate voxel is tested by counting probe hits.
// The prober function counts the number of hits and records the distance to a lit voxel in each probe direction.
// If the hit count is >= insideHitLimit, the candidate voxel is deemed to be "inside", and the filler function is invoked.
//
// The filler function lights voxels within a spherical neighbourhood of the inside voxel (x0,y0,z0), subject to:
// (a) The radius of the sphere is the length of the shortest probe that hits
// (b) Each voxel (x,y,z) within the sphere is tested for a measure of "insideness", using info from the probes.
//     To be considered "inside", the unit vector of the offset of (x,y,z) from (x0,y0,z0) must have a projection
//     onto one of the normalized hitting probes that exceeds a threshold value  projLimit (currently about 0.86?)
//  These criteria enable filling of surface pits without overfilling.
//  All inside voxels are lit (V(x,y,z) = 255), and marked as lit in the candidate list.
//----------------------------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int niter;
	int x, y, z, xx, yy, zz, dirmin, nhits, sweep, n, nadded, nhlim, d, err, iter;	// d2min, 
//	int d2hit[NPROBES];
	double d2min, d2hit[NPROBES];
	int flist[10000][3];
//	int nflist;
	int nsweeps = 4;
	double voxelsize_xy, voxelsize_z;
//	bool nearedge;
	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	probeSetup();

	if (argc != 8) {
		printf("Usage: fill input_tiff filled_tiff voxelsize_xy voxelsize_z probeLen insideHitLimit niter\n");
		printf("   where:   voxelsize_xy is the x-y voxel size (um)\n");
		printf("            voxelsize_z is the z voxel size (um)\n");
		printf("            probeLen is the probe length parameter (um)\n");
		printf("            insideHitLimit is the algorithm decision parameter\n");
		printf("            niter is the number of hole-filling iterations\n");
		return 1;
	}

	printf("Input image stack file: %s\n",argv[1]);
	printf("Output image stack file: %s\n",argv[2]);
	sscanf(argv[3],"%lf",&voxelsize_xy);
	sscanf(argv[4],"%lf",&voxelsize_z);
	sscanf(argv[5],"%lf",&probeLen);
	sscanf(argv[6],"%d",&insideHitLimit);
	sscanf(argv[7],"%d",&niter);
	vsize[0] = voxelsize_xy;
	vsize[1] = voxelsize_xy;
	vsize[2] = voxelsize_z;
	printf("Hole-filling parameters: probeLen: %lf insideHitLimit: %d iterations: %d\n",probeLen,insideHitLimit,niter);
	if (niter < 0) {
		printf("niter < 0: turning off Healing sweep\n");
		HEAL = false;
		niter = -niter;
	}
	candHitLimit = 0.5*insideHitLimit;
	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

	reader->SetFileName(argv[1]);
	try
	{
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
	std::cout << "Image loaded";

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	p = (unsigned char *)(im->GetBufferPointer());

	despeckle();

	maxcand = 0.3*width*height*depth;
//	maxcand = 500000;
	candidate = (CPOINT *)malloc(maxcand*sizeof(CPOINT));
	printf("Allocated candidate\n");
	SetupSlicePts();

	for (iter=0; iter<niter; iter++)
	{
		err = setupCandidateList();
		if (err != 0) {
			printf("Failed\n");
			return 4;
		}
		printf("Number of inside candidates: %d\n",ncand);

		nhlim = insideHitLimit;
		for (sweep=0; sweep<nsweeps; sweep++)
		{
			printf("sweep: %d\n",sweep);
			nadded = 0;
			for (int icand=0; icand<ncand; icand++)
			{
				if (candidate[icand].val != 0)
					continue;
				x = candidate[icand].x;
				y = candidate[icand].y;
				z = candidate[icand].z;
				if (V(x,y,z) != 0) {
					candidate[icand].val = 255;
					continue;
				}
				prober(x,y,z,nhlim,&nhits,&dirmin,&d2min,d2hit);
				if (x < 8 || y < 8 || z < 8 || x > width-1-8 || y > height-1-8 || z > depth-1-8)
				{
					int ne = 6;		// range of edge adjustment
					d = ne;
					if (x < ne) 
						d = MIN(d,x);
					else if (x > width-1-ne)
						d = MIN(d,width-1-x);
					if (y < ne) 
						d = MIN(d,y);
					else if (y > height-1-ne)
						d = MIN(d,height-1-y);
					if (z < ne) 
						d = MIN(d,z);
					else if (z > depth-1-ne)
						d = MIN(d,depth-1-z);
					d = MIN(d,ne);
					nhlim =  insideHitLimit - (ne-d);
				}
				else
				{
					nhlim = insideHitLimit;
				}
				if (nhits >= nhlim)
				{
					n = filler(x,y,z,d2min,d2hit,flist);
					candidate[icand].val = 255;
					nadded += n;
					/*
//					printf("filled: %d %d %d added: %d\n",x,y,z,n);
					for (int k=0; k<n; k++)
					{
						x = flist[k][0];
						y = flist[k][1];
						z = flist[k][2];
						// Look at the neighbours of this newly lit voxel
						for (int i=0; i<26; i++)
						{
							xx = x + probe[i][0];
							if (xx < 0 || xx >= width) continue;
							yy = y + probe[i][1];
							if (yy < 0 || yy >= height) continue;
							zz = z + probe[i][2];
							if (zz < 0 || zz >= depth) continue;
							if (V(xx,yy,zz) == 0)
							{
								prober(xx,yy,zz,nhlim,&nhits,&dirmin,&d2min,d2hit);
								if (nhits >= candHitLimit)
								{
									candidate[ncand].x = xx;
									candidate[ncand].y = yy;
									candidate[ncand].z = zz;
									candidate[ncand].nhits = 0;
									candidate[ncand].val = 0;
									ncand++;
									printf("added: %d %d %d  ncand: %d\n",xx,yy,zz,ncand);
									if (ncand == maxcand)
									{
										printf("Error: candidate list size exceeded\n");
										return 1;
									}
								}
							}
						}
					}
					*/
				}
			}
			printf("Added: %d\n",nadded);
		}
		if (HEAL)
		{
			printf("Healing...\n");
			for (sweep=0; sweep<1; sweep++)
			{
				printf("Healing sweep: %d\n",sweep);
				// Healing code
				nadded = 0;
				for (int icand=0; icand<ncand; icand++)
				{
					if (icand%100000 == 0)
						printf("icand: %d z: %d  nadded: %d\n",icand,candidate[icand].z,nadded);
					if (candidate[icand].val != 0)
						continue;
					x = candidate[icand].x;
					y = candidate[icand].y;
					z = candidate[icand].z;
					if (V(x,y,z) != 0) continue;
					int kmax = 0;
					for (int ip = 0; ip < NPLANES; ip++)
					{
						int knt = 0;
						for (int kf=0; kf<nfore_pts[ip]; kf++)
						{
							int wx = fore_pt[ip][kf][0];
							xx = x+wx;
							if (xx < 0 || xx >= width) continue;
							int wy = fore_pt[ip][kf][1];
							yy = y+wy;
							if (yy < 0 || yy >= height) continue;
							int wz = fore_pt[ip][kf][2];
							zz = z+wz;
							if (zz < 0 || zz >= depth) continue;
							if (V(xx,yy,zz) > 0)
								knt ++;
						}
						if (knt > kmax)
						{
							kmax = knt;
//							best_plane = ip;
						}
					}
//					C(coz,coy,cox) = false;
					if (kmax > 2*epsilonw*KMAX_THRESHOLD)
					{
//						if (queuep.size() == queuep_max)
//						{
//							printf("Error: queuep_max exceeded\n");
//							return 1;
//						}
//						C(coz,coy,cox) = true;
						V(x,y,z) = 255;
						nadded++;
//						data d;
//						d.x = cox;
//						d.y = coy;
//						d.z = coz;
//						d.fco = best_fco;
//						queuep.push_back(d);
					}
				}
				printf("Added: %d\n",nadded);
/*
					// Count lit neighbours, determine best plane only if there are enough
					n = 0;
					for (int i=0; i<26; i++)
					{
						xx = x + probe[i][0];
						if (xx < 0 || xx >= width) continue;
						yy = y + probe[i][1];
						if (yy < 0 || yy >= height) continue;
						zz = z + probe[i][2];
						if (zz < 0 || zz >= depth) continue;
						if (V(xx,yy,zz) != 0)
							n++;
					}
*/
			}
		}
	}

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im);
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
	printf("Created filled image file: %s\n",argv[2]);

	return 0;
}