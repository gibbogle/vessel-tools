/*
 * To floodfill holes in a binary tiff image, working from seed voxels
*/

#include <cstdio>
#include <vector>
/*
#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>
*/
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
//#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkSize.h"
//#include "itkThresholdImageFilter.h"
//#include "itkBinaryThresholdImageFilter.h"
//#include "itkMeanImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

unsigned char *p;

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]

#define MAXFILL 20000

long long width, height, depth, imsize;
typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;

struct xystr
{
	int x, y;
};
struct xyzstr
{
	int x, y, z;
};

typedef xystr XY;
typedef xyzstr XYZ;

XY oldlist[MAXFILL], newlist[MAXFILL], startlist[MAXFILL];
XYZ seedlist[100];
int nold, nseed, nnew, nstart;
int nbr2D[8][2];

void setupNeighbours()
{
	int x, y, k=0;

	for (x=-1; x<=1; x++)
	{
		for (y=-1; y<=1; y++)
		{
			if (x == 0 && y == 0) continue;
			nbr2D[k][0] = x;
			nbr2D[k][1] = y;
			k++;
		}
	}
	printf("Number of 2D neighbours: %d\n",k);
}

// Based on the list of voxels lit in the slice zold (oldlist[], nold), find adjacent unlit voxels
// in next slice (zold+1 or zold-1).
void newseeds(int zold, int znew)
{
	int x, y, z, k;

	z = znew;
	nnew = 0;
	for (k=0; k<nold; k++)
	{
		x = oldlist[k].x;
		y = oldlist[k].y;
		if (z == 41 && y == 777)
			printf("oldlist %d %d  V: %d\n",x,y,V(x,y,z));
		if (V(x,y,z) == 0)
		{
			V(x,y,z) = 255;
			newlist[nnew].x = x;
			newlist[nnew].y = y;
			if (z == 41 && y == 777)
				printf("%d  %d %d\n",nnew,x,y);
			nnew++;
		}
	}
}

void expand(int z)
{
	int i, k, x, y, x0, y0;

	printf("Expand: z: %d  nnew: %d\n",z,nnew);
	for (k=0; k<nnew; k++)
	{
		x0 = newlist[k].x;
		y0 = newlist[k].y;
//		if (z == 28) printf("k: %d seed: %d %d\n",k,x0,y0);
		for (i=0; i<8; i++)
		{
			x = x0 + nbr2D[i][0];
			if (x < 0 || x > width-1) continue;
			y = y0 + nbr2D[i][1];
			if (y < 0 || y > height-1) continue;
			if (V(x,y,z) != 0) continue;
			V(x,y,z) = 255;
			newlist[nnew].x = x;
			newlist[nnew].y = y;
			nnew++;
			if (z == 42 && y == 777) printf("added: %d %d %d  nnew: %d\n",x,y,z,nnew);
			if (nnew == MAXFILL-1)
			{
				printf("Too many seeds at slice z: %d.  Leak?\n",z);
				printf("Trying to add voxel: %d %d %d\n",x,y,z);
				exit(1);
			}
		}
	}
}

void setseeds()
{
	// First, test with seed (xseed,yseed,0)
	nseed = 0;
	seedlist[nseed].x = 444; seedlist[nseed].y = 770; seedlist[nseed].z =  45; nseed++;
	seedlist[nseed].x = 254; seedlist[nseed].y = 820; seedlist[nseed].z =  79; nseed++;
	seedlist[nseed].x = 572; seedlist[nseed].y = 547; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x = 139; seedlist[nseed].y = 155; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x =  29; seedlist[nseed].y = 758; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x = 161; seedlist[nseed].y = 868; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x = 360; seedlist[nseed].y = 828; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x = 359; seedlist[nseed].y = 750; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x = 514; seedlist[nseed].y = 758; seedlist[nseed].z =   0; nseed++;
	seedlist[nseed].x = 606; seedlist[nseed].y = 746; seedlist[nseed].z =   0; nseed++;
}

int main(int argc, char**argv)
{
	int z, k, iseed, zstart, zlast;

	if (argc < 3) {
		printf("Usage: fill input_tiff filled_tiff seedfile\n");
		return 0;
	}

	printf("This program is currently using hard-coded seed locations\n");
	return 1;

	printf("Input image stack file: %s\n",argv[1]);
	printf("Output image stack file: %s\n",argv[2]);

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
		return 1;
	}

	im = reader->GetOutput();
	std::cout << "Image loaded";

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	p = (unsigned char *)(im->GetBufferPointer());

	setupNeighbours();
	// The proposed method cheats by filling zslice by zslice, starting with the seed's z value.
	// The voxels filled in a slice are recorded in a list.  On moving to the next slice,
	// dark voxels in the new slice that are adjacent to those in the list are used as seeds.
	// For this purpose we can use face-adjacency.

	setseeds();
//	if (zstart != 0)
//	{
//		for (k=0; k<nnew; k++)
//			startlist[k] = newlist[k];
//		nstart = nnew;
//	}

	zlast = depth-1;
//	zlast = 61;
	for (iseed=0; iseed<nseed; iseed++)
	{
		printf("Seed: %d  %d %d %d\n",iseed,seedlist[iseed].x,seedlist[iseed].y,seedlist[iseed].z);
		nnew = 1;
		newlist[0].x = seedlist[iseed].x;
		newlist[0].y = seedlist[iseed].y;
		zstart = seedlist[iseed].z;
		for (z=zstart; z<=zlast; z++)
		{
			expand(z);
			if (z == zlast) break;
			for (k=0; k<nnew; k++) {
				oldlist[k] = newlist[k];
				if (z == zstart) 
					startlist[k] = newlist[k];
			}
			if (z == zstart)
				nstart = nnew;
			nold = nnew;
			newseeds(z,z+1);
		}

		if (zstart != 0)
		{
			for (k=0; k<nstart; k++)
				newlist[k] = startlist[k];
			nnew = nstart;
			zstart = seedlist[iseed].z;
			for (z=zstart; z>=0; z--)
			{
				expand(z);
				if (z == 0) break;
				for (k=0; k<nnew; k++)
					oldlist[k] = newlist[k];
				nold = nnew;
				newseeds(z,z-1);
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
		return 1;
	}
	printf("Created filled image file: %s\n",argv[2]);

	return 0;
}