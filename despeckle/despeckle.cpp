/*
 * Despeckle a .tif 
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

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
int width, height, depth, imsize;
unsigned char *p, *despec;
bool dbug;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]
#define D(a,b,c)  despec[(c)*imsize+(b)*width+(a)]

#define USE_COMPRESSION true
#define NBR 2

//----------------------------------------------------------------------------------------------
// Nomenclature:
//   cube = NxNxN block that we are averaging over
//   XYslice = slice of the cube in the xy plane
//   YZslice = slice of the cube in the yz plane
//   ZXslice = slice of the cube in the zx plane
// To keep a cube of side N (N = 2*R+1) completely within the image boundaries,
// R-1 <= x <= width-R, etc.  This defines the "interior".
// For interior voxels:
// for each z:
//    for each y:
//        the N YZslice sums are computed for the starting x value (x=R).  The sum of these is the cube sum.
//        for each x:
//            the slice sum for the leading face YZslice is added to the cube sum, and trailing YZslice is dropped.
//            the YZslice vector is updated.
// In this way at each (x,y,z) to update the sum requires only NxN voxel values to be accessed.
// For the remaining "edge" voxels the sum is computed directly, with no attempt to optimize.
//----------------------------------------------------------------------------------------------
void fastdespec(int nlit, int npar)
{
	int x, y, z, xx, yy, zz, dx, dy, dz, k, sum;
	int z1, z2, kpar,dzpar;
	int xmin,xmax,ymin,ymax,zmin,zmax;

	dbug = false;
	dzpar = depth/npar;
    #pragma omp parallel for num_threads(npar) private(z1, z2, x, y, z, xx, yy, zz, dx, dy, dz, sum)
	for (kpar=0; kpar<npar; kpar++) {
		z1 = kpar*dzpar; 
		if (kpar < npar-1) {
			z2 = z1+dzpar-1;
		} else {
			z2 = depth-1;
		}
		printf("Started thread: %d  z1,z2: %d %d\n",kpar,z1,z2);
		for (z=z1; z<=z2; z++) {
			printf(".");
			for (y=0; y<height; y++) {
				for (x=0; x<width; x++) {
					if (V(x,y,z) != 255) {
						D(x,y,z) = 0;
						continue;
					}
					sum = 0;
					for (dx=-NBR; dx<=NBR; dx++) {
						xx = x + dx;
						if (xx < 0 || xx >= width) continue;
						for (dy=-NBR; dy<=NBR; dy++) {
							yy = y + dy;
							if (yy < 0 || yy >= height) continue;
							for (dz=-NBR; dz<=NBR; dz++) {
								zz = z + dz;
								if (zz < 0 || zz >= depth) continue;
								if (dx == 0 && dy == 0 && dz == 0) continue;
								if (V(xx,yy,zz) == 255) sum++;
							}
						}
					}
					if (sum >= nlit) {
						D(x,y,z) = 255;
					} else {
						D(x,y,z) = 0;
//						V(x,y,z) = 0;
					}
				}
			}
		}
	}
	printf("\n");
}


//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int nlit, npar;
	time_t t1;
	t1 = time(NULL);
	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	if (USE_COMPRESSION)
		printf("Compression on\n");
	else
		printf("Compression off\n");

	if (argc != 5) {
		printf("Usage: despeckle input_tiff output_tiff nlit ncpu\n");
		printf("       where: nlit is the threshold number of lit neighbour voxels (of 26)\n");
		printf("              ncpu is the number of OpenMP threads\n");
		fprintf(fp,"Usage: despeckle input_tiff output_tiff nlit ncpu\n");
		fprintf(fp,"       where: nlit is the threshold number of lit neighbour voxels (of 26)\n");
		fprintf(fp,"              ncpu is the number of OpenMP threads\n");
		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fp);
		return 1;
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	sscanf(argv[3],"%d",&nlit);
	printf("Threshold lit neighbour count nlit: %d\n",nlit);
	sscanf(argv[4],"%d",&npar);
	printf("NCPU: %d\n",npar);

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
		return 2;
	}

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	imsize = width*height;

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	printf("Array size: %d\n",width*height*depth);

	if (depth == 1) {
		printf("Bad image file\n");
		return 2;
	}

	p = (unsigned char *)(im->GetBufferPointer());
	despec = (unsigned char*)malloc(imsize*depth*sizeof(unsigned char));
	fastdespec(nlit, npar);
	memcpy(p,despec,imsize*depth);

	printf("Despeckling completed\n");

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im);
	if (USE_COMPRESSION)
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
		return 3;
	}
	printf("Created despeckled image file: %s\n",argv[2]);
	printf("Elapsed time: %ld seconds \n",(long int)(time(NULL)-t1));

	return 0;
}