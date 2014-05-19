/*
 * To smooth a .tif by simple mean filter 
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
#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkSize.h"
#include "itkThresholdImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMeanImageFilter.h"
//#include "itkBoxMeanImageFilter.h"

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
long long width, height, depth, imsize;
unsigned char *p, *average;
bool dbug;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]
#define A(a,b,c)  average[(c)*imsize+(b)*width+(a)]

#define USE_ITK_FILTER false
#define USE_COMPRESSION true
#define NTHREADS 2

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------
void adder(int xmin, int xmax,int ymin, int ymax, int zmin, int zmax, double *ave)
{
	int xx, yy, zz, cnt;
	double sum;

	sum = 0;
	cnt = 0;
	for (xx=xmin; xx<=xmax; xx++)
	{
		for (yy=ymin; yy<=ymax; yy++)
		{
			for (zz=zmin; zz<=zmax; zz++)
			{
//				if (dbug && yy == 800)
//					printf("%d %d %d  %d\n",xx,yy,zz,V(xx,yy,zz));
				sum += V(xx,yy,zz);
				cnt++;
			}
		}
	}
	*ave = sum/cnt;
//	if (dbug)
//		printf("%d %d  %d %d  %d %d  sum: %f cnt: %d ave: %f\n",xmin,xmax,ymin,ymax,zmin,zmax,sum,cnt,*ave);
}

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
void fastsmth(int R, int npar)
{
	int N, x, y, z, xx, dy, dz, k, sum, sumYZ[64];
	int z1, z2, kpar,dzpar;
	int xmin,xmax,ymin,ymax,zmin,zmax;
	double ave, sumXYZ;

	N = 2*R + 1;
	dbug = false;
	// First handle the "interior" points
	printf("interior: ");
	dzpar = (depth-2*R)/npar;
    #pragma omp parallel for num_threads(npar) private(z1, z2, x, y, z, sumXYZ, xx, sumYZ, dy, dz, sum, k)
	for (kpar=0; kpar<npar; kpar++) 
	{
		z1 = R + kpar*dzpar; 
		if (kpar < npar-1) {
			z2 = z1+dzpar-1;
		} else {
			z2 = depth-R-1;
		}
		printf("Started thread: %d  z1,z2: %d %d\n",kpar,z1,z2);
//	for (z=R; z<=depth-R-1; z++)
	for (z=z1; z<=z2; z++)
	{
		printf(".");
		if (dbug) printf("z: %d\n",z);
		for (y=R; y<=height-R-1; y++)
		{
			if (dbug) printf("y: %d\n",y);
			sumXYZ = 0;				//sumXYZ is the cube sum
			for (xx=0; xx<N; xx++)
			{
				sumYZ[xx] = 0;		// sumYZ[] is the vector of YZslice sums
				for (dy=-R; dy<=R; dy++)
				{
					for (dz=-R; dz<=R; dz++)
					{
						if (dbug) printf("%d %d %d   %lld\n",xx,y+dy,z+dz,(z+dz)*imsize+(y+dy)*width+(xx));
						sumYZ[xx] += V(xx,y+dy,z+dz);
					}
				}
				sumXYZ += sumYZ[xx];
			}
			if (dbug) printf("y: %d sumXYZ: %d\n",y,sumXYZ);
			for (x=R; x<=width-R-1; x++)
			{
				if (dbug) printf("x,R: %d %d\n",x,R);
				if (x > R)
				{
					sum = 0;
					for (dy=-R; dy<=R; dy++)
					{
						for (dz=-R; dz<=R; dz++)
						{
							if (dbug) printf("%d %d %d   %lld\n",x+R,y+dy,z+dz,(z+dz)*imsize+(y+dy)*width+(x+R));
							sum += V(x+R,y+dy,z+dz);
						}
					}
					sumXYZ += sum - sumYZ[0];
					for (k=1; k<N; k++)
						sumYZ[k-1] = sumYZ[k];
					sumYZ[N-1] = sum;
				}
				A(x,y,z) = (unsigned char)(sumXYZ/(N*N*N) + 0.5);
			}
		}
	}
	}
	// Now deal with the "edge" points
	printf("\nedge: ");

	dzpar = depth/npar;
    #pragma omp parallel for num_threads(npar) private(z1, z2, x, y, z, xmin, xmax, ymin, ymax, zmin, zmax, ave)
	for (kpar=0; kpar<npar; kpar++) 
	{
		z1 = kpar*dzpar; 
		if (kpar < npar-1) {
			z2 = z1+dzpar-1;
		} else {
			z2 = depth-1;
		}
		printf("Started thread: %d  z1,z2: %d %d\n",kpar,z1,z2);
//	for (z=0; z<depth; z++)
	for (z=z1; z<=z2; z++)
	{
		printf(".");
		if (z >= R && z <= depth-R-1)
		{
			zmin = z-R;
			zmax = z+R;
			// There are 4 pieces
			for (y=0; y<height; y++)
			{
				ymin = MAX(y-R,0);
				ymax = MIN(y+R,height-1);
				for (x=0; x<R; x++)
				{
					xmin = MAX(x-R,0);
					xmax = x+R;
					adder(xmin,xmax,ymin,ymax,zmin,zmax,&ave);
					A(x,y,z) = (unsigned char)(ave + 0.5);
				}
				for (x=width-R; x<width; x++)
				{
					xmin = x-R;
					xmax = MIN(x+R,width-1);
					adder(xmin,xmax,ymin,ymax,zmin,zmax,&ave);
					A(x,y,z) = (unsigned char)(ave + 0.5);
				}
			}
			for (x=R; x<width-R; x++)
			{
				xmin = x-R;
				xmax = x+R;
				for (y=0; y<R; y++)
				{
					ymin = MAX(y-R,0);
					ymax = y+R;
					adder(xmin,xmax,ymin,ymax,zmin,zmax,&ave);
					A(x,y,z) = (unsigned char)(ave + 0.5);
				}
				for (y=height-R; y<height; y++)
				{
					ymin = y-R;
					ymax = MIN(y+R,height-1);
					adder(xmin,xmax,ymin,ymax,zmin,zmax,&ave);
					A(x,y,z) = (unsigned char)(ave + 0.5);
				}
			}
		}
		else
		{
			if (z < R)
			{
				zmin = 0;
				zmax = z+R;
			}
			else
			{
				zmin = z-R;
				zmax = depth-1;
			}
			for (y=0; y<height; y++)
			{
				ymin = MAX(y-R,0);
				ymax = MIN(y+R,height-1);
				for (x=0; x<width; x++)
				{
					xmin = MAX(x-R,0);
					xmax = MIN(x+R,width-1);
					adder(xmin,xmax,ymin,ymax,zmin,zmax,&ave);
					A(x,y,z) = (unsigned char)(ave + 0.5);
				}
			}
		}
	}
	}
	printf("\n");
}


int main(int argc, char**argv)
{
	int radius, npar;
	long long longsize;
	time_t t1;
	t1 = time(NULL);
	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	if (USE_ITK_FILTER)
		printf("Using itk:::BoxMeanImageFilter with %d threads\n",NTHREADS);
	else
		printf("Using fastsmth\n");
	if (USE_COMPRESSION)
		printf("Compression on\n");
	else
		printf("Compression off\n");

	if (argc != 5) {
		printf("Usage: smooth input_tiff output_tiff radius ncpu\n");
		printf("       where: radius is the radius of the mean smoothing filter\n");
		printf("              ncpu is the number of OpenMP threads\n");
		fprintf(fp,"Usage: smooth input_tiff output_tiff radius ncpu\n");
		fprintf(fp,"       where: radius is the radius of the mean smoothing filter\n");
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
	sscanf(argv[3],"%d",&radius);
	printf("Averaging radius R: %d\n",radius);
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
	longsize = width;
	longsize *= height;
	longsize *= depth;
	printf("Array size: %lld\n",longsize);

	if (!USE_ITK_FILTER)
	{
		p = (unsigned char *)(im->GetBufferPointer());
		average = (unsigned char*)malloc(imsize*depth*sizeof(unsigned char));
		fastsmth(radius, npar);
		memcpy(p,average,imsize*depth);
	} 
	/*
	else
	{
		typedef itk::BoxMeanImageFilter< ImageType, ImageType > S_FilterType;
		S_FilterType::Pointer sfilter = S_FilterType::New();
		ImageType::SizeType indexRadius;
		indexRadius[0] = radius; // radius along x
		indexRadius[1] = radius; // radius along y
		indexRadius[2] = radius; // radius along z
		sfilter->SetRadius( indexRadius );
		sfilter->SetInput( im );
		sfilter->SetNumberOfThreads(NTHREADS);
		try
		{
			sfilter->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cout << e << std::endl;
			return 1;
		}
		im = sfilter->GetOutput();
	}
	*/
	printf("Averaging completed\n");

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
	printf("Created smoothed image file: %s\n",argv[2]);
	printf("Elapsed time: %ld seconds \n",(long int)(time(NULL)-t1));

	return 0;
}