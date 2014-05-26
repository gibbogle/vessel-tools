/*
 * To chop out blobby pieces of a .tif
 * Note: now range indices are 0-based to make them conform with the coordinates in a .am file.
 * This implies a divergence from Irfanview, which uses 1-based indices.
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

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im;
unsigned char *p, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V2(a,b,c)  p2[(c)*xysize+(b)*width+(a)]

long long width, height, depth, xysize;
int stencil_size;
int delta;
double threshold_fraction;
char stencil_shape = 'C';	// C = cube, S = sphere

int check(int x1, int x2, int y1, int y2, int z1, int z2)
{
	int x, y, z;
	int nlimit = threshold_fraction*(x2-x1+1)*(y2-y1+1)*(z2-z1+1);
	bool over = false;
	int count = 0;
	for (x=x1; x<=x2; x++) {
		for (y=y1; y<=y2; y++) {
			for (z=z1; z<=z2; z++) {
				if (V(x,y,z) > 0) {
					count++;
					if (count > nlimit) {
						over = true;
						goto DONE;
					}
				}
			}
		}
	}
DONE: if (!over) return 0;
	for (x=x1; x<=x2; x++) {
		for (y=y1; y<=y2; y++) {
			for (z=z1; z<=z2; z++) {
				V2(x,y,z) = 0;
			}
		}
	}
	return 1;
}


void process(int zstart, int zend)
{
	int radius, x0, y0, z0, x1, y1, z1, x2, y2, z2;
	radius = stencil_size/2;
	x0 = radius;
	y0 = radius;
	z0 = zstart + radius;
	printf("z0: %4d\n",z0);
	for (;;) {
		x1 = x0 - radius;
		y1 = y0 - radius;
		z1 = z0 - radius;
		x2 = x1 + stencil_size + 1;
		y2 = y1 + stencil_size + 1;
		z2 = z1 + stencil_size + 1;
		int res = check(x1,x2,y1,y2,z1,z2);
		if (res != 0) printf("Chopped: x: %4d %4d y: %4d %4d z: %4d %4d\n",x1,x2,y1,y2,z1,z2);
		x0 += delta;
		if (x0-radius+stencil_size+1 > width-1) {
			x0 = radius;
			y0 += delta;
			if (y0-radius+stencil_size+1 > height-1) {
				y0 = radius;
				z0 += delta;
				printf("z0: %4d\n",z0);
				if (z0-radius+stencil_size+1 > zend) break;
			}
		}
	}

}

int main(int argc, char**argv)
{
	int zstart, zend;
	bool use_compression;
	int zparstart[8], zparend[8];
	int kpar, npar = 8;

	if (argc != 7) {
		printf("Usage: deblob input_tiff output_tiff stencil_size threshold_fraction delta ncpu\n");
		printf("       where stencil_size is the cube side length or sphere diameter in voxels\n");
		printf("       threshold_fraction is fraction of voxels that are lit\n");
		printf("       delta is the stepping distance in voxels\n");
		printf("       ncpu is the number of processor threads to use (OpenMP)\n");
		return 0;
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	sscanf(argv[3],"%d",&stencil_size);
	sscanf(argv[4],"%lf",&threshold_fraction);
	sscanf(argv[5],"%d",&delta);
	sscanf(argv[6],"%d",&npar);
	use_compression = true;

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

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	p = (unsigned char *)(im->GetBufferPointer());

	ImageType::Pointer im2 = ImageType::New();
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
	im2->SetRegions(imregion);
	im2->Allocate();
	p2 = (unsigned char *)(im2->GetBufferPointer());

	memcpy(p2, p, width*depth*height);

	int dz = depth/npar;
	printf("dz: %d\n",dz);
	if (dz < stencil_size) {
		npar = depth/stencil_size;
		dz = depth/npar;
		printf("revised npar: %d dz: %d\n",npar, dz);
	}
	for (kpar=0; kpar<npar; kpar++) {
		zparstart[kpar] = kpar*dz;
		zparend[kpar] = (kpar+1)*dz - 1;
		printf("kpar: %d z: %4d %4d\n",kpar,zparstart[kpar],zparend[kpar]);
	}
	zparend[npar-1] = depth-1;

	#pragma omp parallel for num_threads(npar) private(zstart, zend)
	for (kpar=0; kpar<npar; kpar++) {
		zstart = zparstart[kpar];
		zend = zparend[kpar];
		process(zstart, zend);
	}

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im2);
	if (use_compression) {
		writer->UseCompressionOn();
	}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}
	if (use_compression) {
		printf("Created compressed image file: %s\n",argv[2]);
	} else {
		printf("Created uncompressed image file: %s\n",argv[2]);
	}

	return 0;
}