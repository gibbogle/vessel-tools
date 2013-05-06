/*
 * To chop out a piece of a .tif
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
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]

int main(int argc, char**argv)
{
	int x1, x2, y1, y2, z1, z2;
	int x, y, z, xx, yy, zz;
	int width, height, depth, xysize;
	int width2, height2, depth2, xysize2;
	char comp;
	bool use_compression;

	if (argc != 10) {
		printf("Usage: chop input_tiff output_tiff x1 x2 y1 y2 z1 z2 comp\n");
		printf("       where the ranges (x1,x2), (y1,y2), (z1,z2) define the selected region (0-based)\n");
		printf("       (any value < 0 implies use full range for this axis)\n");
		printf("       comp = 'C' to compress image, 'U' otherwise\n");
		return 0;
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	sscanf(argv[3],"%d",&x1);
	sscanf(argv[4],"%d",&x2);
	sscanf(argv[5],"%d",&y1);
	sscanf(argv[6],"%d",&y2);
	sscanf(argv[7],"%d",&z1);
	sscanf(argv[8],"%d",&z2);
	sscanf(argv[9],"%c",&comp);
// Change from 1-based to 0-based indices
//	x1--; x2--; y1--; y2--; z1--; z2--;

//	printf("SPECIAL CASE to convert .nii to .tif\n");

	if (comp == 'C') {
		use_compression = true;
	} else {
		use_compression = false;
	}

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

	if (x1 < 0 || x2 < 0) {
		x1 = 0;
		x2 = width-1;
		printf("Using full range of x\n");
	}
	if (y1 < 0 || y2 < 0) {
		y1 = 0;
		y2 = height-1;
		printf("Using full range of y\n");
	}
	if (z1 < 0 || z2 < 0) {
		z1 = 0;
		z2 = depth-1;
		printf("Using full range of z\n");
	}

	width2 = x2 - x1 + 1;
	height2 = y2 - y1 + 1;
	depth2 = z2 - z1 + 1;
	xysize2 = width2*height2;
	printf("Desired image dimensions: width, height, depth: %d %d %d\n",width2,height2,depth2);

	ImageType::Pointer im2 = ImageType::New();
	ImageType::SizeType imsize; 
	imsize[0] = width2;
	imsize[1] = height2;
	imsize[2] = depth2;
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

	width2 = im2->GetLargestPossibleRegion().GetSize()[0];
	height2 = im2->GetLargestPossibleRegion().GetSize()[1];
	depth2 = im2->GetLargestPossibleRegion().GetSize()[2];
	printf("Cropped image dimensions: width, height, depth: %d %d %d\n",width2,height2,depth2);

	for (xx=0; xx<width2; xx++)
	{
		x = xx + x1;
		for (yy=0; yy<height2; yy++)
		{
			y = yy + y1;
			for (zz=0; zz<depth2; zz++)
			{
				z = zz + z1;
				V2(xx,yy,zz) = V(x,y,z);
			}
		}
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