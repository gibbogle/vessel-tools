/*
 * Make a tiff file to specification for test purposes
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
unsigned char *p;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]

int main(int argc, char**argv)
{
	int n;
	int x, y, z, dx, dy, R, R2;
	int width, height, depth, xysize;
	bool use_compression;

	if (argc != 4) {
		printf("Usage: maketiff output_tiff n R\n");
		printf("       where n is the size (edge of cube)\n");
		printf("             R is the radius of the vessel\n");
		return 0;
	}

	printf("Output image file: %s\n",argv[1]);
	sscanf(argv[2],"%d",&n);
	sscanf(argv[3],"%d",&R);
	use_compression = true;
	R2 = R*R;

	/*
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
	*/
	width = n;
	height = n;
	depth = n;
	xysize = width*height;
	printf("Desired image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	ImageType::Pointer im = ImageType::New();
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
	im->SetRegions(imregion);
	im->Allocate();
	p = (unsigned char *)(im->GetBufferPointer());

//	width = im->GetLargestPossibleRegion().GetSize()[0];
//	height = im->GetLargestPossibleRegion().GetSize()[1];
//	depth = im->GetLargestPossibleRegion().GetSize()[2];
//	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	for (z=0; z<depth; z++)
	{
		for (x=0; x<width; x++)
		{
			dx = x - width/2;
			for (y=0; y<height; y++)
			{
				dy = y - height/2;
				if (dx*dx + dy*dy < R2) 
					V(x,y,z) = 255;
				else
					V(x,y,z) = 0;
			}
		}
	}

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[1]);
	writer->SetInput(im);
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