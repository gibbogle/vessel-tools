//--------------------------------------------------------------------------------------
// To rejoin a tiff file that was cut into two pieces
//--------------------------------------------------------------------------------------
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
ImageType::Pointer im, im1, im2;
unsigned char *p, *p1, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V1(a,b,c)  p1[(c)*xysize1+(b)*width1+(a)]
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


int main(int argc, char**argv)
{
	int x, y, z, xx, yy, zz;
	int width, height, depth, xysize;
	int width1, height1, depth1, xysize1;
	int width2, height2, depth2, xysize2;
	int xstart2, ystart2, zstart2;
	char axis, comp;
	char *LEfile, *GTfile, *joinedfile;
	bool use_compression;

	if (argc != 6) {
		printf("Usage: join input_LE_tiff input_GT_tiff output_tiff axis comp\n");
		printf("       axis = X, Y or Z\n");
		printf("       comp = 'C' to compress image, 'U' otherwise\n");
		return 0;
	}
	
	LEfile = argv[1];
	GTfile = argv[2];
	joinedfile = argv[3];
	printf("Input LE image file: %s\n",LEfile);
	printf("Input GT image file: %s\n",GTfile);
	sscanf(argv[4],"%c",&axis);
	if (axis == 'X') axis = 'x';
	if (axis == 'Y') axis = 'y';
	if (axis == 'Z') axis = 'z';
	printf("Axis chosen: %c\n",axis);
	if (axis != 'x' && axis != 'y' && axis != 'z') {
		printf("Bad axis choice: %c\n",axis);
		return 1;
	}
	sscanf(argv[5],"%c",&comp);
	
	if (comp == 'C') {
		use_compression = true;
	} else {
		use_compression = false;
	}
	
	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader1 = FileReaderType::New();

	reader1->SetFileName(LEfile);
	try
	{
		reader1->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}

	im1 = reader1->GetOutput();
	width1 = im1->GetLargestPossibleRegion().GetSize()[0];
	height1 = im1->GetLargestPossibleRegion().GetSize()[1];
	depth1 = im1->GetLargestPossibleRegion().GetSize()[2];
	xysize1 = width1*height1;
	printf("LE image dimensions: width, height, depth: %d %d %d\n",width1,height1,depth1);
	p1 = (unsigned char *)(im1->GetBufferPointer());

	FileReaderType::Pointer reader2 = FileReaderType::New();
	reader2->SetFileName(GTfile);
	try
	{
		reader2->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}

	im2 = reader2->GetOutput();
	width2 = im2->GetLargestPossibleRegion().GetSize()[0];
	height2 = im2->GetLargestPossibleRegion().GetSize()[1];
	depth2 = im2->GetLargestPossibleRegion().GetSize()[2];
	xysize2 = width2*height2;
	printf("GT image dimensions: width, height, depth: %d %d %d\n",width2,height2,depth2);
	p2 = (unsigned char *)(im2->GetBufferPointer());

	if (axis == 'x' && height1 != height2) {
		printf("Inconsistent heights for X cut: %d %d\n",height1,height2);
		return 2;
	}
	if (axis == 'x' && depth1 != depth2) {
		printf("Inconsistent depths for X cut: %d %d\n",depth1,depth2);
		return 2;
	}
	if (axis == 'y' && width1 != width2) {
		printf("Inconsistent widths for Y cut: %d %d\n",width1,width2);
		return 2;
	}
	if (axis == 'y' && depth1 != depth2) {
		printf("Inconsistent depths for Y cut: %d %d\n",depth1,depth2);
		return 2;
	}
	if (axis == 'z' && width1 != width2) {
		printf("Inconsistent widths for Z cut: %d %d\n",width1,width2);
		return 2;
	}
	if (axis == 'z' && height1 != height2) {
		printf("Inconsistent heights for Z cut: %d %d\n",height1,height2);
		return 2;
	}

	xstart2 = 0;
	ystart2 = 0;
	zstart2 = 0;
	if (axis == 'x') {
		width = width1 + width2;
		height = height1;
		depth = depth1;
		xstart2 = width1;
	}
	if (axis == 'y') {
		width = width1;
		height = height1 + height2;
		depth = depth1;
		ystart2 = height1;
	}
	if (axis == 'z') {
		width = width1;
		height = height1;
		depth = depth1 + depth2;
		zstart2 = depth1;
	}

	ImageType::Pointer im = ImageType::New();
	ImageType::SizeType imsize; 
	ImageType::IndexType imstart; 
	ImageType::RegionType imregion; 

	imsize[0] = width;
	imsize[1] = height;
	imsize[2] = depth;
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	im->SetRegions(imregion);
	im->Allocate();
	p = (unsigned char *)(im->GetBufferPointer());
	xysize = width*height;

	for (x=0; x<width1; x++) {
		for (y=0; y<height1; y++) {
			for (z=0; z<depth1; z++) {
				V(x,y,z) = V1(x,y,z);
			}
		}
	}
	for (x=0; x<width2; x++) {
		xx = xstart2 + x;
		for (y=0; y<height2; y++) {
			yy = ystart2 + y;
			for (z=0; z<depth2; z++) {
				zz = zstart2 + z;
				V(xx,yy,zz) = V2(x,y,z);
			}
		}
	}

	printf("Writing joined file: %s\n   dimensions: width, height, depth: %d %d %d\n",joinedfile,width,height,depth);
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(joinedfile);
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
		return 3;
	}
	printf("Wrote joined file\n");
	return 0;
}