//--------------------------------------------------------------------------------------
// To cut a tiff file into two pieces
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

#include "itkTIFFImageIO.h"
#include "itkRGBAPixel.h"


typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer im, im1, im2;
unsigned char *p, *p1, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V1(a,b,c)  p1[(c)*xysize1+(b)*width1+(a)]
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#ifdef __linux__
void _splitpath(const char* Path, char* Drive, char* Directory, char* Filename, char* Extension);
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	long int x, y, z, xx, yy, zz, pos;
	long int width, height, depth, xysize;
	long int width1, height1, depth1, xysize1;
	long int width2, height2, depth2, xysize2;
	long int xstart2, ystart2, zstart2;
	char axis, comp;
	char LEfile[256], GTfile[256], tag[13];
	char drive[32], dir[256],filename[256], ext[32];
	bool use_compression;

	if (argc != 5) {
		printf("Usage: cut input_tiff axis pos comp\n");
		printf("       axis = X, Y or Z\n");
		printf("       pos = cut position (0-based)\n");
		printf("       comp = 'C' to compress image, 'U' otherwise\n");
		return 1;
	}
	
	printf("Input image file: %s\n",argv[1]);
	sscanf(argv[2],"%c",&axis);
	sscanf(argv[3],"%d",&pos);
	sscanf(argv[4],"%c",&comp);
	if (axis == 'X') axis = 'x';
	if (axis == 'Y') axis = 'y';
	if (axis == 'Z') axis = 'z';
	printf("Axis chosen: %c\n",axis);
	if (axis != 'x' && axis != 'y' && axis != 'z') {
		printf("Bad axis choice: %c\n",axis);
		return 2;
	}
	
	if (comp == 'C') {
		use_compression = true;
	} else {
		use_compression = false;
	}

	_splitpath(argv[1],drive,dir,filename,ext);
	strcpy(LEfile,drive);
	strcat(LEfile,dir);
	strcat(LEfile,filename);
	strcpy(GTfile,LEfile);
	sprintf(tag,"_%cLE%04d.tif",axis,pos);
	strcat(LEfile,tag);
	sprintf(tag,"_%cGT%04d.tif",axis,pos);
	strcat(GTfile,tag);
	printf("LEfile: %s\n",LEfile);
	printf("GTfile: %s\n",GTfile);

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
		return 3;
	}

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	p = (unsigned char *)(im->GetBufferPointer());

	if (axis == 'x' && pos >= width-1) {
		printf("Bad cut position: %c = %d\n",axis,pos);
		return 4;
	}
	if (axis == 'y' && pos >= height-1) {
		printf("Bad cut position: %c = %d\n",axis,pos);
		return 4;
	}
	if (axis == 'z' && pos >= depth-1) {
		printf("Bad cut position: %c = %d\n",axis,pos);
		return 4;
	}

	width1 = width2 = width;
	height1 = height2 = height;
	depth1 = depth2 = depth;
	xstart2 = 0;
	ystart2 = 0;
	zstart2 = 0;

	if (axis == 'x') {
		width1 = pos+1;
		width2 = width - width1;
		xstart2 = width1;
	}
	if (axis == 'y') {
		height1 = pos+1;
		height2 = height - height1;
		ystart2 = height1;
	}
	if (axis == 'z') {
		depth1 = pos+1;
		depth2 = depth - depth1;
		zstart2 = depth1;
	}
	
	printf("depth1: %d  depth2: %d  zstart2: %d\n",depth1,depth2,zstart2);

	ImageType::Pointer im1 = ImageType::New();
	ImageType::SizeType imsize; 
	ImageType::IndexType imstart; 
	ImageType::RegionType imregion; 

	imsize[0] = width1;
	imsize[1] = height1;
	imsize[2] = depth1;
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	im1->SetRegions(imregion);
	im1->Allocate();
	p1 = (unsigned char *)(im1->GetBufferPointer());
	xysize1 = width1*height1;

	for (x=0; x<width1; x++) {
		for (y=0; y<height1; y++) {
			for (z=0; z<depth1; z++) {
				V1(x,y,z) = V(x,y,z);
			}
		}
	}
	printf("Writing LE file: %s  dimensions: width, height, depth: %d %d %d\n",LEfile,width1,height1,depth1);
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
//	writer->SetFileName(argv[2]);
	writer->SetFileName(LEfile);
  typedef  itk::TIFFImageIO TIFFIOType;
//  WriterType::Pointer writer = WriterType::New();
  TIFFIOType::Pointer tiffIO = TIFFIOType::New();
//  tiffIO->SetPixelType(itk::ImageIOBase::GRAYSCALE);
//  writer->SetFileName(outputFilename);
  writer->SetInput(im1);
  writer->SetImageIO(tiffIO);

	writer->SetInput(im1);
	if (use_compression) {
		writer->UseCompressionOn();
		tiffIO->SetCompressionToDeflate();
	}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 5;
	}
	printf("Wrote LE file\n");

	ImageType::Pointer im2 = ImageType::New();
	imsize[0] = width2;
	imsize[1] = height2;
	imsize[2] = depth2;
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	imregion.SetSize(imsize);
	imregion.SetIndex(imstart);
	im2->SetRegions(imregion);
	im2->Allocate();
	p2 = (unsigned char *)(im2->GetBufferPointer());
	xysize2 = width2*height2;

	for (x=0; x<width2; x++) {
		xx = xstart2 + x;
		for (y=0; y<height2; y++) {
			yy = ystart2 + y;
			for (z=0; z<depth2; z++) {
				zz = zstart2 + z;
				V2(x,y,z) = V(xx,yy,zz);
			}
		}
	}

	printf("Writing GT file: %s  dimensions: width, height, depth: %d %d %d\n",GTfile,width2,height2,depth2);

	writer->SetFileName(GTfile);
	writer->SetInput(im2);
	if (use_compression) {
		writer->UseCompressionOn();
		tiffIO->SetCompressionToLZW();
	}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 5;
	}
	printf("Wrote GT file\n");
	
	return 0;
}