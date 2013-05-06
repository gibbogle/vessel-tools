//--------------------------------------------------------------------------------------
// Generate the negative of a tiff image.
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
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int x, y, z;
	int width, height, depth, xysize;
	char comp;
	char *invertfile;
	char drive[32], dir[256],filename[256], ext[32];
	bool use_compression;

	if (argc != 4) {
		printf("Usage: cut input_tiff invert_tiff comp\n");
		printf("       comp = 'C' to compress image, 'U' otherwise\n");
		return 1;
	}
	
	printf("Input image file: %s\n",argv[1]);
	invertfile = argv[2];
	sscanf(argv[3],"%c",&comp);
	
	if (comp == 'C') {
		use_compression = true;
	} else {
		use_compression = false;
	}
	/*
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
	*/
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

	for (x=0;x<width; x++) {
		for (y=0; y<height; y++) {
			for (z=0; z<depth; z++) {
				V(x,y,z) = 255 - V(x,y,z);
			}
		}
	}
	
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(invertfile);
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
		return 4;
	}
	printf("Wrote inverted file\n");

	return 0;
}