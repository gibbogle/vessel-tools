/*
 * Make a tiff file by calling Fortran DLL to acquire the data
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
	int nx, ny, nz, nbytes;
//	int n, x, y, z, dx, dy;
	int width, height, depth, xysize;
	char *tempFile, *outputFile;
	bool use_compression = true;
	FILE *fpdata;

	if (argc != 3) {
		printf("Usage: maketiff tempfile output_tiff\n");
		return 1;
	}

	tempFile = argv[1];
	outputFile = argv[2];
	printf("Output image file: %s\n",outputFile);
	//sscanf(argv[2],"%d",&nx);
	//sscanf(argv[3],"%d",&ny);
	//sscanf(argv[4],"%d",&nz);


	fpdata = fopen(tempFile,"rb");
	if (fpdata == NULL) {
		return 2;
	}
	fread(&nx,4,1,fpdata);
	fread(&ny,4,1,fpdata);
	fread(&nz,4,1,fpdata);

	width = nx;
	height = ny;
	depth = nz;
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

	//strcpy(cmdline,"main_test.exe");
	//strcat(cmdline," ");
	//strcat(cmdline,tempfile);
	//strcat(cmdline," ");
	//strcat(cmdline,argv[2]);
	//strcat(cmdline," ");
	//strcat(cmdline,argv[3]);
	//strcat(cmdline," ");
	//strcat(cmdline,argv[4]);
	//err = system(cmdline);

	nbytes = nx*ny*nz;
	fread(p,nbytes,1,fpdata);

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(outputFile);
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
	if (use_compression) {
		printf("Created compressed image file: %s\n",outputFile);
	} else {
		printf("Created uncompressed image file: %s\n",outputFile);
	}

	return 0;
}