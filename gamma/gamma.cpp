/*
 * To apply gamma correction to a greyscale image
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
unsigned char *p;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define V(a,b,c)  p[(c)*imsize+(b)*width+(a)]

#define USE_COMPRESSION true

int correct(double gamma)
{
	double e, v0, v1;
	int ix, iy, iz;

	printf("correct\n");
	e = 1/gamma;
	for (ix=0; ix<width; ix++) {
		for (iy=0; iy<height; iy++) {
			for (iz=0; iz<depth; iz++) {
				v0 = V(ix,iy,iz);
				v0 = v0/255.0;
				if (v0 != 0)
					v1 = pow(v0,e);
				else
					v1 = 0;
				V(ix,iy,iz) = (unsigned char)(255*v1);
			}
		}
	}
	return 0;
}

int main(int argc, char**argv)
{
//	int npar;
	double gamma;
	time_t t1;
	t1 = time(NULL);
	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	if (USE_COMPRESSION)
		printf("Compression on\n");
	else
		printf("Compression off\n");

	if (argc != 4) {
		printf("Usage: gamma input_tiff output_tiff gamma\n");
		printf("       where: gamma is the gamma correction factor (> 1)\n");
//		printf("              ncpu is the number of OpenMP threads\n");
		fprintf(fp,"Usage: gamma input_tiff output_tiff gamma\n");
		fprintf(fp,"       where: gamma is the gamma correction factor (> 1)\n");
		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}
		fclose(fp);
		return 1;
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	sscanf(argv[3],"%lf",&gamma);
	printf("Gamma correction: %f\n",gamma);
//	sscanf(argv[4],"%d",&npar);
//	printf("NCPU: %d\n",npar);

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
	p = (unsigned char *)(im->GetBufferPointer());

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
	printf("Array size: %d\n",width*height*depth);

	correct(gamma);

	printf("Gamma correction completed\n");

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
	printf("Created gamma-corrected image file: %s\n",argv[2]);

	return 0;
}