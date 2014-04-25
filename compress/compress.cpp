/*
 * To compress a .tif
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
#include "itkTIFFImageIO.h"

typedef itk::Image<unsigned char,3> ImageType; 
ImageType::Pointer im;
int width, height, depth;
unsigned char *p;

int main(int argc, char**argv)
{
	char method;
	FILE *fp;
	char errfile[] = "error.log";
	fp = fopen(errfile,"w");

	if (argc != 4) {
		printf("Usage: compress input_tiff output_tiff method\n");
		printf("compression method: P = Packbits\n");
		printf("                    L = LZW\n");
		printf("                    Z = Zip/Deflate\n");
		fprintf(fp,"Usage: compress input_tiff output_tiff method\n");
		fprintf(fp,"compression method: P = Packbits\n");
		fprintf(fp,"                    L = LZW\n");
		fprintf(fp,"                    Z = Zip/Deflate\n");
		fprintf(fp,"Submitted command line: argc: %d\n",argc);
		for (int i=0; i<argc; i++) {
			fprintf(fp,"argv: %d: %s\n",i,argv[i]);
		}		fclose(fp);
		return 1;	// Wrong command line
	}

	printf("Input image file: %s\n",argv[1]);
	printf("Output image file: %s\n",argv[2]);
	printf("Compression method: %s\n",argv[3]);

	method = *argv[3];

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
		return 2;	// Read error on input file
	}

	im = reader->GetOutput();

	width = im->GetLargestPossibleRegion().GetSize()[0];
	height = im->GetLargestPossibleRegion().GetSize()[1];
	depth = im->GetLargestPossibleRegion().GetSize()[2];

	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(argv[2]);
	writer->SetInput(im);

	typedef  itk::TIFFImageIO TIFFIOType;
	TIFFIOType::Pointer tiffIO = TIFFIOType::New();
	writer->SetImageIO(tiffIO);
	writer->UseCompressionOn();
	if (method == 'P') {
		tiffIO->SetCompressionToPackBits();
		printf("Compressing with Packbits\n");
	} else if (method == 'L') {
		tiffIO->SetCompressionToLZW();
		printf("Compressing with LZW\n");
	} else if (method == 'Z') {
		tiffIO->SetCompressionToDeflate();
		printf("Compressing with Zip/Deflate\n");
	}
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		fprintf(fp,"Write error on output file\n");
		fclose(fp);
		return 3;	// Write error on output file
	}
	printf("Created compressed image file: %s\n",argv[2]);

	return 0;
}