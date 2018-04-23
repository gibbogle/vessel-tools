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
#include <itkBinaryDilateImageFilter.h>
#include "itkFlatStructuringElement.h"
#include "itkBinaryBallStructuringElement.h"

typedef itk::Image<unsigned char,3> ImageType;
ImageType::Pointer inputImage;

unsigned char *p, *p2;

#define V(a,b,c)  p[(c)*xysize+(b)*width+(a)]
#define V2(a,b,c)  p2[(c)*xysize2+(b)*width2+(a)]

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int ImWriter(ImageType::Pointer im, char *filename)
{
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(filename);
	writer->SetInput(im);
	writer->UseCompressionOn();
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}
	return 0;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
int main(int argc, char**argv)
{
	int width, height, depth, xysize, err;
	char *inputFile, *outputFile;
	float ball_radius;
	bool compressdata = true;

	if (argc != 4) {
		printf("Usage: dilate input_tiff output_tiff ball_radius\n");
		return 0;
	}
	inputFile = argv[1];
	outputFile = argv[2];
	sscanf(argv[3],"%f",&ball_radius);
	compressdata = true;
	printf("Input image file: %s\n",inputFile);
	printf("Output image file: %s\n",outputFile);
	printf("Ball radius: %f\n",ball_radius);

	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();

	reader->SetFileName(inputFile);
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
		return 1;
	}

	inputImage = reader->GetOutput();

	width = inputImage->GetLargestPossibleRegion().GetSize()[0];
	height = inputImage->GetLargestPossibleRegion().GetSize()[1];
	depth = inputImage->GetLargestPossibleRegion().GetSize()[2];
	xysize = width*height;
	printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);

//	const unsigned int radiusValue = 2;
	typedef itk::FlatStructuringElement< 3 > StructuringElementType;
	StructuringElementType::RadiusType radius;
	radius.Fill( ball_radius );
	StructuringElementType structuringElement = StructuringElementType::Ball( radius );

	typedef itk::BinaryDilateImageFilter< ImageType, ImageType, StructuringElementType > BinaryDilateImageFilterType;

	BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetInput( reader->GetOutput() );
	dilateFilter->SetKernel( structuringElement );

	printf("Writing dilated tiff\n");
	err = ImWriter(dilateFilter->GetOutput(),outputFile);
	if (err != 0) {
		printf("ImWriter error on output file\n");
		return 1;
	}

	return 0;
}