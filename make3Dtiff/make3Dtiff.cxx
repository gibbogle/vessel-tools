/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ImageSeriesReadWrite.cxx,v $
  Language:  C++
  Date:      $Date: 2009-03-17 20:36:50 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

//  Software Guide : BeginLatex
//
//  This example illustrates how to read a series of 2D slices from independent
//  files in order to compose a volume. The class \doxygen{ImageSeriesReader}
//  is used for this purpose. This class works in combination with a generator
//  of filenames that will provide a list of files to be read. In this
//  particular example we use the \doxygen{NumericSeriesFileNames} class as
//  filename generator. This generator uses a \code{printf} style of string format
//  with a ``\code{\%d}'' field that will be successively replaced by a number specified
//  by the user. Here we will use a format like ``\code{file\%03d.png}'' for reading 
//  PNG files named file001.png, file002.png, file003.png... and so on.
//
//  This requires the following headers as shown.
//
//  \index{itk::ImageSeriesReader!header}
//  \index{itk::NumericSeriesFileNames!header}
//
//  Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkImage.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkTIFFImageIO.h"
#include "itkRescaleIntensityImageFilter.h"

#define V3D(a,b,c)  p3D[(c)*width*height+(b)*width+(a)]
#define V2D_8(a,b)  p2D_8[(b)*width+(a)]
#define V2D_16(a,b)  p2D_16[(b)*width+(a)]
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

bool fexists(const char *filename)
{
    return GetFileAttributes(filename) != INVALID_FILE_ATTRIBUTES;
}

int main( int argc, char ** argv )
{
	bool compress;
  // Verify the number of parameters in the command line
  if( argc < 7 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " 2DTemplateFile bits firstSliceValue lastSliceValue compressFlag output3DImageFile " << std::endl;
	std::cerr << "E.g.: " << argv[0] << " frame%04d.tif 8 0 99 1 result3D.tif " << std::endl;
	std::cerr << "where the 8-bit 2D images are: frame0000.tif, frame0001.tif,..,frame0099.tif " << std::endl;
    return 1;
    }

	// We start by defining the {PixelType} and {ImageType}.
	typedef unsigned char              PixelType8;
	typedef unsigned short             PixelType16;
	const unsigned int Dimension = 3;

	typedef itk::Image< PixelType8, 2 >  InputImageType8;
	typedef itk::Image< PixelType16, 2 >  InputImageType16;
	typedef itk::Image< PixelType8, Dimension >  OutputImageType;
  
	typedef itk::ImageFileReader< InputImageType8 >  ReaderType8;
	typedef itk::ImageFileReader< InputImageType16 >  ReaderType16;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	ReaderType8::Pointer reader8;		//= ReaderType8::New();
  	ReaderType16::Pointer reader16;		//= ReaderType16::New();
	WriterType::Pointer writer = WriterType::New();

	const char * inTemplateFile = argv[1];

	const unsigned int bits = atoi( argv[2] );
	const unsigned int first = atoi( argv[3] );
	const unsigned int last  = atoi( argv[4] );
	const unsigned int compressFlag  = atoi( argv[5] );
	compress = (compressFlag == 1);
	const char * outputFilename = argv[6];
	char infile[256];
	long long width, height, depth;
	double range;
	int maxmax;
	unsigned char *p2D_8;
	unsigned short *p2D_16;
	unsigned char *p3D;

  // Testing fexists()
  //if (fexists("zzz.txt")) {
	 // printf("file exists\n");
  //} else {
	 // printf("file does not exist\n");
  //}
  //return 0;

	if (bits == 8) {
	  	reader8 = ReaderType8::New();
	} else {
		range = pow(2,16);		// 12-bit data stored in 16 bits
  		reader16 = ReaderType16::New();
	}

  	OutputImageType::Pointer image3D;
	image3D = OutputImageType::New();
	OutputImageType::SizeType imsize; 
	OutputImageType::IndexType imstart; 
	imstart[0] = 0;
	imstart[1] = 0;
	imstart[2] = 0;
	OutputImageType::RegionType imregion; 

	// Count files
	depth = 0;
	for (int i=first; i<=last; i++) {
		sprintf(infile, inTemplateFile, i);
		// Check for existence of file
		if (fexists(infile)) depth++;
	}
	
	maxmax = 0;
	int z = -1;
	for (int i=first; i<=last; i++) {
//		int z = i - first;
		sprintf(infile, inTemplateFile, i);
		// Check for existence of file
		if (!fexists(infile)) {
//			printf("Missing file: %s\n",infile);
			continue;
		}
		z++;
//		printf("z: %d  infile: %s\n",z,infile);
		if (bits == 8) {
			reader8->SetFileName(infile);
			try
			{
				reader8->Update();
			}
			catch (itk::ExceptionObject &e)
			{
				std::cout << e << std::endl;
				return 2;
			}
			InputImageType8 *im2D = reader8->GetOutput();
			width = im2D->GetLargestPossibleRegion().GetSize()[0];
			height = im2D->GetLargestPossibleRegion().GetSize()[1];
 			p2D_8 = (unsigned char *)(im2D->GetBufferPointer());
		} else {
			reader16->SetFileName(infile);
			try
			{
				reader16->Update();
			}
			catch (itk::ExceptionObject &e)
			{
				std::cout << e << std::endl;
				return 2;
			}
			InputImageType16 *im2D = reader16->GetOutput();
			width = im2D->GetLargestPossibleRegion().GetSize()[0];
			height = im2D->GetLargestPossibleRegion().GetSize()[1];
 			p2D_16 = (unsigned short *)(im2D->GetBufferPointer());
		}

		if (i == first) {
//			depth = last - first + 1;
			printf("Image dimensions: width, height, depth: %d %d %d\n",width,height,depth);
			imsize[0] = width;
			imsize[1] = height;
			imsize[2] = depth;
			imregion.SetSize(imsize);
			imregion.SetIndex(imstart);
			image3D->SetRegions(imregion);
			image3D->Allocate();
			p3D = (unsigned char *)(image3D->GetBufferPointer());
		}
		int maxpix = 0;
		int max3D = 0;
		for (int x=0; x<width; x++) {
			for (int y=0; y<height; y++) {
				if (bits == 8) {
					maxpix = MAX(maxpix,V2D_8(x,y));
					V3D(x,y,z) = V2D_8(x,y);
				} else {
					maxpix = MAX(maxpix,V2D_16(x,y));
					double r = V2D_16(x,y)/range;
					V3D(x,y,z) = int(MIN(255.,(r*255. + 0.5)));
					max3D = MAX(V3D(x,y,z),max3D);
					//if (z == 0 && x <100 && y == 10) {
					//	printf("%d %d %f  V2D: %d  V3D: %d\n",x,y,r,V2D_16(x,y),V3D(x,y,z));
					//}
				}
			}
		}
		maxmax = MAX(maxpix,maxmax);
		printf("slice i: %3d maxpix: %6d max3D: %6d maxmax: %6d\n",i,maxpix,max3D,maxmax);
	}

    writer->SetFileName( outputFilename );
	writer->SetInput(image3D);
	if (compress) {
		writer->UseCompressionOn();
	} else {
		writer->UseCompressionOff();
	}
	printf("Writing 3D image file\n");
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cout << e << std::endl;
//		fprintf(fp,"Write error on output file\n");
//		fclose(fp);
		return 3;	// Write error on output file
	}
	printf("Created 3D 8-bit image file: %s\n",outputFilename);

  return 0;
}
