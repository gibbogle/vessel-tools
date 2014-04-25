// tiff.cpp
#include "mainwindow.h"

#include <cstdio>
#include <vector>

#include <algorithm>
#include <math.h>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <ctime>

#include "itkImportImageFilter.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int MainWindow::readTiff(const char *tifffile)
{
    printf("readTiff: %s\n",tifffile);
    typedef itk::ImageFileReader<ImageType_u8> FileReaderType_u8;
    FileReaderType_u8::Pointer reader = FileReaderType_u8::New();
    reader->SetFileName(tifffile);
    try
    {
        reader->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        std::cout << e << std::endl;
        return 1;
    }
    im_u8 = reader->GetOutput();
    width = im_u8->GetLargestPossibleRegion().GetSize()[0];
    height = im_u8->GetLargestPossibleRegion().GetSize()[1];
    depth = im_u8->GetLargestPossibleRegion().GetSize()[2];
    xysize = width*height;
    p_raw = (unsigned char *)(im_u8->GetBufferPointer());
    printf("\nImage dimensions: width, height, depth: %d %d %d\n",width,height,depth);
    tiff_read = true;
    return 0;
}
//    ImageType_u8::Pointer im_u8;
//    typedef itk::ImageFileReader<ImageType_u8> FileReaderType_u8;

//--------------------------------------------------------------------------------------
// The data in buffer is used to create the image.
//--------------------------------------------------------------------------------------
int MainWindow::createTiff(const char *tifffile, unsigned char *buffer)
{
    bool use_compression = true;
/*
    unsigned char *p_u8;
    typedef itk::Image<unsigned char,3> ImageType_u8;
    ImageType_u8::Pointer im_u8 = ImageType_u8::New();
    ImageType_u8::SizeType imsize;
    ImageType_u8::IndexType imstart;
    ImageType_u8::RegionType imregion;

    imsize[0] = width;
    imsize[1] = height;
    imsize[2] = depth;
    imstart[0] = 0;
    imstart[1] = 0;
    imstart[2] = 0;
    imregion.SetSize(imsize);
    imregion.SetIndex(imstart);
    im_u8->SetRegions(imregion);
    im_u8->Allocate();
    p_u8 = (unsigned char *)(im_u8->GetBufferPointer());
    // copy buffer to p_u8
    memcpy(p_u8,buffer,(size_t)(width*height*depth));
*/


    typedef unsigned char   PixelType;
    const unsigned int Dimension = 3;
    typedef itk::Image< PixelType, Dimension > ImageType;
    typedef itk::ImportImageFilter< PixelType, Dimension >   ImportFilterType;

    printf("createTiff\n");
    fflush(stdout);
/*
    ImportFilterType::Pointer importFilter = ImportFilterType::New();
    ImportFilterType::SizeType  size;
    size[0]  = width;  // size along X
    size[1]  = height;  // size along Y
    size[2]  = depth;  // size along Z

    ImportFilterType::IndexType start;
    start.Fill( 0 );

    ImportFilterType::RegionType region;
    region.SetIndex( start );
    region.SetSize(  size  );

    importFilter->SetRegion( region );

    double origin[ Dimension ];
    origin[0] = 0.0;    // X coordinate
    origin[1] = 0.0;    // Y coordinate
    origin[2] = 0.0;    // Z coordinate

    importFilter->SetOrigin( origin );

    double spacing[ Dimension ];
    spacing[0] = 1.0;    // along X direction
    spacing[1] = 1.0;    // along Y direction
    spacing[2] = 1.0;    // along Z direction

    importFilter->SetSpacing( spacing );

    const unsigned int numberOfPixels =  size[0] * size[1] * size[2];
    const bool importImageFilterWillOwnTheBuffer = true;
    importFilter->SetImportPointer( buffer, numberOfPixels, importImageFilterWillOwnTheBuffer );

    */

//    printf("Writing 8-bit file: %s  dimensions: width, height, depth: %d %d %d\n",argv[2],width,height,depth);
    typedef itk::ImageFileWriter<ImageType_u8> FileWriterType;
    FileWriterType::Pointer writer = FileWriterType::New();
    writer->SetFileName(tifffile);
    writer->SetInput(im_u8);
//    writer->SetInput(  importFilter->GetOutput()  );
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
//        free(p_u8);
        return 1;
    }
//    free(p_u8);
    printf("Wrote tiff file\n");
    fflush(stdout);

    return 0;
}
