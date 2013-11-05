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

//xysize = width*height;
//#define V_u8(a,b,c)  p_u8[(c)*xysize+(b)*width+(a)]

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

int MainWindow::readTiff(const char *tifffile, int *width, int *height, int *depth)
{
//    typedef itk::Image<unsigned char,3> ImageType_u8;
//    ImageType_u8::Pointer im_u8;
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
    *width = im_u8->GetLargestPossibleRegion().GetSize()[0];
    *height = im_u8->GetLargestPossibleRegion().GetSize()[1];
    *depth = im_u8->GetLargestPossibleRegion().GetSize()[2];
    p_im = (unsigned char *)(im_u8->GetBufferPointer());
    printf("Image dimensions: width, height, depth: %d %d %d\n",*width,*height,*depth);
    tiff_read = true;
    return 0;
}

int MainWindow::createTiff(const char *tifffile, unsigned char *buffer, int width, int height, int depth)
{
    bool use_compression = true;
    unsigned char *p_u8;
    typedef itk::Image<unsigned char,3> ImageType_u8;
//    ImageType_u8::Pointer im_u8;
    typedef itk::ImageFileReader<ImageType_u8> FileReaderType_u8;

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

//    for (x=0; x<width; x++) {
//        for (y=0; y<height; y++) {
//            for (z=0; z<depth; z++) {
//                if (issigned)
//                    V_u8(x,y,z) = V_s16(x,y,z)/256;
//                else
//                    V_u8(x,y,z) = V_u16(x,y,z)/256;
//            }
//        }
//    }

//    printf("Writing 8-bit file: %s  dimensions: width, height, depth: %d %d %d\n",argv[2],width,height,depth);
    typedef itk::ImageFileWriter<ImageType_u8> FileWriterType;
    FileWriterType::Pointer writer = FileWriterType::New();
    writer->SetFileName(tifffile);
    writer->SetInput(im_u8);
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

    return 0;
}
