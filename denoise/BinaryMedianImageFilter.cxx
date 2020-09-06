//  The filter reduces noise both in the background and foreground of
//  the image, as well as smoothing the contours of the regions.

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBinaryMedianImageFilter.h"

int main( int argc, char * argv[] )
{
  if ( argc < 3 ) {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImageFile outputImageFile radius" << std::endl;
    return 1;
  }

  typedef   unsigned char  InputPixelType;
  typedef   unsigned char  OutputPixelType;

  typedef itk::Image< InputPixelType,  3 >   InputImageType;
  typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

  typedef itk::BinaryMedianImageFilter<
               InputImageType, OutputImageType >  FilterType;

  FilterType::Pointer filter = FilterType::New();

//  The size of the neighborhood is defined along every dimension by
//  passing a SizeType object with the corresponding values. The
//  value on each dimension is used as the semi-size of a rectangular
//  box. For example, in a size of \(1,2\) will result in a 3x5 neighborhood.

  const unsigned int radius = atoi( argv[3] );
//  const unsigned int radiusY = atoi( argv[4] );

  InputImageType::SizeType indexRadius;

  indexRadius[0] = radius; // radius along x
  indexRadius[1] = radius; // radius along y
  indexRadius[2] = radius; // radius along z

  filter->SetRadius( indexRadius );

//  The input to the filter can be taken from any other filter, for example
//  a reader. The output can be passed down the pipeline to other filters,
//  for example, a writer. An update call on any downstream filter will
//  trigger the execution of the median filter.

  filter->SetInput( reader->GetOutput() );
  writer->SetInput( filter->GetOutput() );
  writer->UseCompressionOn();
  writer->Update();
  return 0;
}
