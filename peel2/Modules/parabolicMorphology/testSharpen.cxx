#include <iomanip>
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"

#include "plotutils.h"

#include "itkGrayscaleDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"


#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkMorphologicalSharpeningImageFilter.h"
#include "itkTimeProbe.h"
#include "itkMultiThreader.h"


int main(int argc, char * argv[])
{

  int iterations = 1;

  if (argc != 4)
    {
    std::cerr << "Usage: " << argv[0] << " iterations outputimage outputprofile" << std::endl;
    return(EXIT_FAILURE);
    }

  iterations = atoi(argv[1]);

  //itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);
  const int dim = 2;
  
  typedef unsigned char PType;
  typedef itk::Image< PType, dim > IType;
  typedef itk::Image< float, dim > FType;

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();

  // create the input image - we will blur a dot, threshold then blur again
  IType::Pointer input = IType::New();
  IType::SizeType size;
  IType::IndexType index;
  IType::RegionType region;
  IType::SpacingType spacing;
  
  size.Fill(100);
  spacing.Fill(1);
  region.SetSize(size);
  
  index.Fill(50);

  input->SetRegions(region);
  input->SetSpacing(spacing);
  input->Allocate();

  input->FillBuffer(0);
  input->SetPixel(index, 255);

  typedef itk::BinaryBallStructuringElement< PType, dim > SRType;
  typedef itk::GrayscaleDilateImageFilter<IType, IType, SRType> DilateType;
  DilateType::Pointer smallDilate = DilateType::New();
  SRType smallkernel;
  SRType::RadiusType smallrad = smallkernel.GetRadius();
  smallrad.Fill(7);
  smallkernel.SetRadius(smallrad);
  smallkernel.CreateStructuringElement();
  smallDilate->SetKernel(smallkernel);

  typedef itk::SmoothingRecursiveGaussianImageFilter<IType, IType> SmootherType;

  SmootherType::Pointer smoother = SmootherType::New();
  
  smallDilate->SetInput(input);

  smoother->SetInput(smallDilate->GetOutput());
  smoother->SetSigma(3);

  writer->SetInput( smallDilate->GetOutput());
  writer->SetFileName( "input.tif" );
  writer->Update();

  writer->SetInput( smoother->GetOutput());
  writer->SetFileName( "blurrredinput.tif" );
  writer->Update();


  // now to apply the sharpening
  
  typedef itk::MorphologicalSharpeningImageFilter< IType, FType > FilterType;

  FilterType::Pointer filter = FilterType::New();

  filter->SetInput( smoother->GetOutput());
  filter->SetScale(1);
  filter->SetIterations(iterations);
  filter->Update();

  typedef itk::ImageFileWriter< FType > FlWriterType;
  FlWriterType::Pointer flwriter = FlWriterType::New();

  flwriter->SetInput(filter->GetOutput());
  flwriter->SetFileName(argv[2]);
  flwriter->Update();

  // write out profiles
  // input 
  IType::IndexType first, last;
  first[0]=50;
  first[1]=0;
  last[0]=50;
  last[1]=99;

  extractProfile<IType>(smallDilate->GetOutput(), first, last, "inputprof.txt");
  extractProfile<IType>(smoother->GetOutput(), first, last, "blurredprof.txt");
  extractProfile<FType>(filter->GetOutput(), first, last, argv[3]);

  return EXIT_SUCCESS;
}

