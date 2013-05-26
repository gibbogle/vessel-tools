// Segmentation of a tissue section stack that is imaged using a stain
// that shows up blood vessels and faint tissue edges. The faint
// tissue edges aren't interesting and need to be got rid of.
//
// The idea is to create a tissue mask, erode it a bit and then mask
// the original.
//
// Problem from Gib Bogle, via the ITK list.

#include <iostream>
#include "tclap/CmdLine.h"
#include "ioutils.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkBinaryShapeOpeningImageFilter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include <itkMorphologicalWatershedFromMarkersImageFilter.h>

//#if ITK_VERSION_MAJOR >= 4
#include "itkMultiplyImageFilter.h"
//#else
//#include "itkMultiplyByConstantImageFilter.h"
//#endif
#include "itkImageRegionConstIterator.h"

#include <itkSmartPointer.h>
namespace itk
{
    template <typename T>
    class Instance : public T::Pointer {
    public:
        Instance() : SmartPointer<T>( T::New() ) {}
    };
}

#include "morphutils.h"

typedef class CmdLineType
{
public:
  std::string InputIm, OutputImPrefix, suffix;
  float xspacing, yspacing, zspacing;	// not used
  float closingsize;					// box radius(?)
  float sigma;
  float SizeOpening;
  float FirstThresh;
  float peel;
} CmdLineType;

void ParseCmdLine(int argc, char* argv[],
                  CmdLineType &CmdLineObj
                  )
{
  using namespace TCLAP;
  try
    {
    // Define the command line object.
    CmdLine cmd("peelTissueMask ", ' ', "0.9");

    ValueArg<std::string> inArg("","input","input image",true,"result","string");
    cmd.add( inArg );

    ValueArg<std::string> outArg("","outputprefix","output image prefix",true,"result","string");
    cmd.add( outArg );

    ValueArg<std::string> sufArg("","outputsuffix","output image suffix",false,".tif","string");
    cmd.add( sufArg );

    ValueArg<float> xArg("","xspacing","input image xspacing",false, 1.0 ,"float");
    cmd.add( xArg );

    ValueArg<float> yArg("","yspacing","input image xspacing",false, 1.0 ,"float");
    cmd.add( yArg );

    ValueArg<float> zArg("","zspacing","input image xspacing",false, 1.0 ,"float");
    cmd.add( zArg );

    ValueArg<float> sigArg("","sigma","smoothing sigma - can be zero",false, 0.0 ,"float");
    cmd.add( sigArg );

    ValueArg<float> threshArg("","thresh","initial threshold",false, 1.0 ,"float");
    cmd.add( threshArg );

    ValueArg<float> clArg("","closingsize","size of closing filter for marker selection",false, 1.0 ,"string");
    cmd.add( clArg );

    ValueArg<float> pArg("","peelsize","size of erosion applied to mask",false, 1.0 ,"string");
    cmd.add( pArg );

    // Parse the args.
    cmd.parse( argc, argv );

    CmdLineObj.InputIm = inArg.getValue();
    CmdLineObj.OutputImPrefix = outArg.getValue();
    CmdLineObj.suffix = sufArg.getValue();
    CmdLineObj.xspacing = xArg.getValue();
    CmdLineObj.yspacing = yArg.getValue();
    CmdLineObj.zspacing = zArg.getValue();
    CmdLineObj.sigma = sigArg.getValue();
    CmdLineObj.closingsize = clArg.getValue();
    CmdLineObj.SizeOpening = 1000;
    CmdLineObj.FirstThresh = threshArg.getValue();
    CmdLineObj.peel = pArg.getValue();
    }
  catch (ArgException &e)  // catch any exceptions
    {
    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}

template <class PixType, int dim>
void doSeg(const CmdLineType &CmdLineObj)
{
  typedef typename itk::Image<PixType, dim> RawImType;
  typedef typename RawImType::Pointer PRawImType;

  typedef typename itk::Image<unsigned char, dim> MaskImType;
  typedef typename MaskImType::Pointer PMaskImType;

  typedef typename itk::Image<unsigned int, dim> LabImType;
  typedef typename LabImType::Pointer PLabImType;

  // load input = the raw image
  PRawImType input = readIm<RawImType>(CmdLineObj.InputIm);
  printf("Read input image file\n");
  // stuff to do with spacing information should go here - deal with
  // this later.

  // begin filtering to produce a foreground and background marker -
  // this uses world units for filter parameters. There is nothing
  // special about this order of operations.

  // Binarise first as the data seems clean
  // filter, and it probably wouldn't make much difference. I've
  // chosen to apply the threshold estimation step to the stage which
  // will have the most non zero voxels.
  itk::Instance<itk::BinaryThresholdImageFilter<MaskImType, MaskImType> > Thresh;
  Thresh->SetInput(input);
  Thresh->SetUpperThreshold(CmdLineObj.FirstThresh);
  Thresh->SetOutsideValue(1);
  Thresh->SetInsideValue(0);
  printf("closing (1)\n");
  PRawImType closed = doClosingMM<RawImType>(Thresh->GetOutput(), CmdLineObj.closingsize);

  //Gib
  typedef itk::MultiplyImageFilter<MaskImType, MaskImType, RawImType> MultiplyImageFilterType;
  MultiplyImageFilterType::Pointer multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(closed);
  multiplyImageFilter->SetConstant(255);
  writeIm<RawImType>(multiplyImageFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_closed" + CmdLineObj.suffix);

  printf("eroding\n");
  // now do a basic erosion, say by half the closing size. We could do
  // this more efficiently, but don't worry for now
  PMaskImType eroded = doErodeMM<MaskImType>(closed, CmdLineObj.closingsize/4);

  // only keep big objects
  itk::Instance<itk::BinaryShapeOpeningImageFilter<MaskImType> > KeepBig;
  KeepBig->SetInput(eroded);
  KeepBig->SetForegroundValue(1);
  KeepBig->SetLambda(CmdLineObj.SizeOpening);
  KeepBig->SetAttribute("PhysicalSize");

  printf("keeping big\n");
  //Gib
  multiplyImageFilter->SetInput(KeepBig->GetOutput());
  writeIm<RawImType>(multiplyImageFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_eroded_big" + CmdLineObj.suffix);

  printf("dilating\n");
  PMaskImType dilated = doDilateMM<MaskImType>(closed, CmdLineObj.closingsize/4);

  // invert to create background marker
  printf("inverting\n");
  itk::Instance<itk::BinaryThresholdImageFilter<MaskImType, MaskImType> > Invert;
  Invert->SetInput(dilated);
  Invert->SetUpperThreshold(1);
  Invert->SetLowerThreshold(1);
  Invert->SetInsideValue(0);
  Invert->SetOutsideValue(2);

  //Gib
  multiplyImageFilter->SetInput(Invert->GetOutput());
  writeIm<RawImType>(multiplyImageFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_inverted" + CmdLineObj.suffix);

  itk::Instance<itk::MaximumImageFilter<MaskImType, MaskImType, MaskImType> > Comb;
  Comb->SetInput(Invert->GetOutput());
  Comb->SetInput2(KeepBig->GetOutput());

  printf("combining\n");
  writeIm<MaskImType>(Comb->GetOutput(), CmdLineObj.OutputImPrefix + "_marker" + CmdLineObj.suffix);

  //Gib
  multiplyImageFilter->SetInput(Comb->GetOutput());
  writeIm<RawImType>(multiplyImageFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_combined" + CmdLineObj.suffix);

  // now for watershed
  itk::Instance< itk::SmoothingRecursiveGaussianImageFilter <RawImType, RawImType> > Smoother;
  Smoother->SetInput(input);
  Smoother->SetSigma(CmdLineObj.sigma);

  itk::Instance< itk::MorphologicalWatershedFromMarkersImageFilter<RawImType, MaskImType> > WS;
  if (CmdLineObj.sigma > 0)
    {
    WS->SetInput(Smoother->GetOutput());
    }
  else
    {
    WS->SetInput(input);
    }

  WS->SetMarkerImage(Comb->GetOutput());
  // throw away the background label
  itk::Instance<itk::BinaryThresholdImageFilter<MaskImType, MaskImType> > Selector;
  Selector->SetInput(WS->GetOutput());
  Selector->SetUpperThreshold(1);
  Selector->SetLowerThreshold(1);
  Selector->SetInsideValue(1);
  Selector->SetOutsideValue(0);

  writeIm<MaskImType>(Selector->GetOutput(), CmdLineObj.OutputImPrefix + "_wsseg" + CmdLineObj.suffix);

  PMaskImType peeled = doErodeMM<MaskImType>(Selector->GetOutput(), CmdLineObj.peel);

  writeIm<MaskImType>(peeled, CmdLineObj.OutputImPrefix + "_peeled" + CmdLineObj.suffix);

  //Gib
  multiplyImageFilter->SetInput(peeled);
  writeIm<RawImType>(multiplyImageFilter->GetOutput(), CmdLineObj.OutputImPrefix + "_peeled255" + CmdLineObj.suffix);

  itk::Instance<itk::MaskImageFilter<RawImType, MaskImType> > masker;
  masker->SetInput(input);
  masker->SetInput2(peeled);

  writeIm<RawImType>(masker->GetOutput(), CmdLineObj.OutputImPrefix + "_masked" + CmdLineObj.suffix);

}

int main(int argc, char * argv[])
{
  //itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

  CmdLineType CmdLineObj;
  ParseCmdLine(argc, argv, CmdLineObj);
  const int dimension = 3;

  doSeg<unsigned char, dimension>(CmdLineObj);

  return(EXIT_SUCCESS);
}
