#ifndef __itkParabolicErodeDilateImageFilter_h
#define __itkParabolicErodeDilateImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"

namespace itk
{
/**
 * \class ParabolicErodeDilateImageFilter
 * \brief Parent class for morphological operations with parabolic
 * structuring elements.
 *
 * Parabolic structuring elements are the morphological counterpart of
 * the gaussian function in linear processing. Parabolic structuring
 * elements are dimensionally decomposable and fast algorithms are
 * available for computing erosions and dilations along lines.
 * This class implements the "point of contact" algorithm and the
 * "intersection" algorithm. The contact point algorithm is faster at
 * small scales. The intersection algorithm was rediscovered by
 * Felzenszwalb & Huttenlocher, but was actually described and tested
 * by van den Boomgaard many years earlier. The intersection algorithm
 * is faster for scales > 1, and independent of scale. It also seems
 * to perform significantly better in the distance transform
 * application. The intersection algorithm is also used in the IJ
 * article on generalised distance transforms.
 * 
 * Parabolic structuring functions can be used as a fast alternative
 * to the "rolling ball" structuring element classically used in
 * background estimation, for example in ImageJ, have applications in
 * image sharpening and distance transform computation.
 * 
 * This class uses an internal buffer of RealType pixels for each
 * line. This line is cast to the output pixel type when written back
 * to the output image. Since the filter uses dimensional
 * decomposition this approach could result in inaccuracy as pixels
 * are cast back and forth between low and high precision types. Use a
 * high precision output type and cast manually if this is a problem.
 *
 * Boomgaard, R. van den and Dorst, L. and Makram-Ebeid, L.S. and
 * Schavemaker, J. Quadratic structuring functions in mathematical
 * morphology. Mathematical Morphology and its Applications to Image
 * and Signal Processing. 
 *
 * Felzenszwalb, P.F. & Huttenlocher, D.P. Distance Transforms of Sampled Functions.
 * Techreport: Cornell Computing and Information Science, 2004.
 *
 * This filter is threaded. Threading mechanism derived from
 * SignedMaurerDistanceMap extensions by Gaetan Lehman
 *
 * \author Richard Beare, Department of Medicine, Monash University,
 * Australia.  <Richard.Beare@med.monash.edu.au>
 *
**/


template <typename TInputImage,
	  bool doDilate,
          typename TOutputImage= TInputImage >
class ITK_EXPORT ParabolicErodeDilateImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>
{

public:
  /** Standard class typedefs. */
  typedef ParabolicErodeDilateImageFilter  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ParabolicErodeDilateImageFilter, ImageToImageFilter);

  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef TOutputImage                                   OutputImageType;
  typedef typename TInputImage::PixelType                PixelType;
  typedef typename NumericTraits<PixelType>::RealType    RealType;
  typedef typename NumericTraits<PixelType>::ScalarRealType ScalarRealType;
  typedef typename TOutputImage::PixelType  OutputPixelType;

  /** Smart pointer typedef support.  */
  typedef typename TInputImage::Pointer  InputImagePointer;
  typedef typename TInputImage::ConstPointer  InputImageConstPointer;
  typedef typename TInputImage::SizeType    InputSizeType;
  typedef typename TOutputImage::SizeType   OutputSizeType;

  typedef typename OutputImageType::IndexType       OutputIndexType;

  /** a type to represent the "kernel radius" */
  typedef typename itk::FixedArray<ScalarRealType, TInputImage::ImageDimension> RadiusType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);


  typedef typename OutputImageType::RegionType OutputImageRegionType;
  /** Define the image type for internal computations
      RealType is usually 'double' in NumericTraits.
      Here we prefer float in order to save memory.  */

  typedef typename NumericTraits< PixelType >::FloatType   InternalRealType;
  //typedef typename Image<InternalRealType, itkGetStaticConstMacro(ImageDimension) >   RealImageType;

  // set all of the scales the same
  void SetScale(ScalarRealType scale);
  itkSetMacro(Scale, RadiusType);
  itkGetConstReferenceMacro(Scale, RadiusType);

  enum ParabolicAlgorithm {
    NOCHOICE = 0,         // decices based on scale - experimental
    CONTACTPOINT = 1, // sometimes faster at low scale
    INTERSECTION = 2  // default
  };
  /** 
   * Set/Get the method used. Choices are contact point or
   * intersection. Intersection is the default. Contact point can be 
   * faster at small scales.
   */

  itkSetMacro(ParabolicAlgorithm, int);
  itkGetConstReferenceMacro(ParabolicAlgorithm, int);

  /**
   * Set/Get whether the scale refers to pixels or world units -
   * default is false
   */
  itkSetMacro(UseImageSpacing, bool);
  itkGetConstReferenceMacro(UseImageSpacing, bool);
  itkBooleanMacro(UseImageSpacing);
  /** Image related typedefs. */
 
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimension,
                  (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),itkGetStaticConstMacro(OutputImageDimension)>));

  itkConceptMacro(Comparable,
		  (Concept::Comparable<PixelType>));

  /** End concept checking */
#endif

protected:
  ParabolicErodeDilateImageFilter();
  virtual ~ParabolicErodeDilateImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  /** Generate Data */
  void GenerateData( void );
  unsigned int SplitRequestedRegion(unsigned int i, unsigned int num,
    OutputImageRegionType & splitRegion);

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId );

  void GenerateInputRequestedRegion() throw(InvalidRequestedRegionError);
  // Override since the filter produces the entire dataset.
  void EnlargeOutputRequestedRegion(DataObject *output);

  bool m_UseImageSpacing;
  int m_ParabolicAlgorithm;

private:
  ParabolicErodeDilateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  RadiusType m_Scale;
  typename TInputImage::PixelType m_Extreme;

  int m_MagnitudeSign;
  int m_CurrentDimension;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParabolicErodeDilateImageFilter.txx"
#endif


#endif
