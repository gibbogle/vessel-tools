#ifndef __itkBinaryCloseParaImageFilter_h
#define __itkBinaryCloseParaImageFilter_h

#include "itkParabolicErodeImageFilter.h"
#include "itkParabolicDilateImageFilter.h"
#include "itkGreaterEqualValImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

namespace itk
{

/**
 * \class BinaryCloseParaImageFilter
 * \brief Class for binary morphological opening operation.
 *
 * This class uses the parabolic morphology operations to do very
 * efficient erosions by circles/spheres. The operations are efficient
 * because the underlying parabolic operations are separable and
 * the operations are implicitly short circuited in comparison to a
 * full distance transform approach.
 *
 * The basic idea is that a binary erosion or dilation by a circle/sphere
 * can be carried out by thresholding a distance transform. By using
 * the parabolic filters we can avoid computing the entire distance
 * transform and instead only compute the subset we are interested
 * in.
 *
 * Note that the circles and spheres may not be quite what you
 * expect, because this class doesn't explicitly use Bresenham circles
 * as most of the others do. A voxel's centre needs to be less than or
 * equal to the circle radius, rather than any part of the voxel
 * inside the circle.
 *
 * This filter was developed as a result of discussions with
 * M.Starring on the ITK mailing list.
 * 
 * \sa itkParabolicErodeImageFilter
 *
 * \author Richard Beare, Department of Medicine, Monash University,
 * Australia.  <Richard.Beare@med.monash.edu.au>
**/


template <typename TInputImage,
          typename TOutputImage= TInputImage >
class ITK_EXPORT BinaryCloseParaImageFilter:
    public ImageToImageFilter<TInputImage,TOutputImage>

{

public:
  /** Standard class typedefs. */
  typedef BinaryCloseParaImageFilter  Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                   Pointer;
  typedef SmartPointer<const Self>        ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(BinaryCloseParaImageFilter, ImageToImageFilter);


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

  typedef typename NumericTraits< PixelType >::FloatType   InternalRealType;
  // perhaps a bit dodgy, change to int if you want to do enormous
  // binary operations
  typedef  short     InternalIntType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  typedef typename itk::FixedArray<ScalarRealType, TInputImage::ImageDimension> RadiusType;

  void SetRadius(ScalarRealType radius);
  itkSetMacro(Radius, RadiusType);
  itkGetConstReferenceMacro(Radius, RadiusType);

  void SetUseImageSpacing(bool g)
  {
    m_RectErode->SetUseImageSpacing(g);
    m_RectDilate->SetUseImageSpacing(g);
    m_CircErode->SetUseImageSpacing(g);
    m_CircDilate->SetUseImageSpacing(g);
  }
  /**
   * Set/Get whether the erosion is circular/rectangular -
   * default is true (circular)
   */
  itkSetMacro(Circular, bool);
  itkGetConstReferenceMacro(Circular, bool);
  itkBooleanMacro(Circular);

 /** A safe border is added to input image to avoid borders effects
   * and remove it once the closing is done */
  itkSetMacro(SafeBorder, bool);
  itkGetConstReferenceMacro(SafeBorder, bool);
  itkBooleanMacro(SafeBorder);

  /** Image related typedefs. */

  /* add in the traits here */

protected:
  void GenerateData( void );

  BinaryCloseParaImageFilter();
  virtual ~BinaryCloseParaImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  typedef typename itk::Image<InternalRealType, InputImageType::ImageDimension> InternalRealImageType;
  typedef typename itk::Image<InternalIntType, InputImageType::ImageDimension> InternalIntImageType;
  typedef typename itk::ParabolicErodeImageFilter<TInputImage, InternalRealImageType> CircErodeType;
  typedef typename itk::ParabolicErodeImageFilter<TInputImage, InternalIntImageType> RectErodeType;
  typedef typename itk::ParabolicDilateImageFilter<OutputImageType, InternalRealImageType> CircDilateType;
  typedef typename itk::ParabolicDilateImageFilter<OutputImageType, InternalRealImageType> RectDilateType;

  typedef typename itk::GreaterEqualValImageFilter<InternalRealImageType, OutputImageType> CCastTypeA;
  typedef typename itk::BinaryThresholdImageFilter<InternalRealImageType, OutputImageType> CCastTypeB;

  typedef typename itk::GreaterEqualValImageFilter<InternalIntImageType, OutputImageType> RCastTypeA;
  typedef typename itk::BinaryThresholdImageFilter<InternalRealImageType, OutputImageType> RCastTypeB;
  
private:
  BinaryCloseParaImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  RadiusType m_Radius;
  bool m_Circular;
  bool m_SafeBorder;

  typename CircErodeType::Pointer m_CircErode;
  typename CircDilateType::Pointer m_CircDilate;
  typename CCastTypeA::Pointer m_CircCastA;
  typename CCastTypeB::Pointer m_CircCastB;

  typename RectErodeType::Pointer m_RectErode;
  typename RectDilateType::Pointer m_RectDilate;
  typename RCastTypeA::Pointer m_RectCastA;
  typename RCastTypeB::Pointer m_RectCastB;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBinaryCloseParaImageFilter.txx"
#endif


#endif //__itkBinaryCloseParaImageFilter_h
