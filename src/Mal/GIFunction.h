#ifndef MAL_G_IFUNCTION_H
#define MAL_G_IFUNCTION_H

#include <Mal/Config.h>
#include <Mal/GVec.h> //TEMPORAL for GTranslationGF
#include <Mal/GTransform.h> //TEMPORAL for GTranslationGF

namespace mal
{

/*! Generic Function:

     F : I^n -> O^m

  Where I,O are the input domain and the output image sets.

  Can be applied to single values or arrays.
  
  Allows the implementation of generic Coordinate Transforms between
  DSH and Geoms
*/
template <typename InputT, typename OutputT>
class GIFunction
{
public:
    typedef InputT input_type;
    typedef OutputT output_type;
    
public:
    GIFunction() {}
    virtual ~GIFunction() {}
    
    //! \name Explicit (stateless) protocol
    //@{
    virtual void F( const input_type &x, output_type &fx ) const = 0;
    virtual void InvF( input_type &x, const output_type &fx ) const = 0;
    virtual void F( const input_type *vec_x, output_type *vec_fx, int count ) const = 0;
    virtual void InvF( input_type *vec_x, const output_type *vec_fx, int count ) const = 0;
    //@}
};

/*! Interface for evaluation of a function with input/output
    internally bound to specific memory locations.
*/
class IBoundFunction
{
public:
    IBoundFunction() {}
    virtual ~IBoundFunction() {}

    //! \name Reflection
    //@{
    virtual TypeID GetInputTypeID() const = 0;
    virtual TypeID GetOutputTypeID() const = 0;
    //@}
    
    //! \name Eval Interface
    //@{
    virtual void Eval() const = 0;
    virtual void InvEval() const = 0;
    //@}
};

/*! GIFunction bound to specific input/output memory locations  
    - Can be used as a FunctionT, as it inherits from FunctionT and
      forwards required typedefs.
      
    \note DON'T forget to re-bind count if input/output arrays change
    in size!
      
    \todo MAY require independent non-zero strides for input/output
    arrays. This would be most-efficiently implemented in 2 separate
    specialized classes:
    - GBoundFunctionVariableStride => Strides as parameters
    - GBoundFunctionConstantStride => Strides as template arguments
*/
template <typename FunctionT>
class GBoundFunction: public FunctionT, public IBoundFunction
{
public:
    typedef FunctionT function_type;
    typedef typename FunctionT::input_type input_type;
    typedef typename FunctionT::output_type output_type;
    
public:
    GBoundFunction() : m_pIn(0), m_pOut(0), m_Count(0) {}
    virtual ~GBoundFunction() {}
    
    //! \name Bind/Eval Protocol
    //@{
    finline bool Bind( input_type *p_x, output_type *p_fx, int count = 1 ) { m_pIn=p_x; m_pOut=p_x; m_Count=count; return true; }
    finline bool BindIn( input_type *p_x ) { m_pIn=p_x; return true; }
    finline bool BindOut( output_type *p_fx ) { m_pOut=p_fx; return true; }
    finline bool BindCount( int count ) { m_Count = count; return true; }
    finline input_type *GetInput() const { return m_pIn; }
    finline output_type *GetOutput() const { return m_pOut; }
    finline unsigned int GetCount() const { return m_Count; }
    finline void _Eval() const { if( 1==m_Count ) F(*m_pIn,*m_pOut); else F(m_pIn,m_pOut,m_Count); }
    finline void _InvEval() const { if( 1==m_Count ) InvF(*m_pIn,*m_pOut); else InvF(m_pIn,m_pOut,m_Count); }
    //@}

    //! \name IBoundFunction implementation
    //@{
    TypeID GetInputTypeID() const { return pla_type_id<input_type>::value; }
    TypeID GetOutputTypeID() const { return pla_type_id<output_type>::value; }
    void Eval() const { _Eval(); }
    void InvEval() const { _InvEval(); }
    //@}
    
private:
    input_type *m_pIn;
    output_type *m_pOut;
    int m_Count;
};


/*! GIFunction composition F1(F2(x))
    - Can be used as a GIFunction itself, as it forwards required typedefs    
    - Uses an intermediate local value to store F1(x) and InvF2(x) and perform composition.
    \todo Check FunctionT1::output_type == FunctionT2::input_type (Sake uses TypeId<>...)
*/
template <typename FunctionT1, typename FunctionT2>
class GComposedFunction: public GIFunction< typename FunctionT1::input_type, typename FunctionT2::output_type >
{
    MAL_STATIC_ASSERT( sizeof(typename FunctionT1::output_type) == sizeof(typename FunctionT2::input_type) );
                       
public:
    typedef FunctionT1 function1_type;
    typedef FunctionT2 function2_type;
    typedef typename FunctionT1::input_type input_type;
    typedef typename FunctionT2::output_type output_type;
    
public:
    GComposedFunction() {}
    virtual ~GComposedFunction() {}

    //! \name Composed sub-functions retrieval
    //@{
    function1_type &GetF1() { return m_F1; }
    function2_type &GetF2() { return m_F2; }
    const function1_type &GetF1() const { return m_F1; }
    const function2_type &GetF2() const { return m_F2; }
    //@}

    //! \name Function evaluation
    //{@
    finline void _F( const input_type &x, output_type &fx ) const { m_F1.F(x,m_TmpValue); m_F2.F(m_TmpValue,fx); }
    finline void _InvF( input_type &x, const output_type &fx ) const { m_F2.InvF(m_TmpValue,fx); m_F1.InvF(x,m_TmpValue); }
    //@}
    
    //! \name GIFunction implementation
    //@{
    void F( const input_type &x, output_type &fx ) const { _F(x,fx); }
    void InvF( input_type &x, const output_type &fx ) const { _InvF(x,fx); }
    void F( const input_type *vec_x, output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _F(vec_x[i],vec_fx[i]); }
    void InvF( input_type *vec_x, const output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _InvF(vec_x[i],vec_fx[i]); }
    //@}
    
private:
    function1_type m_F1;
    function2_type m_F2;
    mutable typename FunctionT1::output_type m_TmpValue;
};

//----------------------------------------------------------------------
//---- Specific GIFunctions
//----------------------------------------------------------------------

/*! N-dimensional translation */
template <typename T, unsigned N>
class GTranslationGF: public GIFunction< GVec<T,N>, GVec<T,N> >
{
public:
    static const unsigned int cDimension = N;
    typedef GVec<T,N> vec_type;
    typedef GIFunction< GVec<T,N>, GVec<T,N> > parent_type;
    typedef typename parent_type::input_type input_type;
    typedef typename parent_type::output_type output_type;
    
public:
    GTranslationGF() {}

    //! \name Function evaluation
    //{@
    void SetParams( const vec_type &translation ) { m_Translation = translation; }
    finline void _F( const input_type &x, output_type &fx ) const { fx = x + m_Translation; }
    finline void _InvF( input_type &x, const output_type &fx ) const { x = fx - m_Translation; }
    //@}
    
    //! \name GIFunction implementation
    //@{
    void F( const input_type &x, output_type &fx ) const { _F(x,fx); }
    void InvF( input_type &x, const output_type &fx ) const { _InvF(x,fx); }
    void F( const input_type *vec_x, output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _F(vec_x[i],vec_fx[i]); }
    void InvF( input_type *vec_x, const output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _InvF(vec_x[i],vec_fx[i]); }
    //@}
        
private:
    vec_type m_Translation;
};

/*! N-dimensional rigid body transform */
template <typename T, unsigned N>
class GTransformGF: public GIFunction< GTransform<T,N>, GTransform<T,N> >
{
public:
    static const unsigned int cDimension = N;
    typedef GTransform<T,N> transform_type;
    typedef GIFunction< GTransform<T,N>, GTransform<T,N> > parent_type;
    typedef typename parent_type::input_type input_type;
    typedef typename parent_type::output_type output_type;
    
public:
    GTransformGF() {}

    //! \name Function evaluation
    //{@
    void SetParams( const transform_type &transform ) { m_Transform = transform; m_InvTransform = transform.Inverse(); }
    finline void _F( const input_type &x, output_type &fx ) const { fx = m_Transform*x; }
    finline void _InvF( input_type &x, const output_type &fx ) const { x = m_InvTransform*fx; }
    //@}
    
    //! \name GIFunction implementation
    //@{
    void F( const input_type &x, output_type &fx ) const { _F(x,fx); }
    void InvF( input_type &x, const output_type &fx ) const { _InvF(x,fx); }
    void F( const input_type *vec_x, output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _F(vec_x[i],vec_fx[i]); }
    void InvF( input_type *vec_x, const output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _InvF(vec_x[i],vec_fx[i]); }
    //@}

private:
    transform_type m_Transform;
    transform_type m_InvTransform;
};

/*! Generic Copy inputs to outputs (with automatic cast) */
template <typename InputT, typename OutputT>
class GCopyGF: public GIFunction< InputT, OutputT >
{
public:
    typedef GIFunction< InputT, OutputT > parent_type;
    typedef typename parent_type::input_type input_type;
    typedef typename parent_type::output_type output_type;
    
public:
    GCopyGF() {}

    //! \name Function evaluation
    //{@
    finline void _F( const input_type &x, output_type &fx ) const { fx = output_type(x); }
    finline void _InvF( input_type &x, const output_type &fx ) const { x = input_type(fx); }
    //@}
    
    //! \name GIFunction implementation
    //@{
    void F( const input_type &x, output_type &fx ) const { _F(x,fx); }
    void InvF( input_type &x, const output_type &fx ) const { _InvF(x,fx); }
    void F( const input_type *vec_x, output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _F(vec_x[i],vec_fx[i]); }
    void InvF( input_type *vec_x, const output_type *vec_fx, int count ) const { for(int i=0;i<count;i++) _InvF(vec_x[i],vec_fx[i]); }
    //@}
};

//--- ALTRES coord transforms possibles
//CartesianToPolar2
//CartesianToSpherical3
//CartesianToCylindric3
//CopyPositionsAndComputeTangents : vec3 -> pair<vec3,vec3>

//\name Common specific template instantiations
//@{
typedef GTranslationGF<float,2> TranslationGF2f;
typedef GTranslationGF<float,3> TranslationGF3f;
typedef GBoundFunction< TranslationGF2f > TranslationGBF2f;
typedef GBoundFunction< TranslationGF3f > TranslationGBF3f;

typedef GTransformGF<float,2> TransformGF2f;
typedef GTransformGF<float,3> TransformGF3f;
typedef GBoundFunction< TransformGF2f > TransformGBF2f;
typedef GBoundFunction< TransformGF3f > TransformGBF3f;
//@}

/* Tests
static TranslationGF2f m_T2f;
static TransformGBF3f m_TF3f;
static GComposedFunction< TransformGBF3f, GComposedFunction<TransformGF3f,TransformGBF3f> > m_KKKGBF3f;
static GCopyGF< GVec<double,3> ,GVec<float,3> > m_CI2F;
*/

} // namespace mal

#endif //MAL_G_IFUNCTION_H
