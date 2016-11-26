#ifndef GEO_IOBJECT_H
#define GEO_IOBJECT_H

#include <Geo/Config.h>
#include <Geo/shape/IShape.h>
#include <string.h> //TEMP: for memcpy in templated funcs

namespace geo {

namespace bv { class IBoundingVolume; template<unsigned D> class GBoundingVolumeD; }

class IObject;
class EmbeddingData;

/*! Geometric Object Interface. */
class IObject
{
public:
    IObject() {}
    virtual ~IObject() {}

    virtual unsigned int GetDimension() const = 0;

    virtual IShape* GetShapeInterface() = 0;
    virtual const IShape* GetShapeInterface() const = 0;

    virtual void ResetDOF() = 0;
    virtual Real* GetVecDOF_WriteOnly() = 0; //\note Write-Only access, Reads may be out of date!
    virtual const Real* GetVecDOF() const = 0;

    virtual void ComputeBV( bv::IBoundingVolume& bv ) const = 0; //\todo This should NOT be part of the interface... could use an external method to do this, IBV does not make sense, dimension-specific BVD does, move it to IObject2/3!!

    //!\name Embedding
    //@{
    //! SIMPLEST embedding interface, automatically embed using current configuration and default embedding type, non-shared embedding definition. Pass 0-deformer to unembed
    virtual void Embed( const IObject* p_deformer, EEmbeddingMethod em = eEM_Barycentric ) = 0;

    /* Alternative API, similar to ObjectFactory but using an EmbeddingFactory/Library, instead of a ShapeFactory/Library
    virtual bool Embed( const IObject* p_deformer, const EmbeddingDef& embedding_def ) = 0; // non-shared
    virtual bool Embed( const IObject* p_deformer, const EmbeddingID& embedding_id ) = 0; // shared
    //virtual bool Embed( const IObject* p_deformer, const IEmbedding* p_embedding ) = 0; // shared
    */
    //@}

    //!\name BVH
    //@{
    //!
    //virtual const IBVH* CreateBVH() {}; //\todo several BVH types... should be pure virtual but NOT HERE, move to IObject2/3
    //@}
};

/*! Dimension-independent stuff */
class IObjectBase: public IObject
{
public:
    finline IObjectBase() : m_TimeStamp(0), m_pEmbeddingData(0) {}
    ~IObjectBase();

    //!\name Embedding
    //@{
    void Embed( const IObject* p_deformer, EEmbeddingMethod em = eEM_Barycentric );
    //@}

protected:
    finline uint32 GetTimeStamp() const { return m_TimeStamp; }
    finline void Touch() { m_TimeStamp++; }

    finline const EmbeddingData* GetEmbeddingData() const { return m_pEmbeddingData; }
    finline bool IsEmbedded() { return 0 != m_pEmbeddingData; }

    bool SyncEmbedding() const;

protected:
    uint32 m_TimeStamp;
    mutable EmbeddingData* m_pEmbeddingData; //\todo Embedding with lazy evaluation requires this to be mutable... is it a good idea?



    /*\todo INSTEAD of adding EMB and BVH to ALL IObjects, consider
     * having a single Annotations structure that is only allocated if
     * any Annotation is added (EMB,BVH,...) and stores all their
     * pointers and lifetime data.
     */
};

/*! IGeom dimension-specific interface
    - Manages transform
*/
template <unsigned D>
class GIObjectD: public IObjectBase
{
public:
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;

public:
    GIObjectD() : m_Transform( transform_type::Identity() ) {}
    ~GIObjectD() {}

    unsigned int GetDimension() const { return D; }

    //!\name Transform access
    //@{
    void SetTransform( const transform_type& transform ) { if( !IObjectBase::IsEmbedded() ) { m_Transform = transform; IObjectBase::Touch(); }
                                                           else GEO_LOG_WARNING("GObjectSDOF is embedded, cannot SetTransform() directly"); }
    //!\todo Ugly hack to allow settransform from IEmbedding::Apply() without accessing m_Transform directly...
    void SetTransform_Embedded( const transform_type& transform ) { if( IObjectBase::IsEmbedded() ) { m_Transform = transform; IObjectBase::Touch(); }
                                                                    else GEO_LOG_WARNING("GObjectSDOF is NOT embedded, cannot set SetTransform_Embedded()"); }
    const transform_type& GetTransform() const { IObjectBase::SyncEmbedding(); return m_Transform; }
    //@}

    virtual void ComputeBVD( bv::GBoundingVolumeD<D>& bv ) const = 0;

    void ComputeBV( bv::IBoundingVolume& bv ) const { ComputeBVD( static_cast< bv::GBoundingVolumeD<D>& >(bv) ); }

    /*BVH: If there is no "dimensionless" IBVH concept, move it to GIObjectD<D> and use IBVH2, IBVH3...
    //bool SyncBVH() const { GEO_ASSERT(0!=m_pBVH); GetShape()->RefitBVH( m_pBVH, transform, vec_sdof ); }
    //IBVH* GetBVH() const { SyncBVH(); return m_pBVH; }

     The BVH should be created AUTOMATICALLY IFF the Shape has a valid
     BVH (IShape::GetDefaultBVH()??), the BVH geometry DOF should
     allocated in the IObject. IObject timestamp should be used to
     refit it lazily.
    */

protected:
    transform_type m_Transform;
    /*\todo Should be in m_pAnnotations->m_pBVH
    mutable IBVH* m_pBVH; //\todo BVH with lazy evaluation requires this to be mutable... is it a good idea?
    */
};

//! \name Usual dimensions
//@{
typedef GIObjectD<2> IObject2;
typedef GIObjectD<3> IObject3;
//@}

/*! Structured-DOF specific
  - Manages DOF/SDOF allocation
*/
template <unsigned D, typename SDOF_T>
class GObjectSDOF: public GIObjectD<D>
{
public:
    static const unsigned int cDimension = D;
    typedef SDOF_T sdof_type;

public:
    finline GObjectSDOF() : m_NumAllocSDOF(0), m_vecSDOF(0) {}
    finline virtual ~GObjectSDOF() { DeallocSDOF(); }

    //!\name Structured DOF access
    //@{
    sdof_type& GetSDOF( int sdof_id ) { IObjectBase::SyncEmbedding(); IObjectBase::Touch(); return m_vecSDOF[sdof_id]; }
    const sdof_type& GetSDOF( int sdof_id ) const { IObjectBase::SyncEmbedding(); return m_vecSDOF[sdof_id]; }
    finline sdof_type* GetVecSDOF() { IObjectBase::SyncEmbedding(); IObjectBase::Touch(); return m_vecSDOF; }
    finline const sdof_type* GetVecSDOF() const { IObjectBase::SyncEmbedding(); return m_vecSDOF; }
    //@}

    //!\name Real DOF access
    //@{
    finline Real& GetDOF( int dof_id ) { IObjectBase::SyncEmbedding(); IObjectBase::Touch(); return GetVecDOF_WriteOnly()[dof_id]; }
    finline const Real &GetDOF( int dof_id ) const { IObjectBase::SyncEmbedding(); return GetVecDOF()[dof_id]; }
    finline Real* GetVecDOF_WriteOnly() { IObjectBase::Touch(); return reinterpret_cast<Real*>(&m_vecSDOF[0]); }
    finline const Real *GetVecDOF() const { IObjectBase::SyncEmbedding(); return reinterpret_cast<const Real*>(&m_vecSDOF[0]); }
    //@}

protected:
    inline void DeallocSDOF() { if( m_vecSDOF ) delete [] m_vecSDOF; m_NumAllocSDOF = 0; m_vecSDOF = 0; }
    inline void ReallocSDOF( unsigned int num_sdof ) //!< Handles 0-DOF properly
        {
            if( m_NumAllocSDOF < num_sdof )
            {
                DeallocSDOF();
                if( num_sdof > 0 )
                {
                    m_vecSDOF = new sdof_type[num_sdof];
                    m_NumAllocSDOF = num_sdof;
                }
            }
            //Otherwise, there's already enough alloc mem
        }
    inline void SetSDOF( const sdof_type* vec_sdof, unsigned int num_sdof ) //!< Handles 0-DOF and 0-Default-DOF properly
        {
            GEO_ASSERT( num_sdof <= m_NumAllocSDOF );
            if( !IObjectBase::IsEmbedded() )
            {
                if( 0 != vec_sdof )
                    memcpy( m_vecSDOF, vec_sdof, num_sdof*sizeof(sdof_type) );
            }
            else
                GEO_LOG_WARNING("GObjectSDOF is embedded, cannot SetSDOF directly");
        }

protected:
    unsigned int m_NumAllocSDOF; //\todo This is redundant if always equal to shape SDOF count, we'll lose fast Realloc in some cases, but we'll avoid storing m_NumAllocSDOF for EVERY SINGLE IObject (spheres, etc...)
    sdof_type *m_vecSDOF;
};


/*! GObjectES (ES = Exclusive Shape) class template models an instantiated ShapeT in space.
  - Allocates ShapeT statically
  - The Shape instance cannot be changed
  - Allocates ShapeT instance data
    - transform
    - dof vector (internal dynamic allocation)
  \todo m_vecSDOF allocation could be optimized if another template
  parameter NDOF was passed. A dof_type[NDOF] array could be
  statically allocated and used if ShapeT::GetNumSDOF() <= NDOF
  (otherwise new/delete)
*/
template <typename ShapeT>
class GObjectES: public GObjectSDOF<ShapeT::cDimension,typename ShapeT::sdof_type>
{
public:
    static const unsigned int cDimension = ShapeT::cDimension;
    typedef ShapeT shape_type;
    typedef typename ShapeT::transform_type transform_type;
    typedef typename ShapeT::sdof_type sdof_type;
    typedef GObjectSDOF<ShapeT::cDimension,typename ShapeT::sdof_type> osdof_type;

public:
    finline GObjectES() { ResetDOF(); IObjectBase::Touch(); } //\todo MUST Realloc and Set DOF to shape default, not just osdof_type::ReallocSDOF( m_Shape.GetNumSDOF() ); }
    virtual ~GObjectES() {} //delete [] m_vecSDOF done in base class

    finline ShapeT* GetShape() { IObjectBase::Touch(); return &m_Shape; }
    finline const ShapeT* GetShape() const { return &m_Shape; }

    //! \name GIObjectD implementation
    //@{
    IShape* GetShapeInterface() { IObjectBase::Touch(); return GetShape(); }
    const IShape* GetShapeInterface() const { return GetShape(); }
    void ResetDOF()
        {
            if( !IObjectBase::IsEmbedded() )
            {
                // Realloc and Copy SDOF. Handles 0-DOF and 0-Default-DOF shapes properly
                osdof_type::ReallocSDOF( m_Shape.GetNumSDOF() );
                osdof_type::SetSDOF( m_Shape.GetVecDefaultSDOF(), m_Shape.GetNumSDOF() );
                IObjectBase::Touch();
            }
            else
                GEO_LOG_WARNING("GObjectES is embedded, cannot ResetDOF");
        }
    void ComputeBVD( bv::GBoundingVolumeD<cDimension>& bv ) const { m_Shape.ComputeBVD( bv, osdof_type::m_Transform, osdof_type::m_vecSDOF ); }
    //@}

private:
    ShapeT m_Shape;
};

/*! GObjectSS (SS=SharedShape) class template models an instantiated *shared* ShapeT in space.
  - Holds a pointer to an externally allocated ShapeT definition.
  - The Shape instance can be changed
  - Allocates ShapeT instance data:
    - transform
    - dof vector (internal dynamic allocation)
  \todo If the Shape is SHARED, its acces should be STRICTLY CONST!! changes should require a Clone()!!
  \todo m_vecSDOF allocation could be optimized if another template
  parameter NDOF was passed. A dof_type[NDOF] array could be
  statically allocated and used if ShapeT::GetNumSDOF() <= NDOF
  (otherwise new/delete)
*/
template <typename ShapeT>
class GObjectSS: public GObjectSDOF<ShapeT::cDimension,typename ShapeT::sdof_type>
{
public:
    static const unsigned int cDimension = ShapeT::cDimension;
    typedef ShapeT shape_type;
    typedef typename ShapeT::transform_type transform_type;
    typedef typename ShapeT::sdof_type sdof_type;
    typedef GObjectSDOF<ShapeT::cDimension,typename ShapeT::sdof_type> osdof_type;

public:
    finline GObjectSS() : m_pShape(0) {}
    virtual ~GObjectSS() { RemoveShape(); }

    finline void SetShape( const ShapeT* p_shape )
    {
        GEO_ASSERT( m_pShape == 0 && osdof_type::m_vecSDOF == 0 );
        m_pShape = p_shape;
        ResetDOF(); //\note MUST Realloc and Set DOF, not just ReallocSDOF( m_pShape->GetNumSDOF() );
        IObjectBase::Touch();
    }
    finline void RemoveShape()
    {
        m_pShape = 0;
        osdof_type::DeallocSDOF();
        IObjectBase::Touch();
    }

    finline void ResetShape( const ShapeT* p_shape )
    {
        if( IObjectBase::IsEmbedded() ) IObjectBase::Embed(0); //\todo Unembed by now, consider re-embedding
        RemoveShape();
        SetShape(p_shape);
        IObjectBase::Touch();
    }

    //finline ShapeT* GetShape() { return m_pShape; } \note Read-only shape access
    finline const ShapeT* GetShape() const { return m_pShape; }

    //!\name GIObjectD implementation
    //@{
    IShape* GetShapeInterface() { return 0; } //\note Read-only shape access makes ANY non-const IObject to return 0, that is, to pick the non-const overload, even if the returned IShape* is only used as const... ugly, thus, IObject must be made const to retrieve shape
    const IShape* GetShapeInterface() const { return GetShape(); }
    void ResetDOF()
        {
            if( !IObjectBase::IsEmbedded() )
            {
                // Realloc and Copy SDOF. Handles 0-DOF and 0-Default-DOF shapes properly
                osdof_type::ReallocSDOF( m_pShape->GetNumSDOF() );
                this->SetSDOF( m_pShape->GetVecDefaultSDOF(), m_pShape->GetNumSDOF() );
                IObjectBase::Touch();
            }
            else
                GEO_LOG_WARNING("GObjectSS is embedded, cannot ResetDOF");
        }
    void ComputeBVD( bv::GBoundingVolumeD<cDimension>& bv ) const { m_pShape->ComputeBVD( bv, osdof_type::m_Transform, osdof_type::m_vecSDOF ); }
    //@}

private:
    const ShapeT* m_pShape;
};

} //namespace geo

#endif // GEO_IOBJECT_H
