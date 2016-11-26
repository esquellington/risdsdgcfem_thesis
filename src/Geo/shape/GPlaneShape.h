#ifndef GEO_SHAPE_GPLANESHAPE_H
#define GEO_SHAPE_GPLANESHAPE_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()

namespace geo {

/*! 2/3d plane shape. */
template <unsigned D>
class GPlaneShape: public GIShapeD<D>
{
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;
    typedef Real sdof_type; //!< Useless but cannot be void
    typedef mal::GVec<Real,D> vec_type;
    //@}

    typedef GPlaneShape<D> self_type;

public:
    inline GPlaneShape() : m_Normal(vec_type::Zero()), m_CoeffD(0), m_bIsHalfSpace(false) {}
    inline GPlaneShape( const vec_type &normal, Real coeff_d, bool b_half_space )
    : m_Normal(normal), m_CoeffD(coeff_d), m_bIsHalfSpace(b_half_space) {}
    inline ~GPlaneShape() {}

    //!\name Shape params
    //@{
    finline void Init( const vec_type &normal, Real coeff_d, bool b_half_space )
        {
            m_Normal=normal;
            m_CoeffD=coeff_d;
            m_bIsHalfSpace =  b_half_space;
        }
    finline void SetNormalAndCoeffD( const vec_type &normal, Real coeff_d ) { m_Normal=normal; m_CoeffD=coeff_d; }
    finline const vec_type &GetNormal() const { return m_Normal; }
    finline Real GetCoeffD() const { return m_CoeffD; }
    finline void SetHalfSpace( bool b_half_space ) { m_bIsHalfSpace = b_half_space; }
    finline bool IsHalfSpace() const { return m_bIsHalfSpace; }
    //@}

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return shape_type_of<self_type>(); }
    unsigned int GetNumSDOF() const { return 0; }
    const sdof_type *GetVecDefaultSDOF() const { return 0; }
    const Real *GetVecDefaultDOF() const { return 0; }

    void ComputeBVD( bv::GBoundingVolumeD<D> &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
        {
            switch( bv.GetType() )
            {
            case bv::eBV_Sphere2:
            case bv::eBV_Sphere3:
                bv.template As< bv::GSphere<D> >().SetPosRadius( vec_type::Zero(), Real(1000.0f) );
                break;
            case bv::eBV_AABB2:
            case bv::eBV_AABB3:
                bv.template As< bv::GAABB<D> >().SetPosHalfSizes( vec_type::Zero(), vec_type(1000.0f) ); //\todo If axis-aligned, compute proper AABB!
                break;
            case bv::eBV_LSS2:
            case bv::eBV_LSS3:
                bv.template As< bv::GLineSweptSphere<D> >().SetPosRadius( vec_type::Zero(), vec_type::Zero(), Real(1000.0f) );
                break;
            case bv::eBV_Void: break;
            case bv::eBV_Infinite: break;
            default:
                GEO_ASSERT( false ); //wrong type or dimension
                break;
            }
        }

    GIDomainSamplerD<D> *CreateDomainSampler() const { return 0; }
    //@}

protected:
    vec_type m_Normal;
    Real m_CoeffD;
    bool m_bIsHalfSpace;
};

} //namespace geo

#endif // GEO_SHAPE_GPLANESHAPE_H
