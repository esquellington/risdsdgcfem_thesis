#ifndef GEO_SHAPE_GCAPSULESHAPE_H
#define GEO_SHAPE_GCAPSULESHAPE_H

#include <Geo/shape/IShape.h>

// BV types for ComputeBV()
#include <Geo/bv/bv.h>

//#include <Mal/GMatUtils.h>

namespace geo {

/*! 2/3d capsule shape. */
template <unsigned D>
class GCapsuleShape: public GIShapeD<D>
{
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;
    typedef Real sdof_type; //!< Useless but cannot be void
    typedef mal::GVec<Real,D> vec_type;
    //@}

    typedef GCapsuleShape<D> self_type;

public:
    GCapsuleShape( Real radius = Real(1), Real half_height = Real(0.5) ) : m_Radius(radius), m_HalfHeight(half_height) {}
    ~GCapsuleShape() {}

    //!\name Shape params
    //@{
    finline void SetRadiusAndHalfHeight( Real radius, Real half_height ) { m_Radius = radius; m_HalfHeight = half_height; }
    finline void SetRadius( Real radius ) { m_Radius = radius; }
    finline Real GetRadius() const { return m_Radius; }
    finline void SetHalfHeight( Real half_height ) { m_HalfHeight = half_height; }
    finline Real GetHalfHeight() const { return m_HalfHeight; }
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
                bv.template As< bv::GSphere<D> >().SetPosRadius( transform.m_Pos, m_HalfHeight+m_Radius );
                break;
            case bv::eBV_AABB2:
            case bv::eBV_AABB3:
                bv.template As< bv::GAABB<D> >().SetPosHalfSizes( transform.m_Pos, vec_type(m_HalfHeight+m_Radius) ); //\todo could be smaller!!
                break;
            case bv::eBV_LSS2:
            case bv::eBV_LSS3:
                bv.template As< bv::GLineSweptSphere<D> >().SetPosRadius( transform.m_Pos - m_HalfHeight*mal::Column(1,transform.m_Rot),
                                                                          transform.m_Pos + m_HalfHeight*mal::Column(1,transform.m_Rot),
                                                                          m_Radius );
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
    Real m_Radius;
    Real m_HalfHeight;
};

} //namespace geo

#endif // GEO_SHAPE_GCAPSULESHAPE_H
