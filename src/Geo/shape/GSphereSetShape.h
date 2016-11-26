#ifndef GEO_SHAPE_GSPHERESETSHAPE_H
#define GEO_SHAPE_GSPHERESETSHAPE_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()

namespace geo {

/*! 2/3d sphere shape. */
template <unsigned D>
class GSphereSetShape: public GIShapeD<D>
{
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;
    typedef mal::GVec<Real,D> sdof_type;
    typedef mal::GVec<Real,D> vec_type;
    //@}

    typedef GSphereSetShape<D> self_type;

public:
    GSphereSetShape( unsigned int count = 1, Real radius = Real(1) )
    : m_NumSpheres(count), m_Radius(radius)
        { GEO_ASSERT(count>0); }
    ~GSphereSetShape() {}

    //!\name Shape params
    //@{
    finline void SetCount( unsigned int count ) { m_NumSpheres = count; }
    finline unsigned int GetCount() const { return m_NumSpheres; }

    finline void SetRadius( Real radius ) { m_Radius = radius; }
    finline Real GetRadius() const { return m_Radius; }
    //@}

    //!\name Dynamic IShape implementation
    //@{
    EShapeType GetType() const { return shape_type_of<self_type>(); }
    unsigned int GetNumSDOF() const { return m_NumSpheres; }
    const sdof_type *GetVecDefaultSDOF() const { return 0; } //\todo #SDOF > 0 but there are NO DEFAULT SDOF
    const Real *GetVecDefaultDOF() const { return 0; }
    GIDomainSamplerD<D> *CreateDomainSampler() const { return 0; }

    //@}

    void ComputeBVD( bv::GBoundingVolumeD<D> &bv, const transform_type &transform, const sdof_type *vec_sdof ) const
        {
            switch( bv.GetType() )
            {
            case bv::eBV_Sphere2:
            case bv::eBV_Sphere3:
                GEO_ASSERT( false ); //Not yet implemented
                //bv.template As<bv::Sphere2>().SetPosRadius( transform.m_Pos, m_Radius );
                break;
            case bv::eBV_AABB2:
            case bv::eBV_AABB3:
                {
                    bv::GAABB<D> aabb( vec_sdof[0], vec_type(m_Radius) );
                    for( unsigned int i=1; i < m_NumSpheres; i++ )
                        aabb.Merge( vec_sdof[i], m_Radius );
                    bv.template As< bv::GAABB<D> >() = aabb;
                }
                break;
            case bv::eBV_LSS2:
            case bv::eBV_LSS3:
                GEO_ASSERT( false ); //Not yet implemented
                //bv.As<bv::LSS2>().SetPosRadius( transform.m_Pos, transform.m_Pos, m_Radius );
                break;
            case bv::eBV_Void: break;
            case bv::eBV_Infinite: break;
            default:
                GEO_ASSERT( false ); //wrong type or dimension
                break;
            }
        }

protected:
    unsigned int m_NumSpheres;
    Real m_Radius;
};

} //namespace geo

#endif // GEO_SHAPE_GSPHERESETSHAPE_H
