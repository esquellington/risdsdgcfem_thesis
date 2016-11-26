#ifndef GEO_SHAPE_GBOXSHAPE_H
#define GEO_SHAPE_GBOXSHAPE_H

#include <Geo/shape/IShape.h>
#include <Geo/bv/bv.h> // BV types for ComputeBV()

namespace geo {

/*! 2/3d box shape. */
template <unsigned D>
class GBoxShape: public GIShapeD<D>
{
public:
    //! \name Mandatory generic interface for IShape
    //@{
    static const unsigned int cDimension = D;
    typedef mal::GTransform<Real,D> transform_type;
    typedef Real sdof_type; //!< Useless but cannot be void
    typedef mal::GVec<Real,D> vec_type;
    //@}

    typedef GBoxShape<D> self_type;

public:
    GBoxShape( const vec_type &half_sizes = vec_type(Real(1)) ) : m_HalfSizes(half_sizes) {}
    ~GBoxShape() {}


    //!\name Shape params
    //@{
    finline void SetHalfSizes( const vec_type &half_sizes ) { m_HalfSizes = half_sizes; }
    finline const vec_type &GetHalfSizes() const { return m_HalfSizes; }
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
                bv.template As< bv::GSphere<D> >().SetPosRadius( transform.m_Pos, m_HalfSizes.Norm() );
                break;
            case bv::eBV_AABB2:
            case bv::eBV_AABB3:
                {
                    /* 2d case
                    mal::GVec<Real,2> aabb_half_sizes( mal::Abs(transform.m_Rot(0,0)*m_HalfSizes[0])
                                                       + mal::Abs(transform.m_Rot(0,1)*m_HalfSizes[1]),
                                                       mal::Abs(transform.m_Rot(1,0)*m_HalfSizes[0])
                                                       + mal::Abs(transform.m_Rot(1,1)*m_HalfSizes[1]) );
                    bv.As<bv::AABB2>().SetPosHalfSizes( transform.m_Pos, aabb_half_sizes );
                    */
                    //TEMP: compact but inefficient... vec_type aabb_half_sizes( mal::Abs(transform) * mal::Abs(m_HalfSizes) );
                    vec_type aabb_half_sizes(0);
                    for( unsigned int i=0; i<D; i++ )
                        for( unsigned int j=0; j<D; j++ )
                            aabb_half_sizes[i] += mal::Abs( transform.m_Rot(i,j) * m_HalfSizes[i] );
                    bv.template As<bv::GAABB<D> >().SetPosHalfSizes( transform.m_Pos, aabb_half_sizes );
                }
                break;
            case bv::eBV_LSS2:
            case bv::eBV_LSS3:
                GEO_ASSERT( false ); //suboptimal LSS
                bv.template As< bv::GLineSweptSphere<D> >().SetPosRadius( transform.m_Pos, transform.m_Pos, m_HalfSizes.Norm() );
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
    vec_type m_HalfSizes;
};

} //namespace geo

#endif // GEO_SHAPE_GBOXSHAPE_H
