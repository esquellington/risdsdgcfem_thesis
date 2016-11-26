#ifndef GEO_BV_GLINESWEPTSPHERE_H
#define GEO_BV_GLINESWEPTSPHERE_H

#include <Geo/bv/BoundingVolume.h>
#include <Geo/np/ClosestPoints.h>
#include <Geo/np/RayCast.h>

namespace geo {
namespace bv {

template <unsigned D>
class GLineSweptSphere: public GBoundingVolumeD<D>
{
public:
    enum EConstants { cDimension = D };
    typedef mal::GVec<Real,D> vec_type;

public:
    finline GLineSweptSphere() : GBoundingVolumeD<D>((D==2)?eBV_LSS2:eBV_LSS3) {}
    finline GLineSweptSphere( const vec_type &pos0, const vec_type &pos1, Real radius )
    : GBoundingVolumeD<D>((D==2)?eBV_LSS2:eBV_LSS3)
    , m_Pos0(pos0), m_Pos1(pos1), m_Radius(radius) {}
    finline ~GLineSweptSphere() {}

    finline void SetPos0( const vec_type &pos0 ) { m_Pos0 = pos0; IBoundingVolume::Touch(); }
    finline void SetPos1( const vec_type &pos1 ) { m_Pos1 = pos1; IBoundingVolume::Touch(); }
    finline void SetRadius( Real radius ) { m_Radius = radius; IBoundingVolume::Touch(); }
    finline void SetPos( const vec_type &pos0, const vec_type &pos1 ) { m_Pos0 = pos0; m_Pos1 = pos1; IBoundingVolume::Touch(); }
    finline void SetPosRadius( const vec_type &pos0, const vec_type &pos1, Real radius )
        { m_Pos0 = pos0; m_Pos1 = pos1; m_Radius = radius; IBoundingVolume::Touch(); }

    finline const vec_type &GetPos0() const { return m_Pos0; }
    finline const vec_type &GetPos1() const { return m_Pos1; }
    finline Real GetRadius() const { return m_Radius; }

    finline bool IsDegenerate( Real epsilon_length = mal::Epsilon<Real>() ) const { return (m_Pos1-m_Pos0).NormSq() < mal::Sq(epsilon_length); }

    finline void Extend( Real thickness ) { m_Radius += thickness; }

    finline void Merge( const vec_type &point ) { GEO_ASSERT(false); IBoundingVolume::Touch(); }
    finline void Merge( const vec_type &point, Real radius ) { GEO_ASSERT(false); IBoundingVolume::Touch(); }
    finline void Merge( const GLineSweptSphere &lss ) { GEO_ASSERT(false); IBoundingVolume::Touch(); } //!\todo Consider completely-enclosed case

public:
    vec_type m_Pos0;
    vec_type m_Pos1;
    Real m_Radius;
};

//! \name Usual dimensions
//@{
typedef GLineSweptSphere<2> LSS2;
typedef GLineSweptSphere<3> LSS3;

template <> struct bv_type_of<LSS2> { enum { value = eBV_LSS2 }; };
template <> struct bv_type_of<LSS3> { enum { value = eBV_LSS3 }; };
//@}

//Maximum distance between 2 points on the BV
template <unsigned D>
finline Real ComputeLargestExtent( const GLineSweptSphere<D>& lss )
{
    return mal::Norm(lss.m_Pos1-lss.m_Pos0) + 2*lss.m_Radius;
}

//\todo Ugly way to support D==2 and D==3...
template <unsigned D>
finline Real ComputeVolume( const GLineSweptSphere<D>& lss )
{
    Real height( mal::Norm(lss.m_Pos1-lss.m_Pos0) );
    if( D == 2 )
        return mal::Pi<Real>()*mal::Sq(lss.m_Radius) + height*lss.m_Radius; //sphere+rectangle
    else /*D==3*/
        return mal::Pi<Real>()*mal::Sq(lss.m_Radius)*(Real(4.0f/3.0f)*lss.m_Radius + height); //sphere+cylinder
}

template <unsigned D>
finline Real ComputeVolume_Fast( const GLineSweptSphere<D>& lss ) { return ComputeVolume(lss); }

template <unsigned D>
finline bool TestOverlap( const GLineSweptSphere<D>& lss1, const GLineSweptSphere<D>& lss2 )
{
    typedef typename GLineSweptSphere<D>::vec_type vec_type;
    /* If none degenerate, compute closest points between segments,
       if one is degenerate, compute closes point (degenerate) to segment (non-degenerate)
       if both degenerate, test sphere vs sphere.
    */
    bool bIsDegenerate1( lss1.IsDegenerate() );
    bool bIsDegenerate2( lss2.IsDegenerate() );
    if( !bIsDegenerate1 && !bIsDegenerate2 )
    {
        //GEO_LOG_WARNING("LSS");
        Real lambda_a,lambda_b;
        np::GClosestPoints_Segment_Segment<D>( lss1.m_Pos0, lss1.m_Pos1, lss2.m_Pos0, lss2.m_Pos1, lambda_a, lambda_b );
        vec_type pos_a( lss1.m_Pos0 + lambda_a*(lss1.m_Pos1-lss1.m_Pos0) );
        vec_type pos_b( lss2.m_Pos0 + lambda_b*(lss2.m_Pos1-lss2.m_Pos0) );
        return (pos_a-pos_b).NormSq() < mal::Sq( lss1.m_Radius + lss2.m_Radius );
    }
    else if( !bIsDegenerate1 )
    {
        //GEO_LOG_WARNING("LSS:D2");
        vec_type closest_point1( np::GClosestPoint_Point_Segment<D>( lss2.m_Pos0, lss1.m_Pos0, lss1.m_Pos1 ) );
        return ( mal::NormSq( lss2.m_Pos0 - closest_point1 ) <= mal::Sq( lss1.m_Radius + lss2.m_Radius ) );
    }
    else if( !bIsDegenerate2 )
    {
        //GEO_LOG_WARNING("LSS:D1");
        vec_type closest_point2( np::GClosestPoint_Point_Segment<D>( lss1.m_Pos0, lss2.m_Pos0, lss2.m_Pos1 ) );
        return ( mal::NormSq( lss1.m_Pos0 - closest_point2 ) <= mal::Sq( lss1.m_Radius + lss2.m_Radius ) );
    }
    else //( bIsDegenerate1 && bIsDegenerate1 )
    {
        //GEO_LOG_WARNING("LSS:D22")
        return ( lss1.m_Pos0-lss2.m_Pos0 ).NormSq() < mal::Sq( lss1.m_Radius + lss2.m_Radius );
    }
}


}} // namespace geo::bv

#endif // GEO_BV_GLINESWEPTSPHERE_H
