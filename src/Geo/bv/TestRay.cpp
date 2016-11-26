#include "TestRay.h"
#include "bv.h"
#include <Geo/np/RayCast.h>

namespace geo {
namespace bv {

bool TestRay( const BoundingVolume2* p_bv, const np::Ray2& ray )
{
    switch( p_bv->GetType() )
    {
    case eBV_Void2: return false; break;
    case eBV_Infinite2: return true; break;
    case eBV_Sphere2:
        {
            const Sphere2& sphere2( p_bv->As<Sphere2>() );
            Vec2 closest_point( np::GClosestPoint_Point_Segment<2>( sphere2.GetPos(),
                                                                    ray.m_Pos + ray.m_Interval.Min()*ray.m_Dir,
                                                                    ray.m_Pos + ray.m_Interval.Max()*ray.m_Dir ) );
            return mal::NormSq( sphere2.GetPos() - closest_point ) <= mal::Sq( sphere2.GetRadius() );
        }
        break;
    case eBV_AABB2:
        {
            np::RayHit2 tmp_rh;
            return np::GRayCast_CenteredAABB<2>( ray.m_Pos - p_bv->As<AABB2>().GetPos(), ray.m_Dir, ray.m_Interval,
                                                 p_bv->As<AABB2>().GetHalfSizes(),
                                                 tmp_rh );
        }
        break;
    case eBV_LSS2:
        {
            /* We convert ray to lss and reuse LSSvsLSS overlap test
               (supports 0-length degeneracy) \todo This is not
               optimal and MAY have precision problems for long
               rays... */
            LSS2 ray_lss( ray.m_Pos + ray.m_Interval.Min()*ray.m_Dir, ray.m_Pos + ray.m_Interval.Max()*ray.m_Dir, Real(0) );
            return TestOverlap( ray_lss, p_bv->As<LSS2>() );

            /* TEMP: This didn't work, I think
               ///Convert ray into a Segment and perform raycast with interval [0,1] relative to the LSS2
            GEO_LOG_WARNING( "geo::bv::TestRay2( LSS2 ) UNTESTED, may not work!!" );
            const LSS2& lss(  );
            Vec2 ray_pos_rel( ray.m_Pos + ray.m_Dir*ray.m_Interval.Min() - lss.GetPos0() );
            Vec2 ray_dir_rel( ray.m_Dir*ray.m_Interval.Max() - (lss.GetPos1()-lss.GetPos0()) );
            return np::GRayCast_Sphere<2>( ray_pos_rel, ray_dir_rel, Interval(0,1),
                                           Vec2(0,0), lss.GetRadius(),
                                           rh );
            //\note rh.m_Interval is not coherent with ray parametrization here...
            */
        }
        break;
    case eBV_DOP2_K4:
        {
            const bv::DOP2_K4& dop2k4( p_bv->As<bv::DOP2_K4>() );
            Vec2 pos( dop2k4.GetInterval(0).Mid(),
                      dop2k4.GetInterval(1).Mid() );
            Vec2 half_sizes( Real(0.5)*dop2k4.GetInterval(0).Length(),
                             Real(0.5)*dop2k4.GetInterval(1).Length() );
            np::RayHit2 tmp_rh;
            return np::GRayCast_CenteredAABB<2>( ray.m_Pos - pos, ray.m_Dir, ray.m_Interval,
                                                 half_sizes,
                                                 tmp_rh );
        }
        break;
    case eBV_DOP2_K8:
        {
            /*\todo !!
            const bv::DOP2_K4& dop2k4( p_bv->As<bv::DOP2_K4>() );
            Vec2 pos( dop2k4.GetInterval(0).Mid(),
                      dop2k4.GetInterval(1).Mid() );
            Vec2 half_sizes( Real(0.5)*dop2k4.GetInterval(0).Length(),
                             Real(0.5)*dop2k4.GetInterval(1).Length() );
            return np::GRayCast_CenteredAABB<2>( ray.m_Pos - pos, ray.m_Dir, ray.m_Interval,
                                                 half_sizes,
                                                 rh );
            */
            // Clip ray against all slabs, early out if empty cliped interval (see centered aabb sub-case)
            return false;
        }
        break;
    default:
        GEO_LOG_WARNING( "geo::bv::TestRay: BV type %d not yet implemented", p_bv->GetType() );
        return true; //BV default response is TRUE, so that an unimplemented BV does not discard a valid hit
        break;
    }
}


bool TestRay( const BoundingVolume3* p_bv, const np::Ray3& ray )
{
    switch( p_bv->GetType() )
    {
    case eBV_Void3: return false; break;
    case eBV_Infinite3: return true; break;
    case eBV_Sphere3:
        {
            const Sphere3& sphere3( p_bv->As<Sphere3>() );
            Vec3 closest_point( np::GClosestPoint_Point_Segment<3>( sphere3.GetPos(),
                                                                    ray.m_Pos + ray.m_Interval.Min()*ray.m_Dir,
                                                                    ray.m_Pos + ray.m_Interval.Max()*ray.m_Dir ) );
            return mal::NormSq( sphere3.GetPos() - closest_point ) <= mal::Sq( sphere3.GetRadius() );
        }
        break;
    case eBV_AABB3:
        {
            np::RayHit3 tmp_rh;
            return np::GRayCast_CenteredAABB<3>( ray.m_Pos - p_bv->As<AABB3>().GetPos(), ray.m_Dir, ray.m_Interval,
                                                 p_bv->As<AABB3>().GetHalfSizes(),
                                                 tmp_rh );
        }
        break;
    case eBV_LSS3:
        {
            /* We convert ray to lss and reuse LSSvsLSS overlap test
               (supports 0-length degeneracy) \todo This is not
               optimal and MAY have precision problems for long
               rays... */
            LSS3 ray_lss( ray.m_Pos + ray.m_Interval.Min()*ray.m_Dir, ray.m_Pos + ray.m_Interval.Max()*ray.m_Dir, Real(0) );
            return TestOverlap( ray_lss, p_bv->As<LSS3>() );

            /* TEMP: This didn't work, I think
               ///Convert ray into a Segment and perform raycast with interval [0,1] relative to the LSS3
            GEO_LOG_WARNING( "geo::bv::TestRay3( LSS3 ) UNTESTED, may not work!!" );
            const LSS3& lss(  );
            Vec3 ray_pos_rel( ray.m_Pos + ray.m_Dir*ray.m_Interval.Min() - lss.GetPos0() );
            Vec3 ray_dir_rel( ray.m_Dir*ray.m_Interval.Max() - (lss.GetPos1()-lss.GetPos0()) );
            return np::GRayCast_Sphere<3>( ray_pos_rel, ray_dir_rel, Interval(0,1),
                                           Vec3(0,0), lss.GetRadius(),
                                           rh );
            //\note rh.m_Interval is not coherent with ray parametrization here...
            */
        }
        break;
    default:
        GEO_LOG_WARNING( "geo::bv::TestRay: BV type %d not yet implemented", p_bv->GetType() );
        return true; //BV default response is TRUE, so that an unimplemented BV does not discard a valid hit
        break;
    }
}

}} //namespace geo::bv
