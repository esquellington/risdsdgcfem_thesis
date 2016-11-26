#ifndef GEO_NP_RAYCAST_H
#define GEO_NP_RAYCAST_H

#include <Geo/Config.h>
#include "Context.h"
#include <Geo/np/stats.h>

//Fwd declarations
namespace geo { class TriSurfaceShape3; class DCR_MeshSolidShape2; class DCR_TetSolidShape3; } //\todo MAYBE this RC tests should be in a diferent namespace, they're not truly "NP" as they use higher level structures TSS,MSS,etc...

//\todo SPLIT this file into Ray.h, RayCast.h and RayCast_X_Y_Analytic, BVH, etc as done with Overlap and Contact

namespace geo {
namespace np {

//---- Ray Types
template< unsigned D >
struct GRay
{
    finline GRay( const mal::GVec<Real,D> &pos, const mal::GVec<Real,D> &dir, const Interval &interval, Real thickness )
    : m_Pos(pos), m_Dir(dir), m_Interval(interval), m_Thickness(thickness) {}
    mal::GVec<Real,D> m_Pos;
    mal::GVec<Real,D> m_Dir;
    Interval m_Interval;
    Real m_Thickness;
};
typedef GRay<2> Ray2;
typedef GRay<3> Ray3;

template< unsigned D >
struct GRayHit
{
    Interval m_Interval;
    mal::GVec<Real,D> m_Point;
    mal::GVec<Real,D> m_Normal;
    feature_id m_FeatureId;
    // Extra...
    Vec4 m_Extra_BarycentricCoords; //Hit Vertex/Segment/Triangle/Tetrahedron barycentric coords

    //\todo m_Extra could contain hit-object local coords, lambda, barycentric coords of the hit-point, add as required \todo Using a Union is problematic for serialization...
    //union { Real m_Lambda; } m_Extra;
};
typedef GRayHit<2> RayHit2;
typedef GRayHit<3> RayHit3;

template< unsigned D >
struct GRayCache
{
    /* \todo */
};
typedef GRayCache<2> RayCache2;
typedef GRayCache<3> RayCache3;

//! Compute overlap interval between ray and plane
template< unsigned D >
inline bool GRayCast_Plane( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                            const mal::GVec<Real,D> &plane_normal, Real plane_coeff_d,
                            GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                            const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_GPlane );
    Real n_dot_dir( mal::Dot(plane_normal,ray_dir) );
    Real n_dot_p_plus_d( mal::Dot(plane_normal,ray_pos) + plane_coeff_d );
    if( mal::Abs(n_dot_dir) > g_pDefaultContext->m_Epsilon_Dir ) // non-coplanar ray
    {
        Real lambda( -n_dot_p_plus_d / n_dot_dir );
        if( lambda >= interval.Min() && lambda <= interval.Max() )
        {
            rh.m_Interval = Interval(lambda,lambda);
            rh.m_Point = ray_pos + ray_dir * rh.m_Interval.Min();
            rh.m_Normal = plane_normal;
            rh.m_FeatureId = feature_id();
            return true;
        }
        else
            return false;
    }
    else if( mal::Abs(n_dot_p_plus_d) < g_pDefaultContext->m_Epsilon_Length ) // coplanar included ray, \lambda = any
    {
        rh.m_Interval = interval;
        rh.m_Point = ray_pos + ray_dir*rh.m_Interval.Min();
        rh.m_Normal = plane_normal;
        rh.m_FeatureId = feature_id();
        return true;
    }
    else
        return false; // coplanar non-included ray
}

//! Compute overlap interval between ray and halfspace
template< unsigned D >
inline bool GRayCast_HalfSpace( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                                const mal::GVec<Real,D> &plane_normal, Real plane_coeff_d,
                                GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                                const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_GHalfSpace );
    Real n_dot_p_plus_d( mal::Dot(plane_normal,ray_pos) + plane_coeff_d );
    if( n_dot_p_plus_d > 0 ) //above halfspace
    {
        if( GRayCast_Plane<D>(ray_pos,ray_dir,interval,plane_normal,plane_coeff_d,rh,p_rc) )
        {
            //\note rh.m_Interval.Min() == lambda;
            rh.m_Interval.Max() = interval.Max();
            return true;
        }
        else
            return false;
    }
    else // below halfspace
    {
        if( GRayCast_Plane<D>(ray_pos,ray_dir,interval,-plane_normal,-plane_coeff_d,rh,p_rc) )
        {
            rh.m_Interval.Min() = interval.Min();
            //\note rh.m_Interval.Max() == lambda;
        }
        else
            rh.m_Interval = interval;
        rh.m_Point = ray_pos;
        rh.m_Normal = plane_normal;
        return true;
    }
}

//! Compute overlap interval between ray and sphere
template< unsigned D >
inline bool GRayCast_Sphere( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                             const mal::GVec<Real,D> &sphere_pos, Real sphere_radius,
                             GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                             const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_GSphere );
    const mal::GVec<Real,D> ray_pos_rel( ray_pos - sphere_pos );
    Real dot_pos_rel_dir( mal::Dot(ray_pos_rel,ray_dir) );
    Real dir_norm_sq( ray_dir.NormSq() );
    Real sqrt_arg( mal::Sq(dot_pos_rel_dir) - dir_norm_sq * ( ray_pos_rel.NormSq() - mal::Sq(sphere_radius) ) );
    if( sqrt_arg >= Real(0) )
    {
        Real root( mal::Sqrt(sqrt_arg) );
        Real rcp_dir_norm_sq( mal::Rcp(dir_norm_sq) );
        Real lambda0 = (-dot_pos_rel_dir-root)*rcp_dir_norm_sq;
        Real lambda1 = (-dot_pos_rel_dir+root)*rcp_dir_norm_sq;
        if( lambda0 <= interval.Max() && lambda1 >= interval.Min() ) //implicitly lambda0 < lambda1
        {
            rh.m_Interval.Min() = mal::Max(lambda0,Real(interval.Min()));
            rh.m_Interval.Max() = mal::Min(lambda1,interval.Max());
            rh.m_Point = ray_pos + ray_dir * rh.m_Interval.Min();
            rh.m_Normal = mal::SafeNormalized( rh.m_Point - sphere_pos, mal::Epsilon<Real>() ); //Avoids potential div-by-0
            rh.m_FeatureId = feature_id();
            return true;
        }
    }
    return false;
}

/*! Inspired in RTCD 5.3, pg 181.
  Returns true if the ray overlaps the centered AABB along its
  parameter interval, and computes normal from first hit axis.
*/
template< unsigned D >
inline bool GRayCast_CenteredAABB( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                                   const mal::GVec<Real,D> &aabb_half_sizes,
                                   GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                                   const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_GCenteredAABB );
    rh.m_Interval = interval;
    int first_hit_axis(0);
    for( int it_axis=0; it_axis<int(D); it_axis++ )
    {
        // If parallel to slab, either overlaps for any lambda or for none.
        if( mal::Abs(ray_dir[it_axis]) < g_pDefaultContext->m_Epsilon_Dir )
        {
            if( mal::Abs(ray_pos[it_axis]) > aabb_half_sizes[it_axis] )
            {
                //Empty interval and return false
                rh.m_Interval = Interval::Empty();
                return false;
            }
            // Otherwise, current axis does NOT clip the interval, and
            // other axis must be checked as usual.
        }
        else
        {
            Real inv_divisor( mal::Rcp(ray_dir[it_axis]) );
            Real lambda0( ( -aabb_half_sizes[it_axis] - ray_pos[it_axis] ) * inv_divisor );
            Real lambda1( (  aabb_half_sizes[it_axis] - ray_pos[it_axis] ) * inv_divisor );
            if( lambda1 < lambda0 )
            {
                Real tmp( lambda0 );
                lambda0 = lambda1;
                lambda1 = tmp;
            }
            // Clip lambda-interval and update first-axis
            if( rh.m_Interval.Min() < lambda0 ) { rh.m_Interval.Min() = lambda0; first_hit_axis = it_axis; }
            if( rh.m_Interval.Max() > lambda1 ) rh.m_Interval.Max() = lambda1;
            if( rh.m_Interval.IsEmpty() )
                return false;
        }
    }
    // Compute Point and Normal
    rh.m_Point  = ray_pos + rh.m_Interval.Min() * ray_dir;
    rh.m_Normal = mal::GVec<Real,D>::Zero();
    rh.m_Normal[ first_hit_axis ] = ( ray_pos[first_hit_axis] < Real(0) ) ? Real(-1) : Real(1);
    rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E/F sub-features
    return true;
}

//! Compute overlap interval between ray and oriented box
template< unsigned D >
inline bool GRayCast_Box( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                          const mal::GTransform<Real,D> &box_transform, const mal::GVec<Real,D> &box_half_sizes,
                          GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                          const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_GBox );
    mal::GTransform<Real,D> inv_box_transform( box_transform.Inverse() ); //\todo Could be avoided if GTransform.InvMul(p) and InvRot(v) were available
    mal::GVec<Real,D> ray_pos_rel( inv_box_transform * ray_pos );
    mal::GVec<Real,D> ray_dir_rel( inv_box_transform.m_Rot * ray_dir );
    if( GRayCast_CenteredAABB<D>( ray_pos_rel, ray_dir_rel, interval, box_half_sizes, rh, p_rc ) )
    {
        rh.m_Point  = ray_pos + rh.m_Interval.Min() * ray_dir; //== box_transform * rh.m_Point;
        rh.m_Normal = box_transform.m_Rot * rh.m_Normal; //\todo It's just a local axis in AABB, could be optimized
        rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E/F sub-features
        return true;
    }
    else
        return false;
}

template< unsigned D >
inline bool GRayCast_Cylinder( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                               const mal::GVec<Real,D> &cylinder_pos, const mal::GVec<Real,D> &cylinder_axis,
                               Real cylinder_radius, Real cylinder_half_height,
                               GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                               const Context* p_context = g_pDefaultContext )
{
    GEO_LOG_WARNING( "geo::np::GRayCast_Cylinder<D> generic implementation is empty" );
    return false;
}

template<>
inline bool GRayCast_Cylinder<2>( const mal::GVec<Real,2> &ray_pos, const mal::GVec<Real,2> &ray_dir, const Interval &interval,
                                  const mal::GVec<Real,2> &cylinder_pos, const mal::GVec<Real,2> &cylinder_axis,
                                  Real cylinder_radius, Real cylinder_half_height,
                                  GRayHit<2> &rh, GRayCache<2> *p_rc,
                                  const Context* p_context )
{
    GEO_NP_STAT_INC( raycast.m_Cylinder2 );
    Transform2f box_tr;
    box_tr.m_Pos = cylinder_pos;
    mal::GVec<float,2> perpendicular_axis( mal::PerpendicularCW(cylinder_axis) );
    box_tr.m_Rot(0,0) = perpendicular_axis[0]; box_tr.m_Rot(0,1) = cylinder_axis[0];
    box_tr.m_Rot(1,0) = perpendicular_axis[1]; box_tr.m_Rot(1,1) = cylinder_axis[1];
    return GRayCast_Box<2>( ray_pos, ray_dir, interval, box_tr, Vec2( cylinder_radius, cylinder_half_height ), rh, p_rc );
}

/* \todo This algorithm would be generic, but suboptimal for 2d, where a Cylinder is just a Box
template< unsigned D >
inline bool GRayCast_Cylinder( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                               const mal::GVec<Real,D> &cylinder_pos, const mal::GVec<Real,D> &cylinder_axis, Real cylinder_radius, Real cylinder_half_height,
                               GRayHit<D> &rh, GRayCache<D> *p_rc = 0 )
{
    // Init with an Empty reversed interval (max,min) so that (min,max) can be computed incrementally comparing with extreme values (max,min)
    rh.m_Interval = Interval( interval.Max(), interval.Min() );
    np::RayHit2 tmp_rh;
    // Infinite Cylinder
    if( GRayCast_InfinityCylinder( ray_pos, ray_dir, cylinder_pos, cylinder_axis, cylinder_radius, tmp_rh ) )
    {
        // Update first hit if closest
        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
        {
            rh.m_Interval.Min() = tmp_rh.m_Interval.Min();
            rh.m_Point = tmp_rh.m_Point;
            rh.m_Normal = tmp_rh.m_Normal;
            rh.m_FeatureId = tmp_rh.m_FeatureId;
            rh.m_Extra = tmp_rh.m_Extra;
        }
        // Update interval max if farthest
        if( tmp_rh.m_Interval.Max() > rh.m_Interval.Max() ) rh.m_Interval.Max() = tmp_rh.m_Interval.Max();
    }
    // End caps = slabs
    if( GRayCast_Slabs( ray_pos, ray_dir, cylinder_pos, cylinder_axis, cylinder_half_height, tmp_rh ) )
    {
        // Update first hit if closest
        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
        {
            rh.m_Interval.Min() = tmp_rh.m_Interval.Min();
            rh.m_Point = tmp_rh.m_Point;
            rh.m_Normal = tmp_rh.m_Normal;
            rh.m_FeatureId = tmp_rh.m_FeatureId;
            rh.m_Extra = tmp_rh.m_Extra;
        }
        // Update interval max if farthest
        if( tmp_rh.m_Interval.Max() > rh.m_Interval.Max() ) rh.m_Interval.Max() = tmp_rh.m_Interval.Max();
    }
    // If any sub-test has hit, the interval won't be empty
    return !rh.m_Interval.IsEmpty();
}
*/

//! Compute overlap interval between ray and capsule \todo Current implementation is "brute force", could be optimized!!
template< unsigned D >
inline bool GRayCast_Capsule( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, const Interval &interval,
                              const mal::GVec<Real,D> &capsule_pos, const mal::GVec<Real,D> &capsule_axis, Real capsule_radius, Real capsule_half_height,
                              GRayHit<D> &rh, GRayCache<D> *p_rc = 0,
                              const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_GCapsule );
    // Init with an Empty reversed interval (max,min) so that (min,max) can be computed incrementally comparing with extreme values (max,min)
    rh.m_Interval = Interval( interval.Max(), interval.Min() );
    np::RayHit2 tmp_rh;
    if( GRayCast_Cylinder<D>( ray_pos, ray_dir, interval, capsule_pos, capsule_axis, capsule_radius, capsule_half_height, tmp_rh ) )
    {
        // Update first hit if closest
        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
        {
            rh.m_Interval.Min() = tmp_rh.m_Interval.Min();
            rh.m_Point = tmp_rh.m_Point;
            rh.m_Normal = tmp_rh.m_Normal;
            rh.m_FeatureId = tmp_rh.m_FeatureId;
            rh.m_Extra_BarycentricCoords = tmp_rh.m_Extra_BarycentricCoords;
        }
        // Update interval max if farthest
        if( tmp_rh.m_Interval.Max() > rh.m_Interval.Max() ) rh.m_Interval.Max() = tmp_rh.m_Interval.Max();
    }
    if( GRayCast_Sphere<D>( ray_pos, ray_dir, interval, capsule_pos - capsule_half_height*capsule_axis, capsule_radius, tmp_rh ) )
    {
        // Update first hit if closest
        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
        {
            rh.m_Interval.Min() = tmp_rh.m_Interval.Min();
            rh.m_Point = tmp_rh.m_Point;
            rh.m_Normal = tmp_rh.m_Normal;
            rh.m_FeatureId = tmp_rh.m_FeatureId;
            rh.m_Extra_BarycentricCoords = tmp_rh.m_Extra_BarycentricCoords;
        }
        // Update interval max if farthest
        if( tmp_rh.m_Interval.Max() > rh.m_Interval.Max() ) rh.m_Interval.Max() = tmp_rh.m_Interval.Max();
    }
    if( GRayCast_Sphere<D>( ray_pos, ray_dir, interval, capsule_pos + capsule_half_height*capsule_axis, capsule_radius, tmp_rh ) )
    {
        // Update first hit if closest
        if( tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
        {
            rh.m_Interval.Min() = tmp_rh.m_Interval.Min();
            rh.m_Point = tmp_rh.m_Point;
            rh.m_Normal = tmp_rh.m_Normal;
            rh.m_FeatureId = tmp_rh.m_FeatureId;
            rh.m_Extra_BarycentricCoords = tmp_rh.m_Extra_BarycentricCoords;
        }
        // Update interval max if farthest
        if( tmp_rh.m_Interval.Max() > rh.m_Interval.Max() ) rh.m_Interval.Max() = tmp_rh.m_Interval.Max();
    }
    // If any sub-test has hit, the interval won't be empty
    if( rh.m_Interval.Max() < rh.m_Interval.Min() ) rh.m_Interval = Interval(); //Force empty if max<min
    return !rh.m_Interval.IsEmpty();
}

//! Ray vs Segment 2D. Computed Normal is oriented "outwards" assuming Segment==Edge on a CCW polygon
inline bool RayCast_Segment2_SingleSided( const Vec2 &ray_pos, const Vec2 &ray_dir, const Interval &interval,
                                          const Vec2 &segment_p0, const Vec2 &segment_p1,
                                          RayHit2 &rh, RayCache2 *p_rc = 0,
                                          const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Segment2_SS );
    Vec2 segment_dir( segment_p1 - segment_p0 );
    Real det( ray_dir[1]*segment_dir[0] - ray_dir[0]*segment_dir[1] ); //== Dot(ray_dir,PerpendicularCW(segment_dir)) > 0
    if( det > Real(1e-5) ) //\todo USE appropiate g_pDefaultContext->m_Epsilon_Rcp //epsilon for 1/x computation
    {
        Vec2 diff_ab = segment_p0 - ray_pos;
        Real inv_det = mal::Rcp(det);
        Real ray_lambda = inv_det * (diff_ab[1]*segment_dir[0] - diff_ab[0]*segment_dir[1]);
        Real segment_lambda = inv_det * (diff_ab[1]*ray_dir[0] - diff_ab[0]*ray_dir[1]);
        if( ray_lambda >= interval.Min() && ray_lambda <= interval.Max()
            && segment_lambda >= Real(0) && segment_lambda <= Real(1) )
        {
            rh.m_Interval = Interval(ray_lambda,ray_lambda);
            rh.m_Point = ray_pos + ray_lambda * ray_dir;
            rh.m_Normal = mal::PerpendicularCCW( segment_dir.Normalized() );
            rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E sub-features
            rh.m_Extra_BarycentricCoords = Vec4( Real(1) - segment_lambda, segment_lambda, 0, 0 );
            //GEO_LOG_WARNING( "Hit at lambda %f", segment_lambda );
            return true;
        }
        else
            return false;
    }
    else // \todo Parallel, check ray_pos and lambda range
    {
        return false;
    }
}

//! Ray vs Segment 2D. Computed Normal is oriented "outwards" assuming Segment==Edge on a CCW polygon... MAYBE should be actually double-sided normal??
inline bool RayCast_Segment2_DoubleSided( const Vec2 &ray_pos, const Vec2 &ray_dir, const Interval &interval,
                                          const Vec2 &segment_p0, const Vec2 &segment_p1,
                                          RayHit2 &rh, RayCache2 *p_rc = 0,
                                          const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Segment2_DS );
    Vec2 segment_dir( segment_p1 - segment_p0 );
    Real det( ray_dir[1]*segment_dir[0] - ray_dir[0]*segment_dir[1] );
    if( mal::Abs(det) > Real(1e-5) ) //\todo USE appropiate g_pDefaultContext->m_Epsilon_Rcp //epsilon for 1/x computation
    {
        Vec2 diff_ab = segment_p0 - ray_pos;
        Real inv_det = mal::Rcp(det);
        Real ray_lambda = inv_det * (diff_ab[1]*segment_dir[0] - diff_ab[0]*segment_dir[1]);
        Real segment_lambda = inv_det * (diff_ab[1]*ray_dir[0] - diff_ab[0]*ray_dir[1]);
        if( ray_lambda >= interval.Min() && ray_lambda <= interval.Max()
            && segment_lambda >= Real(0) && segment_lambda <= Real(1) )
        {
            rh.m_Interval = Interval(ray_lambda,ray_lambda);
            rh.m_Point = ray_pos + ray_lambda * ray_dir;
            rh.m_Normal = mal::PerpendicularCCW( segment_dir.Normalized() );
            rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E sub-features
            rh.m_Extra_BarycentricCoords = Vec4( Real(1) - segment_lambda, segment_lambda, 0, 0 );
            //GEO_LOG_WARNING( "Hit at lambda %f", segment_lambda );
            return true;
        }
        else
            return false;
    }
    else // \todo Parallel, check ray_pos and lambda range
    {
        return false;
    }
}

//\note Single-sided test
inline bool RayCast_Triangle3_SingleSided( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                           const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2,
                                           RayHit3& rh, RayCache3* p_rc = 0,
                                           const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Triangle3_SS );
    //\note Code from RTCD 5.3.6, pg 191
    Vec3 segment_p( ray_pos + ray_dir * interval.Min() );
    Vec3 segment_q( ray_pos + ray_dir * interval.Max() );
    Vec3 ab = triangle_p1 - triangle_p0;
    Vec3 ac = triangle_p2 - triangle_p0;
    Vec3 qp = segment_p - segment_q;
    Vec3 n = mal::Cross( ab, ac );
    Real d = mal::Dot( qp, n );
    if( d <= Real(0) ) return false;
    Vec3 ap = segment_p - triangle_p0;
    Real t = mal::Dot( ap, n );
    if( t < Real(0) ) return false;
    if( t > d ) return false;
    Vec3 e = mal::Cross( qp, ap );
    Real v = mal::Dot( ac, e );
    if( v < Real(0) || v > d ) return false;
    Real w = -mal::Dot( ab, e );
    if( w < Real(0) || v + w > d ) return false;
    // Hit found, fill rh and return true
    t /= d;
    rh.m_Interval = Interval( interval.Min() + t * interval.Length() );
    rh.m_Point = ray_pos + rh.m_Interval.Min() * ray_dir;
    rh.m_Normal = mal::Normalized( n );
    rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E/F sub-features
    rh.m_Extra_BarycentricCoords = Vec4::Zero(); //\todo Save barycentric coords according to FT
    return true;
}

inline bool RayCast_Triangle3_DoubleSided( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                           const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2,
                                           RayHit3& rh, RayCache3* p_rc = 0,
                                           const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Triangle3_DS );
    //\note Code from RTCD 5.3.6, pg 191
    Vec3 segment_p( ray_pos + ray_dir * interval.Min() );
    Vec3 segment_q( ray_pos + ray_dir * interval.Max() );
    Vec3 ab = triangle_p1 - triangle_p0;
    Vec3 ac = triangle_p2 - triangle_p0;
    Vec3 qp = segment_p - segment_q;
    Vec3 n = mal::Cross( ab, ac );
    Real d = mal::Dot( qp, n );
    // Coplanar?
    //if( mal::Abs(d) < g_pDefaultContext->m_Epsilon_Dir ) return false;
    if( d <= Real(0) )
    {
        // Invert triangle orientation
        n = -n;
        ab = triangle_p2 - triangle_p0;
        ac = triangle_p1 - triangle_p0;
        d = -d;
    }
    Vec3 ap = segment_p - triangle_p0;
    Real t = mal::Dot( ap, n );
    if( t < Real(0) ) return false;
    if( t > d ) return false;
    Vec3 e = mal::Cross( qp, ap );
    Real v = mal::Dot( ac, e );
    if( v < Real(0) || v > d ) return false;
    Real w = -mal::Dot( ab, e );
    if( w < Real(0) || v + w > d ) return false;
    // Hit found, fill rh and return true
    t /= d;
    rh.m_Interval = Interval( interval.Min() + t * interval.Length() );
    rh.m_Point = ray_pos + rh.m_Interval.Min() * ray_dir;
    rh.m_Normal = mal::Normalized( n );
    //\todo THIS MAY NEED to consider inverted triangles (_DoubleSided)
    rh.m_FeatureId = feature_id(); //\todo Consider reporting V/E/F sub-features
    rh.m_Extra_BarycentricCoords = Vec4::Zero(); //\todo Save barycentric coords according to FT
    return true;
}

/* Specialization of RayCast_Triangle3_DoubleSided for
   - Ray == Segment
   - Triangle normal passed explicitly
   - Only rh lambda and point are computed
 */
inline bool RayCast_Triangle3_DoubleSided_FAST_SEGMENT( const Vec3& segment_p, const Vec3& segment_q,
                                                        const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2, const Vec3& triangle_normal_non_unitary,
                                                        Vec3& point, Real& lambda,
                                                        const Context* p_context = g_pDefaultContext )
{
    GEO_NP_STAT_INC( raycast.m_Triangle3_DS );
    Vec3 ab = triangle_p1 - triangle_p0;
    Vec3 ac = triangle_p2 - triangle_p0;
    Vec3 qp = segment_p - segment_q;
    Vec3 n = triangle_normal_non_unitary;
    Real d = mal::Dot( qp, n );
    // Coplanar?
    //if( mal::Abs(d) < g_pDefaultContext->m_Epsilon_Dir ) return false;
    if( d <= Real(0) )
    {
        // Invert triangle orientation
        n = -n;
        ab = triangle_p2 - triangle_p0;
        ac = triangle_p1 - triangle_p0;
        d = -d;
    }
    Vec3 ap = segment_p - triangle_p0;
    Real t = mal::Dot( ap, n );
    if( t < Real(0) ) return false;
    if( t > d ) return false;
    Vec3 e = mal::Cross( qp, ap );
    Real v = mal::Dot( ac, e );
    if( v < Real(0) || v > d ) return false;
    Real w = -mal::Dot( ab, e );
    if( w < Real(0) || v + w > d ) return false;
    // Hit found, fill data and return true
    lambda = t/d;
    point = segment_p - lambda * qp; //== p + lambda * pq
    return true;
}

/* Coherent double-sided RC guarantees that the SAME rh will be
   generated for any permutation of triangle vertices. This is
   REQUIRED by Clip_Triangle3_Tetrahedron3() when used inside
   Clip_TriSurfaceShape3_TetSolidShape3() in order to produce the SAME
   intersection points for a given Tri against all TetSolid
   tetrahedrons, most of which share triangular faces with their
   neighbours.

   IMPORTANT: This is MUCH SLOWER (1.5x) than non-coherent version,
   use with care.
*/
bool RayCast_Triangle3_DoubleSided_COHERENT( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                             const Vec3& triangle_p0, const Vec3& triangle_p1, const Vec3& triangle_p2,
                                             RayHit3& rh, RayCache3* p_rc = 0,
                                             const Context* p_context = g_pDefaultContext );

/* Highe level raycasts
   \todo MAYBE this RC tests should be in a diferent namespace,
   they're not truly "NP" as they use higher level structures
   TSS,MSS,etc...
*/
//@{
bool RayCast_TriSurfaceShape3( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                               const TriSurfaceShape3* p_surface, const Transform3& surface_tr, const Real* p_surface_dof,
                               bool b_double_sided,
                               RayHit3& rh, RayCache3* p_rc = 0,
                               const Context* p_context = g_pDefaultContext );

//\note This method just forwards to BruteForce or BDT variand according to np::Contecxt::m_RCDCR_Method
bool RayCast_DCR_E_DoubleSided( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                const DCR_MeshSolidShape2* p_dcr, uint32 eid,
                                const Mat3x3& Bs, const Mat3x3& invBs, const Mat3x3& Bs_invBm,
                                RayHit2& rh, RayCache2 *p_rc = 0,
                                const Context* p_context = g_pDefaultContext );

bool RayCast_DCR_E_DoubleSided_BruteForce( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                           const DCR_MeshSolidShape2* p_dcr, uint32 eid,
                                           const Mat3x3& Bs, const Mat3x3& invBs, const Mat3x3& Bs_invBm,
                                           RayHit2& rh, RayCache2 *p_rc = 0,
                                           const Context* p_context = g_pDefaultContext );

bool RayCast_DCR_E_DoubleSided_BDT( const Vec2& ray_pos, const Vec2& ray_dir, const Interval& interval,
                                    const DCR_MeshSolidShape2* p_dcr, uint32 eid,
                                    const Mat3x3& Bs, const Mat3x3& invBs, const Mat3x3& Bs_invBm,
                                    RayHit2& rh, RayCache2 *p_rc = 0,
                                    const Context* p_context = g_pDefaultContext );

//\note This method just forwards to BruteForce or BDT variand according to np::Contecxt::m_RCDCR_Method
bool RayCast_DCR_E_DoubleSided( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                const DCR_TetSolidShape3* p_dcr, uint32 eid,
                                const Mat4x4& Bs, const Mat4x4& invBs, const Mat4x4& Bs_invBm,
                                RayHit3& rh, RayCache3 *p_rc = 0,
                                const Context* p_context = g_pDefaultContext );

bool RayCast_DCR_E_DoubleSided_BruteForce( const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                           const DCR_TetSolidShape3* p_dcr, uint32 eid,
                                           const Mat4x4& Bs, const Mat4x4& invBs, const Mat4x4& Bs_invBm,
                                           RayHit3& rh, RayCache3 *p_rc = 0,
                                           const Context* p_context = g_pDefaultContext );

bool RayCast_DCR_E_DoubleSided_BDT(const Vec3& ray_pos, const Vec3& ray_dir, const Interval& interval,
                                   const DCR_TetSolidShape3* p_dcr, uint32 eid,
                                   const Mat4x4& Bs, const Mat4x4& invBs, const Mat4x4& Bs_invBm,
                                   RayHit3& rh, RayCache3 *p_rc = 0,
                                   const Context* p_context = g_pDefaultContext );

//@}

}} //namespace geo::np

#endif // GEO_NP_RAYCAST_H
