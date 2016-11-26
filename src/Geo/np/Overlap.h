#ifndef GEO_NP_OVERLAP_H
#define GEO_NP_OVERLAP_H

#include <Geo/Config.h>
#include <Geo/np/PairwiseCache.h> //\todo While OverlapData not available
#include "Context.h"
#include <Geo/bv/GBDOP.h>

/*! This file contains the prototypes for all default Overlap_X_Y functionality.
  Each method may internally use any available implementation transparently (_Analytic, _SAT, _GJK, _Stochastic,...)
*/

//Fwd declarations
namespace geo { class PolygonalShape2; class TriSurfaceShape3; }

namespace geo {
namespace np {

//----------------------------------------------------------------
// Point Vs X
//----------------------------------------------------------------
bool Overlap_Point2_Triangle2( const Vec2& point,
                               const Vec2& tri0, const Vec2& tri1, const Vec2& tri2,
                               OverlapCache2* p_oc = 0,
                               const Context* p_context = g_pDefaultContext );
bool Overlap_Point2_Triangle2( const Vec2& point,
                               const Mat3x3& tri_invBs, //triangle inverse barycentric transform
                               OverlapCache2* p_oc = 0,
                               const Context* p_context = g_pDefaultContext );

bool Overlap_Point2_Polygonal2( const Vec2& point,
                                const PolygonalShape2* p_polygonal, const Transform2& polygonal_tr, const Real* p_polygonal_dof,
                                OverlapCache2* p_oc = 0,
                                const Context* p_context = g_pDefaultContext );

bool Overlap_Point3_TriSurface3( const Vec3& point,
                                 const TriSurfaceShape3* p_surface, const Transform3& surface_tr, const Real* p_surface_dof,
                                 OverlapCache3* p_oc = 0,
                                 const Context* p_context = g_pDefaultContext );

bool Overlap_Point3_Tetrahedron3( const Vec3& point,
                                  const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3,
                                  OverlapCache3* p_oc = 0,
                                  const Context* p_context = g_pDefaultContext );

//----------------------------------------------------------------
// Segment Vs X
//----------------------------------------------------------------
bool Overlap_Segment2_Triangle2( const Vec2& s0, const Vec2& s1,
                                 const Vec2& tri0, const Vec2& tri1, const Vec2& tri2,
                                 OverlapCache2* p_oc = 0,
                                 const Context* p_context = g_pDefaultContext );
bool Overlap_Segment2_Triangle2( const Vec2& s0, const Vec2& s1,
                                 const Mat3x3& tri_invBs, //triangle inverse barycentric transform
                                 OverlapCache2* p_oc = 0,
                                 const Context* p_context = g_pDefaultContext );
bool Overlap_Segment2_Triangle2_BDOP( const Vec2& s0, const Vec2& s1,
                                      const Mat3x3& tri_invBs, //triangle inverse barycentric transform
                                      const bv::BDOP3& bdop,
                                      OverlapCache2* p_oc = 0,
                                      const Context* p_context = g_pDefaultContext );
bool Overlap_Segment3_Tetrahedron3_BDOP( const Vec3& s0, const Vec3& s1,
                                         const Mat4x4& tet_invBs, //tetrahedron inverse barycentric transform
                                         const bv::BDOP4& bdop,
                                         OverlapCache3* p_oc = 0,
                                         const Context* p_context = g_pDefaultContext );

//----------------------------------------------------------------
// Triangle Vs X
//----------------------------------------------------------------
bool Overlap_Triangle3_Triangle3( const Vec3& p0, const Vec3& p1, const Vec3& p2,
                                  const Vec3& q0, const Vec3& q1, const Vec3& q2,
                                  OverlapCache3* p_oc = 0,
                                  const Context* p_context = g_pDefaultContext );

bool Overlap_Triangle3_Tetrahedron3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                     const Vec3& tet0, const Vec3& tet1, const Vec3& tet2, const Vec3& tet3,
                                     OverlapCache3* p_oc = 0,
                                     const Context* p_context = g_pDefaultContext );

bool Overlap_Triangle3_TriSurface3( const Vec3& tri0, const Vec3& tri1, const Vec3& tri2,
                                    const TriSurfaceShape3* p_surface, const Transform3& surface_tr, const Real* p_surface_dof,
                                    OverlapCache3* p_oc = 0,
                                    const Context* p_context = g_pDefaultContext );

/*
bool Overlap_Point2_Box2( const Vec2& point,
                          const Vec2& box_half_sizes, const Transform2& box_transform,
                          OverlapCache2* p_oc,
                          const Context *p_context = g_pDefaultContext )
{
}
*/

/* \todo NO, all ray-related stuff is in np/RayCast.h, this file should only contain shape-pair overlap tests

   Compute overlap interval between ray and sphere
template< unsigned D >
inline bool GTestOverlap_Ray_Sphere( const mal::GVec<Real,D> &ray_pos, const mal::GVec<Real,D> &ray_dir, Real lambda_max,
                                     const mal::GVec<Real,D> &sphere_pos, Real sphere_radius,
                                     Real &lambda0, Real &lambda1 )
{
    const mal::GVec<Real,D> ray_pos_rel( ray_pos - sphere_pos );
    Real dot_pos_rel_dir( mal::Dot(ray_pos_rel,ray_dir) );
    Real dir_norm_sq( ray_dir.NormSq() );
    Real sqrt_arg( mal::Sq(dot_pos_rel_dir) - dir_norm_sq * ( ray_pos_rel.NormSq() - mal::Sq(sphere_radius) ) );
    if( sqrt_arg >= Real(0) )
    {
        Real root( mal::Sqrt(sqrt_arg) );
        Real rcp_dir_norm_sq( mal::Rcp(dir_norm_sq) );
        lambda0 = (-dot_pos_rel_dir-root)*rcp_dir_norm_sq;
        lambda1 = (-dot_pos_rel_dir+root)*rcp_dir_norm_sq;
        if( lambda0 <= lambda_max && lambda1 >= Real(0) ) //Implicitly lambda0 < lambda1
        {
            lambda0 = mal::Max(lambda0,Real(0));
            lambda1 = mal::Min(lambda1,lambda_max);
            return true;
        }
    }
    return false;
}
*/

}} //namespace geo::np

#endif // GEO_NP_OVERLAP_H
