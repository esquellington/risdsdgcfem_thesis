#include "Contact.h"
#include <Geo/shape/shape.h>

// Include implementations actually used
#include "Contact_X_Y_Analytic.h"
#include "Contact_X_Y_BruteForce.h"
#include "Contact_X_Y_Stochastic.h"
#include "Contact_X_Y_BVH.h"
#include "Contact_X_DCR3_BVH.h"

#include "stats.h"

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// Sphere2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_Sphere2_Sphere2( const Vec2& pos1, Real radius1,
                              const Vec2& pos2, Real radius2,
                              ContactData2& cd, ContactCache2* p_cc,
                              const Context* p_context )
{
    return GContact_Sphere_Sphere_Analytic<2>( pos1, radius1, pos2, radius2, cd, p_cc, p_context );
}

bool Contact_Sphere2_Plane2( const Vec2& sphere_pos, Real sphere_radius,
                             const Vec2& plane_normal, Real plane_coeff_d,
                             ContactData2& cd, ContactCache2* p_cc,
                             const Context* p_context )
{
    return GContact_Sphere_Plane_Analytic<2>( sphere_pos, sphere_radius, plane_normal, plane_coeff_d, cd, p_cc, p_context );
}

bool Contact_Sphere2_MeshSolid2( const Vec2& sphere_pos, Real sphere_radius,
                                 const MeshSolidShape2* p_mesh, const Transform2& mesh_tr, const Real* p_mesh_dof,
                                 ContactData2& cd, ContactCache2* p_cc,
                                 const Context* p_context )
{
    return false; //\todo Contact_Sphere2_MeshSolid2_BruteForce( sphere_pos, sphere_radius, p_mesh, mesh_tr, p_mesh_dof, cd, p_cc, p_context );
}

//--------------------------------------------------------------------------------------------------------------------------------
// Sphere3 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_Sphere3_Sphere3( const Vec3& pos1, Real radius1,
                              const Vec3& pos2, Real radius2,
                              ContactData3& cd, ContactCache3* p_cc,
                              const Context* p_context )
{
    return GContact_Sphere_Sphere_Analytic<3>( pos1, radius1, pos2, radius2, cd, p_cc, p_context );
}

bool Contact_Sphere3_Plane3( const Vec3& sphere_pos, Real sphere_radius,
                             const Vec3& plane_normal, Real plane_coeff_d,
                             ContactData3& cd, ContactCache3* p_cc,
                             const Context* p_context )
{
    return GContact_Sphere_Plane_Analytic<3>( sphere_pos, sphere_radius, plane_normal, plane_coeff_d, cd, p_cc, p_context );
}

//--------------------------------------------------------------------------------------------------------------------------------
// MeshSolid2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_MeshSolid2_Plane2( const MeshSolidShape2* p_mesh, const Transform2& mesh_transform, const Real* p_mesh_dof,
                                const Vec2& plane_normal, Real plane_coeff_d,
                                ContactData2& cd, ContactCache2* p_cc,
                                const Context* p_context )
{
#define __USE_DCR2PLANE
#ifdef __USE_DCR2PLANE
    return Contact_DCRMS2_Plane2_BVH_MSS_Lazy_DCR( p_mesh, mesh_transform, p_mesh_dof, plane_normal, plane_coeff_d, cd, p_cc, p_context );
#else
    return Contact_MeshSolid2_Plane2_BruteForce( p_mesh, mesh_transform, p_mesh_dof, plane_normal, plane_coeff_d, cd, p_cc, p_context );
#endif
}

bool Contact_MeshSolid2_MeshSolid2( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* p_mesh_dof1,
                                    const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* p_mesh_dof2,
                                    ContactData2& cd, ContactCache2* p_cc,
                                    const Context* p_context )
{
    GEO_NP_BEGIN_PROFILE_BLOCK_MS( contact.m_Profile_MSS2_MSS2_ms );
    /*\todo While developing, we'll choose among:
      eMSS2MSS_Bruteforce_MSS
      eMSS2MSS_Stochastic_MSS
      eMSS2MSS_Bruteforce_DCR
      eMSS2MSS_BVH_DCR
      eMSS2MSS_Bruteforce_MSS_Lazy_DCR
      eMSS2MSS_BVH_MSS_Lazy_DCR
    */
    bool bResult(false);
    switch( p_context->m_SS2SS_Method )
    {
    case Context::eSS2SS_Bruteforce_SS: bResult = Contact_MeshSolid2_MeshSolid2_BruteForce( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    case Context::eSS2SS_Stochastic_SS: bResult = Contact_IShape2_IShape2_Stochastic( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    case Context::eSS2SS_Bruteforce_DCR: bResult = Contact_DCRMS2_DCRMS2_BruteForce_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
        //case Context::eSS2SS_BVH_DCR: bResult = Contact_DCRMS2_DCRMS2_BVH_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    case Context::eSS2SS_Bruteforce_SS_Lazy_DCR: bResult = Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    case Context::eSS2SS_BVH_SS_Lazy_DCR: bResult = Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    default: bResult = Contact_MeshSolid2_MeshSolid2_BruteForce( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    }
    GEO_NP_END_PROFILE_BLOCK_MS( contact.m_Profile_MSS2_MSS2_ms );

    return bResult;
}

//--------------------------------------------------------------------------------------------------------------------------------
// TetSolid3 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_TetSolid3_Plane3( const TetSolidShape3* p_solid, const Transform3& solid_tr, const Real* p_solid_dof,
                               const Vec3& plane_normal, Real plane_coeff_d,
                               ContactData3& cd, ContactCache3* p_cc,
                               const Context* p_context )
{
    // GEO_LOG_ERROR("Contact_TetSolid3_Plane3() not yet implemented");
    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_TSS3_Plane_ms );
    bool bResult = Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR( p_solid, solid_tr, p_solid_dof, plane_normal, plane_coeff_d, cd, p_cc, p_context );
    //bool bResult = Contact_TetSolid3_Plane3_BruteForce( p_solid, solid_tr, p_solid_dof, plane_normal, plane_coeff_d, cd, p_cc, p_context );
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_TSS3_Plane_ms );
    return bResult;
}

bool Contact_TetSolid3_TetSolid3( const TetSolidShape3* p_solid1, const Transform3& solid_tr1, const Real* p_solid_dof1,
                                  const TetSolidShape3* p_solid2, const Transform3& solid_tr2, const Real* p_solid_dof2,
                                  ContactData3& cd, ContactCache3* p_cc,
                                  const Context* p_context )
{
    GEO_NP_BEGIN_PROFILE_BLOCK_MS( contact.m_Profile_TSS3_TSS3_ms );
    GEO_NP_BEGIN_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_TSS3_TSS3_ms );
    bool bResult(false);
    switch( p_context->m_SS2SS_Method )
    {
    // case Context::eSS2SS_Bruteforce_MSS: bResult = Contact_MeshSolid2_MeshSolid2_BruteForce( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    // case Context::eSS2SS_Stochastic_MSS: bResult = Contact_IShape2_IShape2_Stochastic( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    // case Context::eSS2SS_Bruteforce_DCR: bResult = Contact_DCRMS2_DCRMS2_BruteForce_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    //     //case Context::eSS2SS_BVH_DCR: bResult = Contact_DCRMS2_DCRMS2_BVH_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    // case Context::eSS2SS_Bruteforce_MSS_Lazy_DCR: bResult = Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR( p_mesh1, mesh_tr1, p_mesh_dof1, p_mesh2, mesh_tr2, p_mesh_dof2, cd, p_cc, p_context ); break;
    case Context::eSS2SS_BVH_SS_Lazy_DCR:
    default:
        bResult = Contact_DCRTS3_DCRTS3_BVH_TSS_Lazy_DCR( p_solid1, solid_tr1, p_solid_dof1,
                                                          p_solid2, solid_tr2, p_solid_dof2,
                                                          cd, p_cc, p_context );
        break;
    }
    GEO_NP_END_PROFILE_BLOCK_MS( contact.m_Profile_TSS3_TSS3_ms );
    GEO_NP_END_PROFILE_BLOCK_ACC_MS( mig2015.m_Profile_TSS3_TSS3_ms );

    return bResult;
}

}} //namespace geo::np
