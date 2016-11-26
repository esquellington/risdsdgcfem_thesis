#ifndef GEO_NP_CONTACT_X_Y_BRUTEFORCE_H
#define GEO_NP_CONTACT_X_Y_BRUTEFORCE_H

#include "ContactData.h"
#include "Context.h"

namespace geo { class MeshSolidShape2; class TetSolidShape3; }

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// MeshSolid2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_MeshSolid2_Plane2_BruteForce( const MeshSolidShape2* p_mesh, const Transform2& mesh_tr, const Real *vec_dof,
                                           const Vec2& plane_normal, Real plane_coeff_d,
                                           ContactData2& cd, ContactCache2* p_cc,
                                           const np::Context* p_context = g_pDefaultContext );

bool Contact_MeshSolid2_MeshSolid2_BruteForce( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real *vec_dof1,
                                               const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real *vec_dof2,
                                               ContactData2& cd, ContactCache2* p_cc,
                                               const np::Context* p_context = g_pDefaultContext );

bool Contact_DCRMS2_DCRMS2_BruteForce_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real *vec_dof1,
                                           const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real *vec_dof2,
                                           ContactData2& cd, ContactCache2* p_cc,
                                           const np::Context* p_context = g_pDefaultContext );

bool Contact_DCRMS2_DCRMS2_BruteForce_MSS_Lazy_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
                                                    const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
                                                    ContactData2& cd, ContactCache2* p_cc,
                                                    const np::Context* p_context = g_pDefaultContext );

//--------------------------------------------------------------------------------------------------------------------------------
// TetSolid3 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_TetSolid3_Plane3_BruteForce( const TetSolidShape3* p_solid, const Transform3& mesh_tr, const Real* vec_dof,
                                          const Vec3& plane_normal, Real plane_coeff_d,
                                          ContactData3& cd, ContactCache3* p_cc,
                                          const np::Context* p_context );

}} //namespace geo::np

#endif // GEO_NP_CONTACT_X_Y_BRUTEFORCE_H
