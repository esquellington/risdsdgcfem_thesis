#ifndef GEO_NP_CONTACT_H
#define GEO_NP_CONTACT_H

#include "ContactData.h"
#include "Context.h"

/*! This file contains the prototypes for all default Contact_X_Y functionality.
  Each method may internally use any available implementation transparently (_Analytic, _Bruteforce, _SAT, _GJK, _Stochastic, ...)
*/

//Fwd declarations
namespace geo { class MeshSolidShape2; class TetSolidShape3; }

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// Sphere2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_Sphere2_Sphere2( const Vec2& pos1, Real radius1,
                              const Vec2& pos2, Real radius2,
                              ContactData2& cd, ContactCache2* p_cc,
                              const Context* p_context = g_pDefaultContext );

//Contact_Sphere2_Capsule2
//Contact_Sphere2_Box2
//Contact_Sphere2_Triangle2
//Contact_Sphere2_Path2

bool Contact_Sphere2_Plane2( const Vec2& sphere_pos, Real sphere_radius,
                             const Vec2& plane_normal, Real plane_coeff_d,
                             ContactData2& cd, ContactCache2* p_cc,
                             const Context* p_context = g_pDefaultContext );

bool Contact_Sphere2_MeshSolid2( const Vec2& sphere_pos, Real sphere_radius,
                                 const MeshSolidShape2* p_mesh, const Transform2& mesh_tr, const Real* p_mesh_dof,
                                 ContactData2& cd, ContactCache2* p_cc,
                                 const Context* p_context = g_pDefaultContext );

//--------------------------------------------------------------------------------------------------------------------------------
// Sphere3 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_Sphere3_Sphere3( const Vec3& pos1, Real radius1,
                              const Vec3& pos2, Real radius2,
                              ContactData3& cd, ContactCache3* p_cc,
                              const Context* p_context = g_pDefaultContext );

bool Contact_Sphere3_Plane3( const Vec3& sphere_pos, Real sphere_radius,
                             const Vec3& plane_normal, Real plane_coeff_d,
                             ContactData3& cd, ContactCache3* p_cc,
                             const Context* p_context = g_pDefaultContext );


//--------------------------------------------------------------------------------------------------------------------------------
// MeshSolid2 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_MeshSolid2_Plane2( const MeshSolidShape2* p_mesh, const Transform2& mesh_transform, const Real* p_mesh_dof,
                                const Vec2& plane_normal, Real plane_coeff_d,
                                ContactData2& cd, ContactCache2* p_cc,
                                const Context* p_context = g_pDefaultContext );

bool Contact_MeshSolid2_MeshSolid2( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* p_mesh_dof1,
                                    const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* p_mesh_dof2,
                                    ContactData2& cd, ContactCache2* p_cc,
                                    const Context* p_context = g_pDefaultContext );

//--------------------------------------------------------------------------------------------------------------------------------
// TetSolid3 Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_TetSolid3_Plane3( const TetSolidShape3* p_solid, const Transform3& solid_tr, const Real* p_solid_dof,
                               const Vec3& plane_normal, Real plane_coeff_d,
                               ContactData3& cd, ContactCache3* p_cc,
                               const Context* p_context = g_pDefaultContext );

bool Contact_TetSolid3_TetSolid3( const TetSolidShape3* p_solid1, const Transform3& solid_tr1, const Real* p_solid_dof1,
                                  const TetSolidShape3* p_solid2, const Transform3& solid_tr2, const Real* p_solid_dof2,
                                  ContactData3& cd, ContactCache3* p_cc,
                                  const Context* p_context = g_pDefaultContext );

}} //namespace geo::np

#endif // GEO_NP_CONTACT_H
