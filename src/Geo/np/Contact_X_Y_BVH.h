#ifndef GEO_NP_CONTACT_X_Y_BVH_H
#define GEO_NP_CONTACT_X_Y_BVH_H

#include "ContactData.h"
#include "Context.h"

namespace geo { class MeshSolidShape2; }

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// MeshSolid2 DCR Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_DCRMS2_Plane2_BVH_MSS_Lazy_DCR( const MeshSolidShape2* p_mesh, const Transform2& mesh_tr, const Real* vec_dof,
                                             const Vec2& plane_normal, Real plane_coeff_d,
                                             ContactData2& cd, ContactCache2* p_cc,
                                             const Context* p_context = g_pDefaultContext );

// bool Contact_DCRMS2_DCRMS2_BVH_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
//                                     const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
//                                     ContactData2& cd, ContactCache2* p_cc,
//                                     const np::Context* p_context = g_pDefaultContext );

bool Contact_DCRMS2_DCRMS2_BVH_MSS_Lazy_DCR( const MeshSolidShape2* p_mesh1, const Transform2& mesh_tr1, const Real* vec_dof1,
                                             const MeshSolidShape2* p_mesh2, const Transform2& mesh_tr2, const Real* vec_dof2,
                                             ContactData2& cd, ContactCache2* p_cc,
                                             const Context* p_context = g_pDefaultContext );

}} //namespace geo::np

#endif // GEO_NP_CONTACT_X_Y_BVH_H
