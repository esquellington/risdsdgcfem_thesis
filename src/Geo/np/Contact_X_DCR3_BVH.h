#ifndef GEO_NP_CONTACT_X_DCR3_BVH_H
#define GEO_NP_CONTACT_X_DCR3_BVH_H

#include "ContactData.h"
#include "Context.h"

namespace geo { class TetSolidShape3; }

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// TetSolidShape3/DCR Vs Y
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_DCRTS3_Plane3_BVH_TSS_Lazy_DCR( const TetSolidShape3* p_solid, const Transform3& solid_tr, const Real* vec_dof,
                                             const Vec3& plane_normal, Real plane_coeff_d,
                                             ContactData3& cd, ContactCache3* p_cc,
                                             const Context* p_context = g_pDefaultContext );

bool Contact_DCRTS3_DCRTS3_BVH_TSS_Lazy_DCR( const TetSolidShape3* p_solid1, const Transform3& solid_tr1, const Real* vec_dof1,
                                             const TetSolidShape3* p_solid2, const Transform3& solid_tr2, const Real* vec_dof2,
                                             ContactData3& cd, ContactCache3* p_cc,
                                             const Context* p_context = g_pDefaultContext );

}} //namespace geo::np

#endif // GEO_NP_CONTACT_X_DCR3_BVH_H
