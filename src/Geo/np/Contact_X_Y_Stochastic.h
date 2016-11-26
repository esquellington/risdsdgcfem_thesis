#ifndef GEO_NP_CONTACT_X_Y_STOCHASTIC_H
#define GEO_NP_CONTACT_X_Y_STOCHASTIC_H

#include "ContactData.h"
#include "Context.h"
#include <Geo/shape/IShape.h>

namespace geo {
namespace np {

//--------------------------------------------------------------------------------------------------------------------------------
// IShape Vs IShape
//--------------------------------------------------------------------------------------------------------------------------------

bool Contact_IShape2_IShape2_Stochastic( const IShape2* p_shape1, const Transform2& tr1, const Real* vec_dof1,
                                         const IShape2* p_shape2, const Transform2& tr2, const Real* vec_dof2,
                                         ContactData2& cd, ContactCache2* p_cc,
                                         const Context* p_context = g_pDefaultContext );

}} //namespace geo::np

#endif // GEO_NP_CONTACT_X_Y_STOCHASTIC_H
