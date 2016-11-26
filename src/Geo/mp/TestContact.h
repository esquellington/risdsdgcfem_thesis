#ifndef GEO_MP_TEST_CONTACT_H
#define GEO_MP_TEST_CONTACT_H

#include <Geo/IObject.h>
#include <Geo/np/Context.h>
#include <Geo/np/ContactData.h>

namespace geo {
namespace mp {

/*\note Pairwise-specific methods follow common static prototype to be
  called from GTestContact<ShapeT1,ShapeT2> and use the
  pairwise-default strategy (Analytic, Bruteforce, Stochastic...), as
  defined by the specific np::Contact_X_Y function.
*/

//---- Object Vs Object
bool TestContact( const IObject2* p_obj1, const IObject2* p_obj2, np::ContactData2& cd, np::ContactCache2* p_cc );
bool TestContact( const IObject3* p_obj1, const IObject3* p_obj2, np::ContactData3& cd, np::ContactCache3* p_cc );

//---- Sphere2 vs X
bool TestContact( const SphereShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const SphereShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );
bool TestContact( const SphereShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const PlaneShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );
bool TestContact( const SphereShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const MeshSolidShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );

//---- Sphere3 vs X
bool TestContact( const SphereShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const SphereShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );
bool TestContact( const SphereShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const PlaneShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );

//---- MeshSolid2 vs X
bool TestContact( const MeshSolidShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const PlaneShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );
bool TestContact( const MeshSolidShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const MeshSolidShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );

//---- TetSolid3 vs X
bool TestContact( const TetSolidShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const PlaneShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );
bool TestContact( const TetSolidShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const TetSolidShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context = np::g_pDefaultContext );

}} //namespace geo::mp

#endif // GEO_MP_TEST_CONTACT_H
