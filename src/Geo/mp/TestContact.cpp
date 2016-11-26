#include "TestContact.h"
#include <Geo/shape/shape.h>
#include <Geo/np/Contact.h> //Default Contact_X_Y() methods
#include <Geo/np/Contact_X_Y_Stochastic.h> //TEMP

#define __USE_DOUBLE_DISPATCHER
#ifdef __USE_DOUBLE_DISPATCHER
#include "DoubleDispatcher.h"
#endif

namespace geo {
namespace mp {

//---- Object Vs Object
bool TestContact( const IObject2* p_obj1, const IObject2* p_obj2, np::ContactData2& cd, np::ContactCache2* p_cc )
{
#ifdef __USE_DOUBLE_DISPATCHER
    GEO_ASSERT( p_obj1->GetDimension() == p_obj2->GetDimension() && p_obj1->GetDimension() == 2 );
    const IShape2* pShape1( static_cast<const IShape2*>(p_obj1->GetShapeInterface()) );
    const IShape2* pShape2( static_cast<const IShape2*>(p_obj2->GetShapeInterface()) );
    const IObject2* pObj1( static_cast<const IObject2*>(p_obj1) );
    const IObject2* pObj2( static_cast<const IObject2*>(p_obj2) );
    return g_pDefaultDoubleDispatcher->TestContact( pShape1, pObj1->GetTransform(), p_obj1->GetVecDOF(),
                                                    pShape2, pObj2->GetTransform(), p_obj2->GetVecDOF(),
                                                    cd, p_cc );
#else
    GEO_ASSERT( p_obj1->GetDimension() == p_obj2->GetDimension() && p_obj1->GetDimension() == 2 );
    const IShape2* pShape1( static_cast<const IShape2*>(p_obj1->GetShapeInterface()) );
    const IShape2* pShape2( static_cast<const IShape2*>(p_obj2->GetShapeInterface()) );
    const IObject2* pObj1( static_cast<const IObject2*>(p_obj1) );
    const IObject2* pObj2( static_cast<const IObject2*>(p_obj2) );
    if( eShape_Sphere2 == pShape1->GetType() && eShape_Sphere2 == pShape2->GetType() )
        return TestContact( static_cast<const SphereShape2*>(pShape1), pObj1->GetTransform(), 0,
                            static_cast<const SphereShape2*>(pShape2), pObj2->GetTransform(), 0,
                            cd, p_cc );
    else if( eShape_Sphere2 == pShape1->GetType() && eShape_Plane2 == pShape2->GetType() )
        return TestContact( static_cast<const SphereShape2*>(pShape1), pObj1->GetTransform(), 0,
                            static_cast<const PlaneShape2*>(pShape2), pObj2->GetTransform(), 0,
                            cd, p_cc );
    else
        GEO_LOG_WARNING( "geo::mp::TestContact: Shape type pair (%d,%d) not yet implemented",
                         pShape1->GetType(), pShape2->GetType() );
    return false;
#endif
}

bool TestContact( const IObject3* p_obj1, const IObject3* p_obj2, np::ContactData3& cd, np::ContactCache3* p_cc )
{
#ifdef __USE_DOUBLE_DISPATCHER
    GEO_ASSERT( p_obj1->GetDimension() == p_obj2->GetDimension() && p_obj1->GetDimension() == 3 );
    const IShape3* pShape1( static_cast<const IShape3*>(p_obj1->GetShapeInterface()) );
    const IShape3* pShape2( static_cast<const IShape3*>(p_obj2->GetShapeInterface()) );
    const IObject3* pObj1( static_cast<const IObject3*>(p_obj1) );
    const IObject3* pObj2( static_cast<const IObject3*>(p_obj2) );
    return g_pDefaultDoubleDispatcher->TestContact( pShape1, pObj1->GetTransform(), p_obj1->GetVecDOF(),
                                                    pShape2, pObj2->GetTransform(), p_obj2->GetVecDOF(),
                                                    cd, p_cc );
#else
    GEO_ASSERT( p_obj1->GetDimension() == p_obj2->GetDimension() && p_obj1->GetDimension() == 3 );
    const IShape3* pShape1( static_cast<const IShape3*>(p_obj1->GetShapeInterface()) );
    const IShape3* pShape2( static_cast<const IShape3*>(p_obj2->GetShapeInterface()) );
    const IObject3* pObj1( static_cast<const IObject3*>(p_obj1) );
    const IObject3* pObj2( static_cast<const IObject3*>(p_obj2) );
    if( eShape_Sphere3 == pShape1->GetType() && eShape_Sphere3 == pShape2->GetType() )
        return TestContact( static_cast<const SphereShape3*>(pShape1), pObj1->GetTransform(), 0,
                            static_cast<const SphereShape3*>(pShape2), pObj2->GetTransform(), 0,
                            cd, p_cc );
    else if( eShape_Sphere3 == pShape1->GetType() && eShape_Plane3 == pShape2->GetType() )
        return TestContact( static_cast<const SphereShape3*>(pShape1), pObj1->GetTransform(), 0,
                            static_cast<const PlaneShape3*>(pShape2), pObj2->GetTransform(), 0,
                            cd, p_cc );
    else
        GEO_LOG_WARNING( "geo::mp::TestContact: Shape type pair (%d,%d) not yet implemented",
                         pShape1->GetType(), pShape2->GetType() );
    return false;
#endif
}

//---- Sphere2 vs X
bool TestContact( const SphereShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const SphereShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_Sphere2_Sphere2( tr1.m_Pos, p_shape1->GetRadius(),
                                        tr2.m_Pos, p_shape2->GetRadius(),
                                        cd, p_cc, p_context );
}

bool TestContact( const SphereShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const PlaneShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_Sphere2_Plane2( tr1.m_Pos, p_shape1->GetRadius(),
                                       p_shape2->GetNormal(), p_shape2->GetCoeffD(),
                                       cd, p_cc, p_context );
}

bool TestContact( const SphereShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const MeshSolidShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context )
{
    //TEMP: Force stochastic as default by now, testing...
    return np::Contact_IShape2_IShape2_Stochastic( p_shape1, tr1, p_dof1,
                                                   p_shape2, tr2, p_dof2,
                                                   cd, p_cc, p_context );
    /*
    return np::Contact_Sphere2_MeshSolid2( tr1.m_Pos, p_shape1->GetRadius(),
                                           p_shape2, tr2, p_dof2,
                                           cd, p_cc, p_context );
    */
}


bool TestContact( const SphereShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const SphereShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_Sphere3_Sphere3( tr1.m_Pos, p_shape1->GetRadius(),
                                        tr2.m_Pos, p_shape2->GetRadius(),
                                        cd, p_cc, p_context );
}

bool TestContact( const SphereShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const PlaneShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_Sphere3_Plane3( tr1.m_Pos, p_shape1->GetRadius(),
                                       p_shape2->GetNormal(), p_shape2->GetCoeffD(),
                                       cd, p_cc, p_context );
}

//---- MeshSolid2 vs X
bool TestContact( const MeshSolidShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const PlaneShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_MeshSolid2_Plane2( p_shape1, tr1, p_dof1,
                                          p_shape2->GetNormal(), p_shape2->GetCoeffD(),
                                          cd, p_cc, p_context );
}

bool TestContact( const MeshSolidShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                  const MeshSolidShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                  np::ContactData2& cd, np::ContactCache2* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_MeshSolid2_MeshSolid2( p_shape1, tr1, p_dof1,
                                              p_shape2, tr2, p_dof2,
                                              cd, p_cc, p_context );
}

//---- TetSolid3 vs X
bool TestContact( const TetSolidShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const PlaneShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_TetSolid3_Plane3( p_shape1, tr1, p_dof1,
                                         p_shape2->GetNormal(), p_shape2->GetCoeffD(),
                                         cd, p_cc, p_context );
}

bool TestContact( const TetSolidShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                  const TetSolidShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                  np::ContactData3& cd, np::ContactCache3* p_cc,
                  const np::Context* p_context )
{
    return np::Contact_TetSolid3_TetSolid3( p_shape1, tr1, p_dof1,
                                            p_shape2, tr2, p_dof2,
                                            cd, p_cc, p_context );
}

}} //namespace geo::mp
