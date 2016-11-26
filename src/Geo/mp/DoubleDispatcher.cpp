#include "DoubleDispatcher.h"
#include <Geo/shape/shape.h>
#include <Geo/np/Contact.h>
#include "TestContact.h"

namespace geo {
namespace mp {

//! Stub for unimplemented type pairs
static bool TestContact_NotImplemented( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                                        const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                                        np::ContactData2& cd, np::ContactCache2* p_cc,
                                        const np::Context* p_context )
{
    GEO_LOG_ERROR( "geo::mp::TestContact() Shape type pair (%d,%d) not yet implemented",
                   (int)p_shape1->GetType(), (int)p_shape2->GetType() );
    return false;
}

static bool TestContact_NotImplemented( const IShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                                        const IShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                                        np::ContactData3& cd, np::ContactCache3* p_cc,
                                        const np::Context* p_context )
{
    GEO_LOG_ERROR( "geo::mp::TestContact() Shape type pair (%d,%d) not yet implemented",
                   (int)p_shape1->GetType(), (int)p_shape2->GetType() );
    return false;
}

//! Symmetric Contact specific call generator.
template <class ShapeT1, class ShapeT2>
bool GTestContact( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                   const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                   np::ContactData2& cd, np::ContactCache2* p_cc,
                   const np::Context* p_context )
{
    GEO_ASSERT( p_shape1->GetType() == shape_type_of<ShapeT1>() && p_shape2->GetType() == shape_type_of<ShapeT2>() );
    return TestContact( static_cast<const ShapeT1*>(p_shape1), tr1, p_dof1,
                        static_cast<const ShapeT2*>(p_shape2), tr2, p_dof2,
                        cd, p_cc );
}

//! Symmetric Contact specific call generator.
template <class ShapeT1, class ShapeT2>
bool GTestContact( const IShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                   const IShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                   np::ContactData3& cd, np::ContactCache3* p_cc,
                   const np::Context* p_context )
{
    GEO_ASSERT( p_shape1->GetType() == shape_type_of<ShapeT1>() && p_shape2->GetType() == shape_type_of<ShapeT2>() );
    return TestContact( static_cast<const ShapeT1*>(p_shape1), tr1, p_dof1,
                        static_cast<const ShapeT2*>(p_shape2), tr2, p_dof2,
                        cd, p_cc );
}

//--------------------------------------------------------------------------------------------------------------------------------
DoubleDispatcher::DoubleDispatcher()
{
    MakeDefault();
}
DoubleDispatcher::~DoubleDispatcher()
{
}

void DoubleDispatcher::MakeDefault()
{
    MakeDefault_Contact();
    //MakeDefault_Overlap();
    //MakeDefault_Proximity();
}

void DoubleDispatcher::MakeDefault_Contact()
{
    //---- 2D
    for( int i=0; i<cNumShapeTypes; i++ )
        for( int j=0; j<cNumShapeTypes; j++ )
            m_DDFST2D[i][j].SetFunction_Contact( TestContact_NotImplemented, false );
    // eShape_Sphere2
    m_DDFST2D[eShape_Sphere2][eShape_Sphere2].SetFunction_Contact( GTestContact<SphereShape2,SphereShape2>, false );
    m_DDFST2D[eShape_Sphere2][eShape_Plane2].SetFunction_Contact( GTestContact<SphereShape2,PlaneShape2>, false );
    m_DDFST2D[eShape_Plane2][eShape_Sphere2].SetFunction_Contact( GTestContact<SphereShape2,PlaneShape2>, true );
    m_DDFST2D[eShape_Sphere2][eShape_MeshSolid2].SetFunction_Contact( GTestContact<SphereShape2,MeshSolidShape2>, false );
    m_DDFST2D[eShape_MeshSolid2][eShape_Sphere2].SetFunction_Contact( GTestContact<SphereShape2,MeshSolidShape2>, true );

    // eShape_MeshSolid2
    m_DDFST2D[eShape_MeshSolid2][eShape_MeshSolid2].SetFunction_Contact( GTestContact<MeshSolidShape2,MeshSolidShape2>, false );
    m_DDFST2D[eShape_MeshSolid2][eShape_Plane2].SetFunction_Contact( GTestContact<MeshSolidShape2,PlaneShape2>, false );
    m_DDFST2D[eShape_Plane2][eShape_MeshSolid2].SetFunction_Contact( GTestContact<MeshSolidShape2,PlaneShape2>, true );

    //---- 3D
    for( int i=0; i<cNumShapeTypes; i++ )
        for( int j=0; j<cNumShapeTypes; j++ )
            m_DDFST3D[i][j].SetFunction_Contact( TestContact_NotImplemented, false );
    // eShape_Sphere3
    m_DDFST3D[eShape_Sphere3][eShape_Sphere3].SetFunction_Contact( GTestContact<SphereShape3,SphereShape3>, false );
    m_DDFST3D[eShape_Sphere3][eShape_Plane3].SetFunction_Contact( GTestContact<SphereShape3,PlaneShape3>, false );
    m_DDFST3D[eShape_Plane3][eShape_Sphere3].SetFunction_Contact( GTestContact<SphereShape3,PlaneShape3>, true );

    // eShape_TetSolid3
    m_DDFST3D[eShape_TetSolid3][eShape_TetSolid3].SetFunction_Contact( GTestContact<TetSolidShape3,TetSolidShape3>, false );
    m_DDFST3D[eShape_TetSolid3][eShape_Plane3].SetFunction_Contact( GTestContact<TetSolidShape3,PlaneShape3>, false );
    m_DDFST3D[eShape_Plane3][eShape_TetSolid3].SetFunction_Contact( GTestContact<TetSolidShape3,PlaneShape3>, true );
}

void DoubleDispatcher::SetFunction_Contact( EShapeType st1, EShapeType st2, test_contact_pair_function_2_type tcpf )
{
    m_DDFST2D[st1][st2].SetFunction_Contact( tcpf, false );
    m_DDFST2D[st2][st1].SetFunction_Contact( tcpf, true );
}

void DoubleDispatcher::SetFunction_Contact( EShapeType st1, EShapeType st2, test_contact_pair_function_3_type tcpf )
{
    m_DDFST3D[st1][st2].SetFunction_Contact( tcpf, false );
    m_DDFST3D[st2][st1].SetFunction_Contact( tcpf, true );
}

bool DoubleDispatcher::TestContact( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                                    const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                                    np::ContactData2& cd, np::ContactCache2* p_cc )
{
    int st1( p_shape1->GetType() );
    int st2( p_shape2->GetType() );
    if( m_DDFST2D[st1][st2].m_Flags.Test(FunctionSet2::eFlipContact) )
    {
        if( m_DDFST2D[st1][st2].m_Contact( p_shape2, tr2, p_dof2,
                                           p_shape1, tr1, p_dof1,
                                           cd, p_cc,
                                           &m_Context ) )
        {
            cd.Flip();
            return true;
        }
        else
            return false;
    }
    else
        return m_DDFST2D[st1][st2].m_Contact( p_shape1, tr1, p_dof1,
                                              p_shape2, tr2, p_dof2,
                                              cd, p_cc,
                                              &m_Context );
}

bool DoubleDispatcher::TestContact( const IShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                                    const IShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                                    np::ContactData3& cd, np::ContactCache3* p_cc )
{
    int st1( p_shape1->GetType() );
    int st2( p_shape2->GetType() );
    if( m_DDFST3D[st1][st2].m_Flags.Test(FunctionSet3::eFlipContact) )
    {
        if( m_DDFST3D[st1][st2].m_Contact( p_shape2, tr2, p_dof2,
                                           p_shape1, tr1, p_dof1,
                                           cd, p_cc,
                                           &m_Context ) )
        {
            cd.Flip();
            return true;
        }
        else
            return false;
    }
    else
        return m_DDFST3D[st1][st2].m_Contact( p_shape1, tr1, p_dof1,
                                              p_shape2, tr2, p_dof2,
                                              cd, p_cc,
                                              &m_Context );
}

}} //namespace geo::mp
