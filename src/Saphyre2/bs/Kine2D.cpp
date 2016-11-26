#include <Saphyre2/bs/Kine2D.h>
#include <Saphyre2/bs/Universe.h>
#include <Saphyre2/ds/Commands.h>
#include <Mal/GConversion.h>

namespace S2
{

// Constructor/Destructor and Methods with controlled scope
Kine2D::Kine2D( ISyncEntity *p_parent )
: ISyncEntity(p_parent)
, m_ShapeId(geo::cInvalidShapeId)
, m_NumDOF(0), m_vecDOF(0)
, m_pShape(0)
{
    m_Transform = Transform2::Identity();
}

Kine2D::~Kine2D()
{
    if( m_vecDOF ) delete[] m_vecDOF;
    if( m_pShape ) delete m_pShape; //\todo Assumes exclusive shape!
}

bool Kine2D::Update( Real dt )
{
    if( !IsValid() ) return false;

    BS_BCMD("Kine2D::Update()");

    GetChannel()->BeginCommand( ds::eCmd_Update );
    {
        GetChannel()->Write( "eid", GetEID() );
        GetChannel()->Write( "dt", (Real)dt );
    }
    GetChannel()->EndCommand();

    // Imperative sync
    bool bResult = Sync_Internal();
    BS_ECMD( bResult );
    return bResult;
}

//---- ISyncEntity internal protocol

//! \pre !IsCreated()
void Kine2D::BeginDef_Internal()
{
    // Default params already set
}

//! \pre IsBeingDefined()
bool Kine2D::EndDef_Internal()
{
    BS_BCMD("Kine2D::EndDef()");

    //---- Send creation command
    GetChannel()->BeginCommand( ds::eCmd_Create );
    {
        GetChannel()->Write( "uid", GetUID() );
        GetChannel()->Write( "peid", GetParent()->GetEID() );
        GetChannel()->Write( "entity_type", (uint32)ds::eEntity_Geom );
        GetChannel()->Write( "geom_type", (uint32)ds::eGeom_Simple );

        // Kine Init
        GetChannel()->Write( "flags", m_Params.m_Flags );
        GetChannel()->Write( "transform", m_Transform );

        // Shape Init
        BS_ASSERT( BSG::GetShapeLibrary().Lookup(m_ShapeId).IsValid() );
        if( false ) //\todo accept m_InternalShapeId and pass it here... geo::cInvalidShapeId != m_ShapeId )
            GetChannel()->Write( "shape_id", (uint32)m_ShapeId );
        else
            GetChannel()->WriteItem( "shape_def", BSG::GetShapeLibrary().Lookup(m_ShapeId) );
    }
    GetChannel()->EndCommand();

    bool bResult(true);
    BS_ECMD( bResult );
    return bResult;
}

//! \pre IsTouched() and was IsLocked()
void Kine2D::Unlock_Internal()
{
    BS_BCMD("Kine2D::EndDef()");

    //---- Send edit command
    GetChannel()->BeginCommand( ds::eCmd_Edit );
    {
        GetChannel()->Write( "eid", GetEID() );
        if( IsTouched(eTouchedTransform) ) GetChannel()->Write( "transform", m_Transform );
        if( IsTouched(eTouchedDOF) ) GetChannel()->WriteArray( "vec_dof", m_vecDOF, m_NumDOF );
    }
    GetChannel()->EndCommand();

    BS_ECMD(true);
}

//---- Kine2D Initialization/Edition Methods
bool Kine2D::SetParams( const Params &params )
{
    if( !IsBeingDefined() ) return false;
    m_Params = params;
    return true;
}

bool Kine2D::AttachShape( geo::ShapeID shape_id )
{
    if( !IsBeingDefined() ) return false;
    BS_ASSERT( BSG::GetShapeLibrary().Lookup(shape_id).IsValid() );
    m_ShapeId = shape_id;
    return true;
}

bool Kine2D::AttachShape( const geo::IShape2 &r_shape )
{
    if( !IsBeingDefined() ) return false;
    m_ShapeId = BSG::GetShapeLibrary().Register(&r_shape);
    return true;
}

const geo::IShape2* Kine2D::GetShape() const
{
    //\todo determine when can it be called
    BS_ASSERT( m_ShapeId != geo::cInvalidShapeId );
    if( !m_pShape ) m_pShape = static_cast<geo::IShape2*>( S2::BSG::GetShapeFactory().CreateExclusive( m_ShapeId ) );
    return m_pShape;
}

geo::IShape2* Kine2D::GetShape()
{
    BS_ASSERT( m_ShapeId != geo::cInvalidShapeId );
    if( !m_pShape ) m_pShape = static_cast<geo::IShape2*>( S2::BSG::GetShapeFactory().CreateExclusive( m_ShapeId ) );
    return m_pShape;
}

bool Kine2D::SetPos( const Point2 &pos )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Transform.m_Pos = pos;
    Touch( eTouchedTransform );
    return true;
}

bool Kine2D::SetRot( const Mat2x2 &rot )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Transform.m_Rot = rot;
    Touch( eTouchedTransform );
    return true;
}

bool Kine2D::SetAngle( Real angle )
{
    return SetRot( mal::GRotation2x2_From(angle) );
}

bool Kine2D::SetTransform( const Transform2 &transform )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    m_Transform = transform;
    Touch( eTouchedTransform );
    return true;
}

bool Kine2D::SetVecDOF( const Real *vec_dof )
{
    if( !IsLocked() && !IsBeingDefined() ) return false;
    for( uint32 i=0; i<m_NumDOF; i++ ) m_vecDOF[i] = vec_dof[i];
    Touch( eTouchedDOF );
    return true;
}

//---- Internal ISyncEntity methods
bool Kine2D::ProcessUpdate( const ds::ReturnIt &rit )
{
    BS_ASSERT(false); //\todo Kines are NEVER updated in DS by now, this message should NOT be received
    m_Transform = rit.Find("transform").Get<Transform2>();

    // Create DOF vector if it's first update (aka Init())
    ds::ReturnIt rit_dof( rit.Find("vec_dof") );
    if( 0 == m_vecDOF )
    {
        m_NumDOF = rit_dof.GetArrayCount();
        m_vecDOF = new Real[m_NumDOF];
    }
    const Real *vec_dof = rit_dof.GetArrayPtr<Real>();
    for( uint32 i=0; i<m_NumDOF; i++ ) m_vecDOF[i] = vec_dof[i];
    return true;
}

} // namespace S2
