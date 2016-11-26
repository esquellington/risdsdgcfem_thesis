#include "ShapeLibrary.h"
#include <Geo/shape/shape.h> //shape types

namespace geo {


ShapeLibrary::ShapeLibrary() : m_pSDS(0), m_bIsExternalSDS(false), m_LastShapeId(0) {}
// ShapeLibrary::ShapeLibrary( ShapeDefStream *p_sds ) : m_pSDS(p_sds), m_bIsExternalSDS(true), m_LastShapeId(0) {}  \todo Disabled to avoid realloc conundrums

ShapeLibrary::~ShapeLibrary()
{
    if(!m_bIsExternalSDS && 0 != m_pSDS) delete m_pSDS;
}

bool ShapeLibrary::Reserve( uint32 sds_size, Flags32 realloc_flags )
{
    if(m_bIsExternalSDS || 0 != m_pSDS) return false;
    if( 0 != m_pSDS ) delete m_pSDS;
    m_pSDS = new ShapeDefStream( sds_size, 1<<10, realloc_flags );
    return 0 != m_pSDS;
}

//---- Generic shape definition
ShapeID ShapeLibrary::Register( const IShape &shape, const char* name ) { return Register( &shape, name ); }
ShapeID ShapeLibrary::Register( const IShape *p_shape, const char* name )
{
    switch( p_shape->GetType() )
    {
    case eShape_Plane2: return Register_Plane2( static_cast<const PlaneShape2*>(p_shape)->GetNormal(),
                                                static_cast<const PlaneShape2*>(p_shape)->GetCoeffD(),
                                                static_cast<const PlaneShape2*>(p_shape)->IsHalfSpace(), name ); break;
    case eShape_Plane3: return Register_Plane3( static_cast<const PlaneShape3*>(p_shape)->GetNormal(),
                                                static_cast<const PlaneShape3*>(p_shape)->GetCoeffD(),
                                                static_cast<const PlaneShape3*>(p_shape)->IsHalfSpace(), name ); break;
    case eShape_Sphere2: return Register_Sphere2( static_cast<const SphereShape2*>(p_shape)->GetRadius(), name ); break;
    case eShape_Sphere3: return Register_Sphere3( static_cast<const SphereShape3*>(p_shape)->GetRadius(), name ); break;
    case eShape_Capsule2: return Register_Capsule2( static_cast<const CapsuleShape2*>(p_shape)->GetRadius(),
                                                    static_cast<const CapsuleShape2*>(p_shape)->GetHalfHeight(), name ); break;
    case eShape_Capsule3: return Register_Capsule3( static_cast<const CapsuleShape3*>(p_shape)->GetRadius(),
                                                    static_cast<const CapsuleShape3*>(p_shape)->GetHalfHeight(), name ); break;
    case eShape_Box2: return Register_Box2( static_cast<const BoxShape2*>(p_shape)->GetHalfSizes(), name ); break;
    case eShape_Box3: return Register_Box3( static_cast<const BoxShape3*>(p_shape)->GetHalfSizes(), name ); break;
    case eShape_Polygonal2:
        {
            const PolygonalShape2 *pPS2( static_cast<const PolygonalShape2*>(p_shape) );
            return Register_Polygonal2( pPS2, name );
        }
        break;
    case eShape_Polygonal3:
        {
            const PolygonalShape3 *pPS3( static_cast<const PolygonalShape3*>(p_shape) );
            return Register_Polygonal3( pPS3, name );
        }
        break;
    case eShape_Path2:
        {
            const PathShape2 *pSS2( static_cast<const PathShape2*>(p_shape) );
            return Register_PathShape2( pSS2, name );
        }
        break;
    case eShape_MeshSolid2:
        {
            const MeshSolidShape2 *pMSS2( static_cast<const MeshSolidShape2*>(p_shape) );
            return Register_MeshSolidShape2( pMSS2, name );
        }
        break;
    case eShape_TetSolid3:
        {
            const TetSolidShape3 *pTSS3( static_cast<const TetSolidShape3*>(p_shape) );
            return Register_TetSolidShape3( pTSS3, name );
        }
        break;
    case eShape_TriSurface3:
        {
            const TriSurfaceShape3 *pTSS3( static_cast<const TriSurfaceShape3*>(p_shape) );
            return Register_TriSurfaceShape3( pTSS3, name );
        }
        break;
    default: return cInvalidShapeId; break;
    }
}

ShapeID ShapeLibrary::Register( const ShapeDef &shape_def, const char* name )
{
    if( eType_Property_Group == shape_def.GetType() ) { GEO_LOG_WARNING( "shape_def is of eType_Property_Group, assuming correct... (loaded from TXT?)" ); }
    else GEO_ASSERT( eType_ShapeDef == shape_def.GetType() );
    if( 0 == name ) return m_pSDS->WriteItem( ++m_LastShapeId, shape_def ).GetId();
    else return m_pSDS->WriteItem( name, shape_def ).GetId();
}

//---- Simple shape definition
ShapeID ShapeLibrary::Register_Plane2( const Vec2 &normal, Real coeff_d, bool b_half_space, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Plane2 );
        m_pSDS->Write( "normal", normal );
        m_pSDS->Write( "coeff_d", coeff_d );
        m_pSDS->Write( "is_half_space", b_half_space );
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_Plane3( const Vec3 &normal, Real coeff_d, bool b_half_space, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Plane3 );
        m_pSDS->Write( "normal", normal );
        m_pSDS->Write( "coeff_d", coeff_d );
        m_pSDS->Write( "is_half_space", b_half_space );
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_Sphere2( Real radius, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Sphere2 );
        m_pSDS->Write( "radius", radius );
    }
    return m_pSDS->EndComplex().GetId();
}
ShapeID ShapeLibrary::Register_Sphere3( Real radius, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Sphere3 );
        m_pSDS->Write( "radius", radius );
    }
    return m_pSDS->EndComplex().GetId();
}
ShapeID ShapeLibrary::Register_Capsule2( Real radius, Real half_height, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Capsule2 );
        m_pSDS->Write( "radius", radius );
        m_pSDS->Write( "half_height", half_height );
    }
    return m_pSDS->EndComplex().GetId();
}
ShapeID ShapeLibrary::Register_Capsule3( Real radius, Real half_height, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Capsule3 );
        m_pSDS->Write( "radius", radius );
        m_pSDS->Write( "half_height", half_height );
    }
    return m_pSDS->EndComplex().GetId();
}
ShapeID ShapeLibrary::Register_Box2( const Vec2 &half_sizes, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Box2 );
        m_pSDS->Write( "half_sizes", half_sizes );
    }
    return m_pSDS->EndComplex().GetId();
}
ShapeID ShapeLibrary::Register_Box3( const Vec3 &half_sizes, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Box3 );
        m_pSDS->Write( "half_sizes", half_sizes );
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_Polygonal2( const PolygonalShape2 *p_ps2, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        uint32 num_vertices = p_ps2->GetNumVertices();
        m_pSDS->Write( "shape_type", (uint32)eShape_Polygonal2 );
        m_pSDS->Write( "version", (uint32)PolygonalShape2::cVersion );
        m_pSDS->Write( "num_vertices", num_vertices );
        m_pSDS->Write( "is_closed", p_ps2->IsClosed() );
        m_pSDS->Write( "radius", p_ps2->GetRadius() );
        if( 0 != p_ps2->GetPointsVDL() )
            m_pSDS->WriteArray( "vec_points",
                                reinterpret_cast<const uint8*>(p_ps2->GetPointsVDL()),
                                num_vertices*sizeof(Vec2) );
        if( 0 != p_ps2->GetNormalsVDL() )
            m_pSDS->WriteArray( "vec_normals",
                                reinterpret_cast<const uint8*>(p_ps2->GetNormalsVDL()),
                                num_vertices*sizeof(Vec2) );
        if( 0 != p_ps2->GetTangentsVDL() )
            m_pSDS->WriteArray( "vec_tangents",
                                reinterpret_cast<const uint8*>(p_ps2->GetTangentsVDL()),
                                num_vertices*sizeof(Vec2) );
        if( 0 != p_ps2->GetLambdasVDL() )
            m_pSDS->WriteArray( "vec_lambdas",
                                reinterpret_cast<const uint8*>(p_ps2->GetLambdasVDL()),
                                num_vertices*sizeof(Real) );
        if( 0 != p_ps2->GetFlagsVDL() )
            m_pSDS->WriteArray( "vec_flags",
                                reinterpret_cast<const uint8*>(p_ps2->GetFlagsVDL()),
                                num_vertices*sizeof(Flags32) );
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_Polygonal3( const PolygonalShape3 *p_ps3, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        uint32 num_vertices = p_ps3->GetNumVertices();
        m_pSDS->Write( "shape_type", (uint32)eShape_Polygonal3 );
        m_pSDS->Write( "version", (uint32)PolygonalShape3::cVersion );
        m_pSDS->Write( "num_vertices", num_vertices );
        m_pSDS->Write( "is_closed", p_ps3->IsClosed() );
        m_pSDS->Write( "radius", p_ps3->GetRadius() );
        if( 0 != p_ps3->GetPointsVDL() )
            m_pSDS->WriteArray( "vec_points",
                                reinterpret_cast<const uint8*>(p_ps3->GetPointsVDL()),
                                num_vertices*sizeof(Vec3) );
        if( 0 != p_ps3->GetNormalsVDL() )
            m_pSDS->WriteArray( "vec_normals",
                                reinterpret_cast<const uint8*>(p_ps3->GetNormalsVDL()),
                                num_vertices*sizeof(Vec3) );
        if( 0 != p_ps3->GetTangentsVDL() )
            m_pSDS->WriteArray( "vec_tangents",
                                reinterpret_cast<const uint8*>(p_ps3->GetTangentsVDL()),
                                num_vertices*sizeof(Vec3) );
        if( 0 != p_ps3->GetLambdasVDL() )
            m_pSDS->WriteArray( "vec_lambdas",
                                reinterpret_cast<const uint8*>(p_ps3->GetLambdasVDL()),
                                num_vertices*sizeof(Real) );
        if( 0 != p_ps3->GetFlagsVDL() )
            m_pSDS->WriteArray( "vec_flags",
                                reinterpret_cast<const uint8*>(p_ps3->GetFlagsVDL()),
                                num_vertices*sizeof(Flags32) );
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_PathShape2( const PathShape2 *p_ps2, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_Path2 );
        m_pSDS->Write( "version", (uint32)PathShape2::cVersion );
        m_pSDS->Write( "num_v", p_ps2->GetNumV() );
        m_pSDS->Write( "num_e", p_ps2->GetNumE() );
        m_pSDS->Write( "is_closed", p_ps2->IsClosed() );
        //Write internal type arrays as byte blocks
        //m_pSDS->WriteArray( "vec_points", p_ps2->GetVecPoints(), p_ps2->GetNumV() );
        m_pSDS->WriteArray( "vec_points",
                            reinterpret_cast<const uint8*>(p_ps2->GetVecPoints()),
                            p_ps2->GetNumV()*sizeof(Vec2) );
        m_pSDS->WriteArray( "vec_e",
                            reinterpret_cast<const uint8*>(p_ps2->GetVecE()),
                            p_ps2->GetNumAllocE()*sizeof(PathShape2::edge_type) );
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_MeshSolidShape2( const MeshSolidShape2 *p_mss2, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_MeshSolid2 );
        m_pSDS->Write( "version", (uint32)MeshSolidShape2::cVersion );
        m_pSDS->Write( "num_v", p_mss2->GetNumV() );
        m_pSDS->Write( "num_p", p_mss2->GetNumP() );
        m_pSDS->Write( "num_he", p_mss2->GetNumHE() );
        m_pSDS->Write( "num_boundary_p", p_mss2->GetNumBoundaryP() );
        m_pSDS->Write( "num_boundary_he", p_mss2->GetNumBoundaryHE() );
        m_pSDS->Write( "num_l", p_mss2->GetNumL() );
        // Write ALL type arrays as byte blocks
        m_pSDS->WriteArray( "vec_points",
                            reinterpret_cast<const uint8*>(p_mss2->GetVecPoints()),
                            p_mss2->GetNumV()*sizeof(Vec2) );
        m_pSDS->WriteArray( "vec_v",
                            reinterpret_cast<const uint8*>(p_mss2->GetVecV()),
                            p_mss2->GetNumAllocV()*sizeof(MeshSolidShape2::vertex_type) );
        m_pSDS->WriteArray( "vec_p",
                            reinterpret_cast<const uint8*>(p_mss2->GetVecP()),
                            p_mss2->GetNumAllocP()*sizeof(MeshSolidShape2::polygon_type) );
        m_pSDS->WriteArray( "vec_he",
                            reinterpret_cast<const uint8*>(p_mss2->GetVecHE()),
                            p_mss2->GetNumAllocHE()*sizeof(MeshSolidShape2::half_edge_type) );
        m_pSDS->WriteArray( "vec_l",
                            reinterpret_cast<const uint8*>(p_mss2->GetVecL()),
                            p_mss2->GetNumL()*sizeof(MeshSolidShape2::polygon_layer_type) );
        // Optional annotations
        if( p_mss2->GetDCR() != 0 )
        {
            const DCR_MeshSolidShape2* pDCR( p_mss2->GetDCR() );
            m_pSDS->BeginComplex( "dcr", eType_Property_Object );
            {
                m_pSDS->Write( "version", (uint32)DCR_MeshSolidShape2::cVersion );
                m_pSDS->Write( "num_e", pDCR->m_NumElements );
                m_pSDS->Write( "num_p", pDCR->m_NumPatches );
                m_pSDS->Write( "num_s", pDCR->m_NumSegments );
                m_pSDS->Write( "num_v", pDCR->m_NumVertices );
                m_pSDS->WriteArray( "vec_e",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecE),
                                    pDCR->m_NumElements*sizeof(DCR_MeshSolidShape2::Element) );
                m_pSDS->WriteArray( "vec_p",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecP),
                                    pDCR->m_NumPatches*sizeof(DCR_MeshSolidShape2::Patch) );
                m_pSDS->WriteArray( "vec_s",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecS),
                                    pDCR->m_NumSegments*sizeof(DCR_MeshSolidShape2::Segment) );
                m_pSDS->WriteArray( "vec_v",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecV),
                                    pDCR->m_NumVertices*sizeof(Vec2) );
            }
            m_pSDS->EndComplex();
        }
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_TetSolidShape3( const TetSolidShape3 *p_tss3, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_TetSolid3 );
        m_pSDS->Write( "version", (uint32)TetSolidShape3::cVersion );
        m_pSDS->Write( "num_v", p_tss3->GetNumV() );
        m_pSDS->Write( "num_t", p_tss3->GetNumT() );
        m_pSDS->Write( "num_bf", p_tss3->GetNumBF() );
        m_pSDS->Write( "num_l", p_tss3->GetNumL() );
        // Write ALL type arrays as byte blocks
        m_pSDS->WriteArray( "vec_points",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecPoints()),
                            p_tss3->GetNumV()*sizeof(Vec3) );
        m_pSDS->WriteArray( "vec_v",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecV()),
                            p_tss3->GetNumV()*sizeof(TetSolidShape3::vertex_type) );
        m_pSDS->WriteArray( "vec_t",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecT()),
                            p_tss3->GetNumT()*sizeof(TetSolidShape3::tetrahedron_type) );
        m_pSDS->WriteArray( "vec_bf",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecBF()),
                            p_tss3->GetNumBF()*sizeof(TetSolidShape3::boundary_face_type) );
        m_pSDS->WriteArray( "vec_l",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecL()),
                            p_tss3->GetNumL()*sizeof(TetSolidShape3::tetrahedron_layer_type) );
        // Optional annotations
        if( p_tss3->GetDCR() != 0 )
        {
            const DCR_TetSolidShape3* pDCR( p_tss3->GetDCR() );
            m_pSDS->BeginComplex( "dcr", eType_Property_Object );
            {
                m_pSDS->Write( "version", (uint32)DCR_TetSolidShape3::cVersion );
                m_pSDS->Write( "num_e", pDCR->m_NumElements );
                m_pSDS->Write( "num_p", pDCR->m_NumPatches );
                m_pSDS->Write( "num_t", pDCR->m_NumTriangles );
                m_pSDS->Write( "num_v", pDCR->m_NumVertices );
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
                m_pSDS->Write( "num_npid", pDCR->m_NumNPID );
#endif
                m_pSDS->WriteArray( "vec_e",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecE),
                                    pDCR->m_NumElements*sizeof(DCR_TetSolidShape3::Element) );
                m_pSDS->WriteArray( "vec_p",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecP),
                                    pDCR->m_NumPatches*sizeof(DCR_TetSolidShape3::Patch) );
                m_pSDS->WriteArray( "vec_t",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecT),
                                    pDCR->m_NumTriangles*sizeof(DCR_TetSolidShape3::Triangle) );
                m_pSDS->WriteArray( "vec_v",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecV),
                                    pDCR->m_NumVertices*sizeof(Vec3) );
#ifdef __ENABLE_TETSS_DCR_PATCH_TOPOLOGY_VAR_LENGTH
                m_pSDS->WriteArray( "vec_npid",
                                    reinterpret_cast<const uint8*>(pDCR->m_vecNPID),
                                    pDCR->m_NumNPID*sizeof(uint32) );
#endif
            }
            m_pSDS->EndComplex();
        }
    }
    return m_pSDS->EndComplex().GetId();
}

ShapeID ShapeLibrary::Register_TriSurfaceShape3( const TriSurfaceShape3 *p_tss3, const char* name )
{
    if( 0 == name ) m_pSDS->BeginComplex( ++m_LastShapeId, eType_ShapeDef );
    else m_pSDS->BeginComplex( name, eType_ShapeDef );
    {
        m_pSDS->Write( "shape_type", (uint32)eShape_TriSurface3 );
        m_pSDS->Write( "version", (uint32)TriSurfaceShape3::cVersion );
        m_pSDS->Write( "flags", (int32)p_tss3->GetFlags() ); //\todo Support GFlags<T> serialization!?
        m_pSDS->Write( "num_v", p_tss3->GetNumV() );
        m_pSDS->Write( "num_t", p_tss3->GetNumT() );
        // Write ALL type arrays as byte blocks
        m_pSDS->WriteArray( "vec_points",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecPoints()),
                            p_tss3->GetNumV()*sizeof(Vec3) );
        m_pSDS->WriteArray( "vec_v",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecV()),
                            p_tss3->GetNumV()*sizeof(TriSurfaceShape3::vertex_type) );
        m_pSDS->WriteArray( "vec_t",
                            reinterpret_cast<const uint8*>(p_tss3->GetVecT()),
                            p_tss3->GetNumT()*sizeof(TriSurfaceShape3::triangle_type) );
    }
    return m_pSDS->EndComplex().GetId();
}

//...

ShapeDef ShapeLibrary::Lookup( ShapeID sid ) const
{
    return m_pSDS->Find(sid);
}

ShapeID ShapeLibrary::GetShapeIdByName( const char *name ) const
{
    ShapeDef shape_def = m_pSDS->Find(name);
    if( shape_def.IsValid() ) return shape_def.GetId();
    else return cInvalidShapeId;
}

ShapeDef ShapeLibrary::GetFirstShapeDef() const
{
    return m_pSDS->Begin();
}

void ShapeLibrary::Clear()
{
    m_pSDS->Clear();
    m_LastShapeId = 0;
}

bool ShapeLibrary::Save( const char *file_name, bool b_binary ) const
{
    if( b_binary ) return m_pSDS->SaveBin( file_name );
    else return m_pSDS->SaveTxt( file_name );
}

/* We allow realloc during Load() because it's *NOT* incremental, but
   afterwards we restore the original realloc policy, set in Reserve(),
   that allows Identifier realloc but disallows Data realloc to avoid
   shape data pointer invalidation (it's happened!)

   \note Binary load does not require incremental realloc because it
         uses exact memory size saved on the .bin resource.
*/
bool ShapeLibrary::Load( const char *file_name, bool b_binary )
{
    Flags32 rf( m_pSDS->GetReallocFlags() );
    Clear();
    m_pSDS->SetReallocFlags( util::ItemStream::eRealloc_Everything );
    bool bResult( b_binary ? m_pSDS->LoadBin( file_name ) : m_pSDS->LoadTxt( file_name ) );
    m_pSDS->SetReallocFlags( rf );
    return bResult;
}

} //namespace geo
