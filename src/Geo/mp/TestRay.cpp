#include "TestRay.h"
#include <Geo/shape/shape.h>
#include <Geo/np/RayCast.h>

namespace geo {
namespace mp {

bool TestRay( const IObject2* p_obj, const np::Ray2& ray, np::RayHit2& rh, np::RayCache2* p_rc )
{
    GEO_ASSERT( 2 == p_obj->GetDimension() );
    const IObject2* pObj( static_cast<const IObject2*>(p_obj) );
    const IShape* pShape( p_obj->GetShapeInterface() );
    switch( pShape->GetType() )
    {
    case eShape_Plane2:
        {
            const PlaneShape2* pPlane = static_cast<const PlaneShape2*>(pShape);
            //\todo if( ray.m_Thickness > 0 )
            if( pPlane->IsHalfSpace() )
                return np::GRayCast_HalfSpace<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                                  pPlane->GetNormal(), pPlane->GetCoeffD(),
                                                  rh, p_rc );
            else
                return np::GRayCast_Plane<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                              pPlane->GetNormal(), pPlane->GetCoeffD(),
                                              rh, p_rc );
        }
        break;
    case eShape_Sphere2:
        {
            const SphereShape2* pSphere = static_cast<const SphereShape2*>(pShape);
            /* \todo
            if( ray.m_Thickness > 0 )
                return np::GThickRayCast_Sphere<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval, ray.m_Thickness,
                                                    pObj->GetTransform().m_Pos, pSphere->GetRadius(),
                                                    rh, p_rc );
            else
            */
                return np::GRayCast_Sphere<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                               pObj->GetTransform().m_Pos, pSphere->GetRadius(),
                                               rh, p_rc );
        }
        break;
    case eShape_Capsule2:
        {
            const CapsuleShape2* pCapsule = static_cast<const CapsuleShape2*>(pShape);
            Transform2 tr2( pObj->GetTransform() );
            return np::GRayCast_Capsule<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                            tr2.m_Pos, mal::Column(1,tr2.m_Rot), pCapsule->GetRadius(), pCapsule->GetHalfHeight(),
                                            rh, p_rc );
        }
        break;
    case eShape_Box2:
        {
            const BoxShape2* pBox = static_cast<const BoxShape2*>(pShape);
            // \todo if( ray.m_Thickness > 0 )
            return np::GRayCast_Box<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                        pObj->GetTransform(), pBox->GetHalfSizes(),
                                        rh, p_rc );
        }
        break;
        /*
    case eShape_SphereSet2:
        {
            const BoxShape2* pBox = static_cast<const BoxShape2*>(pShape);
            // \todo if( ray.m_Thickness > 0 )
            return np::GRayCast_Box<2>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                        pObj->GetTransform(), pBox->GetHalfSizes(),
                                        rh, p_rc );
        }
        break;
        */
    case eShape_Polygonal2: //\todo Works but it's slow, should use some kind of BVH
        {
            const PolygonalShape2* pPS2( static_cast<const PolygonalShape2*>(pShape) );
            const GObjectSDOF<2,Vec2>* pPSGO( static_cast<const GObjectSDOF<2,Vec2>* >(pObj) );
            Transform2 tr2( static_cast<const IObject2*>(pObj)->GetTransform() );
            int num_edges( pPS2->IsClosed() ? pPS2->GetNumVertices() : pPS2->GetNumVertices()-1 );
            const Vec2* actual_sdof( pPSGO->GetVecSDOF() );
            unsigned int num_hits(0);
            for( int it_edge=0; it_edge < num_edges; it_edge++ )
            {
                Vec2 p0( tr2 * pPS2->V_Pos( it_edge, actual_sdof ) );
                int next_vid( (it_edge+1) % pPS2->GetNumVertices() );
                Vec2 p1( tr2 * pPS2->V_Pos( next_vid, actual_sdof ) );
                np::RayHit2 tmp_rh;
                if( np::RayCast_Segment2_SingleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, p0, p1, tmp_rh, p_rc ) )
                {
                    if( num_hits == 0 || tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
                    {
                        rh = tmp_rh;
                        //Compute rh.m_FeatureId from barycentric coords computed by RayCast_Segment2
                        if( rh.m_Extra_BarycentricCoords[1] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, (feature_index_type)it_edge );
                        else if( rh.m_Extra_BarycentricCoords[0] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, (feature_index_type)next_vid );
                        else rh.m_FeatureId.Set( eFT_Segment, (feature_index_type)it_edge );
                        // \todo OK: Normal points outwards if CCW
                    }
                    num_hits++;
                }
            }
            return num_hits > 0;
        }
        break;
    case eShape_Path2: //\todo Works but it's slow, should use some kind of BVH
        {
            const PathShape2* pPS2( static_cast<const PathShape2*>(pShape) );
            const GObjectSDOF<2,Vec2>* pPSGO( static_cast<const GObjectSDOF<2,Vec2>* >(pObj) );
            Transform2 tr2( static_cast<const IObject2*>(pObj)->GetTransform() );
            const Vec2* actual_sdof( pPSGO->GetVecSDOF() );
            unsigned int num_hits(0);
            for( unsigned int it_edge=0; it_edge < pPS2->GetNumE(); it_edge++ )
            {
                uint32 vid0( pPS2->E_OriginVID(it_edge) );
                uint32 vid1( pPS2->E_FinalVID(it_edge) );
                Vec2 p0( tr2 * pPS2->V_Pos( vid0, actual_sdof ) );
                Vec2 p1( tr2 * pPS2->V_Pos( vid1, actual_sdof ) );
                np::RayHit2 tmp_rh;
                // Test according to edge curve type
                switch( pPS2->E_CurveType(it_edge) )
                {
                case PathShape2::edge_type::eCT_Line:
                    {
                        if( np::RayCast_Segment2_SingleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, p0, p1, tmp_rh, p_rc ) )
                        {
                            if( num_hits == 0 || tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
                            {
                                rh = tmp_rh;
                                //Compute rh.m_FeatureId from barycentric coords computed by RayCast_Segment2
                                if( rh.m_Extra_BarycentricCoords[1] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, (feature_index_type)vid0 );
                                else if( rh.m_Extra_BarycentricCoords[0] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, (feature_index_type)vid1 );
                                else rh.m_FeatureId.Set( eFT_Segment, (feature_index_type)it_edge );
                                // \todo OK: Normal points outwards if CCW
                            }
                            num_hits++;
                        }
                    }
                    break;
                case PathShape2::edge_type::eCT_Bezier3:
                    {
                        // Bezier curve
                        Vec2 bezier3_a( tr2 * pPS2->E_Data(it_edge).m_ParamA ); //\todo If params are also SDOF, use actual_sdof here!!
                        Vec2 bezier3_b( tr2 * pPS2->E_Data(it_edge).m_ParamB );
                        Vec2 b0( p0 );
                        const unsigned int cNumSamples(10);
                        for( unsigned int i=1; i<cNumSamples; i++ )
                        {
                            Vec2 b1 = Eval_Bezier3( p0, p1, bezier3_a, bezier3_b, Real(i) / (cNumSamples-1) );
                            if( np::RayCast_Segment2_SingleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, b0, b1, tmp_rh, p_rc ) )
                            {
                                if( num_hits == 0 || tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
                                {
                                    rh = tmp_rh;
                                    //Compute rh.m_FeatureId from barycentric coords computed by RayCast_Segment2
                                    if( i == 0 && rh.m_Extra_BarycentricCoords[1] < gs_EpsilonVertexLambda )
                                        rh.m_FeatureId.Set( eFT_Vertex, (feature_index_type)vid0 );
                                    else if( i == cNumSamples-1 && rh.m_Extra_BarycentricCoords[0] < gs_EpsilonVertexLambda )
                                        rh.m_FeatureId.Set( eFT_Vertex, (feature_index_type)vid1 );
                                    else
                                        rh.m_FeatureId.Set( eFT_Segment, (feature_index_type)it_edge );
                                    // \todo OK: Normal points outwards if CCW
                                }
                                num_hits++;
                            }
                            b0 = b1;
                        }
                    }
                    break;
                default: break;
                }
            }
            return num_hits > 0;
        }
        break;
    case eShape_MeshSolid2: //\todo Works but it's slow, should use some kind of BVH and optional DCR
        {
            //IMPORTANT: MSS boundary is CW, not CCW, therefore we'll pass reversed p1->p0 segments to RayCast_Segment2_SingleSided()
            const MeshSolidShape2* pMSS2( static_cast<const MeshSolidShape2*>(pShape) );
            const GObjectSDOF<2,Vec2>* pPSGO( static_cast<const GObjectSDOF<2,Vec2>* >(pObj) );
            Transform2 tr2( static_cast<const IObject2*>(pObj)->GetTransform() );
            const Vec2* actual_sdof( pPSGO->GetVecSDOF() );
            unsigned int num_hits(0);
            np::RayHit2 tmp_rh;
            for( unsigned int it_bp=0; it_bp<pMSS2->GetNumBoundaryP(); it_bp++ )
            {
                unsigned int it_he( pMSS2->BP_FirstHEID(it_bp) );
                Vec2 p0( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) );
                do
                {
                    it_he = pMSS2->HE_Next(it_he);
                    Vec2 p1( tr2 * pMSS2->V_Pos( pMSS2->HE_OriginVID(it_he), actual_sdof ) );
                    if( np::RayCast_Segment2_SingleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, p1, p0, tmp_rh, p_rc ) ) //p1->p0 to force CCW
                    {
                        if( num_hits == 0 || tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
                        {
                            rh = tmp_rh;
                            //Compute rh.m_FeatureId from barycentric coords computed by RayCast_Segment2 \note p0,p1 reversed!
                            if( rh.m_Extra_BarycentricCoords[0] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, pMSS2->HE_OriginVID(it_he) );
                            else if( rh.m_Extra_BarycentricCoords[1] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, pMSS2->HE_FinalVID(it_he) );
                            else rh.m_FeatureId.Set( eFT_Segment, (feature_index_type)it_he );
                        }
                        num_hits++;
                    }
                    p0 = p1;
                } while ( it_he != pMSS2->BP_FirstHEID(it_bp) );
            }
            return num_hits > 0;
        }
        break;
    default:
        GEO_LOG_WARNING( "geo::mp::TestRay: Shape type %d not yet implemented", pShape->GetType() );
        return false;
        break;
    }
}

bool TestRay( const IObject3* p_obj, const np::Ray3& ray, np::RayHit3& rh, np::RayCache3* p_rc )
{
    GEO_ASSERT( 3 == p_obj->GetDimension() );
    const IObject3* pObj( static_cast<const IObject3*>(p_obj) );
    const IShape* pShape( p_obj->GetShapeInterface() );
    switch( pShape->GetType() )
    {
    case eShape_Plane3:
        {
            const PlaneShape3* pPlane = static_cast<const PlaneShape3*>(pShape);
            //\todo if( ray.m_Thickness > 0 )
            if( pPlane->IsHalfSpace() )
                return np::GRayCast_HalfSpace<3>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                                  pPlane->GetNormal(), pPlane->GetCoeffD(),
                                                  rh, p_rc );
            else
                return np::GRayCast_Plane<3>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                              pPlane->GetNormal(), pPlane->GetCoeffD(),
                                              rh, p_rc );
        }
        break;
    case eShape_Sphere3:
        {
            const SphereShape3* pSphere = static_cast<const SphereShape3*>(pShape);
            /* \todo
            if( ray.m_Thickness > 0 )
                return np::GThickRayCast_Sphere<3>( ray.m_Pos, ray.m_Dir, ray.m_Interval, ray.m_Thickness,
                                                    pObj->GetTransform().m_Pos, pSphere->GetRadius(),
                                                    rh, p_rc );
            else
            */
            return np::GRayCast_Sphere<3>( ray.m_Pos, ray.m_Dir, ray.m_Interval,
                                           pObj->GetTransform().m_Pos, pSphere->GetRadius(),
                                           rh, p_rc );
        }
        break;
    case eShape_TetSolid3: //\todo Works but it's slow, should use some kind of BVH and optional DCR
        {
            //IMPORTANT: MSS boundary is CW, not CCW, therefore we'll pass reversed p1->p0 segments to RayCast_Segment2_SingleSided()
            const TetSolidShape3* pTSS3( static_cast<const TetSolidShape3*>(pShape) );
            const GObjectSDOF<3,Vec3>* pPSGO( static_cast<const GObjectSDOF<3,Vec3>* >(pObj) );
            Transform3 tr3( static_cast<const IObject3*>(pObj)->GetTransform() );
            const Vec3* actual_sdof( pPSGO->GetVecSDOF() );
            unsigned int num_hits(0);
            np::RayHit3 tmp_rh;
            for( unsigned int it_bf=0; it_bf<pTSS3->GetNumBF(); it_bf++ )
            {
                Vec3 p1( tr3 * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,0), actual_sdof ) );
                Vec3 p2( tr3 * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,1), actual_sdof ) );
                Vec3 p3( tr3 * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,2), actual_sdof ) );
                if( np::RayCast_Triangle3_SingleSided( ray.m_Pos, ray.m_Dir, ray.m_Interval, p1, p2, p3, tmp_rh, p_rc ) ) //p1->p0 to force CCW
                {
                    if( num_hits == 0 || tmp_rh.m_Interval.Min() < rh.m_Interval.Min() )
                    {
                        rh = tmp_rh;
                        //\todo Compute rh.m_FeatureId from barycentric coords computed by RayCast_Segment2 \note p0,p1 reversed!
                        // if( rh.m_Extra_BarycentricCoords[0] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, pMSS2->HE_OriginVID(it_he) );
                        // else if( rh.m_Extra_BarycentricCoords[1] < gs_EpsilonVertexLambda ) rh.m_FeatureId.Set( eFT_Vertex, pMSS2->HE_FinalVID(it_he) );
                        // else rh.m_FeatureId.Set( eFT_Segment, (feature_index_type)it_he );
                        rh.m_FeatureId.Set( eFT_Triangle, it_bf );
                    }
                    num_hits++;
                }
            }
            return num_hits > 0;
        }
        break;
    default:
        GEO_LOG_WARNING( "geo::mp::TestRay: Shape type %d not yet implemented", pShape->GetType() );
        return false;
        break;
    }
}

}} //namespace geo::mp
