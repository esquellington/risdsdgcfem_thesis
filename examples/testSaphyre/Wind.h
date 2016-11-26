#ifndef TEST_SAPHYRE_WIND_H
#define TEST_SAPHYRE_WIND_H

#include "Scene.h"
#include <Saphyre2/bs/RayCastQuery.h>

#define __ENABLE_WIND_FRAGMENTS_TO_NODES

//!\sa http://en.wikipedia.org/wiki/Dynamic_pressure
inline float ComputeWindPressureFromSpeed( float speed )
{
    const float cAirDensity( 1.25f / 1000.0f ); //1.25kg/m^3 (dry air at sea level)
    return 0.5f*cAirDensity*mal::Sq(speed); // 1/2 \rho v^2
}

/* Unidirectional wind with constant speed */
class Wind2D
{
public:
    inline Wind2D( Scene& scene )
    : m_rScene(scene)
    , m_Pos(0,0), m_Dir(1,0), m_Size(1), m_Length(100)
    , m_Speed(10) //m/s
    , m_Pressure( ComputeWindPressureFromSpeed(m_Speed) ) //N/m^2
    , m_Resolution(10)
    , m_bufferG(0), m_bufferW(0)
    , m_pRCQ(0)
        {
            m_pRCQ = m_rScene.GetUniverse()->CreateRayCastQuery2D();
        }
    inline ~Wind2D()
        {
            if( 0 != m_bufferG ) delete[] m_bufferG;
            if( 0 != m_bufferW ) delete[] m_bufferW;
        }

    void Init( const Vec2f& pos, const Vec2f& dir, float size, float length,
               float speed,
               uint32 resolution )
        {
            m_Pos = pos;
            m_Dir = dir;
            m_Size = size;
            m_Length = length;
            m_Speed = speed;
            m_Pressure = ComputeWindPressureFromSpeed( speed );
            m_Resolution = resolution;

            // Setup transforms, by now Identity rot and dir along X axis
            m_TransformW2G = Transform2f( m_Pos, mal::GMat2x2_From_Columns(m_Dir,mal::PerpendicularCW(m_Dir)) );
            m_TransformG2W = mal::Inverse( m_TransformW2G );

            //---- Create buffers
            if( 0 != m_bufferG ) delete[] m_bufferG;
            if( 0 != m_bufferW ) delete[] m_bufferW;
            // Create G-buffer
            m_bufferG = new FragmentG[m_Resolution];
            memset( &m_bufferG[0], 0, m_Resolution*sizeof(FragmentG) );
            // Create W-buffer
            m_bufferW = new FragmentW[m_Resolution];
            memset( &m_bufferW[0], 0, m_Resolution*sizeof(FragmentW) );
        }

    void Update( double dt )
        {
            //---- (1) Generate G-buffer: Throw rays, gather new g-fragments
            uint32 num_hits(0);
            for( uint32 it_gf=0; it_gf<m_Resolution; it_gf++ )
            {
                FragmentG& gf( m_bufferG[it_gf] );
                m_pRCQ->Run( m_TransformW2G * GetFragmentLocalPos(it_gf),
                             m_Dir,
                             Intervalf(0,m_Length),
                             0.0f,
                             S2::eRCQF_Default );
                S2::RayCastQuery2D::Hit2D h2d;
                if( m_pRCQ->Closest(h2d)
                    && mal::Dot( h2d.m_RayHit.m_Normal, m_Dir ) < 0 ) //backface culling
                {
                    gf.m_pEntity = h2d.m_pEntity;
                    gf.m_FeatureId = h2d.m_RayHit.m_FeatureId;
                    gf.m_BarycentricCoords = mal::GRange<0,1>(h2d.m_RayHit.m_Extra_BarycentricCoords);
                    gf.m_Depth = h2d.m_RayHit.m_Interval.Min();
                    gf.m_Pos = m_TransformG2W * h2d.m_RayHit.m_Point;
                    gf.m_Normal = m_TransformG2W.Rot() * h2d.m_RayHit.m_Normal;
                    //GEO_LOG_WARNING( "Wind hits entity %s", gf.m_pEntity->GetName() );
                    num_hits++;
                }
                else
                {
                    gf.m_pEntity = 0;
                    gf.m_FeatureId = geo::feature_id();//invalid
                    gf.m_BarycentricCoords = Vec2f::Zero();
                    gf.m_Depth = 0;
                    gf.m_Pos = Vec2f::Zero();
                    gf.m_Normal = Vec2f::Zero();
                }
            }
            //APP_LOG( "Num FragmentG = %d", num_hits );
            //---- (2) \todo Generate W-buffer simulating wind \todo Maybe using previous W-buffer
            float fragment_area( m_Size / m_Resolution );
            for( uint32 it_wf=0; it_wf<m_Resolution; it_wf++ )
            {
                FragmentW& wf( m_bufferW[it_wf] );
                FragmentG& gf( m_bufferG[it_wf] );
                wf.m_Fn = fragment_area * m_Pressure * mal::Dot(Vec2f(1,0),gf.m_Normal) * gf.m_Normal; //wind dir local = (1,0), normal force
                wf.m_Ft = Vec2f::Zero();//\todo Compute tangential force => friction
            }

            //---- (3) Apply W-buffer forces to simulation nodes
#ifdef __ENABLE_WIND_FRAGMENTS_TO_NODES
            //-- (a) Per-fragment apply on nodes
            for( uint32 it_wf=0; it_wf<m_Resolution; it_wf++ )
            {
                FragmentG& gf( m_bufferG[it_wf] );
                if( 0 != gf.m_pEntity ) // occluded
                {
                    FragmentW& wf( m_bufferW[it_wf] );
                    Vec2f force( m_TransformW2G.Rot() * (wf.m_Fn + wf.m_Ft) );
                    switch( gf.m_pEntity->GetType() )
                    {
                    case S2::eEntity_Solid2D:
                        {
                            S2::Solid2D* pS2D = static_cast<S2::Solid2D*>(gf.m_pEntity);
                            const geo::MeshSolidShape2* pMSS2 = pS2D->GetMeshGO()->GetShape();
                            pS2D->Lock();
                            {
                                //apply force to feature_id
                                if( gf.m_FeatureId.IsVertex() )
                                {
                                    // Node
                                    geo::feature_index_type vid = gf.m_FeatureId.AsVertex();
                                    pS2D->ApplyForce( vid, force );
                                }
                                else if( gf.m_FeatureId.IsSegment() )
                                {
                                    // Boundary edge
                                    geo::feature_index_type heid = gf.m_FeatureId.AsSegment();
                                    geo::feature_index_type vid0 = pMSS2->HE_OriginVID( heid );
                                    geo::feature_index_type vid1 = pMSS2->HE_FinalVID( heid );
                                    pS2D->ApplyForce( vid0, gf.m_BarycentricCoords[0]*force );
                                    pS2D->ApplyForce( vid1, gf.m_BarycentricCoords[1]*force );
                                }
                            }
                            pS2D->Unlock();
                        }
                        break;
                    default: break;
                    }
                }
            }
#else
            //-- (b) Per-node read from W-buffer
            for( Scene::EntityIterator it=m_rScene.GetEntityIterator(); it.IsValid(); ++it )
            {
                S2::IEntity* pEntity( *it );
                switch( pEntity->GetType() )
                {
                case S2::eEntity_Solid2D:
                    {
                        S2::Solid2D* pS2D = static_cast<S2::Solid2D*>(pEntity);
                        const geo::MeshSolidShape2* pMSS2 = pS2D->GetMeshGO()->GetShape();
                        if( true ) //\todo if( overlaps wind region ) Use Wind BV and MeshGO BV
                        {
                            pS2D->Lock();
                            {
                                // For each boundary polygon \todo ACTUALLY, only for the EXTERNAL ONE...
                                for( unsigned int it_bp=0; it_bp < pMSS2->GetNumBoundaryP(); it_bp++ )
                                {
                                    // For each edge
                                    unsigned int it_he( pMSS2->BP_FirstHEID(it_bp) );
                                    do
                                    {
                                        unsigned int nid( pMSS2->HE_OriginVID(it_he) );
                                        Vec2f pos( pS2D->GetMeshGO()->GetVecSDOF()[nid] ); //\todo Use edge midpoint instead?
                                        Vec2f pos_w( m_TransformG2W*pos );
                                        if( pos_w.x() > 0
                                            && pos_w.x() < m_Length
                                            && pos_w.y() > -0.5f*m_Size
                                            && pos_w.y() < 0.5f*m_Size ) //\todo if( overlaps wind region )
                                        {
                                            uint32 gf_index = GetFragmentIndexAtLocalPos(pos_w);
                                            FragmentG& gf( m_bufferG[ gf_index ] );
                                            FragmentW& wf( m_bufferW[ gf_index ] );
                                            //TEMP: We compare depth with epsilon to avoid self-occlusion...
                                            if( pos_w.x() <= gf.m_Depth + 0.1f ) //\todo if( not occluded in G-buffer )
                                            {
                                                Vec2f Fn,Ft,p;
                                                if( GetFragmentW(gf_index,p,Fn,Ft) )
                                                    pS2D->ApplyForce( nid, Fn + Ft ); //\todo Apply W-buffer force
                                            }
                                        }
                                        it_he = pMSS2->HE_Next(it_he);
                                    } while ( it_he != pMSS2->BP_FirstHEID(it_bp) );
                                }

                                /*OLD: for each node (slow, includes internal nodes!)
                                for( uint32 it_node=0;
                                     it_node < pS2D->GetMeshGO()->GetShape()->GetNumV();
                                     it_node++ )
                                {
                                    Vec2f pos( pS2D->GetMeshGO()->GetVecSDOF()[it_node] );
                                    Vec2f pos_w( m_TransformG2W*pos );
                                    if( pos_w.x() > 0
                                        && pos_w.x() < m_Length
                                        && pos_w.y() > -0.5f*m_Size
                                        && pos_w.y() < 0.5f*m_Size ) //\todo if( overlaps wind region )
                                    {
                                        uint32 gf_index = GetFragmentIndexAtLocalPos(pos_w);
                                        FragmentG& gf( m_bufferG[ gf_index ] );
                                        FragmentW& wf( m_bufferW[ gf_index ] );
                                        //TEMP: We compare depth with epsilon to avoid self-occlusion...
                                        if( pos_w.x() <= gf.m_Depth + 0.1f ) //\todo if( not occluded in G-buffer )
                                        {
                                            Vec2f Fn,Ft,p;
                                            if( GetFragmentW(gf_index,p,Fn,Ft) )
                                                pS2D->ApplyForce( it_node, Fn + Ft ); //\todo Apply W-buffer force
                                        }
                                    }
                                }
                                */
                            }
                            pS2D->Unlock();
                        }
                    }
                    break;
                default: break;
                }
            }
#endif //__ENABLE_WIND_FRAGMENTS_TO_NODES
        }

    inline uint32 GetResolution() const { return m_Resolution; }
    inline Vec2f GetFragmentLocalPos( int gf_id ) const
        {
            return Vec2f( 0, (float(gf_id)/m_Resolution - 0.5f) * m_Size );
        }
    inline int GetFragmentIndexAtLocalPos( const Vec2f& pos_w ) const
        {
            return mal::Clamp<float>( ( (pos_w.y() + 0.5f*m_Size) / m_Size ) * m_Resolution, 0, m_Resolution-1 ); //\todo may be inaccurate
        }

    inline bool GetFragmentG( int gf_id, Vec2f& pos, Vec2f& normal ) const
        {
            if( 0 != m_bufferG[gf_id].m_pEntity )
            {
                pos = m_TransformW2G * m_bufferG[gf_id].m_Pos;
                normal = m_TransformW2G.Rot() * m_bufferG[gf_id].m_Normal;
                return true;
            }
            else
                return false;
        }

    inline bool GetFragmentW( int wf_id, Vec2f& pos, Vec2f& fn, Vec2f& ft ) const
        {
            if( 0 != m_bufferG[wf_id].m_pEntity )
            {
                pos = m_TransformW2G * m_bufferG[wf_id].m_Pos;
                fn = m_TransformW2G.Rot() * m_bufferW[wf_id].m_Fn;
                ft = m_TransformW2G.Rot() * m_bufferW[wf_id].m_Ft;
                return true;
            }
            else
                return false;
        }


public: //TEMP: For AppRenderer by now...

    Scene& m_rScene;

    Vec2f m_Pos;
    Vec2f m_Dir;
    float m_Size;
    float m_Length;
    float m_Speed;
    float m_Pressure;
    uint32 m_Resolution;

    Transform2f m_TransformW2G; //Wind to Global transform
    Transform2f m_TransformG2W; //Global to Wind transform

    //\note All fragment magnitudes are stored in wind-relative coords
    struct FragmentG
    {
        S2::ISyncEntity* m_pEntity; //Entity receiving the wind
        geo::feature_id m_FeatureId;
        Vec2f m_BarycentricCoords; //u,v
        float m_Depth;
        Vec2f m_Pos; //\todo Can be derived from m_FragmentPos (XY) and m_Depth (Z)
        Vec2f m_Normal;
    };
    FragmentG* m_bufferG;

    //\note All fragment magnitudes are stored in wind-relative coords
    struct FragmentW
    {
        Vec2f m_Fn;
        Vec2f m_Ft;

        Vec2f m_vecAccForce[2]; //AccPressure/FragmentSurface distrib accross 2 vtx

        //\todo unnecessary by now uint32 m_vecVID[2];
    };
    FragmentW* m_bufferW;

    S2::RayCastQuery2D* m_pRCQ;
};

#endif //TEST_SAPHYRE_WIND_H
