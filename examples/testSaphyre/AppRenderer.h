#ifndef TEST_SAPHYRE_APP_RENDERER_H
#define TEST_SAPHYRE_APP_RENDERER_H

#include "AppTask.h"
#include "SceneRenderer.h"
#include <Safra/Safra.h>
#include "Params.h"

/*! App-specific Scene Renderer with volatile data (raycast hits, intersections, etc...) */
class AppRenderer: public SceneRenderer
{
public:
    AppRenderer( const AppTask &app_task, sfr::gfx::IRenderer *p_renderer )
    : SceneRenderer( app_task.GetParams(), app_task.GetScene(), p_renderer )
    , m_rAppTask( app_task )
    , m_rParams( app_task.GetParams() )
        {}

    bool Render()
    {
        SceneRenderer::Render();

        if( m_rAppTask.GetDimension() == 2 )
        {
            // Draw RayCast stuff
            if( m_rAppTask.m_bIsEnabledRayCast )
            {
                //DrawPoint2( m_rAppTask.m_DragPoint2D, sfr::gfx::Style(sfr::gfx::Color(1,1,1),10) );
                DrawPoint2( m_rAppTask.m_RayCastPoint, sfr::gfx::Style(sfr::gfx::Color(1,1,0,1),5) );
                DrawSegment2( m_rAppTask.m_RayCastPoint,
                              m_rAppTask.m_RayCastPoint + m_rAppTask.m_RayCastDir,
                              sfr::gfx::Style(sfr::gfx::Color(1,1,0,1),2) );
                for( unsigned int it_rh=0; it_rh < m_rAppTask.m_vecRH2D.size(); it_rh++ )
                {
                    const geo::np::RayHit2 rh( m_rAppTask.m_vecRH2D[it_rh].second );
                    if( !rh.m_Interval.IsEmpty() )
                    {
                        // Point
                        DrawPoint2( rh.m_Point, sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),5) );
                        // Normal
                        DrawSegment2( rh.m_Point, rh.m_Point + rh.m_Normal, sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),1) );
                        // Interval
                        DrawSegment2( rh.m_Point,
                                      sfr::Vec2( m_rAppTask.m_RayCastPoint.x() + rh.m_Interval.Max()*m_rAppTask.m_RayCastDir.x(),
                                                 m_rAppTask.m_RayCastPoint.y() + rh.m_Interval.Max()*m_rAppTask.m_RayCastDir.y() ),
                                      sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0.5,1),1) );
                        // FeatureId
                        char str[1024];
                        sprintf( str, "%s (%d,%d)", m_rAppTask.m_vecRH2D[it_rh].first->GetName(), rh.m_FeatureId.m_Type, rh.m_FeatureId.m_Index );
                        DrawLabel( sfr::Vec3(rh.m_Point.x(),rh.m_Point.y(),0),
                                   str, sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0.5,1),1) );
                    }
                }
            }
            if( m_rAppTask.m_bIsEnabledDrag )
                DrawPoint( m_rAppTask.m_DragNodePoint3D, sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),15) );
            if( m_rAppTask.m_bIsEnabledMove )
                DrawPoint2( m_rAppTask.m_MoveOrRotatePoint2D, sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),5) );

#ifdef __ENABLE_WIND
            if( m_rAppTask.m_bIsEnabledWind )
            {
                const Wind2D& wind( *m_rAppTask.m_pWind );
                // Draw Wind transform
                DrawRefSys2( wind.m_TransformW2G.Pos(), wind.m_TransformW2G.Rot(), 1.0f, sfr::gfx::Style() );
                // Draw Wind area
                DrawBox2( wind.m_TransformW2G.Pos() + 0.5f*wind.m_Dir*wind.m_Length,
                          wind.m_TransformW2G.Rot(),
                          0.5f*Vec2f(wind.m_Length,wind.m_Size),
                          sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,0.5,1),0.5) );
                // Draw G-buffer
                for( uint32 it_gf=0; it_gf<wind.m_Resolution; it_gf++ )
                {
                    Wind2D::FragmentG& gf( wind.m_bufferG[it_gf] );
                    /* Redundant
                       Vec2f fragment_pos( wind.m_TransformW2G * wind.GetFragmentLocalPos(it_gf) );
                       DrawSegment2( fragment_pos,
                       fragment_pos + 0.5f*wind.m_Dir,
                       sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,1),2) );
                    */
                    Vec2f p,n;
                    if( wind.GetFragmentG(it_gf,p,n) )
                    {
                        DrawSegment2( p, p+0.25f*n, sfr::gfx::Style(sfr::gfx::Color(0.75,0.75,1,1),1) );

                        Vec2f fragment_pos( wind.m_TransformW2G * wind.GetFragmentLocalPos(it_gf) );
                        DrawSegment2( fragment_pos,
                                      fragment_pos + gf.m_Depth*wind.m_Dir,
                                      sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,1,0.5),1) );


                        /*TEMP: Disabled, noisy
                        // FeatureId
                        char str[1024];
                        //sprintf( str, "%s (%d,%d)", gf.m_pEntity->GetName(), gf.m_FeatureId.m_Type, gf.m_FeatureId.m_Index );
                        //sprintf( str, "%s (%d,%d) => [%f,%f]", gf.m_pEntity->GetName(), gf.m_FeatureId.m_Type, gf.m_FeatureId.m_Index, gf.m_BarycentricCoords[0], gf.m_BarycentricCoords[1] );
                        //sprintf( str, "%llx (%d,%d)", (machine_uint_type)gf.m_pEntity, gf.m_FeatureId.m_Type, gf.m_FeatureId.m_Index );
                        //sprintf( str, "%d, (%f,%f)", wind.GetFragmentIndexAtLocalPos(gf.m_Pos), gf.m_Pos.x(), gf.m_Pos.y() );
                        DrawLabel( sfr::Vec3(p.x(),p.y(),0),
                        str, sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0.5),1) );
                        */
                    }
                    else //Unoccluded wind
                    {
                        Vec2f fragment_pos( wind.m_TransformW2G * wind.GetFragmentLocalPos(it_gf) );
                        DrawSegment2( fragment_pos,
                                      fragment_pos + wind.m_Length*wind.m_Dir,
                                      sfr::gfx::Style(sfr::gfx::Color(0.25,1,0.25,0.5),1) );

                    }
                }
                /* TEMP: Disabled, noisy
                // Draw F-buffer
                for( uint32 it_wf=0; it_wf<wind.m_Resolution; it_wf++ )
                {
                Vec2f p,Fn,Ft;
                if( wind.GetFragmentW(it_wf,p,Fn,Ft) )
                {
                Vec2f force( Fn + Ft );
                DrawSegment2( p, p+force, sfr::gfx::Style(sfr::gfx::Color(1.0,0.0,0),1) );
                }
                }
                */
            }
#endif

            // Draw Contact stuff
            if( m_rAppTask.m_bIsEnabledContact2D )
            {
                // Selected objects
                if( 0 != m_rAppTask.m_pContactKine2D_1 )
                    DrawPoint2( m_rAppTask.m_pContactKine2D_1->GetPos(), sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),5) );
                if( 0 != m_rAppTask.m_pContactKine2D_2 )
                    DrawPoint2( m_rAppTask.m_pContactKine2D_2->GetPos(), sfr::gfx::Style(sfr::gfx::Color(0,0,1,1),5) );
                // ContactData: should have its own Draw_ContactData2() method somewhere, i guess
                const geo::np::ContactData2 &cd( m_rAppTask.m_ContactData2 );
                if( cd.Size() > 0 )
                {
                    for( unsigned int it_cp=0; it_cp<cd.Size(); it_cp++ )
                    {
                        // Pos
                        DrawPoint2( cd.GetCP(it_cp).m_Pos1, sfr::gfx::Style(sfr::gfx::Color(0,1,1,1),5) );
                        DrawPoint2( cd.GetCP(it_cp).m_Pos2, sfr::gfx::Style(sfr::gfx::Color(1,0,1,1),5) );
                        // Normal
                        DrawSegment2( cd.GetCP(it_cp).m_Pos1,
                                      cd.GetCP(it_cp).m_Pos1 + cd.GetCP(it_cp).m_Normal * cd.GetCP(it_cp).m_Depth,
                                      sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),2) );
                    }
                }
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                for( unsigned int it_ip=0; it_ip<cd.m_VD.m_vecIP.size(); it_ip++ )
                    DrawPoint2( 0.5f*(cd.m_VD.m_vecIP[it_ip].first + cd.m_VD.m_vecIP[it_ip].second), sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),6) );
                for( unsigned int it_cp=0; it_cp<cd.m_VD.m_vecCP.size(); it_cp++ )
                {
                    DrawPoint2( cd.m_VD.m_vecCP[it_cp].first, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0,1),4) );
                    DrawPoint2( cd.m_VD.m_vecCP[it_cp].second, sfr::gfx::Style(sfr::gfx::Color(1,0.5,0,1),4) );
                }
                for( unsigned int it_np=0; it_np<cd.m_VD.m_vecNP.size(); it_np++ )
                {
                    DrawPoint2( cd.m_VD.m_vecNP[it_np].first, sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),4) );
                    DrawPoint2( cd.m_VD.m_vecNP[it_np].second, sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),4) );
                }
                for( unsigned int it_fp=0; it_fp<cd.m_VD.m_vecFP.size(); it_fp++ )
                {
                    DrawPoint2( cd.m_VD.m_vecFP[it_fp].first, sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,1,1),2) );
                    DrawPoint2( cd.m_VD.m_vecFP[it_fp].second, sfr::gfx::Style(sfr::gfx::Color(0.25,0.25,1,1),2) );
                }
                for( unsigned int it_rnp=0; it_rnp<cd.m_VD.m_vecRNP.size(); it_rnp++ )
                {
                    DrawPoint2( cd.m_VD.m_vecRNP[it_rnp].first, sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0.5,1),2) );
                    DrawPoint2( cd.m_VD.m_vecRNP[it_rnp].second, sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),2) );
                }
                for( unsigned int it_pca=0; it_pca<cd.m_VD.m_vecPCA.size(); it_pca++ )
                    DrawSegment2( cd.m_VD.m_vecPCA[it_pca].first,
                                  cd.m_VD.m_vecPCA[it_pca].first + cd.m_VD.m_vecPCA[it_pca].second,
                                  sfr::gfx::Style(sfr::gfx::Color(1,1,0,1),2) );
                //IM
                for( unsigned int it_im_p=0; it_im_p<cd.m_VD.m_vecIM_XPoint.size(); it_im_p++ )
                {
                    DrawPoint2( cd.m_VD.m_vecIM_XPoint[it_im_p],
                                sfr::gfx::Style( sfr::gfx::Color(1,0,1,1), 8 ) );
                    char str[128];
                    sprintf( str, "( %f, %f )", cd.m_VD.m_vecIM_XLength[it_im_p].first, cd.m_VD.m_vecIM_XLength[it_im_p].second );
                    DrawLabel2( cd.m_VD.m_vecIM_XPoint[it_im_p], str,
                                sfr::gfx::Style(sfr::gfx::Color(0.5,0.5,0.5,1),1) );
                }
                for( unsigned int it_im_s=0; it_im_s<cd.m_VD.m_vecIM_XSegment.size(); it_im_s++ )
                    DrawSegment2( cd.m_VD.m_vecIM_XSegment[it_im_s].first,
                                  cd.m_VD.m_vecIM_XSegment[it_im_s].second,
                                  sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 4 ) );
                for( unsigned int it_ib=0; it_ib<cd.m_VD.m_vecIM_IB1_Points.size(); it_ib++ )
                {
                    sfr::Vec2 p( cd.m_VD.m_vecIM_IB1_Points[it_ib] );
                    DrawPoint2( p, sfr::gfx::Style( sfr::gfx::Color(0,1,0,1), 4 ) );
                    sfr::Vec2 n( cd.m_VD.m_vecIM_IB1_Normals[it_ib] );
                    sfr::Real r( cd.m_VD.m_vecIM_IB1_Radii[it_ib] );
                    sfr::Vec2 t( mal::PerpendicularCW( n ) );
                    const float cOffset( mal::RandomF<float>(0,0.1f) );
                    DrawSegment2( p + cOffset*n - t*r, p + cOffset*n + t*r, sfr::gfx::Style( sfr::gfx::Color(0.5,1,0,1), 2 ) );
                }
                for( unsigned int it_ib=0; it_ib<cd.m_VD.m_vecIM_IB2_Points.size(); it_ib++ )
                {
                    sfr::Vec2 p( cd.m_VD.m_vecIM_IB2_Points[it_ib] );
                    DrawPoint2( p, sfr::gfx::Style( sfr::gfx::Color(0,1,0,1), 4 ) );
                    sfr::Vec2 n( cd.m_VD.m_vecIM_IB2_Normals[it_ib] );
                    sfr::Real r( cd.m_VD.m_vecIM_IB2_Radii[it_ib] );
                    sfr::Vec2 t( mal::PerpendicularCW( n ) );
                    const float cOffset( mal::RandomF<float>(0,0.1f) );
                    DrawSegment2( p + cOffset*n - t*r, p + cOffset*n + t*r, sfr::gfx::Style( sfr::gfx::Color(0.5,1,0,1), 2 ) );
                }
                for( unsigned int it_ib=0; it_ib<cd.m_VD.m_vecIM_IB_Discarded.size(); it_ib++ )
                    DrawPoint2( cd.m_VD.m_vecIM_IB_Discarded[it_ib],
                                sfr::gfx::Style( sfr::gfx::Color(1,1,1,1), 4 ) );
                for( unsigned int it_im_m=0; it_im_m<cd.m_VD.m_vecIM_Map.size(); it_im_m++ )
                    DrawSegment2( cd.m_VD.m_vecIM_Map[it_im_m].first,
                                  cd.m_VD.m_vecIM_Map[it_im_m].second,
                                  sfr::gfx::Style( sfr::gfx::Color(0,1,0,1), 2 ) );
                for( unsigned int it_im_m=0; it_im_m<cd.m_VD.m_vecIM_RMap.size(); it_im_m++ )
                    DrawSegment2( cd.m_VD.m_vecIM_RMap[it_im_m].first + sfr::Vec2(0.1,0.1),
                                  cd.m_VD.m_vecIM_RMap[it_im_m].second + sfr::Vec2(0.1,0.1),
                                  sfr::gfx::Style( sfr::gfx::Color(0,0,1,1), 1 ) );
#endif
            }

            // Draw Pressure stuff
            if( m_rAppTask.m_bIsEnabledPressure )
                DrawDisk( sfr::Vec3(m_rAppTask.m_PressurePoint.x(),m_rAppTask.m_PressurePoint.y(),0),
                          sfr::Mat3x3::Identity(),
                          m_rAppTask.m_PressureRadius,
                          sfr::gfx::Style(sfr::gfx::Color(0,0,1,1),5) );
        }
        else //dim == 3
        {
            if( m_rAppTask.m_bIsEnabledDrag )
                DrawSphere( m_rAppTask.m_DragNodePoint3D, sfr::Mat3x3::Identity(), 0.1, sfr::gfx::Style(sfr::gfx::Color(0,1,0,1), 1, sfr::gfx::Style::eSolid) );
            if( m_rAppTask.m_bIsEnabledMove )
                DrawSphere( m_rAppTask.m_MoveOrRotateHitPoint3D, sfr::Mat3x3::Identity(), 0.1, sfr::gfx::Style(sfr::gfx::Color(1,0,0,0.5), 1, sfr::gfx::Style::eSolid) );
            // Draw Contact stuff
            if( m_rAppTask.m_bIsEnabledContact3D )
            {
                if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_SelectedObjects ) )
                {
                    // Selected objects
                    if( 0 != m_rAppTask.m_pContactKine3D_1 )
                        DrawPoint( m_rAppTask.m_pContactKine3D_1->GetPos(), sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),5) );
                    if( 0 != m_rAppTask.m_pContactKine3D_2 )
                        DrawPoint( m_rAppTask.m_pContactKine3D_2->GetPos(), sfr::gfx::Style(sfr::gfx::Color(0,0,1,1),5) );
                }
                // ContactData: should have its own Draw_ContactData2() method somewhere, i guess
                const geo::np::ContactData3& cd( m_rAppTask.m_ContactData3 );

                if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_ContactPoints | Params::ContactRender::eDF_ContactNormal ) )
                {
                    if( cd.Size() > 0 )
                    {
                        for( unsigned int it_cp=0; it_cp<cd.Size(); it_cp++ )
                        {
                            // Pos
                            if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_ContactPoints ) )
                            {
                                DrawPoint( cd.GetCP(it_cp).m_Pos1, sfr::gfx::Style(sfr::gfx::Color(0,1,1,1),4) );
                                DrawPoint( cd.GetCP(it_cp).m_Pos2, sfr::gfx::Style(sfr::gfx::Color(1,0,1,1),4) );
                            }
                            // Normal
                            if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_ContactNormal ) )
                            {
                                DrawVector( cd.GetCP(it_cp).m_Pos1,
                                            cd.GetCP(it_cp).m_Normal * cd.GetCP(it_cp).m_Depth,
                                            sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),3) );
                                // DrawSegment( cd.GetCP(it_cp).m_Pos1,
                                //              cd.GetCP(it_cp).m_Pos1 + cd.GetCP(it_cp).m_Normal * cd.GetCP(it_cp).m_Depth,
                                //              sfr::gfx::Style(sfr::gfx::Color(0,1,0,1),2) );
                            }
                        }
                    }
                }
#ifdef __ENABLE_GEO_NP_CONTACTDATA_VIZDATA
                //DCR3 vs Primitive
                if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2Primitive_IC ) )
                    for( auto segment : cd.m_VD.m_vecDCR3_Vs_Primitive_vecClippedT_Segments )
                        DrawSegment( segment.first, segment.second, sfr::gfx::Style( Vec4f(0,0.5,0,1), 1 ) );
                //DCR3 vs DCR3
                mal::SetRandomSeed(666);
                if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IC ) )
                    for( auto ic : cd.m_VD.m_vecDCR2DCR_vecIC )
                    {
                        Vec4f color( mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) );
                        // for( auto segment : ic )
                        //     //DrawSegment( segment.first, segment.second, sfr::gfx::Style( color, 1 ) );
                        //     DrawVector( segment.first, segment.second-segment.first, sfr::gfx::Style( color, 1 ) );
                        for( unsigned int it_segment=0; it_segment<ic.size(); it_segment++ )
                        {
                            // color[3] = 1 - (float(it_segment) / ic.size()); //Unnecessary, CCW sufficiently proved
                            DrawSegment( ic[it_segment].first, ic[it_segment].second, sfr::gfx::Style( color, 6 ) );
                            //\todo Nice to se CCW orientation DrawVector( ic[it_segment].first, ic[it_segment].second-ic[it_segment].first, sfr::gfx::Style( color, 1 ) );
                        }
                    }
                // Unconnected CTP
                if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IC_UnconnectedCTP ) )
                    for( auto p : cd.m_VD.m_vecDCR2DCR_vecUnconnectedCTP )
                        DrawPoint( p, sfr::gfx::Style( Vec4f(1,1,1,1), 8 ) );

                /*TEMP: Supeseded by m_vecDCR2DCR_IB_vecInternalT_Side_x_IB
                for( auto it : cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side[0] )
                    DrawTriangle( it.first, it.second, it.third, sfr::gfx::Style( Vec4f(0.5,0,0,1), 1, sfr::gfx::Style::eSolid ) );
                for( auto it : cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side[1] )
                    DrawTriangle( it.first, it.second, it.third, sfr::gfx::Style( Vec4f(0,0,0.5,1), 1, sfr::gfx::Style::eSolid ) );
                */

                for( int it_side=0; it_side<2; it_side++ )
                {
                    // SeedHE
                    if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IC_SeedHE ) )
                        for( auto he : cd.m_VD.m_vecDCR2DCR_IB_vecSeedHE_Side[it_side] )
                        {
                            Vec4f color = it_side == 0
                                          ? Vec4f(0.75,0,0,1)
                                          : Vec4f(0,0,0.75,1);
                            DrawVector( he.first, he.second-he.first, sfr::gfx::Style( color, 2 ) );
                            // DrawPoint( he.second, sfr::gfx::Style( Vec4f(0.75,0,0,1), 5 ) );
                        }

                    // IB.T
                    if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IB_FloodT ) )
                        for( unsigned int it_ib=0; it_ib<cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[it_side].size(); it_ib++ )
                        {
                            float lambda01( 1 - float(it_ib)/cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[it_side].size() );
                            switch( m_rParams.contact_render.m_FloodRM )
                            {
                            case Params::eRM_Wire:
                                {
                                    Vec4f color = it_side == 0
                                                  ? Vec4f( 0.5*lambda01, 0, 0, 1 )
                                                  : Vec4f( 0, 0, 0.5*lambda01, 1 );
                                    for( auto tri : cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[it_side][it_ib] )
                                        DrawTriangle( tri.first, tri.second, tri.third, sfr::gfx::Style( color, 4, sfr::gfx::Style::eWire ) );
                                }
                                break;
                            case Params::eRM_Solid:
                                {
                                    Vec4f color = it_side == 0
                                                  ? Vec4f( 0.25 + 0.5*lambda01, 0, 0, 1 )
                                                  : Vec4f( 0, 0, 0.25 + 0.5*lambda01, 1 );
                                    for( auto tri : cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[it_side][it_ib] )
                                    {
                                        Vec3f p0( tri.first );
                                        Vec3f p1( tri.second );
                                        Vec3f p2( tri.third );
                                        Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
                                        float L( - mal::Dot( normal, m_LightDir1 ) + mal::Dot( normal, m_LightDir2 ) );
                                        float factor( m_AmbientCoeff + m_DirectCoeff*mal::Clamp01(L) );
                                        Vec4f c( factor * color );
                                        c[3] = 1;
                                        DrawTriangle( tri.first, tri.second, tri.third, sfr::gfx::Style( c, 1, sfr::gfx::Style::eSolid ) );
                                    }
                                }
                                break;
                            case Params::eRM_Flat:
                                {
                                    if( it_side == 0 )
                                    {
                                        // Vec4f color = it_side == 0
                                        //               ? Vec4f( 0.5*lambda01, 0, 0, 1 )
                                        //               : Vec4f( 0, 0, 0.5*lambda01, 1 );
                                        // for( auto tri : cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[it_side][it_ib] )
                                        //     DrawTriangle( tri.first, tri.second, tri.third, sfr::gfx::Style( color, 1, sfr::gfx::Style::eWire ) );
                                    }
                                    else
                                    {

                                        Vec4f color = it_side == 0
                                                      ? Vec4f( 0.25 + 0.5*lambda01, 0, 0, 0.2 )
                                                      : Vec4f( 0, 0, 0.25 + 0.5*lambda01, 1 );
                                        for( auto tri : cd.m_VD.m_vecDCR2DCR_IB_vecInternalT_Side_x_IB[it_side][it_ib] )
                                        {
                                            Vec3f p0( tri.first );
                                            Vec3f p1( tri.second );
                                            Vec3f p2( tri.third );
                                            Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
                                            float L( - mal::Dot( normal, m_LightDir1 ) + mal::Dot( normal, m_LightDir2 ) );
                                            float factor( m_AmbientCoeff + m_DirectCoeff*mal::Clamp01(mal::Abs(L)) );
                                            Vec4f c( factor * color );
                                            c[3] = 1;
                                            DrawTriangle( tri.first, tri.second, tri.third, sfr::gfx::Style( c, 1, sfr::gfx::Style::eSolid ) );
                                        }
                                    }
                                }
                                break;
                            default: break;
                            }
                        }
                    // IB and IB.P
                    if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IB_FloodP | Params::ContactRender::eDF_DCR2DCR_IB_PosNormal ) )
                        for( unsigned int it_ib=0; it_ib<cd.m_VD.m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[it_side].size(); it_ib++ )
                        {
                            float lambda01( 1 - float(it_ib)/cd.m_VD.m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[it_side].size() );
                            Vec4f color = it_side == 0
                                          ? Vec4f( 0.5*lambda01, 0, 0, 1 )
                                          : Vec4f( 0, 0, 0.5*lambda01, 1 );
                            // IB pos/normal
                            if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IB_PosNormal ) )
                                DrawVector( cd.m_VD.m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[it_side][it_ib].first,
                                            0.1f*cd.m_VD.m_vecDCR2DCR_IB_PosAndNormal_Side_x_IB[it_side][it_ib].second,
                                            sfr::gfx::Style( color, 5 ) );
                            // IB.P pos/normal
                            if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IB_FloodP ) )
                                for( auto patch : cd.m_VD.m_vecDCR2DCR_IB_vecP_Side_x_IB[it_side][it_ib] )
                                    DrawVector( patch.first, 0.025f*patch.second, sfr::gfx::Style( color, 3 ) );
                        }
                    /* Clipped triangles (whole)
                    {
                        Vec4f color = it_side == 0
                                      ? Vec4f(1,1,0,0.25f)
                                      : Vec4f(0,1,1,0.25f);
                        for( auto ct : cd.m_VD.m_vecDCR2DCR_IB_vecClippedT_Side[it_side] )
                            DrawTriangle( ct.first, ct.second, ct.third, sfr::gfx::Style( color, 1, sfr::gfx::Style::eSolid ) );
                    }
                    */
                    // Clipped Edges (clipped edges CCW or radial from centroid)
                    if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IB_CE ) )
                        for( auto segment : cd.m_VD.m_vecDCR2DCR_IB_vecCE_Side[it_side] )
                        {
                            Vec4f color = it_side == 0
                                          ? Vec4f(0.75,0.75,0,1)
                                          : Vec4f(0,0.75,0.75,1);
                            DrawSegment( segment.first, segment.second, sfr::gfx::Style( color, 2 ) );
                        }
                    // Clipped Triangle aggregate stuff (centroid + normal)
                    if( m_rParams.contact_render.m_DF.Test( Params::ContactRender::eDF_DCR2DCR_IB_CT ) )
                        for( auto ct : cd.m_VD.m_vecDCR2DCR_IB_vecCT_Side[it_side] )
                        {
                            Vec4f color = it_side == 0
                                          ? Vec4f(0.75,0.75,0,1)
                                          : Vec4f(0,0.75,0.75,1);
                            DrawVector( ct.first, ct.second, sfr::gfx::Style( color, 2 ) );
                        }
                }
#endif
            }
        }
        return true;
    }

private:
    const AppTask& m_rAppTask;
    const Params& m_rParams;
};

#endif //TEST_SAPHYRE_APP_RENDERER_H
