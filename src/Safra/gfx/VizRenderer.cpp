#include "VizRenderer.h"
#include <Mal/GConversion.h>

namespace sfr { namespace gfx
{

void VizRenderer::RenderItem( const util::VizItem &it )
{
    if( !it.IsComplex() )
    {
        /*
        // Old viz types would go here
        Vec3f pos(0,0,0);
        gfx::Color color(0.4,0.4,0.4);
        switch( it.GetType() )
        {
        case S2::debug::eType_VizPoint:
            {
                const S2::debug::VizPoint &viz_point = it.Get<S2::debug::VizPoint>();
                pos = viz_point.m_Pos;
                color = *(gfx::Color*)&viz_point.m_Color;
                DrawPoint( viz_point.m_Pos, gfx::Style(color,viz_point.m_Size) );
            }
            break;
        case S2::debug::eType_VizVec:
            {
                const S2::debug::VizVec &viz_vec = it.Get<S2::debug::VizVec>();
                pos = viz_vec.m_Pos;
                color = *(gfx::Color*)&viz_vec.m_Color;
                DrawVector( viz_vec.m_Pos, viz_vec.m_Vec, gfx::Style(color, viz_vec.m_Size) );
            }
            break;
        case S2::debug::eType_VizPosRot:
            {
                const S2::debug::VizPosRot &viz_posrot = it.Get<S2::debug::VizPosRot>();
                pos = viz_posrot.m_Pos;
                color = *(gfx::Color*)&viz_posrot.m_Color;
                DrawRefSys( viz_posrot.m_Pos, mal::GMat3x3_From(viz_posrot.m_Quat), viz_posrot.m_Size, gfx::Style() );
            }
            break;
        case S2::debug::eType_VizSegment:
            {
                const S2::debug::VizSegment &viz_segment = it.Get<S2::debug::VizSegment>();
                pos = 0.5f*(viz_segment.m_Pos0+viz_segment.m_Pos1);
                color = *(gfx::Color*)&viz_segment.m_Color;
                DrawVector( viz_segment.m_Pos0, viz_segment.m_Pos1-viz_segment.m_Pos0,
                            gfx::Style(color,viz_segment.m_Size) );
            }
            break;
        case S2::debug::eType_VizCircle:
            {
                const S2::debug::VizCircle &viz_circle = it.Get<S2::debug::VizCircle>();
                pos = viz_circle.m_Pos;
                color = *(gfx::Color*)&viz_circle.m_Color;
                DrawDisk( viz_circle.m_Pos, mal::GMat3x3_From(viz_circle.m_Quat), viz_circle.m_Radius,
                          gfx::Style(color,1.0f,gfx::Style::eWire) );
            }
            break;
        case S2::debug::eType_VizSphere:
            {
                const S2::debug::VizSphere &viz_sphere = it.Get<S2::debug::VizSphere>();
                pos = viz_sphere.m_Pos;
                color = *(gfx::Color*)&viz_sphere.m_Color;
                DrawSphere( viz_sphere.m_Pos, Mat3x3::Identity(), viz_sphere.m_Radius, gfx::Style(color) );
            }
            break;
        case S2::debug::eType_VizGrid2D:
            {
                const S2::debug::VizGrid2D &viz_grid2d = it.Get<S2::debug::VizGrid2D>();
                pos = viz_grid2d.m_Pos + 1.01*viz_grid2d.m_Dir2;
                color = *(gfx::Color*)&viz_grid2d.m_Color;
                DrawGrid2D( viz_grid2d.m_Pos,
                            viz_grid2d.m_Dir1.Normalized(), viz_grid2d.m_Dir2.Normalized(),
                            viz_grid2d.m_Dir1.Norm(),
                            viz_grid2d.m_Dir2.Norm(),
                            viz_grid2d.m_Dim1, viz_grid2d.m_Dim2,
                            gfx::Style(color,viz_grid2d.m_Size) );
            }
            break;
        default: break;
        }
        if( it.IsNamed() )
            DrawLabel( pos, it.GetName(), gfx::Style(color) );
        */
    }
    else //IsComplex
    {
        //S2::BSG::GetLogStream() << util::BInfo() << "VizItem" << it << util::EInfo();
        util::VizItem sub_it( it.GetSubItem() );
        // Get VizItem style if specified
        gfx::Style style;
        util::VizItem style_it = sub_it.Find("style");
        if( style_it.IsValid() )
        {
            SFR_ASSERT( util::eType_VizStyle == style_it.GetType() );
            style_it = style_it.GetSubItem();
            style.m_Color = gfx::Color( style_it.Find("color").SafeGet<Vec4f>(style.m_Color) );
            style.m_PenSize = style_it.Find("pen_size").SafeGet<float>(style.m_PenSize);
            style.m_Flags = style_it.Find("flags").SafeGet<Flags32>(style.m_Flags);
        }

        switch( it.GetType() )
        {
        case util::eType_VizSet:
            for( util::VizItem set_it = it.GetSubItem(); set_it.IsValid(); ++set_it )
                RenderItem(set_it);
            break;
            // Linear viz
        case util::eType_VizPoint:
            DrawPoint( sub_it.Find("pos").Get<Vec3f>(), style );
            break;
        case util::eType_VizSegment:
            DrawSegment( sub_it.Find("pos1").Get<Vec3f>(), sub_it.Find("pos2").Get<Vec3f>(), style );
            break;
        case util::eType_VizVec:
            DrawVector( sub_it.Find("pos").Get<Vec3f>(), sub_it.Find("vec").Get<Vec3f>(), style );
            break;
        case util::eType_VizRefSys:
            DrawRefSys( sub_it.Find("pos").Get<Vec3f>(),
                        mal::GRotation3x3_From(sub_it.Find("rot").Get<Quatf>()),
                        sub_it.Find("scale").Get<float>(),
                        style );
            break;
        case util::eType_VizTransform3:
            {
                Transform3f tr3( sub_it.Find("transform").Get<Transform3f>() );
                DrawRefSys( tr3.m_Pos, tr3.m_Rot, sub_it.Find("scale").Get<float>(), style );
            }
            break;
        case util::eType_VizTransform2:
            {
                Transform2f tr2( sub_it.Find("transform").Get<Transform2f>() );
                DrawRefSys2( tr2.m_Pos, tr2.m_Rot, sub_it.Find("scale").Get<float>(), style );
            }
            break;

            // Planar viz
        case util::eType_VizTriangle:
            DrawTriangle( sub_it.Find("p1").Get<Vec3f>(),
                          sub_it.Find("p2").Get<Vec3f>(),
                          sub_it.Find("p3").Get<Vec3f>(),
                          style );
            break;
        case util::eType_VizQuad:
            DrawQuad( sub_it.Find("p1").Get<Vec3f>(),
                      sub_it.Find("p2").Get<Vec3f>(),
                      sub_it.Find("p3").Get<Vec3f>(),
                      sub_it.Find("p4").Get<Vec3f>(),
                      style );
            break;
        case util::eType_VizRectangle:
            DrawRectangle( sub_it.Find("pos").Get<Vec3f>(),
                           mal::GRotation3x3_From(sub_it.Find("rot").Get<Quatf>()),
                           sub_it.Find("half_sizes").Get<Vec2f>(),
                           style );
            break;
        case util::eType_VizDisk:
            DrawDisk( sub_it.Find("pos").Get<Vec3f>(),
                      mal::GRotation3x3_From(sub_it.Find("rot").Get<Quatf>()),
                      sub_it.Find("radius").Get<float>(),
                      style );
            break;
            // Volume viz
        case util::eType_VizBox:
            DrawBox( sub_it.Find("pos").Get<Vec3f>(),
                     mal::GRotation3x3_From(sub_it.Find("rot").Get<Quatf>()),
                     sub_it.Find("half_sizes").Get<Vec3f>(),
                     style );
            break;
        case util::eType_VizSphere:
            DrawSphere( sub_it.Find("pos").Get<Vec3f>(),
                        mal::GRotation3x3_From(sub_it.Find("rot").Get<Quatf>()),
                        sub_it.Find("radius").Get<float>(),
                        style );
            break;
        case util::eType_VizCapsule: SFR_ASSERT(false); break;
        case util::eType_VizCylinder: SFR_ASSERT(false); break;
        case util::eType_VizTetrahedron:
            {
                Vec3f p0,p1,p2,p3;
                p0 = sub_it.Find("p0").Get<Vec3f>();
                p1 = sub_it.Find("p1").Get<Vec3f>();
                p2 = sub_it.Find("p2").Get<Vec3f>();
                p3 = sub_it.Find("p3").Get<Vec3f>();
                DrawTriangle( p0, p2, p1, style );
                DrawTriangle( p0, p3, p2, style );
                DrawTriangle( p0, p1, p3, style );
                DrawTriangle( p1, p2, p3, style );
            }
            break;
            // Misc
        case util::eType_VizAABB: SFR_ASSERT(false); break;
        case util::eType_VizGrid2D:
            DrawGrid( sub_it.Find("pos").Get<Vec3f>(),
                      sub_it.Find("dir1").Get<Vec3f>().Normalized(), sub_it.Find("dir2").Get<Vec3f>().Normalized(),
                      sub_it.Find("dir1").Get<Vec3f>().Norm(), sub_it.Find("dir2").Get<Vec3f>().Norm(),
                      sub_it.Find("size1").Get<uint32>(), sub_it.Find("size2").Get<uint32>(),
                      style );
            break;

        default: break;
        }

        if( it.IsNamed() && sub_it.Find("pos").IsValid() )
            DrawLabel( sub_it.Find("pos").Get<Vec3f>(), it.GetName(), style );
    }
}

bool VizRenderer::Render()
{
    // Render propagated renderer
    GetRenderer()->Render();

    // Axis
    if( m_DrawFlags.Test(eDraw_Axis) )
        DrawRefSys( Vec3::Zero(), Mat3x3::Identity(), 25.0f, gfx::Style() );

    // Coordinate-plane grids with 5m spacing
    Real grid_size = 50.0f;
    Real half_grid_size = 0.5f * grid_size;
    if( m_DrawFlags.Test(eDraw_GridX) )
        DrawGrid( Vec3(0,-half_grid_size,-half_grid_size),
                  Vec3(0,1,0), Vec3(0,0,1),
                  grid_size, grid_size,
                  10, 10,
                  gfx::Style(gfx::Color(0.5,0.1,0.1,1),1) );
    if( m_DrawFlags.Test(eDraw_GridY) )
        DrawGrid( Vec3(-half_grid_size,0,-half_grid_size),
                  Vec3(1,0,0), Vec3(0,0,1),
                  grid_size, grid_size,
                  10, 10,
                  gfx::Style(gfx::Color(0.1,0.5,0.1,1),1) );
    if( m_DrawFlags.Test(eDraw_GridZ) )
        DrawGrid( Vec3(-half_grid_size,-half_grid_size,0),
                  Vec3(1,0,0), Vec3(0,1,0),
                  grid_size, grid_size,
                  10, 10,
                  gfx::Style(gfx::Color(0.1,0.1,0.5,1),1) );

    // Render all items in the stream
    for( util::VizItem it = m_pVS->Begin(); it.IsValid(); ++it )
        RenderItem(it);

    return true;
}

} } // namespace sfr::gfx
