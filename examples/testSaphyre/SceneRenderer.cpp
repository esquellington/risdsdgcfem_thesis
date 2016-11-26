#include "SceneRenderer.h"
#include "Scene.h"
#include <Safra/core/gfx/IRenderer.h>
#include <Geo/shape/shape.h>
#include "Params.h"

//#define __ENABLE_DRAW_SOLID3D_DEGENERATE_ELEMENTS
//#define __ENABLE_DRAW_SOLID3D_VERTICES_AND_EDGES
//#define __ENABLE_DRAW_SOLID3D_BOUNDARY
//#define __ENABLE_DRAW_SOLID3D_BOUNDARY_FACES_FLAT_SHADING
//#define __ENABLE_DRAW_SOLID3D_HACKED_GROUND

#define __ENABLE_RENDER_MIG2015
#ifdef __ENABLE_RENDER_MIG2015
#  include <GL/gl.h>
#endif

bool SceneRenderer::Render()
{
    // Render propagated renderer
    GetRenderer()->Render();

    //RenderLabel( Vec3f(0,0,0), "(0,0,0)", sfr::gfx::Color(1,1,1) );

    // TEMP: DCLFEM
    // Axis
    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_Axis) )
        DrawRefSys( Vec3f::Zero(), Mat3x3f::Identity(), 25.0f, sfr::gfx::Style() );

    // Coordinate-plane grids with 5m spacing
    float grid_size = 50.0f;
    float half_grid_size = 0.5f * grid_size;
    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_GridX) )
        DrawGrid( Vec3f(0,-half_grid_size,-half_grid_size),
                  Vec3f(0,1,0), Vec3f(0,0,1),
                  grid_size, grid_size,
                  10, 10,
                  sfr::gfx::Style(sfr::gfx::Color(0.5,0.1,0.1,1),1) );
    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_GridY) )
        DrawGrid( Vec3f(-half_grid_size,0,-half_grid_size),
                  Vec3f(1,0,0), Vec3f(0,0,1),
                  grid_size, grid_size,
                  10, 10,
                  sfr::gfx::Style(sfr::gfx::Color(0.1,0.5,0.1,1),1) );
    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_GridZ) )
        DrawGrid( Vec3f(-half_grid_size,-half_grid_size,0),
                  Vec3f(1,0,0), Vec3f(0,1,0),
                  grid_size, grid_size,
                  10, 10,
                  sfr::gfx::Style(sfr::gfx::Color(0.1,0.1,0.5,1),1) );
    //

    // Draw entities polimorphically
    for( Scene::EntityIterator it=m_rScene.GetEntityIterator(); it.IsValid(); ++it )
        Draw_Entity(*it);

    return true;
}

void SceneRenderer::Draw_Entity( const S2::IEntity *p_entity )
{
    switch( p_entity->GetType() )
    {
    case S2::eEntity_Kine2D: Draw_Kine2D(static_cast<const S2::Kine2D*>(p_entity),0.5,0.5,0.5); break;
    case S2::eEntity_Kine3D: Draw_Kine3D(static_cast<const S2::Kine3D*>(p_entity),0.5,0.5,0.5); break;
    case S2::eEntity_Particle2D: Draw_Particle2D( static_cast<const S2::Particle2D*>(p_entity), 1, 1, 1 ); break;
    case S2::eEntity_Particle3D: Draw_Particle3D( static_cast<const S2::Particle3D*>(p_entity), 1, 1, 1 ); break;
    case S2::eEntity_ParticleSys2D: Draw_PS2D( static_cast<const S2::ParticleSys2D*>(p_entity), 0, 0, 1 ); break;
    case S2::eEntity_Fluid2D: Draw_Fluid2D( static_cast<const S2::Fluid2D*>(p_entity), 0, 0, 1 ); break;
    case S2::eEntity_Solid2D: Draw_Solid2D( static_cast<const S2::Solid2D*>(p_entity), 1, 0, 1 ); break;
    case S2::eEntity_Solid3D: Draw_Solid3D( static_cast<const S2::Solid3D*>(p_entity), 0.8, 0.1, 1 ); break;
    default: break;
    }
}

void SceneRenderer::Draw_Kine2D( const S2::Kine2D *p_kine2d, float r, float g, float b )
{
    Draw_Shape2( p_kine2d->GetShape(),
                 p_kine2d->GetTransform(), p_kine2d->GetVecDOF(),
                 sfr::gfx::Style(sfr::gfx::Color(r,g,b,1),1.0f,sfr::gfx::Style::eWire),
                 m_rParams.scene_render.m_KineRM );
}

void SceneRenderer::Draw_Kine3D( const S2::Kine3D *p_kine3d, float r, float g, float b )
{
    Draw_Shape3( p_kine3d->GetShape(),
                 p_kine3d->GetTransform(), p_kine3d->GetVecDOF(),
                 sfr::gfx::Style(sfr::gfx::Color(r,g,b,1),1.0f,sfr::gfx::Style::eWire),
                 m_rParams.scene_render.m_KineRM );
}

void SceneRenderer::Draw_Particle3D( const S2::Particle3D *p_particle3d, float r, float g, float b )
{
    DrawSphere( p_particle3d->GetPos(), Mat3x3f::Identity(), 0.25*p_particle3d->GetMass(),
                sfr::gfx::Style(sfr::gfx::Color(r,g,b,1)) );
}

void SceneRenderer::Draw_Particle2D( const S2::Particle2D *p_particle2d, float r, float g, float b )
{
    S2::Vec2 pos(p_particle2d->GetPos());
    DrawDisk( Vec3f(pos[0],pos[1],0), Mat3x3f::Identity(), 0.25*p_particle2d->GetMass(),
              sfr::gfx::Style(sfr::gfx::Color(r,g,b,1),1.0f,sfr::gfx::Style::eWire) );
}

/*
  void SceneRenderer::Draw_Rigid( const S2::Rigid *p_rigid, float r, float g, float b )
  {
  RenderRefSys( p_rigid->GetPos(), p_rigid->GetOrientation(), 3.0f );
  RenderBox( p_rigid->GetPos(), p_rigid->GetOrientation(),
  p_rigid->GetBodyInertiaTensor().Norm()*S2::Vec3f( 1.0f/mal::Sqrt( p_rigid->GetBodyInertiaTensor().x() ),
  1.0f/mal::Sqrt( p_rigid->GetBodyInertiaTensor().y() ),
  1.0f/mal::Sqrt( p_rigid->GetBodyInertiaTensor().z() ) ),
  sfr::gfx::Color(r,g,b) );
  }
*/

void SceneRenderer::Draw_PS2D( const S2::ParticleSys2D *p_ps2d, float r, float g, float b )
{
    for( unsigned int it_particles=0; it_particles < p_ps2d->GetNumParticles(); it_particles++ )
    {
        S2::Point2 pos = p_ps2d->GetPos(it_particles);
        DrawPoint2( pos, sfr::gfx::Style(sfr::gfx::Color(r,g,b,1),5) );
    }
}

void SceneRenderer::Draw_Fluid2D( const S2::Fluid2D *p_fluid2d, float r, float g, float b )
{
    for( unsigned int it_particles=0; it_particles < p_fluid2d->GetNumParticles(); it_particles++ )
    {
        S2::Point2 pos = p_fluid2d->GetPos(it_particles);
        DrawPoint2( pos, sfr::gfx::Style(sfr::gfx::Color(r,g,b,1),5) );
    }
}

void SceneRenderer::Draw_Solid2D( const S2::Solid2D *p_solid2d, float r, float g, float b )
{
    if( 0 != p_solid2d->GetEmbeddedGO() )
        Draw_Object_Polygonal2( p_solid2d->GetEmbeddedGO(),
                                sfr::gfx::Style( sfr::gfx::Color(1,1,1,1), 2 ),
                                m_rParams.scene_render.m_EmbRM );

    const S2::Solid2D::geo_object_type *pMSSGO( p_solid2d->GetMeshGO() );
    Draw_Shape_MeshSolid2( pMSSGO->GetShape(), pMSSGO->GetTransform(), pMSSGO->GetVecDOF(),
                           sfr::gfx::Style( sfr::gfx::Color(1,1,1,1), 1, sfr::gfx::Style::eSolid ),
                           m_rParams.scene_render.m_DynRM, m_rParams.scene_render.m_DynDCRRM );

#define __ENABLE_DCLFEM_RENDER
#ifdef __ENABLE_DCLFEM_RENDER
    const geo::MeshSolidShape2* pMSS( pMSSGO->GetShape() );
    Transform2f tr( pMSSGO->GetTransform() );
    const geo::Vec2* actual_sdof( ( 0 != pMSSGO->GetVecSDOF() ) ? pMSSGO->GetVecSDOF() : pMSS->GetVecDefaultSDOF() );
    // Boundary Polygons
    for( unsigned int it_bp=0; it_bp < pMSS->GetNumBoundaryP(); it_bp++ )
    {
        unsigned int it_he( pMSS->BP_FirstHEID(it_bp) );
        int num_neighbours(0);
        Vec2f p0( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) );
        do
        {
            if( pMSS->HE_Sym(it_he) != geo::cInvalidFeatureIndex ) num_neighbours++;
            it_he = pMSS->HE_Next(it_he);
            Vec2f p1( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) );
            //DrawSegment2( p0, p1, sfr::gfx::Style( sfr::gfx::Color(r,g,b,1), 2 ) );
            //DrawPoint2( p0, sfr::gfx::Style( sfr::gfx::Color(r,g,b,1), 4 ) );
            DrawSegment2( p0, p1, sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 4 ) ); //TEMP: For DCLFEM
            // DrawPoint2( p0, sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 4 ) ); //TEMP: For DCLFEM
            p0 = p1;
        } while ( it_he != pMSS->BP_FirstHEID(it_bp) );
    }
    //
    /* Interior degenerate triangles
    for( unsigned int it_p=0; it_p < pMSS->GetNumP(); it_p++ )
    {
        unsigned int it_he( pMSS->P_FirstHEID(it_p) );
        Vec2f p0( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
        Vec2f p1( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
        Vec2f p2( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
        GEO_ASSERT( it_he == pMSS->P_FirstHEID(it_p) ); //Closed triangle
        if( mal::Det( mal::GMat2x2_From_Columns( p1-p0, p2-p0 ) ) <= 0 )
        {
            //DrawTriangle2( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color(1,0,0,1), 2, sfr::gfx::Style::eWire ) );
            DrawTriangle2( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color(1,0,0,1), 5, sfr::gfx::Style::eWire ) ); //TEMP: For DCLFEM
        }
    }
    */
    /* Interior undegenerate triangles (eWire) TEMP: DCLFEM
       for( unsigned int it_p=0; it_p < pMSS->GetNumP(); it_p++ )
       {
       unsigned int it_he( pMSS->P_FirstHEID(it_p) );
       Vec2f p0( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
       Vec2f p1( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
       Vec2f p2( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
       GEO_ASSERT( it_he == pMSS->P_FirstHEID(it_p) ); //Closed triangle
       if( mal::Det( mal::GMat2x2_From_Columns( p1-p0, p2-p0 ) ) > 0 )
       {
       //DrawTriangle2( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color(r,g,b,1), 1, sfr::gfx::Style::eWire ) );
       DrawTriangle2( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color(0.5,0.5,0.5,1), 3, sfr::gfx::Style::eWire ) ); //TEMP: DCLFEM
       }
       }
    */
    // Interior undegenerate triangles (eSolid) TEMP: DCLFEM
    for( unsigned int it_p=0; it_p < pMSS->GetNumP(); it_p++ )
    {
        unsigned int it_he( pMSS->P_FirstHEID(it_p) );
        Vec2f p0( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
        Vec2f p1( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
        Vec2f p2( tr * pMSS->V_Pos( pMSS->HE_OriginVID(it_he), actual_sdof ) ); it_he = pMSS->HE_Next(it_he);
        GEO_ASSERT( it_he == pMSS->P_FirstHEID(it_p) ); //Closed triangle
        if( mal::Det( mal::GMat2x2_From_Columns( p1-p0, p2-p0 ) ) > 0 )
        {
            float lambda01( float(it_p & 0xFFFFFFFE)/(pMSS->GetNumP()-1) ); //Remove last bit of idx to ensure both tris in the same rect have the same color
            DrawTriangle2( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color( 0.1*r + lambda01*0.6*r, 0.7*g, 0.7*b,1), 2, sfr::gfx::Style::eSolid ) );
        }
    }
#endif
}

void SceneRenderer::Draw_Solid3D( const S2::Solid3D *p_solid3d, float r, float g, float b )
{
    const int cNumColors(3);
    static const sfr::gfx::Color vec_colors[cNumColors] = { sfr::gfx::Color(0.9 ,0.3, 0   ,1),
                                                            sfr::gfx::Color(0   ,0.9, 0.3 ,1),
                                                            sfr::gfx::Color(0.3 ,0  , 0.9 ,1) };
    static int s_color_index(0);
    s_color_index++;

    if( 0 != p_solid3d->GetEmbeddedGO() )
        Draw_Object_TriSurface3( p_solid3d->GetEmbeddedGO(),
                                 // sfr::gfx::Style( sfr::gfx::Color(0.9,0.3,0,1), 3, sfr::gfx::Style::eSolid ), //Single color
                                 sfr::gfx::Style( vec_colors[s_color_index%cNumColors], 3, sfr::gfx::Style::eSolid ), //Per-object color (%3) for multi_armadillo1K_x_4.skr scene
                                 m_rParams.scene_render.m_EmbRM );

    const S2::Solid3D::geo_object_type *pTSSGO( p_solid3d->GetMeshGO() );
    Draw_Shape_TetSolid3( pTSSGO->GetShape(), pTSSGO->GetTransform(), pTSSGO->GetVecDOF(),
                          sfr::gfx::Style( sfr::gfx::Color(0.5,0.1,0.6,1), 2, sfr::gfx::Style::eSolid ),
                          m_rParams.scene_render.m_DynRM, m_rParams.scene_render.m_DynDCRRM );

    const geo::TetSolidShape3* pTSS3( pTSSGO->GetShape() );
    Transform3f tr( pTSSGO->GetTransform() );
    const geo::Vec3 *actual_sdof( ( 0 != pTSSGO->GetVecSDOF() ) ? pTSSGO->GetVecSDOF() : pTSS3->GetVecDefaultSDOF() );

#ifdef __ENABLE_DRAW_SOLID3D_DEGENERATE_ELEMENTS
    // Interior degenerate elements
    // Internal | Tetrahedrons => Edges //TEMP: For DCLFEM
    sfr::gfx::Style style( sfr::gfx::Color(1,0,0,1), 5, sfr::gfx::Style::eWire );
    for( unsigned int it_tet=0; it_tet < pTSS3->GetNumT(); it_tet++ )
    {
        // HE, expanded from barycenter to differentiate from symmetrics
        Vec3f p0( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tet,0), actual_sdof ) );
        Vec3f p1( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tet,1), actual_sdof ) );
        Vec3f p2( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tet,2), actual_sdof ) );
        Vec3f p3( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tet,3), actual_sdof ) );
        if( mal::Det( mal::GMat3x3_From_Columns( p1-p0, p2-p0, p3-p0 ) ) < 0 )
        {
            DrawTriangle( p0, p2, p1, style );
            DrawTriangle( p0, p3, p2, style );
            DrawTriangle( p0, p1, p3, style );
            DrawTriangle( p1, p2, p3, style );
        }
    }
#endif
#ifdef __ENABLE_DRAW_SOLID3D_VERTICES_AND_EDGES
    // Boundary | Triangles => Vertices and Edges //TEMP: For DCLFEM
    for( unsigned int it_bf=0; it_bf < pTSS3->GetNumBF(); it_bf++ )
    {
        // HE, expanded from barycenter to differentiate from symmetrics
        Vec3f p0( tr * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,0), actual_sdof ) );
        Vec3f p1( tr * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,1), actual_sdof ) );
        Vec3f p2( tr * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,2), actual_sdof ) );
        DrawTriangle( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 1 ) );
        DrawPoint( p0, sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 5 ) );
        DrawPoint( p1, sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 5 ) );
        DrawPoint( p2, sfr::gfx::Style( sfr::gfx::Color(1,1,0,1), 5 ) );
        /*TEMP: Normals
          Vec3f normal( mal::Normalized( mal::Cross( p1-p0, p2-p0 ) ) );
          Vec3f barycenter( (1.0f/3.0f) * (p0+p1+p2) );
          DrawVector( barycenter, 0.25f*normal, sfr::gfx::Style( sfr::gfx::Color(0.2,0.2,1,1), 5 ) );
        */
    }
#endif //__ENABLE_DRAW_SOLID3D_VERTICES_AND_EDGES

#ifdef __ENABLE_DRAW_SOLID3D_BOUNDARY
    // Boundary | Triangles => Faces //TEMP: For DCLFEM
    for( unsigned int it_bf=0; it_bf < pTSS3->GetNumBF(); it_bf++ )
    {
        // HE, expanded from barycenter to differentiate from symmetrics
        Vec3f p0( tr * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,0), actual_sdof ) );
        Vec3f p1( tr * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,1), actual_sdof ) );
        Vec3f p2( tr * pTSS3->V_Pos( pTSS3->BF_VID(it_bf,2), actual_sdof ) );
#  ifdef __ENABLE_DRAW_SOLID3D_BOUNDARY_FACES_FLAT_SHADING
        Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
        float L1( -mal::Dot( normal, m_LightDir1 ) );
        float L2( -mal::Dot( normal, m_LightDir2 ) );
        float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
        float factor( m_AmbientCoeff + m_DirectCoeff*L );
        Vec4f color( factor*r,
                     factor*g,
                     factor*b,
                     1 );
        DrawTriangle( p0, p1, p2, sfr::gfx::Style( color, 2, sfr::gfx::Style::eSolid ) ); //TEMP: For DCLFEM
#  else
        float lambda01( float(it_bf & 0xFFFFFFFE)/(pTSS3->GetNumT()-1) ); //Remove last bit of idx to ensure both tris in the same rect have the same color
        //DrawTriangle( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color( 0.1*r + lambda01*0.6*r, 0.7*g, 0.7*b,0.75), 2, sfr::gfx::Style::eSolid ) ); //TEMP: For DCLFEM
        DrawTriangle( p0, p1, p2, sfr::gfx::Style( sfr::gfx::Color( 0.1*r + lambda01*0.6*r, 0.7*g, 0.7*b,1), 2, sfr::gfx::Style::eSolid ) ); //TEMP: For DCLFEM
#  endif
    }
#endif //__ENABLE_DRAW_SOLID3D_BOUNDARY

#ifdef __ENABLE_DRAW_SOLID3D_HACKED_GROUND
    // Ground
    const float cHalfSize = 2.0f;
    DrawQuad( Vec3f(-cHalfSize,-1,-cHalfSize),
              Vec3f( cHalfSize,-1,-cHalfSize),
              Vec3f( cHalfSize,-1, cHalfSize),
              Vec3f(-cHalfSize,-1, cHalfSize),
              sfr::gfx::Style( sfr::gfx::Color( 0.5, 0.5, 0.5, 0.75 ), 1, sfr::gfx::Style::eSolid ) ); //TEMP: For DCLFEM
#endif
}

//\todo Low level draw methods, only depend on geo and Safra, could be refactored out into sfr::VizRenderer or similar

void SceneRenderer::Draw_DCR( const geo::DCR_MeshSolidShape2* p_dcr,
                              const geo::MeshSolidShape2* p_mss, const geo::Transform2& tr, const geo::Real* vec_dof,
                              Params::ERenderMethod rm, const Vec4f& dcr_color )
{
    if( rm == Params::eRM_None ) return;

    std::vector<geo::Vec2> alloc_v(p_dcr->m_NumVertices); //\todo Scope-allocation, consider optimizing with a Viz-specific scratchpad!
    geo::Vec2* vec_v = &alloc_v[0]; //\todo p_context->m_ScratchPad.NewPOD<Vec2>(pDCR1->m_NumVertices); //POD => Don't call ctor/dtor
    const geo::Vec2* default_sdof( p_mss->GetVecDefaultSDOF() );
    const geo::Vec2* actual_sdof( ( 0 != vec_dof ) ? reinterpret_cast<const geo::Vec2*>(vec_dof) : default_sdof );

    // Apply baricentric transform and tr to all V
    for( unsigned int it_e=0; it_e<p_dcr->m_NumElements; it_e++ )
    {
        const geo::DCR_MeshSolidShape2::Element& ed( p_dcr->m_vecE[it_e] );
        uint32 vec_nid[3];
        p_mss->P_VecVID( it_e, vec_nid, 3 );
        geo::Mat2x2 B( mal::GMat2x2_From_Columns(actual_sdof[vec_nid[1]]-actual_sdof[vec_nid[0]],
                                                 actual_sdof[vec_nid[2]]-actual_sdof[vec_nid[0]]) );
        geo::Mat2x2 B0( mal::GMat2x2_From_Columns(default_sdof[vec_nid[1]]-default_sdof[vec_nid[0]],
                                                  default_sdof[vec_nid[2]]-default_sdof[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
        geo::Mat2x2 B_InvB0( B * mal::Inverse(B0) );
        // Optimized version
        geo::Transform2 tr_B_InvB0( tr.m_Pos, tr.m_Rot * B_InvB0 );
        geo::Vec2 p0_0( tr.m_Rot * actual_sdof[vec_nid[0]] );
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
            vec_v[ed.m_FirstVID + it_vie] = tr_B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]]) + p0_0;
    }
    // Per-patch segments
    mal::SetRandomSeed(666);
    for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
    {
        Vec4f patch_color( m_rParams.scene_render.m_DCR_UsePatchColor ? mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) : dcr_color );
        const geo::DCR_MeshSolidShape2::Patch& pd( p_dcr->m_vecP[it_patch] );
        for( unsigned int it_sid=0; it_sid < pd.m_NumSegments; it_sid++ )
        {
            unsigned int sid( pd.m_FirstSID + it_sid );
            DrawSegment2( vec_v[p_dcr->m_vecS[sid].GetVID(0)],
                          vec_v[p_dcr->m_vecS[sid].GetVID(1)],
                          sfr::gfx::Style( patch_color, 1, sfr::gfx::Style::eSolid ) );
            //TEMP: show segment direction
            // Vec2f n( mal::SafeNormalized( mal::PerpendicularCW(vec_v[p_dcr->m_vecS[sid].m_vecVID[1]]-vec_v[p_dcr->m_vecS[sid].m_vecVID[0]] ) ) );
            // DrawTriangle2( (Vec2f)vec_v[p_dcr->m_vecS[sid].m_vecVID[0]] - n*0.01,
            //                (Vec2f)vec_v[p_dcr->m_vecS[sid].m_vecVID[0]] + n*0.01,
            //                (Vec2f)vec_v[p_dcr->m_vecS[sid].m_vecVID[1]],
            //                sfr::gfx::Style( color, 2, sfr::gfx::Style::eSolid ) );
        }
    }
}

void SceneRenderer::Draw_DCR( const geo::DCR_TetSolidShape3* p_dcr,
                              const geo::TetSolidShape3* p_tss, const geo::Transform3& tr, const geo::Real* vec_dof,
                              Params::ERenderMethod rm, const Vec4f& dcr_color )
{
    if( rm == Params::eRM_None ) return;

    std::vector<geo::Vec3> alloc_v(p_dcr->m_NumVertices); //\todo Scope-allocation, consider optimizing with a Viz-specific scratchpad!
    geo::Vec3* vec_v = &alloc_v[0]; //\todo p_context->m_ScratchPad.NewPOD<Vec2>(pDCR1->m_NumVertices); //POD => Don't call ctor/dtor
    const geo::Vec3* default_sdof( p_tss->GetVecDefaultSDOF() );
    const geo::Vec3* actual_sdof( ( 0 != vec_dof ) ? reinterpret_cast<const geo::Vec3*>(vec_dof) : default_sdof );

    // Apply baricentric transform and tr to all V
    for( unsigned int it_e=0; it_e<p_dcr->m_NumElements; it_e++ )
    {
        const geo::DCR_TetSolidShape3::Element& ed( p_dcr->m_vecE[it_e] );
        uint32 vec_nid[4] = { p_tss->T_VID(it_e,0), p_tss->T_VID(it_e,1), p_tss->T_VID(it_e,2), p_tss->T_VID(it_e,3) };
        geo::Mat3x3 B( mal::GMat3x3_From_Columns(actual_sdof[vec_nid[1]]-actual_sdof[vec_nid[0]],
                                                 actual_sdof[vec_nid[2]]-actual_sdof[vec_nid[0]],
                                                 actual_sdof[vec_nid[3]]-actual_sdof[vec_nid[0]] ) );
        geo::Mat3x3 B0( mal::GMat3x3_From_Columns(default_sdof[vec_nid[1]]-default_sdof[vec_nid[0]],
                                                  default_sdof[vec_nid[2]]-default_sdof[vec_nid[0]],
                                                  default_sdof[vec_nid[3]]-default_sdof[vec_nid[0]]) ); //\todo Could be precomputed in DCR::ED
        geo::Mat3x3 B_InvB0( B * mal::Inverse(B0) );
        /*\note Direct but slower code
          for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
          vec_v[ed.m_FirstVID + it_vie] = tr * ( B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]])
          + actual_sdof[vec_nid[0]] );
        */
        // Optimized version
        geo::Transform3 tr_B_InvB0( tr.m_Pos, tr.m_Rot * B_InvB0 );
        geo::Vec3 p0_0( tr.m_Rot * actual_sdof[vec_nid[0]] );
        for( unsigned int it_vie=0; it_vie<ed.m_NumVertices; it_vie++ )
            vec_v[ed.m_FirstVID + it_vie] = tr_B_InvB0 * (p_dcr->m_vecV[ed.m_FirstVID + it_vie] - default_sdof[vec_nid[0]]) + p0_0;
    }
    // Per-patch triangles
    mal::SetRandomSeed(666);
    switch( rm )
    {
    case Params::eRM_Wire:
        if( m_rParams.scene_render.m_DCR_BoundaryOnly )
        {
            for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
            {
                Vec4f patch_color( m_rParams.scene_render.m_DCR_UsePatchColor ? mal::RandomV( Vec4f(0,0,0,0.5), Vec4f(1,1,1,0.5) ) : dcr_color );
                const geo::DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_patch] );
                for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                {
                    unsigned int tid( pd.m_FirstTID + it_tid );
                    for( unsigned int it_ntit=0; it_ntit < 3; it_ntit++ )
                    {
                        unsigned int ntid( p_dcr->m_vecT[tid].GetNTID(it_ntit) );
                        if( tid < ntid
                            && !p_dcr->Is_TID_In_PID( ntid, it_patch ) )
                            DrawSegment( vec_v[p_dcr->m_vecT[tid].GetVID(it_ntit)],
                                         vec_v[p_dcr->m_vecT[tid].GetVID((it_ntit+1)%3)],
                                         sfr::gfx::Style( patch_color, 4, sfr::gfx::Style::eWire ) );
                    }
                }
            }
        }
        else
        {
            for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
            {
                Vec4f patch_color( m_rParams.scene_render.m_DCR_UsePatchColor ? mal::RandomV( Vec4f(0,0,0,0.25), Vec4f(1,1,1,0.25) ) : dcr_color );
                const geo::DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_patch] );
                for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                {
                    unsigned int tid( pd.m_FirstTID + it_tid );
                    DrawTriangle( vec_v[p_dcr->m_vecT[tid].GetVID(0)],
                                  vec_v[p_dcr->m_vecT[tid].GetVID(1)],
                                  vec_v[p_dcr->m_vecT[tid].GetVID(2)],
                                  sfr::gfx::Style( patch_color, 1, sfr::gfx::Style::eWire ) );
                }
            }
        }
        break;
    case Params::eRM_Solid:
        {
            for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
            {
                Vec4f patch_color( m_rParams.scene_render.m_DCR_UsePatchColor ? mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) : dcr_color );
                const geo::DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_patch] );
                for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                {
                    unsigned int tid( pd.m_FirstTID + it_tid );
                    DrawTriangle( vec_v[p_dcr->m_vecT[tid].GetVID(0)],
                                  vec_v[p_dcr->m_vecT[tid].GetVID(1)],
                                  vec_v[p_dcr->m_vecT[tid].GetVID(2)],
                                  sfr::gfx::Style( patch_color, 1, sfr::gfx::Style::eSolid ) );
                }
            }
        }
        break;
    case Params::eRM_Flat:
        {
#ifdef __ENABLE_RENDER_MIG2015
            glBegin( GL_TRIANGLES );
            for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
            {
                Vec4f patch_color( m_rParams.scene_render.m_DCR_UsePatchColor ? mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) : dcr_color );
                const geo::DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_patch] );
                for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                {
                    unsigned int tid( pd.m_FirstTID + it_tid );
                    Vec3f p0( vec_v[p_dcr->m_vecT[tid].GetVID(0)] );
                    Vec3f p1( vec_v[p_dcr->m_vecT[tid].GetVID(1)] );
                    Vec3f p2( vec_v[p_dcr->m_vecT[tid].GetVID(2)] );
                    Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
                    float L1( -mal::Dot( normal, m_LightDir1 ) );
                    float L2( -mal::Dot( normal, m_LightDir2 ) );
                    float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
                    float factor( m_AmbientCoeff + m_DirectCoeff*L );
                    Vec4f color( factor*patch_color );
                    color[3] = 1;
                    glColor4fv( color.AsArray() );
                    glVertex3fv( p0.AsArray() );
                    glVertex3fv( p1.AsArray() );
                    glVertex3fv( p2.AsArray() );
                }
            }
            glEnd();
#else
            for( unsigned int it_patch=0; it_patch < p_dcr->m_NumPatches; it_patch++ )
            {
                Vec4f patch_color( m_rParams.scene_render.m_DCR_UsePatchColor ? mal::RandomV( Vec4f(0,0,0,1), Vec4f(1,1,1,1) ) : dcr_color );
                const geo::DCR_TetSolidShape3::Patch& pd( p_dcr->m_vecP[it_patch] );
                for( unsigned int it_tid=0; it_tid < pd.m_NumTriangles; it_tid++ )
                {
                    unsigned int tid( pd.m_FirstTID + it_tid );
                    Vec3f p0( vec_v[p_dcr->m_vecT[tid].GetVID(0)] );
                    Vec3f p1( vec_v[p_dcr->m_vecT[tid].GetVID(1)] );
                    Vec3f p2( vec_v[p_dcr->m_vecT[tid].GetVID(2)] );
                    Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
                    float L1( -mal::Dot( normal, m_LightDir1 ) );
                    float L2( -mal::Dot( normal, m_LightDir2 ) );
                    float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
                    float factor( m_AmbientCoeff + m_DirectCoeff*L );
                    Vec4f color( factor*patch_color );
                    color[3] = 1;
                    DrawTriangle( p0, p1, p2, sfr::gfx::Style( color, 1, sfr::gfx::Style::eSolid ) );
                }
            }
#endif
        }
        break;
    default: break;
    }
}

void SceneRenderer::Draw_BVH( const geo::BVH_MeshSolidShape2* p_bvh )
{
    /*
      for( unsigned int it_node=0; it_node<p_bvh->m_vecNode.size(); it_node++ )
      {
      geo::bv::AABB2 aabb( p_bvh->m_vecNode[it_node].m_Geometry.m_BV );
      DrawBox2( aabb.GetPos(), sfr::Mat2x2::Identity(), aabb.GetHalfSizes(), sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),1) );
      }
    */
    /* AABB2
       p_bvh->Map( [this]( const geo::BVH_MeshSolidShape2& bvh, const typename geo::BVH_MeshSolidShape2::Node& node, uint32 level )
       {
       float height01( float(level) / bvh.GetHeight() );
       const geo::bv::AABB2& aabb( node.m_Geometry.m_BV );
       Vec3f pos( aabb.GetPos().x(), aabb.GetPos().y(), height01 );
       DrawQuad( pos + Vec3f(-aabb.GetHalfSizes().x(),-aabb.GetHalfSizes().y(),0),
       pos + Vec3f( aabb.GetHalfSizes().x(),-aabb.GetHalfSizes().y(),0),
       pos + Vec3f( aabb.GetHalfSizes().x(), aabb.GetHalfSizes().y(),0),
       pos + Vec3f(-aabb.GetHalfSizes().x(), aabb.GetHalfSizes().y(),0),
       sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
       return true;
       }
       );
    */
    // DOP2_K8
    p_bvh->Map( [this]( const geo::BVH_MeshSolidShape2& bvh, const typename geo::BVH_MeshSolidShape2::Node& node, uint32 level )
                {
                    float height01( float(level) / bvh.GetHeight() );
                    // Get KDOP segments
                    std::pair<geo::Vec2,geo::Vec2> vec_segments[8];
                    unsigned int num_segments = geo::bv::DOP2_K8_Segments( node.m_Geometry.m_BV, vec_segments );
                    // Viz clipped segments, regardless of actual order
                    for( unsigned int it_segment=0; it_segment < num_segments; it_segment++ )
                        DrawSegment( Vec3f( vec_segments[it_segment].first[0], vec_segments[it_segment].first[1], height01 ),
                                     Vec3f( vec_segments[it_segment].second[0], vec_segments[it_segment].second[1], height01 ),
                                     sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
                    return true;
                }
        );
    /* Sphere2
       p_bvh->Map( [this]( const geo::BVH_MeshSolidShape2& bvh, const typename geo::BVH_MeshSolidShape2::Node& node, uint32 level )
       {
       float height01( float(level) / bvh.GetHeight() );
       const geo::bv::Sphere2& sphere( node.m_Geometry.m_BV );
       Vec3f pos( sphere.GetPos().x(), sphere.GetPos().y(), height01 );
       DrawDisk( pos, Mat3x3f::Identity(), sphere.m_Radius,
       sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
       return true;
       }
       );
    */
}

void SceneRenderer::Draw_BVH( const geo::BVH_TetSolidShape3* p_bvh )
{
    /* AABB3
       p_bvh->Map( [this]( const geo::BVH_TetSolidShape3& bvh, const typename geo::BVH_TetSolidShape3::Node& node, uint32 level )
       {
       float height01( float(level) / bvh.GetHeight() );
       const geo::bv::AABB3& aabb( node.m_Geometry.m_BV );
       sfr::Vec3 pos( aabb.GetPos() );
       DrawBox( pos, sfr::Mat3x3::Identity(), aabb.GetHalfSizes(), sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
       return true;
       }
       );
    */
    /* DOP3_K6
       p_bvh->Map( [this]( const geo::BVH_TetSolidShape3& bvh, const typename geo::BVH_TetSolidShape3::Node& node, uint32 level )
       {
       float height01( float(level) / bvh.GetHeight() );
       const geo::bv::DOP3_K6& dop3k6( node.m_Geometry.m_BV );
       sfr::Vec3 pos( dop3k6.GetPos() );
       DrawBox( pos, sfr::Mat3x3::Identity(),
       0.5f* sfr::Vec3( dop3k6.GetInterval<0>().Length(),
       dop3k6.GetInterval<1>().Length(),
       dop3k6.GetInterval<2>().Length() ),
       sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
       return true;
       }
       );
    */
    // DOP3_K14 TEMP: Draw only AABB3
    p_bvh->Map( [this]( const geo::BVH_TetSolidShape3& bvh, const typename geo::BVH_TetSolidShape3::Node& node, uint32 level )
                {
                    float height01( float(level) / bvh.GetHeight() );
                    const geo::bv::DOP3_K14& dop3k14( node.m_Geometry.m_BV );
                    sfr::Vec3 pos( dop3k14.GetPos() );
                    DrawBox( pos, sfr::Mat3x3::Identity(),
                             0.5f*sfr::Vec3( dop3k14.GetInterval<0>().Length(),
                                             dop3k14.GetInterval<1>().Length(),
                                             dop3k14.GetInterval<2>().Length() ),
                             sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
                    return true;
                }
        );
    //
    /* Sphere3
       p_bvh->Map( [this]( const geo::BVH_TetSolidShape3& bvh, const typename geo::BVH_TetSolidShape3::Node& node, uint32 level )
       {
       float height01( float(level) / bvh.GetHeight() );
       const geo::bv::Sphere3& sphere( node.m_Geometry.m_BV );
       Vec3f pos( sphere.GetPos() );
       DrawSphere( pos, Mat3x3f::Identity(), sphere.m_Radius, sfr::gfx::Style( sfr::gfx::Color( 1,1,1-height01,0.1 ), 1, sfr::gfx::Style::eWire ) );
       return true;
       }
       );
    */
}

void SceneRenderer::Draw_Shape2( const geo::IShape2 *p_shape, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    switch( p_shape->GetType() )
    {
    case geo::eShape_Plane2: Draw_Shape_Plane2( static_cast<const geo::PlaneShape2*>(p_shape), tr, style, m_rParams.scene_render.m_GroundRM ); break;
    case geo::eShape_Sphere2: Draw_Shape_Sphere2( static_cast<const geo::SphereShape2*>( p_shape ), tr, style, rm ); break;
    case geo::eShape_Box2: Draw_Shape_Box2( static_cast<const geo::BoxShape2*>( p_shape ), tr, style, rm ); break;
    case geo::eShape_Capsule2: Draw_Shape_Capsule2( static_cast<const geo::CapsuleShape2*>( p_shape ), tr, style, rm ); break;
    case geo::eShape_Polygonal2: Draw_Shape_Polygonal2( static_cast<const geo::PolygonalShape2*>( p_shape ), tr, vec_dof, style, rm ); break;
    case geo::eShape_Path2: Draw_Shape_Path2( static_cast<const geo::PathShape2*>( p_shape ), tr, vec_dof, style, rm ); break;
    case geo::eShape_MeshSolid2: Draw_Shape_MeshSolid2( static_cast<const geo::MeshSolidShape2*>( p_shape ), tr, vec_dof, style, rm, m_rParams.scene_render.m_KineDCRRM ); break;
    default: SFR_LOG_ERROR( "SceneRenderer: Unknwon shape type %d", p_shape->GetType() ); break;
    }
    /*\todo MAYBE
      geo::IDomainSampler2 *pDS( static_cast<const geo::IShape2*>( p_shape )->CreateDomainSampler() );
      if( 0 != pDS ) { pDS->UpdateDOF(vec_dof); Draw_DomainSampler2( pDS, tr, style ); delete pDS; }
    */
}

void SceneRenderer::Draw_Shape3( const geo::IShape3 *p_shape, const S2::Transform3 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    switch( p_shape->GetType() )
    {
    case geo::eShape_Plane3: Draw_Shape_Plane3( static_cast<const geo::PlaneShape3*>(p_shape), tr, style, m_rParams.scene_render.m_GroundRM ); break;
    case geo::eShape_Sphere3: Draw_Shape_Sphere3( static_cast<const geo::SphereShape3*>( p_shape ), tr, style, rm ); break;
    case geo::eShape_TetSolid3: Draw_Shape_TetSolid3( static_cast<const geo::TetSolidShape3*>( p_shape ), tr, vec_dof, style, rm, m_rParams.scene_render.m_KineDCRRM ); break;
    default: SFR_LOG_ERROR( "SceneRenderer: Unknwon shape type %d", p_shape->GetType() ); break;
    }
}

void SceneRenderer::Draw_Shape_Plane2( const geo::PlaneShape2 *p_plane2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    sfr::Vec2 n( p_plane2->GetNormal() );
    sfr::Vec2 p(0,0);
    if( n[0] != sfr::Real(0) )
        p[0] = -p_plane2->GetCoeffD() / n[0];
    else if( n[1] != sfr::Real(0) )
        p[1] = -p_plane2->GetCoeffD() / n[1];
    sfr::Real size(100);
    sfr::Vec2 p1( p + size*sfr::Vec2(-n[1],n[0]) );
    sfr::Vec2 p2( p - size*sfr::Vec2(-n[1],n[0]) );
    DrawSegment2( p1, p2, style );
}

void SceneRenderer::Draw_Shape_Plane3( const geo::PlaneShape3 *p_plane3, const S2::Transform3 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    sfr::Vec3 n( p_plane3->GetNormal() );
    sfr::Vec3 p(0,0,0);
    sfr::Vec3 u(0,0,0);
    if( n[0] != sfr::Real(0) )
    {
        p[0] = -p_plane3->GetCoeffD() / n[0];
        u = mal::Normalized( mal::Cross(n,sfr::Vec3(0,1,0) ) );
    }
    else if( n[1] != sfr::Real(0) )
    {
        p[1] = -p_plane3->GetCoeffD() / n[1];
        u = mal::Normalized( mal::Cross(n,sfr::Vec3(0,0,1) ) );
    }
    else if( n[2] != sfr::Real(0) )
    {
        p[2] = -p_plane3->GetCoeffD() / n[2];
        u = mal::Normalized( mal::Cross(n,sfr::Vec3(1,0,0) ) );
    }
    sfr::Vec3 v( mal::Cross(n,u) );
    sfr::Real size(12);
    sfr::Vec3 p1( p - size*u - size*v );
    sfr::Vec3 p2( p + size*u - size*v );
    sfr::Vec3 p3( p + size*u + size*v );
    sfr::Vec3 p4( p - size*u + size*v );

    // Overwrite style with RM
    sfr::gfx::Style s( style );
    switch( rm )
    {
    case Params::eRM_Wire: s.m_Flags = sfr::gfx::Style::eWire; break;
    case Params::eRM_Solid: s.m_Flags = sfr::gfx::Style::eSolid; break;
    case Params::eRM_Flat:
        {
            s.m_Flags = sfr::gfx::Style::eSolid;
            float L1( -mal::Dot( n, m_LightDir1 ) );
            float L2( -mal::Dot( n, m_LightDir2 ) );
            float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
            float factor( m_AmbientCoeff + m_DirectCoeff*L );
            s.m_Color = factor * style.m_Color;
            s.m_Color[3] = style.m_Color[3];
        }
        break;
    default: break;
    }

    DrawQuad( p1, p2, p3, p4, s );
}

void SceneRenderer::Draw_Shape_Sphere2( const geo::SphereShape2 *p_sphere2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    DrawDisk2( tr.m_Pos, p_sphere2->GetRadius(), style, 32 );
    /*TEMP: Works fine, no need to clutter viz
      geo::IDomainSampler2 *pDS( p_sphere2->CreateDomainSampler() );
      if( 0 != pDS ) { pDS->UpdateDOF(vec_dof); Draw_DomainSampler2( pDS, tr, style ); delete pDS; }
    */
}

void SceneRenderer::Draw_Shape_Sphere3( const geo::SphereShape3 *p_sphere3, const S2::Transform3 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    // Overwrite style with RM
    sfr::gfx::Style s( style );
    switch( rm )
    {
    case Params::eRM_Wire: s.m_Flags = sfr::gfx::Style::eWire; break;
    case Params::eRM_Solid: s.m_Flags = sfr::gfx::Style::eSolid; break;
    case Params::eRM_Flat:
        {
            s.m_Flags = sfr::gfx::Style::eSolid;
            //\todo
            // float L1( -mal::Dot( normal, m_LightDir1 ) );
            // float L2( -mal::Dot( normal, m_LightDir2 ) );
            // float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
            // float factor( m_AmbientCoeff + m_DirectCoeff*L );
            // s.m_Color = factor * style.m_Color;
            // s.m_Color[3] = style.m_Color[3];
        }
        break;
    default: break;
    }

    DrawSphere( tr.m_Pos, tr.m_Rot, p_sphere3->GetRadius(), s );
    /*TEMP: Works fine, no need to clutter viz
      geo::IDomainSampler2 *pDS( p_sphere2->CreateDomainSampler() );
      if( 0 != pDS ) { pDS->UpdateDOF(vec_dof); Draw_DomainSampler2( pDS, tr, style ); delete pDS; }
    */
}

void SceneRenderer::Draw_Shape_Box2( const geo::BoxShape2 *p_box2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    DrawBox2( tr.m_Pos, tr.m_Rot, p_box2->GetHalfSizes(), style );
}

void SceneRenderer::Draw_Shape_Capsule2( const geo::CapsuleShape2 *p_capsule2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    DrawDisk2( tr.m_Pos - mal::Column(1,tr.m_Rot) * p_capsule2->GetHalfHeight(), p_capsule2->GetRadius(), style, 32 );
    DrawDisk2( tr.m_Pos + mal::Column(1,tr.m_Rot) * p_capsule2->GetHalfHeight(), p_capsule2->GetRadius(), style, 32 );
    DrawBox2( tr.m_Pos, tr.m_Rot, sfr::Vec2( p_capsule2->GetRadius(), p_capsule2->GetHalfHeight() ), style );
}

void SceneRenderer::Draw_Shape_Polygonal2( const geo::PolygonalShape2 *p_polygonal2, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    const geo::Vec2* actual_sdof( ( 0 != vec_dof ) ? reinterpret_cast<const geo::Vec2*>(vec_dof) : p_polygonal2->GetVecDefaultSDOF() );

    int num_edges( p_polygonal2->IsClosed() ? p_polygonal2->GetNumVertices() : p_polygonal2->GetNumVertices()-1 );
    for( int it_edge=0; it_edge < num_edges; it_edge++ )
    {
        sfr::Vec2 p0( tr * p_polygonal2->V_Pos( it_edge, actual_sdof ) );
        sfr::Vec2 p1( tr * p_polygonal2->V_Pos( (it_edge+1) % p_polygonal2->GetNumVertices(), actual_sdof ) );
        DrawSegment2( p0, p1, style );
    }
}

void SceneRenderer::Draw_Shape_Path2( const geo::PathShape2 *p_path2, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    const geo::Vec2* actual_sdof( ( 0 != vec_dof ) ? reinterpret_cast<const geo::Vec2*>(vec_dof) : p_path2->GetVecDefaultSDOF() );

    for( unsigned int it_edge=0; it_edge < p_path2->GetNumE(); it_edge++ )
    {
        sfr::Vec2 p0( tr * p_path2->V_Pos( p_path2->E_OriginVID(it_edge), actual_sdof ) );
        sfr::Vec2 p1( tr * p_path2->V_Pos( p_path2->E_FinalVID(it_edge), actual_sdof ) );
        switch( p_path2->E_CurveType(it_edge) )
        {
        case geo::PathShape2::edge_type::eCT_Line:
            DrawSegment2( p0, p1, style );
            break;
        case geo::PathShape2::edge_type::eCT_Bezier3:
            {
                // Bezier curve
                sfr::Vec2 bezier3_a( tr * p_path2->E_Data(it_edge).m_ParamA ); //\todo If params are also SDOF, use actual_sdof here!!
                sfr::Vec2 bezier3_b( tr * p_path2->E_Data(it_edge).m_ParamB );
                sfr::Vec2 b0( p0 );
                const unsigned int cNumSamples(10);
                for( unsigned int i=1; i<cNumSamples; i++ )
                {
                    sfr::Vec2 b1 = geo::Eval_Bezier3( p0, p1, bezier3_a, bezier3_b, sfr::Real(i) / (cNumSamples-1) );
                    DrawSegment2( b0, b1, style );
                    b0 = b1;
                }
            }
            break;
        default: break;
        }
    }
}

void SceneRenderer::Draw_Shape_MeshSolid2( const geo::MeshSolidShape2 *p_meshsolid2, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm, Params::ERenderMethod dcr_rm )
{
    if( rm != Params::eRM_None )
    {
        const geo::Vec2 *actual_sdof( 0 != vec_dof ? reinterpret_cast<const geo::Vec2*>(vec_dof) : p_meshsolid2->GetVecDefaultSDOF() );

        // Boundary Polygons
        for( unsigned int it_bp=0; it_bp < p_meshsolid2->GetNumBoundaryP(); it_bp++ )
        {
            unsigned int it_he( p_meshsolid2->BP_FirstHEID(it_bp) );
            int num_neighbours(0);
            sfr::Vec2 p0( tr * p_meshsolid2->V_Pos( p_meshsolid2->HE_OriginVID(it_he), actual_sdof ) );
            do
            {
                if( p_meshsolid2->HE_Sym(it_he) != geo::cInvalidFeatureIndex ) num_neighbours++;
                it_he = p_meshsolid2->HE_Next(it_he);
                sfr::Vec2 p1( tr * p_meshsolid2->V_Pos( p_meshsolid2->HE_OriginVID(it_he), actual_sdof ) );
                DrawSegment2( p0, p1, style );
                p0 = p1;
            } while ( it_he != p_meshsolid2->BP_FirstHEID(it_bp) );
        }

        /*TEMP: Testing backwards BP iteration
          for( unsigned int it_bp=0; it_bp < p_meshsolid2->GetNumBoundaryP(); it_bp++ )
          {
          unsigned int it_he( p_meshsolid2->BP_FirstHEID(it_bp) );
          int num_neighbours(0);
          sfr::Vec2 p0( tr * p_meshsolid2->V_Pos( p_meshsolid2->HE_FinalVID(it_he), actual_sdof ) );
          do
          {
          if( p_meshsolid2->HE_Sym(it_he) != geo::cInvalidFeatureIndex ) num_neighbours++;
          it_he = p_meshsolid2->HE_Prev(it_he);
          sfr::Vec2 p1( tr * p_meshsolid2->V_Pos( p_meshsolid2->HE_FinalVID(it_he), actual_sdof ) );
          DrawSegment2( p0, p1, style );
          p0 = p1;
          } while ( it_he != p_meshsolid2->BP_FirstHEID(it_bp) );
          }
        */

        /*TEMP: Works fine, no need to clutter viz
          geo::IDomainSampler2 *pDS( p_meshsolid2->CreateDomainSampler() );
          if( 0 != pDS ) { pDS->UpdateDOF( reinterpret_cast<const geo::Real*>(actual_sdof) ); Draw_DomainSampler2( pDS, tr, style ); delete pDS; }
        */
    }

    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_DCR)
        && p_meshsolid2->GetDCR() )
        Draw_DCR( p_meshsolid2->GetDCR(),
                  p_meshsolid2, tr, vec_dof, dcr_rm, Vec4f(0,1,1,1) );
    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_BVH)
        && p_meshsolid2->GetBVH() )
        Draw_BVH( p_meshsolid2->GetBVH() );
}

void SceneRenderer::Draw_Shape_TetSolid3( const geo::TetSolidShape3 *p_tetsolid3, const S2::Transform3 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm, Params::ERenderMethod dcr_rm )
{
    if( rm != Params::eRM_None )
    {
        const geo::Vec3 *actual_sdof( 0 != vec_dof ? reinterpret_cast<const geo::Vec3*>(vec_dof) : p_tetsolid3->GetVecDefaultSDOF() );

        // Overwrite style with RM
        sfr::gfx::Style s( style );
        switch( rm )
        {
        case Params::eRM_Wire:
            s.m_Flags = sfr::gfx::Style::eWire;
            s.m_PenSize = 2;
            break;
        case Params::eRM_Solid:
        case Params::eRM_Flat:
            s.m_Flags = sfr::gfx::Style::eSolid; break;
        default: break;
        }

        // Boundary | Triangles => Faces
        for( unsigned int it_bf=0; it_bf < p_tetsolid3->GetNumBF(); it_bf++ )
        {
            Vec3f p0( tr * p_tetsolid3->V_Pos( p_tetsolid3->BF_VID(it_bf,0), actual_sdof ) );
            Vec3f p1( tr * p_tetsolid3->V_Pos( p_tetsolid3->BF_VID(it_bf,1), actual_sdof ) );
            Vec3f p2( tr * p_tetsolid3->V_Pos( p_tetsolid3->BF_VID(it_bf,2), actual_sdof ) );
            Vec4f color( style.m_Color );
            if( rm == Params::eRM_Flat )
            {
                Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
                float L1( -mal::Dot( normal, m_LightDir1 ) );
                float L2( -mal::Dot( normal, m_LightDir2 ) );
                float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
                float factor( m_AmbientCoeff + m_DirectCoeff*L );
                color *= factor;
                color[3] = s.m_Color[3];
            }
            DrawTriangle( p0, p1, p2, sfr::gfx::Style( color, s.m_PenSize, s.m_Flags ) );
        }
    }

    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_Coords) )
    {
        char str[128];
        snprintf( str, 128, "(%f,%f,%f)", tr.Pos().x(), tr.Pos().y(), tr.Pos().z() );
        DrawLabel( tr.Pos(), str, sfr::gfx::Style(sfr::gfx::Color(1,1,1,1),1) );
    }

    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_DCR)
        && p_tetsolid3->GetDCR() )
        Draw_DCR( p_tetsolid3->GetDCR(),
                  p_tetsolid3, tr, vec_dof, dcr_rm, Vec4f(0.75,0.1,0.75,1) );
    if( m_rParams.scene_render.m_DF.Test(Params::SceneRender::eDF_BVH)
        && p_tetsolid3->GetBVH() )
        Draw_BVH( p_tetsolid3->GetBVH() );
}

void SceneRenderer::Draw_DomainSampler2( geo::IDomainSampler2 *p_ds2, const S2::Transform2 &tr, const sfr::gfx::Style &style )
{
    /* Samples
       const unsigned int cNumPOB(100);
       geo::point_on_boundary_id_type vec_pob[cNumPOB];
       p_ds2->CreateRandomPOB( vec_pob, cNumPOB );
       for( unsigned int it_pob=0; it_pob<cNumPOB; it_pob++ )
       {
       S2::Vec2 pos( tr * p_ds2->POB_Position( vec_pob[it_pob] ) );
       S2::Vec2 normal( tr.m_Rot * p_ds2->POB_Normal( vec_pob[it_pob] ) );
       DrawPoint2( pos, sfr::gfx::Style(style.m_Color,3) );
       DrawSegment2( pos, pos + normal, style );
       }
       p_ds2->APOB_Destroy( vec_pob, cNumPOB );
    */
    /* Neighbours
       const unsigned int cNumNeighbours(5);
       geo::point_on_boundary_id_type vec_neighbours[cNumNeighbours];
       geo::point_on_boundary_id_type neighbour_seed_pob = p_ds2->CreateRandomPOB();
       p_ds2->CreateRandomNeighbourPOB( vec_neighbours, cNumNeighbours,
       0.5f, neighbour_seed_pob );
       for( unsigned int it_neighbour=0; it_neighbour<cNumNeighbours; it_neighbour++ )
       {
       S2::Vec2 pos( tr * p_ds2->POB_Position( vec_neighbours[it_neighbour] ) );
       S2::Vec2 normal( tr.m_Rot * p_ds2->POB_Normal( vec_neighbours[it_neighbour] ) );
       DrawPoint2( pos, sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),5) );
       DrawSegment2( pos, pos + normal, sfr::gfx::Style(sfr::gfx::Color(1,0,0,1),5) );
       }
       //IMPORTANT: We DecRef instead of Clear() because p_ds2 may have other POB allocated for CD
       p_ds2->POB_Destroy( neighbour_seed_pob );
       p_ds2->APOB_Destroy( vec_neighbours, cNumNeighbours );
    */

    // Stepping
    const unsigned int cNumSteps(16);
    geo::point_on_boundary_id_type vec_steps[cNumSteps];
    vec_steps[0] = p_ds2->CreateRandomPOB();
    S2::Vec2 p0( tr * p_ds2->POB_Position( vec_steps[0] ) );
    DrawPoint2( p0, sfr::gfx::Style(sfr::gfx::Color(1,0,1,1),6) );
    for( unsigned int it_step=1; it_step<cNumSteps; it_step++ )
    {
        vec_steps[it_step] = p_ds2->StepPOB( vec_steps[it_step-1],
                                             mal::PerpendicularCW( p_ds2->POB_Normal( vec_steps[it_step-1] ) ),
                                             0.5, 0 );
        S2::Vec2 p1( tr * p_ds2->POB_Position( vec_steps[it_step] ) );
        DrawSegment2( p0, p1, sfr::gfx::Style(sfr::gfx::Color(1,0,1,1),3) );
        p0 = p1;
    }
    p_ds2->APOB_Destroy( vec_steps, cNumSteps );
    //

    //p_ds2->Trace();
    //SFR_LOG_WARNING("Clear!");
    //p_ds2->Clear();
    //p_ds2->Trace();
}


void SceneRenderer::Draw_Object_Polygonal2( const geo::IObject* p_go, const sfr::gfx::Style& style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;
    const geo::PolygonalShape2* pPS2( static_cast<const geo::PolygonalShape2*>(p_go->GetShapeInterface()) );
    const geo::GObjectSDOF<2,geo::Vec2> *pPSGO( static_cast<const geo::GObjectSDOF<2,geo::Vec2> *>(p_go) );
    Transform2f tr( static_cast<const geo::IObject2*>(p_go)->GetTransform() );
    int num_edges( pPS2->IsClosed() ? pPS2->GetNumVertices() : pPS2->GetNumVertices()-1 );
    const geo::Vec2 *actual_sdof( ( 0 != pPSGO->GetVecSDOF() ) ? pPSGO->GetVecSDOF() : pPS2->GetVecDefaultSDOF() );
    for( int it_edge=0; it_edge < num_edges; it_edge++ )
    {
        sfr::Vec2 p0( tr * pPS2->V_Pos( it_edge, actual_sdof ) );
        sfr::Vec2 p1( tr * pPS2->V_Pos( (it_edge+1) % pPS2->GetNumVertices(), actual_sdof ) );
        DrawSegment2( p0, p1, style );
    }
}

void SceneRenderer::Draw_Object_TriSurface3( const geo::IObject* p_go, const sfr::gfx::Style& style, Params::ERenderMethod rm )
{
    if( rm == Params::eRM_None ) return;

    // Overwrite style with RM
    sfr::gfx::Style s( style );
    switch( rm )
    {
    case Params::eRM_Wire:
        s.m_Flags = sfr::gfx::Style::eWire;
        s.m_PenSize = 1;
        s.m_Color[3] = 0.1f;
        break;
    case Params::eRM_Solid:
    case Params::eRM_Flat:
        s.m_Flags = sfr::gfx::Style::eSolid; break;
    default: break;
    }

    const geo::TriSurfaceShape3* pTSS3( static_cast<const geo::TriSurfaceShape3*>(p_go->GetShapeInterface()) );
    const geo::GObjectSDOF<3,geo::Vec3> *pTSSGO( static_cast<const geo::GObjectSDOF<3,geo::Vec3> *>(p_go) );
    Transform3f tr( pTSSGO->GetTransform() );
    const geo::Vec3* actual_sdof( ( 0 != pTSSGO->GetVecSDOF() ) ? pTSSGO->GetVecSDOF() : pTSS3->GetVecDefaultSDOF() );

#ifdef __ENABLE_RENDER_MIG2015
    glBegin( GL_TRIANGLES );
    for( unsigned int it_tri=0; it_tri < pTSS3->GetNumT(); it_tri++ )
    {
        Vec3f p0( pTSS3->V_Pos( pTSS3->T_VID(it_tri,0), actual_sdof ) );
        Vec3f p1( pTSS3->V_Pos( pTSS3->T_VID(it_tri,1), actual_sdof ) );
        Vec3f p2( pTSS3->V_Pos( pTSS3->T_VID(it_tri,2), actual_sdof ) );
        Vec4f color( s.m_Color );
        Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
        float L1( -mal::Dot( normal, m_LightDir1 ) );
        float L2( -mal::Dot( normal, m_LightDir2 ) );
        float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
        float factor( m_AmbientCoeff + m_DirectCoeff*L );
        color *= factor;
        color[3] = s.m_Color[3];
        glColor4fv( color.AsArray() );
        glVertex3fv( p0.AsArray() );
        glVertex3fv( p1.AsArray() );
        glVertex3fv( p2.AsArray() );
    }
    glEnd();
#else
    for( unsigned int it_tri=0; it_tri < pTSS3->GetNumT(); it_tri++ )
    {
        Vec3f p0( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tri,0), actual_sdof ) );
        Vec3f p1( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tri,1), actual_sdof ) );
        Vec3f p2( tr * pTSS3->V_Pos( pTSS3->T_VID(it_tri,2), actual_sdof ) );
        Vec4f color( s.m_Color );
        if( rm == Params::eRM_Flat )
        {
            Vec3f normal( mal::SafeNormalized( mal::Cross( p1-p0, p2-p0 ) ) );
            float L1( -mal::Dot( normal, m_LightDir1 ) );
            float L2( -mal::Dot( normal, m_LightDir2 ) );
            float L( mal::Clamp01(mal::Clamp01(L1)+mal::Clamp01(L2)) );
            float factor( m_AmbientCoeff + m_DirectCoeff*L );
            color *= factor;
            color[3] = s.m_Color[3];
        }
        DrawTriangle( p0, p1, p2, sfr::gfx::Style( color, s.m_PenSize, s.m_Flags ) );
    }
#endif
}
