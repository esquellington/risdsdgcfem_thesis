#ifndef TEST_SAPHYRE_SCENE_RENDERER_H
#define TEST_SAPHYRE_SCENE_RENDERER_H

#include "Scene.h"
#include <Safra/core/gfx/IRenderer.h>
#include <Geo/shape/shape.h>
#include "Params.h"

/*! Renders a s2::Scene using a propagated Renderer

  Draws the objects inside a Scene in a standard way, based only on
  the information available in the S2::Object itself.

  Optionally draws the S2::DebugInfo
*/
class SceneRenderer: public sfr::gfx::RendererPropagated
{
public:
    SceneRenderer( const Params& params, const Scene &scene, sfr::gfx::IRenderer *p_renderer )
    : sfr::gfx::RendererPropagated( p_renderer )
    , m_rParams(params)
    , m_rScene( scene )
    , m_AmbientCoeff( 0.1f )//0.2f )
    , m_DirectCoeff( 1 - m_AmbientCoeff )
    , m_LightDir1( mal::Normalized(Vec3f(0 ,-1,-1)) )//0,-1,0 )
    , m_LightDir2( 0.25*mal::Normalized(Vec3f(0 ,1, -1)) )
    {}

    ~SceneRenderer() {}
    bool Render();

protected:
    void Draw_Entity( const S2::IEntity *p_entity );
    void Draw_Kine2D( const S2::Kine2D *p_kine2d, float r, float g, float b );
    void Draw_Kine3D( const S2::Kine3D *p_kine3d, float r, float g, float b );
    void Draw_Particle3D( const S2::Particle3D *p_particle3d, float r, float g, float b );
    void Draw_Particle2D( const S2::Particle2D *p_particle2d, float r, float g, float b );
    void Draw_PS2D( const S2::ParticleSys2D *p_ps2d, float r, float g, float b );
    void Draw_Fluid2D( const S2::Fluid2D *p_fluid2d, float r, float g, float b );
    void Draw_Solid2D( const S2::Solid2D *p_solid2d, float r, float g, float b );
    void Draw_Solid3D( const S2::Solid3D *p_solid3d, float r, float g, float b );

private:
    //\todo Low level draw methods, only depend on geo and Safra, could be refactored out into sfr::VizRenderer or similar
    void Draw_DCR( const geo::DCR_MeshSolidShape2* p_dcr, const geo::MeshSolidShape2* p_mss, const geo::Transform2& tr, const geo::Real* vec_dof, Params::ERenderMethod rm, const Vec4f& dcr_color );
    void Draw_DCR( const geo::DCR_TetSolidShape3* p_dcr, const geo::TetSolidShape3* p_tss, const geo::Transform3& tr, const geo::Real* vec_dof, Params::ERenderMethod rm, const Vec4f& dcr_color );
    void Draw_BVH( const geo::BVH_MeshSolidShape2* p_bvh );
    void Draw_BVH( const geo::BVH_TetSolidShape3* p_bvh );

    void Draw_Shape2( const geo::IShape2 *p_shape, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape3( const geo::IShape3 *p_shape, const S2::Transform3 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm );

    void Draw_Shape_Plane2( const geo::PlaneShape2 *p_plane2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Plane3( const geo::PlaneShape3 *p_plane3, const S2::Transform3 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Sphere2( const geo::SphereShape2 *p_sphere2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Sphere3( const geo::SphereShape3 *p_sphere3, const S2::Transform3 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Box2( const geo::BoxShape2 *p_box2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Capsule2( const geo::CapsuleShape2 *p_capsule2, const S2::Transform2 &tr, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Polygonal2( const geo::PolygonalShape2 *p_polygonal2, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_Path2( const geo::PathShape2 *p_path2, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm );
    void Draw_Shape_MeshSolid2( const geo::MeshSolidShape2 *p_meshsolid2, const S2::Transform2 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm, Params::ERenderMethod dcr_rm );
    void Draw_Shape_TetSolid3( const geo::TetSolidShape3 *p_tetsolid3, const S2::Transform3 &tr, const S2::Real *vec_dof, const sfr::gfx::Style &style, Params::ERenderMethod rm, Params::ERenderMethod dcr_rm );

    void Draw_Object_Polygonal2( const geo::IObject* p_go, const sfr::gfx::Style& style, Params::ERenderMethod rm );
    void Draw_Object_TriSurface3( const geo::IObject* p_go, const sfr::gfx::Style& style, Params::ERenderMethod rm );

    void Draw_DomainSampler2( geo::IDomainSampler2 *p_ds2, const S2::Transform2 &tr, const sfr::gfx::Style &style );

public:
    const Params& m_rParams;
    const Scene& m_rScene;

    float m_AmbientCoeff;
    float m_DirectCoeff;
    Vec3f m_LightDir1;
    Vec3f m_LightDir2;
};

#endif // TEST_SAPHYRE_SCENE_RENDERER_H
