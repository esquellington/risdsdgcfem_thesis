#ifndef SFR_GFX_VIZ_RENDERER_H
#define SFR_GFX_VIZ_RENDERER_H

#include <Safra/core/gfx/IRenderer.h>
#include <util/VizStream.h>

namespace sfr { namespace gfx
{

/*! Renders util::ItemStream Viz data independently from where its
  generated (DS or BS)
*/
class VizRenderer: public sfr::gfx::RendererPropagated
{
public:
    enum EDrawFlags {
        eDraw_Axis   = (1<<0),
        eDraw_GridX  = (1<<1),
        eDraw_GridY  = (1<<2),
        eDraw_GridZ  = (1<<3),
        eDraw_Grids  = (eDraw_GridX | eDraw_GridY | eDraw_GridZ),
        eDraw_Default = (eDraw_Axis),
        eDraw_All    = (eDraw_Axis | eDraw_Grids)
    };
    
public:
    VizRenderer( const util::VizStream *p_vs, sfr::gfx::IRenderer *p_renderer, int32 draw_flags = eDraw_Default )
    : sfr::gfx::RendererPropagated( p_renderer )
    , m_pVS( p_vs )
    , m_DrawFlags( draw_flags )
    {}    
    
    ~VizRenderer() {}

    bool Render();

private:
    void RenderItem( const util::VizItem &it );
    
private:
    const util::VizStream *m_pVS;
    Flags32 m_DrawFlags;
};

} } // namespace sfr::gfx

#endif // SFR_GFX_VIZ_RENDERER_H
