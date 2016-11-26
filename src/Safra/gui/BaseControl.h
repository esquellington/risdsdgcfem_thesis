#ifndef SFR_GUI_BASE_CONTROL_H
#define SFR_GUI_BASE_CONTROL_H

#include <Safra/Config.h>
#include <Safra/core/gfx/types.h>
#include <Safra/core/gfx/IRenderer.h>
#include <Safra/core/gui/IMouseListener.h>
#include <Safra/core/gui/IKeyboardListener.h>
#include <Safra/core/gui/IWidget.h>

namespace sfr { namespace gui
{

enum EDesktopControlFlags {
    eDCF_None                        = 0,
    eDCF_HighlightOnFocus            = (1<<0),
    eDCF_AlwaysConsumeMouseEvents    = (1<<1), //!< Any mouse event inside box will be consumed even if no child does
    eDCF_AlwaysConsumeKeyboardEvents = (1<<2), //!< Any keyboard event inside box will be consumed even if no child does
    eDCF_AutoRemoveAndDeleteChildren = (1<<3), //!< On delete, remove and delete all children automatically
    // Corners
    eDCF_Corner_UL                   = (1<<4),
    eDCF_Corner_UR                   = (1<<5),
    eDCF_Corner_BL                   = (1<<6),
    eDCF_Corner_BR                   = (1<<7),
    // Default
    eDCF_Default                     = (eDCF_None | eDCF_AlwaysConsumeMouseEvents)// | eDCF_HighlightOnFocus)
};

/*! Base Desktop Control interface & implementation.
  Manages:
  - Recursive Update and Draw
  - Input focus, Assigning or Retrieving it from children
    controls:  
    - Keeps track of its Focused immediate child, retrieving its focus
      (maybe recursively) if lost.    
    - Forwards ANY keyboard event to the potential FocusedChild, or
      rejects it if nonexistent.    
    - Forwards ANY mouse event to the potential FocusedChild first. If
      FocusedChild does not exist or rejects the event, it'll be
      iteratively forwarded to other children overlapping mouse position
      until one of them accepts it.
  - Ignores children modification by default

  \note Does NEVER delete children
  automatically. RemoveAndDeleteChildren() must be explicitly called
  from a reimplemented virtual destructor OR from outside.
*/
class BaseControl
{
protected:
    enum EStatusFlags {
        eStatus_None    = 0,
        eStatus_Visible = (1<<0), //!< Draws if set
        eStatus_Active  = (1<<1), //!< Updates if set
        eStatus_Focused = (1<<2), //!< Focused if set
        eStatus_Default = (eStatus_Visible | eStatus_Active)
    };
    
public:
    BaseControl();
    virtual ~BaseControl();

    //\name Definition and Queries
    //@{
    virtual bool SetFlags( Flags32 flags ) { m_Flags = flags; return true; }
    virtual bool SetPosRel( Vec2 pos_rel ) { if( pos_rel != GetPosRel() ) { m_PosRel = pos_rel; m_TimeStamp++; } return true; }
    virtual bool SetSizes( Vec2 sizes ) { if( sizes != GetSizes() ) { m_Sizes = sizes; m_TimeStamp++; } return true; }
    virtual bool SetSizes_Layout( Vec2 sizes ) { m_Sizes = sizes; return true; } //TEMPORAL!!
    virtual void SetColor( gfx::Color color ) { m_Color = color; }
    
    virtual Flags32 GetFlags() const { return m_Flags; }
    virtual Vec2 GetPosRel() const { return m_PosRel; }
    virtual Vec2 GetPosAbs() const { Vec2 parent_pos(0,0); if( GetParent() ) parent_pos = GetParent()->GetPosAbs(); return parent_pos+m_PosRel; }
    virtual Vec2 GetSizes() const { return m_Sizes; }
    virtual gfx::Color GetColor() const { return m_Color; }

    virtual Vec2 GetNaturalSizes() const { return m_Sizes; }
    //@}    

    //\name Focus Handling
    //@{
    virtual bool OnFocusGained() { m_StatusFlags.Enable(eStatus_Focused); return true; }
    virtual void OnFocusLost() { if( m_pFocusedChild ) m_pFocusedChild->OnFocusLost(); m_pFocusedChild = 0; m_StatusFlags.Disable(eStatus_Focused); }
    bool HasFocus() const { return m_StatusFlags.Test(eStatus_Focused); }
    //@}

    //\name Status Handling
    //@{
    inline bool IsVisible() const { return m_StatusFlags.Test(eStatus_Visible); }
    inline void SetVisible( bool b_visible ) { if(b_visible) m_StatusFlags.Enable(eStatus_Visible); else m_StatusFlags.Disable(eStatus_Visible); }
    inline bool IsActive() const { return m_StatusFlags.Test(eStatus_Active); }
    inline void SetActive( bool b_active ) { if(b_active) m_StatusFlags.Enable(eStatus_Active); else m_StatusFlags.Disable(eStatus_Active); }
    //@}

    //\name Control Nesting/Hierarchy
    //@{
    inline bool SetParent( BaseControl *p_ctrl ) { m_pParent = p_ctrl; return true; }
    inline bool SetNextSibling( BaseControl *p_ctrl ) { m_pNextSibling = p_ctrl; return true; }
    bool AddChild( BaseControl *p_ctrl );
    bool RemoveChild( BaseControl *p_ctrl );

    inline BaseControl *GetParent() const { return m_pParent; }
    inline BaseControl *GetNextSibling() const { return m_pNextSibling; }
    inline BaseControl *GetFirstChild() const { return m_pFirstChild; }

    /*\todo This was a good idea but does not compile strictly due to T() cast, specially in 64b platforms
    template<typename T> inline void SetParentDataSlot( T pds ) { m_PDS = reinterpret_cast<void*>(pds); }
    template<typename T> inline T GetParentDataSlot() const { return T(m_PDS); }
    */
    inline void SetParentTimeStamp( uint32 ts ) { m_ParentTimeStamp = ts; }
    inline uint32 GetParentTimeStamp() const { return m_ParentTimeStamp; }
    inline uint32 GetTimeStamp() const { return m_TimeStamp; }
    
    inline bool UnlinkParentAndSibling() { SetParent(0); SetNextSibling(0); return true; }
    void RemoveAndDeleteChildren();  
    //@}

    //!\name Update-cycle methods
    //@{
    virtual bool Draw( gfx::IRenderer *p_renderer );
    virtual bool Update( float dt );
    //@}

public:
    //!\name Input processing
    //@{
    virtual bool OnMouseButtonPressed( EMouseButton mb, Vec2 p );
    virtual bool OnMouseButtonReleased( EMouseButton mb, Vec2 p );
    virtual bool OnMouseMotion( Vec2 p );
    virtual bool OnKeyPressed( EKey key, Vec2 p );
    virtual bool OnKeyReleased( EKey key, Vec2 p );
    //@}

public:
    //!\name IWidget specific API retrieval
    //@{
    virtual IWidget::ITreeAPI *GetAPI_Tree() { return 0; }
    virtual IWidget::IPropertyTreeAPI *GetAPI_PropertyTree() { return 0; }
    virtual IWidget::IPlotAPI *GetAPI_Plot() { return 0; }
    //@}
    
protected:
    inline bool IsPointInside( const Vec2 &point_abs )
        {
            Vec2 point_rel( point_abs - GetPosAbs() );
            return ( point_rel.x() >= 0 && point_rel.x() <= m_Sizes.x()
                     && point_rel.y() >= 0 && point_rel.y() <= m_Sizes.y() );
        }
    
protected:
    Flags32 m_Flags;
    Vec2 m_PosRel;
    Vec2 m_Sizes;
    gfx::Color m_Color;
    double m_Time;
    Flags32 m_StatusFlags;
    uint32 m_TimeStamp;
    
    BaseControl *m_pParent;
    BaseControl *m_pNextSibling;
    BaseControl *m_pFirstChild;
    BaseControl *m_pFocusedChild;
    //void *m_PDS;
    uint32 m_ParentTimeStamp;
};

} } // namespace sfr::gui

#endif // SFR_GUI_BASE_CONTROL_H
