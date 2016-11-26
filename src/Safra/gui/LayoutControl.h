#ifndef SFR_GUI_LAYOUT_CONTROL_H
#define SFR_GUI_LAYOUT_CONTROL_H

#include "BaseControl.h"

namespace sfr { namespace gui
{

inline Vec2 Combine( const Vec2 &sizes0, const Vec2 &pos1, const Vec2 &sizes1 )
{
    return Vec2( mal::Max(sizes0[0],pos1[0]+sizes1[0]), mal::Max(sizes0[1],pos1[1]+sizes1[1]) );
}

/*! Complex Control that offers layout functionality.
  Layout is recomputed whenever:
  - A child is added or removed (=> sets m_bIsDirtyLayout = true)
  - A child has an obsolete parent-observed timestamp (this happens
    when a child changes size or subhierarchy structure, for example)
  \note Layout has no background color by default (transparent)

  \note Layout has NaturalSizes DIFFERENT from Sizes: NaturalSizes are
  the sizes computed in RecomputeLayout() to bound children controls,
  while actual Sizes may have been externally set in order to force
  alignment accross "parallel layouts", This gave enormous problems in
  TreeControl when trying to achieve lazy layout recomputation, and this
  was the cornerstone.
*/
class ILayoutControl: public BaseControl
{
public:
    enum ELayoutMode {
        eLayoutMode_Min,    //Left  | Top
        eLayoutMode_Max,    //Right | Bottom
        eLayoutMode_Center
    };

public:
    ILayoutControl()
    : m_MarginsLU(0,0)
    , m_MarginsRD(0,0)
    , m_MinSeparation(10)
    , m_Mode(eLayoutMode_Min)
    , m_bIsDirtyLayout(true)
    {
        SetColor( gfx::Color(0,0,0,0) ); //\todo Disable drawing
    }
    virtual ~ILayoutControl() {}

    //! Set Layout internal margins (begin,end) and min separation between subctl
    inline void SetMargins( Vec2 margin_lu, Vec2 margin_rd )
    {
        m_TimeStamp++;
        m_bIsDirtyLayout = true;
        m_MarginsLU = margin_lu;
        m_MarginsRD = margin_rd;
    }

    inline void SetMinSeparation( Real min_separation )
    {
        m_TimeStamp++;
        m_bIsDirtyLayout = true;
        m_MinSeparation = min_separation;
    }

    inline void SetMode( ELayoutMode lm ) { m_Mode = lm; }

    inline void SetDirtyLayout( bool b_dirty ) { m_bIsDirtyLayout = b_dirty; }

    //\name Control Nesting/Hierarchy
    //@{
    virtual bool AddChild( BaseControl *p_ctrl )
        {
            bool bResult = BaseControl::AddChild(p_ctrl);
            if( bResult ) p_ctrl->SetParentTimeStamp(0);
            m_bIsDirtyLayout = true;
            return bResult;
        }
    virtual bool RemoveChild( BaseControl *p_ctrl )
        {
            m_bIsDirtyLayout = true;
            return BaseControl::RemoveChild(p_ctrl);
        }
    //@}

    //!< Adapt layout to externally imposed sizes
    virtual bool SetSizes( Vec2 sizes )
    {
        if( sizes != GetSizes() )
        {
            m_bIsDirtyLayout = true;
            switch( m_Mode )
            {
            case eLayoutMode_Min:
                BaseControl::SetSizes(sizes);
                break;
            case eLayoutMode_Max:
                if( m_NaturalSizes != Vec2(0,0) )
                {
                    SetPosRel( GetPosRel() + Vec2(sizes[0]-m_NaturalSizes[0],0) );
                    BaseControl::SetSizes(m_NaturalSizes);
                }
                else
                    BaseControl::SetSizes(sizes);
                break;
            case eLayoutMode_Center: //\todo
                BaseControl::SetSizes(sizes);
                break;
            }
        }
        return true;
    }

    //\name Other
    //@{
    virtual bool Update( float dt )
    {
        bool res = BaseControl::Update(dt);

#define __ENABLE_LAZY_RECOMPUTE_LAYOUT
#ifdef __ENABLE_LAZY_RECOMPUTE_LAYOUT
        // Check if layout is dirty or any child has obsolete parent-observed timestamp
        bool bDirty( m_bIsDirtyLayout );
        for( BaseControl *pChild=GetFirstChild(); pChild != 0 && !bDirty; pChild=pChild->GetNextSibling() )
            bDirty = pChild->GetTimeStamp() != pChild->GetParentTimeStamp();
        // If dirty, recompute layout and set timestamps to sync
        if( bDirty )
        {
            RecomputeLayout();
            for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
                pChild->SetParentTimeStamp( pChild->GetTimeStamp() );
            m_bIsDirtyLayout = false;
        }
#else
        RecomputeLayout();
#endif
        return res;
    }
    //@}

    virtual Vec2 GetNaturalSizes() const { return m_NaturalSizes; }

public:
    //!\name Specific layout computation
    //@{
    //! Compute sizes bottom-up and place sub-controls top-down
    virtual void RecomputeLayout() = 0;
    //@}

protected:
    Vec2 m_MarginsLU;
    Vec2 m_MarginsRD;
    Real m_MinSeparation;
    ELayoutMode m_Mode;
    bool m_bIsDirtyLayout;
    Vec2 m_NaturalSizes;
};


/*! Horizontal Layout Complex Control.
  Concatenates children with a given sparation margin.
*/
class HorizontalLayout: public ILayoutControl
{
public:
    HorizontalLayout() {}
    ~HorizontalLayout() {}

public:
    //----Specific layout computation
    void RecomputeLayout()
    {
        // Compute subctrl horizontal layout positions and aggregated sizes
        Vec2 pos( m_MarginsLU );
        Vec2 sizes(0,0);
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
        {
            pChild->SetPosRel(pos);
            Vec2 natural_sizes( pChild->GetNaturalSizes() );
            natural_sizes[0] += (0!=pChild->GetNextSibling()) ? m_MinSeparation : 0;
            sizes = Combine( sizes, pChild->GetPosRel(), natural_sizes );
            pos[0] += natural_sizes[0];
        }
        // Set common v-size to all subctrl
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            pChild->SetSizes_Layout( Vec2(pChild->GetNaturalSizes()[0],sizes[1]) );
        // Set widget size as aggregated sizes + margins
        SetSizes( sizes + m_MarginsRD ); //m_MarginsLU already accounted for in pos_rel
        m_NaturalSizes = sizes;
    }
};

/*! VerticalLayout Layout Complex Control.
  Stacks children with a given sparation margin.
*/
class VerticalLayout: public ILayoutControl
{
public:
    VerticalLayout() {}
    ~VerticalLayout() {}

public:
    //----Specific layout computation
    void RecomputeLayout()
    {
        // Compute subctrl vertical layout positions and aggregated sizes
        Vec2 pos( m_MarginsLU );
        Vec2 sizes(0,0);
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
        {
            pChild->SetPosRel(pos);
            Vec2 natural_sizes( pChild->GetNaturalSizes() );
            natural_sizes[1] += (0!=pChild->GetNextSibling()) ? m_MinSeparation : 0;
            sizes = Combine( sizes, pChild->GetPosRel(), natural_sizes );
            pos[1] += natural_sizes[1];
        }
        // Set common h-size to all subctrl
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            pChild->SetSizes_Layout( Vec2(sizes[0],pChild->GetNaturalSizes()[1]) );
        // Set widget size as aggregated sizes + margins
        SetSizes( sizes + m_MarginsRD ); //m_MarginsLU already accounted for in pos_rel
        m_NaturalSizes = sizes;
    }
};

/*! TableLayout Layout Complex Control.
  Combines a Vertical layout with Horizontal tabbing for each children
  controls that are HorizontalLayouts
*/

/*\todo It's complicated, because usually the children of a
  VerticalLayout are not directly HorizontalLayouts, and therefore the
  potential table structure is hidden... eg: TreeControl
class TableLayout: public ILayoutControl
{
public:
    TableLayout() {}
    ~TableLayout() {}

public:
    //----Specific layout computation
    void RecomputeLayout()
    {
        // Compute subctrl vertical layout positions and aggregated sizes
        Vec2 pos( m_MarginsLU );
        Vec2 sizes(0,0);
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
        {
            pChild->SetPosRel(pos);
            Vec2 natural_sizes( pChild->GetNaturalSizes() );
            natural_sizes[1] += (0!=pChild->GetNextSibling()) ? m_MinSeparation : 0;
            sizes = Combine( sizes, pChild->GetPosRel(), natural_sizes );
            pos[1] += natural_sizes[1];
        }
        // Set common h-size to all subctrl
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            pChild->SetSizes_Layout( Vec2(sizes[0],pChild->GetNaturalSizes()[1]) );
        // Set widget size as aggregated sizes + margins
        SetSizes( sizes + m_MarginsRD ); //m_MarginsLU already accounted for in pos_rel
        m_NaturalSizes = sizes;
    }
};
*/

/*! Frame Layout.
  Adds a Frame to SINGLE children control... mostly required by
  TreeControl.
*/
class FrameLayout: public ILayoutControl
{
public:
    FrameLayout() : m_CountRL(0) {}
    ~FrameLayout() {}

    virtual bool AddChild( BaseControl *p_control )
    {
        if( 0 != GetFirstChild() )
            SFR_LOG_ERROR("FrameLayout:: Only works well with a SINGLE child!!");
        return ILayoutControl::AddChild(p_control);
    }

public:
    //----Specific layout computation
    void RecomputeLayout()
    {
        m_CountRL++;
        // Compute children bounding box...
        Vec2 pos( m_MarginsLU );
        Vec2 sizes(0,0);
        //\note Actually, only 1 child is expected...
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
        {
            pChild->SetPosRel(pos);
            Vec2 natural_sizes( pChild->GetNaturalSizes() );
            sizes[0] = mal::Max( sizes[0], natural_sizes[0] );
            sizes[1] = mal::Max( sizes[1], natural_sizes[1] );
        }
        // Set common h-size to all subctrl
        for( BaseControl *pChild=GetFirstChild(); pChild != 0; pChild=pChild->GetNextSibling() )
            pChild->SetSizes_Layout( sizes );
        // Set widget size as aggregated sizes + margins
        SetSizes( sizes + m_MarginsRD ); //m_MarginsLU already accounted for in pos_rel
        m_NaturalSizes = sizes;
    }
protected:
    uint32 m_CountRL;
};

}} //namespace sfr::gui

#endif //SFR_GUI_LAYOUT_CONTROL_H
