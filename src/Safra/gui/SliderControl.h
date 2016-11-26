#ifndef SFR_GUI_SLIDER_CONTROL_H
#define SFR_GUI_SLIDER_CONTROL_H

#include "BaseControl.h"
#include <boost/function.hpp>    

namespace sfr { namespace gui
{

/*! Simple text slider that calls a callback when the value changes.

  \note D_OnChangeCall is ONLY called when the Slider pos changes,
  even if mouse events are processed. Otherwise, clicking on the exact
  slider pos or releasing the slider without moving or, even moving
  beyond limits without changing the actual position would force
  unnecessary D_OnChangeCall invocations.
  
  \todo Allow specifying if D_OnChangeCall is invoked OnRelease
  (cheap) or also during OnMove (may generate lots of events...)

  \todo Allow non-linear scales (quadratic, logarithmic...)
  \todo Allow draw style (frame, line, handle shape...)
  \todo Allow horizontal/vertical slider
*/
class SliderControl: public BaseControl
{
public:    
    typedef boost::function<void (void *)> D_OnChangeCall;
    
public:
    SliderControl( uint32 slider_size, void *p_user_data = 0 )
    : m_pUserData(p_user_data)
    , m_bIsPressed(false)
    , m_SliderSize(slider_size)
    , m_HandlePos(0)
    , m_HandleRadius(5)
    , m_HandleColorReleased(0.7,0.7,0.7,1)
    , m_HandleColorPressed(1,1,1,1)
    {
        SetSizes( Vec2(m_SliderSize,m_HandleRadius) );
    }
    ~SliderControl() {}
    
    inline void SetOnChangeCall( D_OnChangeCall on_slider_call ) { m_OnChangeCall = on_slider_call; }
    
    inline void SetUserData( void *p_user_data ) { m_pUserData = p_user_data; }
    inline void *GetUserData() const { return m_pUserData; }   

    inline uint32 GetSliderSize() const { return (uint32)m_SliderSize; }
    inline uint32 GetHandlePos() const { return (uint32)m_HandlePos; }
    // Set handle pos, return true iff changed
    inline bool SetHandlePos( int32 handle_pos )
        {
            int32 new_pos = mal::Clamp(handle_pos,int32(0),m_SliderSize);
            if( new_pos == m_HandlePos ) return false;
            m_HandlePos = new_pos;
            return true;
        }

    Vec2 GetNaturalSizes() const { return Vec2( m_SliderSize, 2*m_HandleRadius ); }
    
    //! Toggle color when pressed
    bool OnMouseButtonPressed( EMouseButton mb, Vec2 p )
    {
        if( !IsPointInside(p) )
            return false;
        else //\todo If pressed on handle, then move it smoothly, otherwise, handle jumps to clicked pos
        {
            if( RecomputeHandlePos( p ) ) OnSliderChanged();
            m_bIsPressed = true;
            //Handle color LabelControl::SetFontColor( m_FontColorPressed );
            return true;
        }
    }

    //! Call delegate and toggle color when released
    bool OnMouseButtonReleased( EMouseButton mb, Vec2 p )
    {
        if( m_bIsPressed )
        {
            if( RecomputeHandlePos( p ) ) OnSliderChanged();
            m_bIsPressed = false; //\todo Hack because sometimes focus was lost but ispressed kept being true...
        }

        if( !IsPointInside(p) )
        {
            if( HasFocus() ) OnFocusLost(); //Explicitly lose focus when released outside
            return false;
        }
        else
            return true;
    }

    // Change value iff handle changes
    bool OnMouseMotion( Vec2 p )
    {
        if( m_bIsPressed )
        {
            if( RecomputeHandlePos( p ) ) OnSliderChanged();
            return true;
        }
        else
            return false;
        //return IsPointInside(p); //\todo this works too, but leaks
        //events to next controller if m_bIsPressed is true and the
        //mouse is outside the slider, which is confusing
    }
    
    //! Toggle color if focus lost
    void OnFocusLost()
    {
        m_bIsPressed = false;
        BaseControl::OnFocusLost();
    }

    bool Draw( gfx::IRenderer *p_renderer )
    {
        Vec2 pos_abs( GetPosAbs() );
        Vec2 natural_sizes( GetNaturalSizes() );
        Vec2 half_natural_sizes( 0.5f*natural_sizes );
        Vec2 half_sizes( 0.5f*GetSizes() );
        Vec2 hcenter_pos_abs( pos_abs[0], pos_abs[1] + half_sizes[1] );
        // Draw frame
        /*
        p_renderer->DrawAABB2( hcenter_pos_abs + Vec2(0,-half_natural_sizes[1]),
                               hcenter_pos_abs + Vec2(natural_sizes[0],half_natural_sizes[1]),
                               gfx::Style( m_HandleColorReleased, 1.0f, gfx::Style::eWire ) );
        */
        // Draw slider rail
        p_renderer->DrawSegment2( hcenter_pos_abs,
                                  hcenter_pos_abs + Vec2(natural_sizes[0],0),
                                  gfx::Style( m_HandleColorReleased, 1.0f, gfx::Style::eWire ) );
        // Draw slider ticks
        // Draw handle
        /* \todo Allow selecting handle style...
        p_renderer->DrawAABB2( hcenter_pos_abs + Vec2(m_HandlePos-m_HandleRadius,-m_HandleRadius),
                               hcenter_pos_abs + Vec2(m_HandlePos+m_HandleRadius,m_HandleRadius),
                               gfx::Style( (m_bIsPressed) ? m_HandleColorPressed : m_HandleColorReleased,
                                           1.0f,
                                           gfx::Style::eSolid ) );                                           
        */
        p_renderer->DrawDisk2( hcenter_pos_abs + Vec2(m_HandlePos,0), m_HandleRadius,
                               gfx::Style( (m_bIsPressed) ? m_HandleColorPressed : m_HandleColorReleased,
                                           1.0f,
                                           gfx::Style::eSolid ),
                               4 );
        return true;
    }

protected:
    // Recompute handle pos return true iff changed
    inline bool RecomputeHandlePos( const Vec2 &p )
        {
            //clamp pos to rectangle
            Vec2 pos_abs( GetPosAbs() );
            Vec2 clamped_pos( mal::Clamp( p, pos_abs, pos_abs+GetNaturalSizes() ) );
            return SetHandlePos( int32(clamped_pos[0] - pos_abs[0]) );
        }
    virtual void OnSliderChanged()
        {
            if( !m_OnChangeCall.empty() ) m_OnChangeCall(m_pUserData);
        }
    
protected:
    void *m_pUserData;
    D_OnChangeCall m_OnChangeCall;
    bool m_bIsPressed;
    int32 m_SliderSize;
    int32 m_HandlePos;
    int32 m_HandleRadius;
    gfx::Color m_HandleColorReleased;
    gfx::Color m_HandleColorPressed;
};

template <typename T>
class GNIR_SliderControl: public SliderControl
{
public:
    GNIR_SliderControl( uint32 slider_size,
                        GProperty_NumberInRange<T> &nir,
                        void *p_user_data = 0 )
    : SliderControl( slider_size, p_user_data )
    , m_rNIR(nir)
        {
            OnValueChanged();
        }

    void OnValueChanged()
        {
            float64 lambda = mal::Clamp01( float64(m_rNIR.m_Value-m_rNIR.m_Min) / (m_rNIR.m_Max - m_rNIR.m_Min) );
            SetHandlePos( int32(lambda * m_SliderSize) );
        }    
    
protected:
    virtual void OnSliderChanged()
        {
            float64 lambda = mal::Clamp01( float64(m_HandlePos) / m_SliderSize );
            m_rNIR.m_Value = m_rNIR.m_Min + T( lambda*(m_rNIR.m_Max-m_rNIR.m_Min) );
            SliderControl::OnSliderChanged();
        }
private:
    GProperty_NumberInRange<T> &m_rNIR;
};

template <typename T>
class GSliderControl: public SliderControl
{
public:
    GSliderControl( uint32 slider_size,
                    T &value, T min, T max,
                    void *p_user_data = 0 )
    : SliderControl( slider_size, p_user_data )
    , m_rValue(value), m_Min(min), m_Max(max)
        {
            //SFR_LOG_WARNING( "GSliderControl init %f, %f, %f", m_rValue, m_Min, m_Max );
            OnValueChanged();
        }

    void OnValueChanged()
        {
            float64 lambda = mal::Clamp01( float64(m_rValue-m_Min) / (m_Max - m_Min) );
            SetHandlePos( int32(lambda * m_SliderSize) );            
        }
    
protected:
    virtual void OnSliderChanged()
        {
            float64 lambda = mal::Clamp01( float64(m_HandlePos) / m_SliderSize );
            m_rValue = m_Min + T( lambda*(m_Max-m_Min) );
            SliderControl::OnSliderChanged();
        }
private:
    T &m_rValue;
    T m_Min;
    T m_Max;
};

}} //namespace sfr::gui

#endif //SFR_GUI_SLIDER_CONTROL_H
