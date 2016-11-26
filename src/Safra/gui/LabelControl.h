#ifndef SFR_GUI_LABEL_CONTROL_H
#define SFR_GUI_LABEL_CONTROL_H

#include "BaseControl.h"
#include <string.h>

namespace sfr { namespace gui
{

/*! Static Label Control
  \note Changing text or fontsize yields a resize()
  \note If set, FieldLength clamps the displayed string.
*/
class LabelControl: public BaseControl
{
public:
    enum EConstants { cMaxLength = 64 };

public:
    LabelControl( const char *str = 0 )
    : m_TextLength(0)
    , m_FieldLength(0)
    , m_FontColor(0,0,0,1)
    , m_FontSizes( 8, 13 ) //( 9, 15 ) //\todo Allow setting font type...
        {
            if( 0 != str ) SetLabel(str);
            else m_Text[0] = 0;
            SetColor( gfx::Color(0,0,0,0) ); //Label has no background by default
        }
    ~LabelControl() {}

    inline void SetFontColor( const gfx::Color &color ) { m_FontColor = color; }
    inline const gfx::Color &GetFontColor() const { return m_FontColor; }
    void SetFontSizes( const Vec2 &font_sizes ) { m_FontSizes = font_sizes; SetSizes( GetNaturalSizes() ); }
    void SetFieldLength( int fixed_length )
        {
            m_FieldLength = mal::Clamp(fixed_length,0,cMaxLength-1);
            // Clamp string
            if( m_FieldLength > 0 && m_FieldLength < m_TextLength )
            {
                m_Text[m_FieldLength-1] = '\\';
                m_Text[m_FieldLength] = 0;
            }
            SetSizes( GetNaturalSizes() );
        }

    void SetLabel( const char *str, int length = 0 )
        {
            m_TextLength = length;
            if( m_TextLength == 0 ) m_TextLength = strlen(str);
            SFR_ASSERT( m_TextLength < cMaxLength );
            strcpy( &m_Text[0], str );
            // Clamp string
            if( m_FieldLength > 0 && m_FieldLength < m_TextLength )
            {
                m_Text[m_FieldLength-1] = '\\';
                m_Text[m_FieldLength] = 0;
            }
            SetSizes( GetNaturalSizes() );
        }
    inline const char *GetLabel() const { return &m_Text[0]; }
    inline int GetLength() const { return m_TextLength; }

    //!\name Specializations
    //@{
    Vec2 GetNaturalSizes() const
        {
            if( m_FieldLength > 0 )
                return Vec2( m_FieldLength*m_FontSizes.x(), m_FontSizes.y() );
            else
                return Vec2( m_TextLength*m_FontSizes.x(), m_FontSizes.y() );
        }
    bool Draw( gfx::IRenderer *p_renderer )
    {
        BaseControl::Draw(p_renderer);
        Vec2 pos( GetPosAbs()
                  + Vec2( 0,
                          m_FontSizes.y()-1) // Font offset +height-1 pixel
            );
        p_renderer->DrawLabel( Vec3(pos.x(), pos.y(), 0),
                               m_Text,
                               gfx::Style( m_FontColor, m_FontSizes.x() ) );
        return true;
    }
    //@}

private:
    int m_TextLength;
    int m_FieldLength;
    char m_Text[cMaxLength];

    gfx::Color m_FontColor;
    Vec2 m_FontSizes;
};

}} //namespace sfr::gfx

#endif //SFR_GUI_LABEL_CONTROL_H
