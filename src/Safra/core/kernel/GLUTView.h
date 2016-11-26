#ifndef SFR_CORE_KERNEL_GLUT_VIEW_H
#define SFR_CORE_KERNEL_GLUT_VIEW_H

#include <Safra/core/IView.h>
#if defined(WIN32)
#  include <windows.h>
#  include <GL/glut.h>
#else // defined(POSIX)
#  ifdef SAFRA_ENABLE_GLEW
//#warning "GLUTKernel with GLEW support enabled"
#    include <GL/glew.h>
#  endif
#  include <GL/glut.h>
#endif

#include <Safra/core/gfx/GLRenderer.h>
#include <Safra/core/gfx/GLCamera.h>

#define __ENABLE_SFR_VIEW_COMPUTE_FPS
#ifdef __ENABLE_SFR_VIEW_COMPUTE_FPS
#  include <stdio.h>
#  include <util/Chronos.h>
#endif

namespace sfr
{

/*! Implements IView using GL/GLUT, including mouse/keyboard */
class GLUTView: public IView
{
public:
    GLUTView( int width, int height, Flags32 flags )
    : m_ViewId(-1)
    , m_pCamera(0)
    , m_pRenderer(0)
    , m_pMouseListener(0)
    , m_pKeyboardListener(0)
    , m_pDesktop(0)
    , m_Flags( flags )
#ifdef __ENABLE_SFR_VIEW_COMPUTE_FPS
    , m_FilteredFPS(0)
#endif
    {
        m_pCamera = new gfx::GLCamera();
        m_pCamera->InitPerspective( Vec3(5.0,7.0,10.0),
                                    Vec3(0,0,0),
                                    Vec3(0,1,0),
                                    45.0, 1.0,
                                    1.0, 500.0,
                                    width, height );
    }

    ~GLUTView()
    {
        if( m_pCamera ) delete m_pCamera;
        if( m_pDesktop ) delete m_pDesktop;
    }

    // IUpdateTarget implementation
    virtual bool Update( double dt )
    {
        glutSetWindow( GetViewId() );
        glutPostRedisplay();
        return true;
    }

    //!\name IView Interface Implementation
    //@{
    virtual void BeginDraw()
    {
        glClearColor(1,1,1,1);
        //glClearColor(0.2,0.2,0.2,1);
        //glClearColor(0,0,0,1);
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glEnable( GL_DEPTH_TEST );
        glEnable( GL_BLEND );
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        glDisable(GL_LIGHTING);
    }

    virtual bool Draw()
    {
        if( GetCamera() )
        {
            GetCamera()->PrepareToShoot();
            if( GetRenderer() ) return GetRenderer()->Render();
        }
        return false;
    }

    virtual void EndDraw()
    {
        if( GetDesktop() && GetDesktop()->IsEnabled() ) GetDesktop()->Draw();
#ifdef __ENABLE_SFR_VIEW_COMPUTE_FPS
        if( m_Flags.Test(sfr::IView::eFlags_EnableFPS)
            && GetCamera()
            && GetRenderer() )
        {
            char str[128];
            float fps( 1.0f/m_Chronos.GetDt() );
            m_FilteredFPS = (m_FilteredFPS != 0)
                            ? 0.9f * m_FilteredFPS + 0.1f * fps
                            : fps;
            sprintf( str, "FPS %d", (int)m_FilteredFPS );
            GetCamera()->PrepareToShootViewport();
            glDisable(GL_DEPTH_TEST);
            int w,h;
            GetViewportSizes( w, h );
            GetRenderer()->DrawLabel2( sfr::Vec2(0.85f*w,0.05f*h),
                                       str, sfr::gfx::Style(sfr::gfx::Color(0,0,0,1),4) );
            glEnable(GL_DEPTH_TEST);
        }
#endif
    }

    void Reshape( int w, int h )
    {
        glViewport(0, 0, w, h);
        if( m_pCamera )
        {
            m_pCamera->SetViewportSizes( w , h );
            m_pCamera->SetAspect( (float) w / (float) h );
        }
    }

    void GetViewportSizes( int &w, int &h ) const
    {
        glutSetWindow( GetViewId() );
        w = glutGet( GLUT_WINDOW_WIDTH );
        h = glutGet( GLUT_WINDOW_HEIGHT );
    }
    //@}

    //!\name Camera Management
    //@{
    gfx::ICamera *GetCamera() { return m_pCamera; }
    void SetCamera( gfx::ICamera *p_camera ) { m_pCamera = p_camera; }
    //@}

    //!\name Renderer Management
    //@{
    virtual void SetRenderer( gfx::IRenderer *p_renderer ) { m_pRenderer = p_renderer; }
    gfx::IRenderer *GetRenderer() const { return m_pRenderer; }
    //@}

    //!\name GUI related methods
    //@{
    void SetDesktop( gui::IDesktop *p_desktop ) { m_pDesktop = p_desktop; }

    Flags32 GetFlags() const { return m_Flags; }
    void SetFlags( Flags32 flags ) { m_Flags = flags; }
    void SetMouseListener( gui::IMouseListener *p_mouse_listener ) { m_pMouseListener = p_mouse_listener; }
    void SetKeyboardListener( gui::IKeyboardListener *p_keyboard_listener ) { m_pKeyboardListener = p_keyboard_listener; }

    gui::IDesktop *GetDesktop() const { return m_pDesktop; }
    gui::IMouseListener *GetMouseListener() const { return m_pMouseListener; }
    gui::IKeyboardListener *GetKeyboardListener() const { return m_pKeyboardListener; }
    //}@

    //!\name GLUT-Specific Methods
    //@{
    void SetViewId( int view_id ) { m_ViewId = view_id; }
    int GetViewId() const { return m_ViewId; }
    //@}

private:
    int m_ViewId;
    gfx::ICamera *m_pCamera;
    gfx::IRenderer *m_pRenderer;
    gui::IMouseListener *m_pMouseListener;
    gui::IKeyboardListener *m_pKeyboardListener;
    gui::IDesktop *m_pDesktop;
    Flags32 m_Flags;

#ifdef __ENABLE_SFR_VIEW_COMPUTE_FPS
    util::Chronos m_Chronos;
    float m_FilteredFPS;
#endif

};

} //sfr

#include <string.h>

namespace sfr {
/*! ITextView implementation
*/
class GLUTTextView: public GLUTView, public ITextView
{
public:
    GLUTTextView( int width, int height,
                  Flags32 flags,
                  int text_rows, int text_columns )
    : GLUTView(width,height,flags)
    , m_FontHeight(24), m_FontWidth(24), m_FontType(GLUT_BITMAP_TIMES_ROMAN_24)
    , m_DefaultColor(0.5,0.5,0.5,1)
    , m_Rows(text_rows), m_Cols(text_columns)
    , m_CursorRow(0), m_CursorCol(0), m_bCursorVisible(false)
    , m_Text(0)
    {
        m_Text = new char[text_rows*text_columns];
        SetFontSize(13);
        Clear();
    }
    ~GLUTTextView() { if(m_Text) delete m_Text; }

    //!\name Camera Management
    //@{
    gfx::ICamera *GetCamera() { return 0; }
    void SetCamera( gfx::ICamera *p_camera ) { SFR_ASSERT(false); }
    //@}

    void SetFontSize( int size_in_pixels )
    {
        m_FontHeight = size_in_pixels;
        m_FontWidth = size_in_pixels;
        switch( size_in_pixels )
        {
        case 13: m_FontType = GLUT_BITMAP_8_BY_13; m_FontWidth = 8; break;
        case 15: m_FontType = GLUT_BITMAP_9_BY_15; m_FontWidth = 9; break;
        case 11: m_FontType = GLUT_BITMAP_TIMES_ROMAN_10; m_FontWidth = 10; break;
        case 24: m_FontType = GLUT_BITMAP_TIMES_ROMAN_24; m_FontWidth = 24; break;
        case 10: m_FontType = GLUT_BITMAP_HELVETICA_10; m_FontWidth = 10; break;
        case 12: m_FontType = GLUT_BITMAP_HELVETICA_12; m_FontWidth = 12; break;
        case 18: m_FontType = GLUT_BITMAP_HELVETICA_18; m_FontWidth = 18; break;
        default:
            m_FontType = GLUT_BITMAP_HELVETICA_18;
        }
    }
    void SetDefaultColor( const gfx::Color &color ) { m_DefaultColor = color; }

    void GetCursor( int &row, int &col ) const { row = m_CursorRow; col = m_CursorCol; }
    void SetCursor( int row, int col ) { m_CursorRow = row; m_CursorCol = col; }
    void MoveCursor( int rows, int cols )
    {
        m_CursorRow += rows;
        m_CursorRow = (m_CursorRow < 0) ? 0 : (m_CursorRow < m_Rows) ? m_CursorRow : m_Rows-1;
        m_CursorCol += cols;
        m_CursorCol = (m_CursorCol < 0) ? 0 : (m_CursorCol < m_Cols) ? m_CursorCol : m_Cols-1;
    }
    void SetCursorVisible( bool b_visible ) { m_bCursorVisible = b_visible; }

    char Read() const { return Read(m_CursorRow,m_CursorCol); }

    int Write( const char *text, int len )
    {
        //Write string ensuring that it fits AND that there is 1 empty
        //char at the end to write '\0'
        int num_free_chars = m_Cols - m_CursorCol - 1;
        int num_chars = (len < num_free_chars) ? len : num_free_chars;
        char *p_buf = &m_Text[m_CursorRow*m_Cols + m_CursorCol];
        for( int i=0; i<num_chars; i++, p_buf++ ) *p_buf = text[i];
        *p_buf = '\0'; //Ending mark that MAY be overwritten by next Write in the same row
        m_CursorCol += num_chars;
        return num_chars;
    }

    int Write( const char *text ) { return Write( text, strlen( text ) ); }
    int Write( char c ) { return Write(&c,1); }

    int Space( int count )
    {
        //Write string ensuring that it fits AND that there is 1 empty
        //char at the end to write '\0'
        int num_free_chars = m_Cols - m_CursorCol - 1;
        int num_chars = (count < num_free_chars) ? count : num_free_chars;
        char *p_buf = &m_Text[m_CursorRow*m_Cols + m_CursorCol];
        for( int i=0; i<num_chars; i++, p_buf++ ) *p_buf = ' ';
        *p_buf = '\0'; //Ending mark that MAY be overwritten by next Write in the same row
        m_CursorCol += num_chars;
        return num_chars;
    }
    int Tab() { return Space(4); }

    int BackSpace()
    {
        if( m_CursorCol > 0 )
        {
            m_CursorCol--;
            Del();
        }
        return 0;
    }

    int Del()
    {
        for( int i=m_CursorCol; i<m_Cols-1; i++ )
            Write( Read(m_CursorRow,i+1), m_CursorRow, i );
        return 0;
    }

    void Endl()
    {
        if( m_CursorCol < m_Cols )
        {
            m_Text[m_CursorRow*m_Cols + m_CursorCol] = '\0';
            m_CursorRow++;
            m_CursorCol = 0;
        }
        if( m_CursorRow == m_Rows )
        {
            m_CursorRow = 0;
            Write("****Wrapped****");
            Endl();
        }
    }

    void Clear()
    {
        m_CursorRow=0;
        m_CursorCol=0;
        memset(m_Text,'\0',m_Rows*m_Cols);
    }
    const char *GetText() const { return m_Text; }
    const char *GetLine() const { return &m_Text[m_CursorRow*m_Cols]; }

    //!\name IView Implementation
    //@{
    bool Draw()
    {
        if( GetRenderer() )
        {
            int w,h;
            GetViewportSizes(w,h);
            for( int i=0; i<m_CursorRow+1; i++ )
                DrawLine( -1.0f,
                          1.0f - 2.0f*( float(i+1)*float(m_FontHeight)/h ),
                          &m_Text[i*m_Cols],
                          m_DefaultColor );
            if( m_bCursorVisible )
                DrawChar( -1.0f + 2.0f*( float(m_CursorCol)*float(m_FontWidth)/w ),
                          1.0f - 2.0f*( float(m_CursorRow+1)*float(m_FontHeight)/h ),
                         '_',
                          gfx::Color(1,0,0,1) );
            return true;
        }
        return false;
    }

    void Reshape( int w, int h )
    {
        GLUTView::Reshape(w,h);
    }
    //@}

    //! \name Disambiguation
    //@{
    bool Update( double dt ) { return GLUTView::Update(dt); }
    void GetViewportSizes( int &w, int &h ) const { GLUTView::GetViewportSizes(w,h); }
    void SetRenderer( gfx::IRenderer *p_renderer ) { GLUTView::SetRenderer(p_renderer); }
    gfx::IRenderer *GetRenderer() const { return GLUTView::GetRenderer(); }
    Flags32 GetFlags() const { return GLUTView::GetFlags(); }
    void SetFlags( Flags32 flags ) { GLUTView::SetFlags(flags); }
    void SetMouseListener( gui::IMouseListener *p_mouse_listener ) { GLUTView::SetMouseListener(p_mouse_listener); }
    void SetKeyboardListener( gui::IKeyboardListener *p_keyboard_listener ) { GLUTView::SetKeyboardListener(p_keyboard_listener); }
    gui::IMouseListener *GetMouseListener() const { return GLUTView::GetMouseListener(); }
    gui::IKeyboardListener *GetKeyboardListener() const { return GLUTView::GetKeyboardListener(); }
    gui::IDesktop *GetDesktop() const { return GLUTView::GetDesktop(); }
    //}@

private:
    inline void Write( char c, int row, int col ) { m_Text[row*m_Cols + col] = c; }
    inline char Read( int row, int col ) const { return m_Text[row*m_Cols + col]; }

    void DrawLine( float x, float y, const char *text, const gfx::Color &color )
    {
        glColor3f(color[0],color[1],color[2]);
        glRasterPos2f( x, y );
        for( unsigned int i=0; i < strlen(text); i++ )
            glutBitmapCharacter( m_FontType, text[i] );
    }

    void DrawChar( float x, float y, char c, const gfx::Color &color )
    {
        glColor3f(color[0],color[1],color[2]);
        glRasterPos2f( x, y );
        glutBitmapCharacter( m_FontType, c );
    }

private:
    int m_FontHeight, m_FontWidth;
    void *m_FontType;
    gfx::Color m_DefaultColor;
    int m_Rows, m_Cols;
    int m_CursorRow, m_CursorCol;
    bool m_bCursorVisible;
    char *m_Text;
};


} // namespace sfr

#endif // SFR_CORE_KERNEL_GLUT_VIEW_H
