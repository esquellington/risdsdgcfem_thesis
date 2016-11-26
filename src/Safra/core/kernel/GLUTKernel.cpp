#define SAFRA_ENABLE_GLEW
#ifdef SAFRA_ENABLE_GLEW
#  include <GL/glew.h>
#endif
#include <Safra/core/kernel/GLUTKernel.h>
#include <Safra/core/kernel/GLUTView.h>
#include <Safra/gui/CameraControllerML.h>

#define __ENABLE_DESKTOP
#ifdef __ENABLE_DESKTOP
#  include <Safra/gui/Desktop.h>
#  include <iostream>
#endif

namespace sfr
{

std::list<ITask *> GLUTKernel::m_listTasks;
std::vector<GLUTView*> GLUTKernel::m_vecViews;
std::map<int,GLUTView*> GLUTKernel::m_mapWindowsToViews;
gfx::IRenderer *GLUTKernel::m_pDefaultRenderer = NULL;

unsigned int GLUTKernel::m_LastTick = 0;
unsigned int GLUTKernel::m_LastCreateWindowEndX = 64;

#define SFR_GLUT_CREATE_WINDOW_SPACING_X 10

GLUTKernel::GLUTKernel()
{
    // Ugly code required for glutInit to work
    int dummy_argc = 0;
    glutInit(&dummy_argc,0);

    // Initialise Renderer
    if( !m_pDefaultRenderer )
        m_pDefaultRenderer = new gfx::GLRenderer();
}

GLUTKernel::~GLUTKernel()
{
}

bool GLUTKernel::AddTask( ITask *p_task )
{
    m_listTasks.push_back( p_task );
    return true;
}

IView *GLUTKernel::CreateView( const char *name, Flags32 flags, int width, int height )
{
    // GLUT Window Init
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutInitWindowSize( width, height );
    glutInitWindowPosition( m_LastCreateWindowEndX, 0 );
    int window_id = glutCreateWindow(name);
    m_LastCreateWindowEndX += width + SFR_GLUT_CREATE_WINDOW_SPACING_X;
    glViewport( 0, 0, width, height );

    // GLUTView Init
    GLUTView *pView = new GLUTView( width, height, flags );
    pView->SetViewId( window_id );
    pView->SetRenderer( GetDefaultRenderer() );
    m_vecViews.push_back( pView );
    m_mapWindowsToViews[window_id] = pView;
    SetGLUTCallbacks();

#ifdef __ENABLE_DESKTOP
    if( flags.Test(sfr::IView::eFlags_EnableDesktop) )
    {
        pView->SetDesktop( new sfr::gui::Desktop(pView) );
        pView->GetDesktop()->Resize( width, height );
        /*TEMPORAL
          pView->GetDesktop()->CreateWidget( name, 0, Vec2(10,10) );
          std::cout << "Creating Widget " << name << std::endl;
        */
    }
#endif

    return pView;
}

ITextView *GLUTKernel::CreateTextView( const char *name, Flags32 flags,
                                       int width, int height,
                                       int text_rows, int text_columns )
{
    // GLUT Window Init
    glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE );
    glutInitWindowSize( width, height );
    glutInitWindowPosition( m_LastCreateWindowEndX, 0 );
    int window_id = glutCreateWindow(name);
    m_LastCreateWindowEndX += width + SFR_GLUT_CREATE_WINDOW_SPACING_X;
    glViewport( 0, 0, width, height );

    // GLUTView Init
    GLUTTextView *pTextView = new GLUTTextView( width, height, flags, text_rows, text_columns );
    pTextView->SetViewId( window_id );
    pTextView->SetRenderer( GetDefaultRenderer() );
    m_vecViews.push_back( pTextView );
    m_mapWindowsToViews[window_id] = pTextView;
    SetGLUTCallbacks();

    return pTextView;
}

sfr::gfx::IRenderer *GLUTKernel::GetDefaultRenderer() const
{
    return m_pDefaultRenderer;
}

void GLUTKernel::Run()
{
    // Check if any GLUT window has been created and, if not, create a
    // default one (GLUT crashes otherwise)
    if( 0 == m_vecViews.size() )
    {
        int window_id = glutCreateWindow("Default");
        m_vecViews.push_back( NULL );
        m_mapWindowsToViews[window_id] = NULL;
        SetGLUTCallbacks();
    }

#ifdef SAFRA_ENABLE_GLEW
    // Needs a valid GL context, created when the first glutCreateWindow is called.
    GLenum res = glewInit();
    if( GLEW_OK != res )
    {
        SFR_LOG_ERROR( "GLEW error: %s", glewGetErrorString(res) );
    }
    else
    {
        SFR_LOG( "GL version %s", glGetString(GL_VERSION) );
        SFR_LOG( "GLEW version %s", glewGetString(GLEW_VERSION) );
    }
#endif

    // Start GLUT Loop (never returns)
    glutMainLoop();
}

void GLUTKernel::Kill()
{
    exit(0);
}

//---- helpers
static gui::EKey TranslateSpecialToKey( int key )
{
    //std::cout << "TSK " << key << std::endl;

    // Check modifiers
    int modifiers = glutGetModifiers();
    if( 0 != (modifiers & GLUT_ACTIVE_CTRL) ) return gui::eKey_Control;
    else if( 0 != (modifiers & GLUT_ACTIVE_ALT) ) return gui::eKey_Alt;
    else if( 0 != (modifiers & GLUT_ACTIVE_SHIFT) ) return gui::eKey_Shift;

    // Check specials
    switch( key )
    {
    case GLUT_KEY_UP: return gui::eKey_Up; break;
    case GLUT_KEY_DOWN: return gui::eKey_Down; break;
    case GLUT_KEY_LEFT: return gui::eKey_Left; break;
    case GLUT_KEY_RIGHT: return gui::eKey_Right; break;
    case GLUT_KEY_F1: return gui::eKey_F1; break;
    case GLUT_KEY_F2: return gui::eKey_F2; break;
    case GLUT_KEY_F3: return gui::eKey_F3; break;
    case GLUT_KEY_F4: return gui::eKey_F4; break;
    case GLUT_KEY_F5: return gui::eKey_F5; break;
    case GLUT_KEY_F6: return gui::eKey_F6; break;
    case GLUT_KEY_F7: return gui::eKey_F7; break;
    case GLUT_KEY_F8: return gui::eKey_F8; break;
    case GLUT_KEY_F9: return gui::eKey_F9; break;
    case GLUT_KEY_F10: return gui::eKey_F10; break;
    case GLUT_KEY_F11: return gui::eKey_F11; break;
    case GLUT_KEY_F12: return gui::eKey_F12; break;
    default: break;
    }
    return gui::EKey(key);
}

//---- static GLUT Callbacks

// Register global GLUTKernel callbacks
void GLUTKernel::SetGLUTCallbacks()
{
    glutIdleFunc( GLUT_CB_OnIdle );
    glutDisplayFunc( GLUT_CB_OnDisplay );
    glutKeyboardFunc( GLUT_CB_OnKeyDown );
    glutKeyboardUpFunc( GLUT_CB_OnKeyUp );
    glutSpecialFunc( GLUT_CB_OnSpecialDown );
    glutSpecialUpFunc( GLUT_CB_OnSpecialUp );
    glutReshapeFunc( GLUT_CB_OnReshape );

    glutMouseFunc( GLUT_CB_OnMouse );
    glutMotionFunc( GLUT_CB_OnMouseMotion );
    glutPassiveMotionFunc( GLUT_CB_OnMousePassiveMotion );

    //glutVisibilityFunc(CallBackVisibilityFunc);
}

void GLUTKernel::GLUT_CB_OnIdle()
{
    for(std::list<ITask *>::iterator it_task = m_listTasks.begin();
        it_task != m_listTasks.end();
        ++it_task)
        (*it_task)->Run(m_LastTick);
    m_LastTick++;
}

void GLUTKernel::GLUT_CB_OnDisplay()
{
    IView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView )
    {
        pView->BeginDraw();
        bool bDraw = pView->Draw();
        pView->EndDraw();
        if( bDraw )
            glutSwapBuffers();
    }
}

void GLUTKernel::GLUT_CB_OnKeyDown( unsigned char key, int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView ) pView->OnKeyPressed( gui::EKey(key), x, y );
}

void GLUTKernel::GLUT_CB_OnKeyUp( unsigned char key, int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView ) pView->OnKeyReleased( gui::EKey(key), x, y );
}

void GLUTKernel::GLUT_CB_OnSpecialDown( int key, int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView ) pView->OnKeyPressed( TranslateSpecialToKey(key), x, y );
}

void GLUTKernel::GLUT_CB_OnSpecialUp( int key, int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView ) pView->OnKeyReleased( TranslateSpecialToKey(key), x, y );
}

void GLUTKernel::GLUT_CB_OnReshape( int w, int h )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView )
    {
        pView->Reshape(w,h);
        if( pView->GetDesktop() )
            pView->GetDesktop()->Resize(w,h);
    }
}

void GLUTKernel::GLUT_CB_OnMouse( int button, int state, int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView )
    {
        gui::EMouseButton mb;
        switch(button)
        {
        case GLUT_LEFT_BUTTON: mb = gui::eMB_Left; break;
        case GLUT_MIDDLE_BUTTON: mb = gui::eMB_Middle; break;
        case GLUT_RIGHT_BUTTON: mb = gui::eMB_Right; break;
        default: return; break;
        }
        if( state == GLUT_DOWN ) pView->OnMouseButtonPressed( mb, x, y );
        else pView->OnMouseButtonReleased( mb, x, y );
    }
}

void GLUTKernel::GLUT_CB_OnMouseMotion( int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView ) pView->OnMouseMotion(x,y);
}

void GLUTKernel::GLUT_CB_OnMousePassiveMotion( int x, int y )
{
    GLUTView *pView = GetViewByGLUTId( glutGetWindow() );
    if( pView ) pView->OnMouseMotion(x,y);
}

GLUTView *GLUTKernel::GetViewByGLUTId( int glut_window_id )
{
    std::map<int,GLUTView*>::iterator it = m_mapWindowsToViews.find( glut_window_id );
    if( it != m_mapWindowsToViews.end() )
        return it->second;
    else
        return NULL;
}


} // namespace sfr
