#ifndef SFR_CORE_KERNEL_GLUT_KERNEL_H
#define SFR_CORE_KERNEL_GLUT_KERNEL_H

#include <Safra/core/IKernel.h>
#include <list>
#include <vector>
#include <map>

namespace sfr
{

class GLUTView;

namespace gfx { class IRenderer; }

class GLUTKernel: public IKernel
{
public:
    GLUTKernel();
    ~GLUTKernel();

    bool AddTask( ITask *p_task );

    IView *CreateView( const char *name, Flags32 flags, int width, int height );
    ITextView *CreateTextView( const char *name, Flags32 flags,
                               int width, int height,
                               int text_rows, int text_columns );

    sfr::gfx::IRenderer *GetDefaultRenderer() const;

    void Run();
    void Kill();

private:
    void SetGLUTCallbacks();
    static void GLUT_CB_OnIdle(void);
    static void GLUT_CB_OnDisplay(void);
    static void GLUT_CB_OnKeyDown( unsigned char key, int x, int y );
    static void GLUT_CB_OnKeyUp( unsigned char key, int x, int y );
    static void GLUT_CB_OnSpecialDown( int key, int x, int y );
    static void GLUT_CB_OnSpecialUp( int key, int x, int y );
    static void GLUT_CB_OnReshape( int w, int h );

    static void GLUT_CB_OnMouse( int button, int state, int x, int y );
    static void GLUT_CB_OnMouseMotion( int x, int y );
    static void GLUT_CB_OnMousePassiveMotion( int x, int y );

    static GLUTView *GetViewByGLUTId( int glut_window_id );

private:

    static unsigned int m_LastTick;

    static std::list<ITask *> m_listTasks;

    static std::vector<GLUTView*> m_vecViews;
    static std::map<int,GLUTView*> m_mapWindowsToViews;
    static unsigned int m_LastCreateWindowEndX; //!< Used to avoid overlaping on creation

    static gfx::IRenderer *m_pDefaultRenderer;
};

} // namespace sfr

#endif // SFR_CORE_KERNEL_GLUT_KERNEL_H
