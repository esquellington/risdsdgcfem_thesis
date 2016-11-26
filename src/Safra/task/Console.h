#ifndef SFR_TASK_CONSOLE_H
#define SFR_TASK_CONSOLE_H

#include <Safra/core/ITask.h>
#include <Safra/core/IView.h>

#include <util/skr.h>

namespace sfr
{

/*! Input Console for a given skr::Language/Interpreter that uses an ITextView window.

  The console adds some default commands to the language and some
  default vars to the interpreter's user_context.
*/
class Console: public ITask, gui::IKeyboardListener
{
public:
    Console( ITextView *p_view_text );
    ~Console();

    ITextView *GetTextView() const { return m_pTextView; }

    //! Extend Language Def with local commands
    void ExtendLanguageDef( skr::LanguageDef &language_def, const char *prefix = "" );

    //! Set Interpreter for console commands
    void SetInterpreter( skr::Interpreter *p_interpreter );
    
    //! IKeyboardListener implementation
    bool OnKeyPressed( sfr::gui::EKey key, int x, int y );
    
    //! ITask implementation
    void Run( unsigned int tick );    
    
private:
    //! \name SFR Console language definition
    //@{
    static bool Method_Quit( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params );
    static bool Method_Cls( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params );
    static bool Method_Load( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params );
    static bool Method_Save( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params );
    //@}
    
private:
    ITextView *m_pTextView;
    skr::Interpreter *m_pInterpreter;    
};

} // namespace sfr

#endif // SFR_TASK_CONSOLE_H
