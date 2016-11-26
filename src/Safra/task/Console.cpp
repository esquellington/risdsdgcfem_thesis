#include "Console.h"

namespace sfr
{


//---- Console Implementation
Console::Console( ITextView *p_view_text )
: m_pTextView(p_view_text)
{
    // install our IKeyboardListener
    m_pTextView->SetKeyboardListener(this);

    // clear the view
    m_pTextView->Clear();
    m_pTextView->SetCursorVisible(true);
}

Console::~Console()
{
}

//! Extend Language Def with local commands
void Console::ExtendLanguageDef( skr::LanguageDef &language_def, const char *prefix )
{
    // Init Console language definition and Interpreter context
    language_def.DefineMethod( prefix, "quit","Quit console",&Method_Quit);
    language_def.DefineMethod( prefix, "exit","== quit",&Method_Quit);
    language_def.DefineMethod( prefix, "cls","Clear console",&Method_Cls);
    language_def.DefineMethod( prefix, "save","(TODO) Save script to a .skr file",&Method_Save);
}

//! Set Interpreter for console commands
void Console::SetInterpreter( skr::Interpreter *p_interpreter )
{
    m_pInterpreter = p_interpreter;
    m_pInterpreter->GetUserContext().Add("p_console",reinterpret_cast<void*>(this));
}

bool Console::OnKeyPressed( sfr::gui::EKey key, int x, int y )
{
    if( 0 == m_pInterpreter )
        return false;

    //std::cout << "Console " << key << std::endl;

    switch( key )
    {
        // Normal (=ASCII code)
    case sfr::gui::eKey_BackSpace: m_pTextView->BackSpace(); break;
    case sfr::gui::eKey_Del: m_pTextView->Del(); break;

        // Specials (>255)
    case sfr::gui::eKey_Up: m_pTextView->MoveCursor(-1,0); break;
    case sfr::gui::eKey_Down: m_pTextView->MoveCursor(1,0); break;
    case sfr::gui::eKey_Left: m_pTextView->MoveCursor(0,-1); break;
    case sfr::gui::eKey_Right: m_pTextView->MoveCursor(0,1); break;

        // Modifiers (>255)
        case sfr::gui::eKey_Control: m_pTextView->Write("CONTROL"); break;
        case sfr::gui::eKey_Alt: m_pTextView->Write("ALT"); break;
        case sfr::gui::eKey_Shift: m_pTextView->Write("SHIFT"); break;

        case sfr::gui::eKey_F1: m_pTextView->Write("F1"); break;
        case sfr::gui::eKey_F2: m_pTextView->Write("F2"); break;
        case sfr::gui::eKey_F3: m_pTextView->Write("F3"); break;
        case sfr::gui::eKey_F4: m_pTextView->Write("F4"); break;
        case sfr::gui::eKey_F5: m_pTextView->Write("F5"); break;
        case sfr::gui::eKey_F6: m_pTextView->Write("F6"); break;
        case sfr::gui::eKey_F7: m_pTextView->Write("F7"); break;
        case sfr::gui::eKey_F8: m_pTextView->Write("F8"); break;
        case sfr::gui::eKey_F9: m_pTextView->Write("F9"); break;
        case sfr::gui::eKey_F10: m_pTextView->Write("F10"); break;
        case sfr::gui::eKey_F11: m_pTextView->Write("F11"); break;
        case sfr::gui::eKey_F12: m_pTextView->Write("F12"); break;

    case sfr::gui::eKey_Return:
        if( !m_pInterpreter->Execute(m_pTextView->GetLine()) )
            std::cout << "Error while executing script" << std::endl;
        m_pTextView->Endl();
        break;
    case sfr::gui::eKey_Escape: m_pTextView->Clear(); break;
    default: m_pTextView->Write(key); break;
    }

    return true;
}

void Console::Run( unsigned int tick )
{
    //Update view (Render)
    m_pTextView->Update( 0 );
}


//---- SFR Console language definition
bool Console::Method_Quit( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
{
    exit(0);
    return true;
}
bool Console::Method_Cls( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
{
    reinterpret_cast<Console*>(user_context.Get<void*>("p_console"))->GetTextView()->Clear();
    return true;
}

bool Console::Method_Save( skr::VarSeq &user_context, skr::VarSeq &results, const skr::VarSeq &params )
{
    return false;
}


} // namespace sfr
