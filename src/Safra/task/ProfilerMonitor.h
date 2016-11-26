#ifndef SFR_TASK_PROFILER_MONITOR_H
#define SFR_TASK_PROFILER_MONITOR_H

#include <Safra/core/ITask.h>
#include <Safra/core/IView.h>
#include <util/ItemStream.h>
#include <util/Profiler.h>
#include <stdio.h> //TEMPORAL: for sprintf...

namespace sfr
{

class ProfilerMonitor: public ITask
{
public:
    ProfilerMonitor( const util::ItemStream *p_pis, ITextView *p_view_text )
    : m_pPIS(p_pis)
    , m_pTextView(p_view_text)
    , m_AccTotalTime(0)
    , m_NumSamples(0)
    {
        m_pTextView->Clear();
    }
    ~ProfilerMonitor() {}

    void Run( unsigned int tick )
    {        
        m_pTextView->Clear();
        
        //Write entries from ProfileItemStream
        char str[1024]; // \todo STACK-SMASH again... static buffers are shit, solve it once (use sNprintf and a SINGLE static global g_sTmpStr
        for( util::ItemStream::ItemIt it = m_pPIS->Begin(); it.IsValid(); ++it )
        {            
            switch( it.GetType() )
            {
            case util::eTypeProfileEntry:
                {
                    if( it.IsArray() )
                    {
                        const util::Profiler::Entry *vec_entry = it.GetArrayPtr<util::Profiler::Entry>();
                        unsigned int num_entry = it.GetArrayCount();
                        double total_elapsed_time = (num_entry>0) ? vec_entry[0].m_ElapsedTime : 0;
                        
                        if( num_entry > 0 )
                        {
                            m_AccTotalTime += total_elapsed_time;
                            m_NumSamples++;
                        }
                        
                        for( unsigned int it_entry=0; it_entry<num_entry; it_entry++ )
                        {
                            int num_chars = 0;
                            
                            sprintf( str, "[%2d]", vec_entry[it_entry].m_Level );
                            num_chars += m_pTextView->Write( str );

                            sprintf( str, " %6.3fms %3.0f%% %1.0f%% ",
                                     vec_entry[it_entry].m_ElapsedTime * 1000.0, //miliseconds
                                     (vec_entry[it_entry].m_ElapsedTime/total_elapsed_time) * 100.0,
                                     (vec_entry[it_entry].m_UnaccountedTime/vec_entry[it_entry].m_ElapsedTime) * 100.0 );
                            num_chars += m_pTextView->Write( str );
                            
                            for( int i=0; i<vec_entry[it_entry].m_Level-1; i++ )
                                num_chars += m_pTextView->Space(4);

                            sprintf( str, "%3.0f%% %s",
                                     vec_entry[it_entry].m_ParentFraction * 100.0,
                                     vec_entry[it_entry].m_Name );
                            //num_chars += m_pTextView->Write( vec_entry[it_entry].m_Name );
                            num_chars += m_pTextView->Write( str );
                            
                            m_pTextView->Endl();
                        }
                    }
                    else
                    {
                        const util::Profiler::Entry &entry = it.Get<util::Profiler::Entry>();
                        sprintf(str,"[%2d] %s = %6.3f", entry.m_Level, entry.m_Name, entry.m_ElapsedTime );
                        for( int i=0; i<entry.m_Level; i++ ) m_pTextView->Tab();
                        m_pTextView->Write( str );
                        m_pTextView->Endl();
                    }
                }
                break;
            default:
                m_pTextView->Write("Unknown entry type!"); m_pTextView->Endl();
                break;
            }
        }

        // Total time avg
        sprintf( str, "Avg Total: %3.3fms", 1000.0*m_AccTotalTime/m_NumSamples );
        //\todo sprintf( str, "Avg Window (%d): %3.3fms", m_WindowSize, 1000.0*m_AccWindotTime/m_WindowSize );
        m_pTextView->Write( str );
        m_pTextView->Endl();
        
        //Update view (Render)
        m_pTextView->Update( 0 );
    }
    
private:
    const util::ItemStream *m_pPIS;
    ITextView *m_pTextView;

    double m_AccTotalTime;
    int m_NumSamples;
};

} // namespace sfr

#endif // SFR_TASK_PROFILER_MONITOR_H
