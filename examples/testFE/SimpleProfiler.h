#ifndef UTIL_SIMPLE_PROFILER_H
#define UTIL_SIMPLE_PROFILER_H

#include <util/Chronos.h>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <assert.h>

namespace util
{

//! Very simple profiler
/*!
 */
class SimpleProfiler
{
public:

    static const unsigned int cMaxEntries = 128;
    static const unsigned int cMaxNameLength = 64;

public:
    struct Entry
    {
        Entry() {}
        Entry( const char *name, int level )
        : m_Level(level)
        , m_UnaccountedTime(0)
        , m_ParentFraction(0)
            {
                strncpy( m_Name, name, cMaxNameLength-1 );
            }
        char m_Name[cMaxNameLength];
        int m_Level;
        double m_ElapsedTime;
        double m_UnaccountedTime;
        float m_ParentFraction;
    };

    struct StackEntry
    {
        StackEntry() {}
        StackEntry( double start_time, int index, int parent_index )
        : m_StartTime(start_time)
        , m_Index(index)
        , m_ParentIndex(parent_index)
        , m_NumChildren(0)
            {}
        double m_StartTime;
        int m_Index;       //!< Index in vecEntry
        int m_ParentIndex;
        int m_NumChildren;
    };

public:
    SimpleProfiler( unsigned int hint_max_entries = 0 )
    : m_Level(0)
        {
            if( hint_max_entries > 0 )
            {
                m_stackOpenEntry.reserve(hint_max_entries);
                m_vecEntry.reserve(hint_max_entries);
            }
            m_Clock.ResetTime();
        }

    ~SimpleProfiler() {}

    void Clear()
        {
            m_Level = 0;
            m_stackOpenEntry.clear();
            m_vecEntry.clear();
        }

    const std::vector<Entry> &GetEntries() const { return m_vecEntry; }

    void Begin( const char *name )
        {
            // Find new entry parent idx
            int parent_idx = (m_stackOpenEntry.size()>0) ? m_stackOpenEntry.back().m_Index : -1;
            // Add new entry to parent count
            if( -1 != parent_idx )
                m_stackOpenEntry.back().m_NumChildren++;
            // Add new entry to both stack and vec
            m_stackOpenEntry.push_back( StackEntry( m_Clock.GetTime(), m_vecEntry.size(), parent_idx ) );
            m_vecEntry.push_back( Entry(name, m_Level) );
            m_Level++;
        }

    void End()
        {
            assert( m_Level > 0 );

            const StackEntry &stack_entry = m_stackOpenEntry.back();
            Entry &entry = m_vecEntry[stack_entry.m_Index];

            // Compute elapsed time
            entry.m_ElapsedTime = m_Clock.GetTime() - stack_entry.m_StartTime;

            // Discount unaccounted time on parent entry. It will be negative until it ends.
            if( -1 != stack_entry.m_ParentIndex )
                m_vecEntry[stack_entry.m_ParentIndex].m_UnaccountedTime -= entry.m_ElapsedTime;

            // Correct unaccounted time if it has children, otherwise set to 0
            if( stack_entry.m_NumChildren > 0 )
                entry.m_UnaccountedTime += entry.m_ElapsedTime;
            else
                entry.m_UnaccountedTime = 0;

            // Compute fraction for all immediate children
            for( unsigned int it_children=stack_entry.m_Index; it_children<m_vecEntry.size(); it_children++ )
                if( m_vecEntry[it_children].m_Level == entry.m_Level+1 )
                    m_vecEntry[it_children].m_ParentFraction = m_vecEntry[it_children].m_ElapsedTime/entry.m_ElapsedTime;

            // Remove from open stack
            m_stackOpenEntry.pop_back();
            m_Level--;
        }

    void Print()
        {
            unsigned int num_entry = m_vecEntry.size();
            double total_elapsed_time = (num_entry>0) ? m_vecEntry[0].m_ElapsedTime : 0;
            printf( "----SimpleProfiler----\n" );
            for( unsigned int it_entry=0; it_entry<num_entry; it_entry++ )
            {
                printf( "[%2d] %6.3fms %3.0f%% %1.0f%% ",
                        m_vecEntry[it_entry].m_Level,
                        m_vecEntry[it_entry].m_ElapsedTime * 1000.0, //ms
                        (m_vecEntry[it_entry].m_ElapsedTime/total_elapsed_time) * 100.0,
                        (m_vecEntry[it_entry].m_UnaccountedTime/m_vecEntry[it_entry].m_ElapsedTime) * 100.0 );
                // Tabulation
                for( int i=0; i<m_vecEntry[it_entry].m_Level-1; i++ )
                    printf("    ");
                printf( "%3.0f%% %s \n",
                        m_vecEntry[it_entry].m_ParentFraction * 100.0,
                        m_vecEntry[it_entry].m_Name );
            }
        }

private:
    std::vector<StackEntry> m_stackOpenEntry;
    int m_Level;
    util::Chronos m_Clock;

    std::vector<Entry> m_vecEntry;
};

} //end namespace util

#endif // UTIL_SIMPLE_PROFILER_H
