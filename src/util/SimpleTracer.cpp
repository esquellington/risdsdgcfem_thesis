#include "SimpleTracer.h"
#include <iostream>
#include <string>

namespace util
{

SimpleTracer* g_pTracer;

bool SimpleTracer::Init()
{
    m_SliceId = 0;
    m_NestedScopeLevel = 0;
    return true;
}

bool SimpleTracer::ShutDown() { return true; }
bool SimpleTracer::BeginSlice( float slice_time )
{
    m_SliceId++;
    m_SliceTime = slice_time;
    return true;
}
bool SimpleTracer::EndSlice( bool b_flush )
{
    return true;
}
bool SimpleTracer::BeginScope( const char* name )
{
    UTIL_ASSERT( m_NestedScopeLevel < cMaxNestedScope );
    m_NestedScopeNameStack[m_NestedScopeLevel++].Set( name );
    return true;
}
bool SimpleTracer::EndScope()
{
    UTIL_ASSERT( m_NestedScopeLevel > 0 );
    m_NestedScopeLevel--;
    return true;
}
bool SimpleTracer::TraceLocal( const char* name, const float& value )
{
    std::string nested_scope_name( "/" );
    for( unsigned int it_scope=0; it_scope<m_NestedScopeLevel; it_scope++ )
        nested_scope_name = nested_scope_name + std::string( m_NestedScopeNameStack[it_scope] ) + '/';
    std::cout << "S["<< m_SliceId << "]" << nested_scope_name << name << " = " << value << std::endl;
    return true;
}
bool SimpleTracer::TraceLocal( const char* name, const float* vec_values, uint32 num_values )
{
    std::string nested_scope_name( "/" );
    for( unsigned int it_scope=0; it_scope<m_NestedScopeLevel; it_scope++ )
        nested_scope_name = nested_scope_name + std::string( m_NestedScopeNameStack[it_scope] ) + '/';
    std::cout << "S["<< m_SliceId << "]" << nested_scope_name << name << " = ";
    for( unsigned int it_value=0; it_value<num_values; it_value++ )
        std::cout << vec_values[it_value] << " ";
    std::cout << std::endl;
    return true;
}
bool SimpleTracer::TraceGlobal( const char* name, const float& value )
{
    std::string nested_scope_name( "/" );
    for( unsigned int it_scope=0; it_scope<m_NestedScopeLevel; it_scope++ )
        nested_scope_name = nested_scope_name + std::string( m_NestedScopeNameStack[it_scope] ) + '/';
    std::cout << "S["<< m_SliceId << "]" << '/' << name << " = " << value << " (set in " << nested_scope_name << ")" << std::endl;
    return true;
}
bool SimpleTracer::TraceGlobal( const char* name, const float* vec_values, uint32 num_values )
{
    std::string nested_scope_name( "/" );
    for( unsigned int it_scope=0; it_scope<m_NestedScopeLevel; it_scope++ )
        nested_scope_name = nested_scope_name + std::string( m_NestedScopeNameStack[it_scope] ) + '/';
    std::cout << "S["<< m_SliceId << "]" << '/' << name << " = ";
    for( unsigned int it_value=0; it_value<num_values; it_value++ )
        std::cout << vec_values[it_value] << " ";
    std::cout << " (set in " << nested_scope_name << ")" << std::endl;
    return true;
}

}
