#include "StringUtils.h"
#include <string.h>

namespace util
{

bool IsLetter( char c )
{
    return ( '_' == c //underscore prefix
             || '%' == c //free param prefix
             || ('a' <= c && c <= 'z')
             || ('A' <= c && c <= 'Z') );
}

bool IsDigit( char c )
{
    return ('0' <= c && c <= '9');
}
bool IsNumericSymbol( char c )
{
    return IsCharInSet(c,"+-.");
}

bool IsCharInSet( char c, const char *set_char )
{
    if( 0 == set_char ) return false;
    const char *it_char=set_char;
    while( *it_char != '\0' )
    {
        if( *it_char != c )
            ++it_char;
        else
            return true;
    }
    return false;
}

//! '\n' is automatically skipped
const char *Skip( const char *str, char c )
{
    const char *it_char;
    for( it_char=str;
         *it_char != '\0' && (*it_char == c || *it_char == '\n');
         ++it_char );
    return (*it_char=='\0') ? 0 : it_char;
}
//! '\n' is automatically skipped
const char *Skip( const char *str, const char *skipped_chars )
{
    if( 0 == skipped_chars ) return str;
    const char *it_char;
    for( it_char=str;
         *it_char != '\0' && (IsCharInSet(*it_char,skipped_chars) || *it_char == '\n');
         ++it_char );
    return (*it_char=='\0') ? 0 : it_char;
}
//! '\n' is NOT automatically skipped
const char *SkipNotNL( const char *str, char c )
{
    const char *it_char;
    for( it_char=str;
         *it_char != '\0' && *it_char == c;
         ++it_char );
    return (*it_char=='\0') ? 0 : it_char;
}
//! '\n' is NOT automatically skipped
const char *SkipNotNL( const char *str, const char *skipped_chars )
{
    if( 0 == skipped_chars ) return str;
    const char *it_char;
    for( it_char=str;
         *it_char != '\0' && IsCharInSet(*it_char,skipped_chars);
         ++it_char );
    return (*it_char=='\0') ? 0 : it_char;
}
const char *Find( const char *str, char c )
{
    const char *it_char;
    for( it_char=str;
         *it_char != '\0' && *it_char != c;
         ++it_char );
    return (*it_char=='\0') ? 0 : it_char;
}
const char *Find( const char *str, const char *find_chars )
{
    return strpbrk(str,find_chars);
}

} //namespace util
