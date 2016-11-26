#ifndef SFR_STRING_UTILS_H
#define SFR_STRING_UTILS_H

namespace util
{

bool IsLetter( char c );
bool IsDigit( char c );
bool IsNumericSymbol( char c );

bool IsCharInSet( char c, const char *set_char );
const char *Skip( const char *str, char c ); //skips '\n'
const char *Skip( const char *str, const char *skipped_chars ); //skips '\n'
const char *SkipNotNL( const char *str, char c ); //skips '\n'
const char *SkipNotNL( const char *str, const char *skipped_chars );
const char *Find( const char *str, char c );
const char *Find( const char *str, const char *find_chars );

} //namespace util

#endif //SFR_STRING_UTILS_H
