#ifndef S2_DS_CONFIG_H
#define S2_DS_CONFIG_H

#include <pla_types.h>
#include <Saphyre2/Config.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define DS_ASSERT( x )
#  define DS_LOG( x, ... )
#  define DS_LOG_WARNING( x, ... )
#  define DS_LOG_ERROR( x, ... )
#  define DS_LOG_ASSERT( c, x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define DS_ASSERT( x ) assert(x)
#  define DS_LOG( x, ... ) { printf("<DS LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define DS_LOG_WARNING( x, ... ) { printf("<DS WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define DS_LOG_ERROR( x, ... ) { printf("<DS ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define DS_LOG_ASSERT( c, x, ... ) { if(!(c)) { printf("<DS ASSERT FAILED> "); printf( x, ##__VA_ARGS__ ); printf("\n"); assert(c); } }
#endif

// Paranoid checks can be enabled/disabled independently from regular ones
#define __S2_DS_ENABLE_PARANOID
#ifdef __S2_DS_ENABLE_PARANOID
#  include <assert.h>
#  include <stdio.h>
#  define DS_PARANOID_ASSERT( x ) assert(x)
#else
#  define DS_PARANOID_ASSERT( x )
#endif

// Enable debug stuff
#define __S2_DS_ENABLE_STATS
#define __S2_DS_ENABLE_PARAMS

#ifdef __S2_DS_ENABLE_STATS
#  define S2_DS_STAT_SET( name, value ) { (name) = (value); }
#  define S2_DS_STAT_INC( name ) { (name)++; }
#  define S2_DS_STAT_DEC( name ) { (name)--; }
#  define S2_DS_STAT_ADD( name, value ) { (name) += (value); }
#else //__S2_DS_ENABLE_STATS
#  define S2_DS_STAT_SET( name, value )
#  define S2_DS_STAT_INC( name )
#  define S2_DS_STAT_DEC( name )
#  define S2_DS_STAT_ADD( name, value )
#endif //__S2_DS_ENABLE_STATS

#endif // S2_DS_CONFIG_H
