#ifndef S2_MS_CONFIG_H
#define S2_MS_CONFIG_H

#include <pla_types.h>
#include <Saphyre2/Config.h>

// Some useful macros
#ifdef PROFILE_FINAL
#  define MS_ASSERT( x )
#  define MS_LOG( x, ... )
#  define MS_LOG_WARNING( x, ... )
#  define MS_LOG_ERROR( x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define MS_ASSERT( x ) assert(x)
#  define MS_LOG( x, ... ) { printf("<MS LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define MS_LOG_WARNING( x, ... ) { printf("<MS WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define MS_LOG_ERROR( x, ... ) { printf("<MS ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#endif

#define __S2_MS_ENABLE_STATS

#ifdef __S2_MS_ENABLE_STATS
#  define S2_MS_STAT_SET( name, value ) { (name) = (value); }
#  define S2_MS_STAT_INC( name ) { (name)++; }
#  define S2_MS_STAT_DEC( name ) { (name)--; }
#  define S2_MS_STAT_ADD( name, value ) { (name) += (value); }
#else //__S2_MS_ENABLE_STATS
#  define S2_MS_STAT_SET( name, value )
#  define S2_MS_STAT_INC( name )
#  define S2_MS_STAT_DEC( name )
#  define S2_MS_STAT_ADD( name, value )
#endif //__S2_MS_ENABLE_STATS

#endif // S2_MS_CONFIG_H
