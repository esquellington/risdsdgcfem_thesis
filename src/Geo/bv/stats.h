#ifndef GEO_BV_STATS_H
#define GEO_BV_STATS_H

#include <Geo/Config.h>

#ifdef __GEO_ENABLE_STATS

namespace geo {
namespace bv {

struct Stats
{
    Stats() { Reset(); }
    inline void Reset() { for( unsigned int i=0; i<sizeof(Stats) ; i++ ) reinterpret_cast<char*>(this)[i] = 0; }
    struct TestOverlapStats
    {
        uint32 m_GAABB_GAABB;
        uint32 m_GSphere_GSphere;

        uint32 m_DOP2K4_DOP2K4;
        uint32 m_DOP2K8_DOP2K8;

        uint32 m_DOP3K6_DOP3K6;
        uint32 m_DOP3K14_DOP3K14;

        uint32 m_GBDOP_GBDOP;
        uint32 m_GBDOP_GBDOP_Axis;
    } testoverlap;
};
extern Stats *g_pStats;

}} //namespace geo::bv

#  define GEO_BV_STAT_SET( name, value ) { if(geo::bv::g_pStats) { (geo::bv::g_pStats->name) = (value); } }
#  define GEO_BV_STAT_INC( name ) { if(geo::bv::g_pStats) { (geo::bv::g_pStats->name)++; } }
#  define GEO_BV_STAT_DEC( name ) { if(geo::bv::g_pStats) { (geo::bv::g_pStats->name)--; } }
#  define GEO_BV_STAT_ADD( name, value ) { if(geo::bv::g_pStats) { (geo::bv::g_pStats->name) += (value); } }

#else //__GEO_ENABLE_STATS

#  define GEO_BV_STAT_SET( name, value )
#  define GEO_BV_STAT_INC( name )
#  define GEO_BV_STAT_DEC( name )
#  define GEO_BV_STAT_ADD( name, value )

#endif //__GEO_ENABLE_STATS

#endif // GEO_BV_STATS_H
