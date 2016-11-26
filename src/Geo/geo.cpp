#include "geo.h"

#include <Geo/Config.h>
#include <Geo/mp/DoubleDispatcher.h>

#include "geo_params.h"
#include "geo_stats.h"

#if defined(__GEO_ENABLE_PARAMS) || defined(__GEO_ENABLE_STATS)
#  define __GEO_ENABLE_ARCHETYPES
#endif

#ifdef __GEO_ENABLE_ARCHETYPES
#  include <util/Archetype.h>
#endif

#ifdef __GEO_ENABLE_PARAMS
#  include "np/Context.h"
#  include <string>
#endif

namespace geo {

//\name Global Stuff
//@{
static bool gs_bIsInitialized(false);
#ifdef __GEO_ENABLE_ARCHETYPES
util::ArchetypeLibrary *gs_pArchetypeLibrary;
#endif
namespace mp { DoubleDispatcher *g_pDefaultDoubleDispatcher(0); }
namespace np { Context *g_pDefaultContext(0); }
//@}

bool Init()
{
    if( IsInialized() ) return true;
    bool bOk(true);
    if( !bOk ) GEO_LOG_ERROR( "geo::Init() error!" );
    gs_bIsInitialized = bOk;
    mp::g_pDefaultDoubleDispatcher = new mp::DoubleDispatcher();
    np::g_pDefaultContext = mp::g_pDefaultDoubleDispatcher->GetContext();

#ifdef __GEO_ENABLE_ARCHETYPES
    gs_pArchetypeLibrary = new util::ArchetypeLibrary();
#  ifdef __GEO_ENABLE_PARAMS
    np::Context::InitArchetype( *gs_pArchetypeLibrary );
#  endif
    Init_Stats();
#endif //__GEO_ENABLE_ARCHETYPES

    return gs_bIsInitialized;
}

bool ShutDown()
{
    if( !IsInialized() ) return false;
    bool bOk(true);
    if( !bOk ) GEO_LOG_ERROR( "geo::ShutDown() error!" );
    gs_bIsInitialized = false;
    delete mp::g_pDefaultDoubleDispatcher;
    mp::g_pDefaultDoubleDispatcher = 0;
    np::g_pDefaultContext = 0; //no delete, it's alloc by g_pDefaultDoubleDispatcher
#ifdef __GEO_ENABLE_ARCHETYPES
    delete gs_pArchetypeLibrary;
    gs_pArchetypeLibrary = 0;
#endif
    ShutDown_Stats();
    return bOk;
}

bool IsInialized() { return gs_bIsInitialized; }

util::ItemStream::ItemItRW QueryParams( util::ItemStream &params_is )
{
#ifdef __GEO_ENABLE_PARAMS
    params_is.BeginComplex( "geo::Params", eType_Property_Group );
    {
        //-- global
        //-- shape
        //-- bv
        //-- bp
        //...
        //-- mp
        //default DDT
        //-- np
        params_is.BeginComplex( "<NP>", eType_Property_Group );
        {
            //default Context
            gs_pArchetypeLibrary->ExportInstance( "Archetype_geo_np_Context", np::g_pDefaultContext, params_is );
        }
        params_is.EndComplex();
        //-- viz
        //DDF... etc...
    }
    return params_is.EndComplex();
#else
    return util::ItemStream::ItemItRW();
#endif
}

void SyncParams( util::ItemStream::ItemItRW params_it )
{
#ifdef __GEO_ENABLE_PARAMS
    GEO_ASSERT( std::string("geo::Params") == params_it.GetName() );
    //-- np
    // default Context
    gs_pArchetypeLibrary->SyncInstance( "Archetype_geo_np_Context",
                                        np::g_pDefaultContext,
                                        params_it.GetSubItem().Find("<NP>").GetSubItem().Find("default_geo_np_Context") );
#endif
}

} //namespace geo
