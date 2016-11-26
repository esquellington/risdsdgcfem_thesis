#ifndef GEO_GEO_STATS_H
#define GEO_GEO_STATS_H

#include <util/ItemStream.h>

namespace geo {

// Init/ShutDown stats subsystem
bool Init_Stats();
bool ShutDown_Stats();

// Begin/End stats collection
void BeginStats();
void EndStats();

// Query/Sync stats report
util::ItemStream::ItemItRW QueryStats( util::ItemStream &stats_is );
void SyncStats( util::ItemStream::ItemItRW stats_it );

}

#endif //GEO_GEO_STATS_H
