#include <Saphyre2/ds/DSG.h>

namespace S2 { namespace ds
{

double DSG::static_DebugTime = 0.0;
util::VizStream* DSG::static_pVizIS = NULL;
util::LogStream DSG::static_LogStream;
util::Profiler DSG::static_Profiler(100);
util::ItemStream* DSG::static_pProfIS = NULL;
util::ItemStream* DSG::static_pStatIS = NULL;
util::ItemStream* DSG::static_pParamIS = NULL;

//---- \name Init/Shutdown methods

/*! Initializes the internal services that all S2 classes can use. */
bool DSG::Init()
{
    static_LogStream.Open("S2_DS.log.xml", "---------- S2DarkSide Log File ------------" );
    return true;
}

/*! Shuts down the internal services that all S2 classes can use. */
void DSG::ShutDown()
{  
    // We close the Log last to allow reporting shutdown messages for
    // other internal services.
    static_LogStream.Close();
}

} } // namespace S2::ds
