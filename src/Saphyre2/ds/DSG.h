#ifndef S2_DS_DSG_H
#define S2_DS_DSG_H

#include <util/LogStream.h>
#include <util/Profiler.h>
#include <util/ItemStream.h>
#include <util/VizStream.h>

namespace S2 { namespace ds {

/*! This class gives access to all services (logging, factories, etc...)
  that may be used anywhere in the S2DarkSide

  It is initialized/finished by the S2::DarkSideST
*/
class DSG
{
public:
    //! \name Init/Shutdown methods
    //@{
    static bool Init();
    static void ShutDown();
    //@}

    //! \name Specific service querying methods
    //@{
    static double GetDebugTime() { return static_DebugTime; }
    static util::VizStream *GetVizIS() { return static_pVizIS; }    
    static util::LogStream &GetLogStream() { return static_LogStream; }
    static util::Profiler &GetProfiler() { return static_Profiler; }
    static util::ItemStream *GetProfIS() { return static_pProfIS; }
    static util::ItemStream *GetStatIS() { return static_pStatIS; }
    static util::ItemStream *GetParamIS() { return static_pParamIS; }
    //@}

    //! \name Debug
    //@{
    static void SetDebugTime( double time ) { static_DebugTime = time; }
    static void SetVizIS( util::VizStream *p_vs ) { static_pVizIS = p_vs; }
    static void SetProfIS( util::ItemStream *p_pis ) { static_pProfIS = p_pis; }
    static void SetStatIS( util::ItemStream *p_pis ) { static_pStatIS = p_pis; }
    static void SetParamIS( util::ItemStream *p_pis ) { static_pParamIS = p_pis; }
    //@}
    
private:
    static double static_DebugTime;
    static util::VizStream *static_pVizIS;
    static util::LogStream static_LogStream;
    static util::Profiler static_Profiler;
    static util::ItemStream *static_pProfIS;
    static util::ItemStream *static_pStatIS;
    static util::ItemStream *static_pParamIS;
};

} } // namespace S2::ds

// Global logging and profiling macros
#define DS_ERROR( x ) ds::DSG::GetLogStream() << util::BError() << x << util::EError()
#define DS_WARNING( x ) ds::DSG::GetLogStream() << util::BWarning() << x << util::EWarning()
#define DS_INFO( x ) ds::DSG::GetLogStream() << util::BInfo() << x << util::EInfo()
#define DS_BCMD( x ) ds::DSG::GetLogStream() << util::BCmd(x)
#define DS_ECMD( b_result ) ds::DSG::GetLogStream() << ((b_result)?"OK":"Error") << util::ECmd()

#define DS_BPROF( x ) ds::DSG::GetProfiler().Begin(x)
#define DS_EPROF() ds::DSG::GetProfiler().End()

#endif // S2_DS_DSG_H
