#ifndef S2_BS_BSG_H
#define S2_BS_BSG_H

#include <Saphyre2/bs/Config.h>

#include <util/ItemStream.h>
#include <util/LogStream.h>

#include <Saphyre2/ds/IDarkSide.h>

#include <Geo/shape/ShapeLibrary.h>
#include <Geo/shape/ShapeFactory.h>
#include <Geo/ObjectFactory.h>

#include <string>
#include <map>

/*! \namespace S2 Saphyre2 BS Namespace.
  \brief This namespace collects all BS related classes and global functions.
*/
namespace S2 {

// Forward declarations
class Universe;

//! BS Entity types
enum EEntityType {
    eEntity_Universe = 0,
    eEntity_Kine2D,
    eEntity_Kine3D,
    eEntity_Particle2D,
    eEntity_Particle3D,
    eEntity_ParticleSys2D,
    eEntity_ParticleSys3D,
    eEntity_Fluid2D,
    eEntity_Solid2D,
    eEntity_Solid3D,

    eNumEntityTypes
};

//! BrightSide Global Services
class BSG
{
public:
    //! \name S2BS Init/Shutdown
    //@{
    static bool Init( const std::string &log_file ); //\todo NO direct output to file, just to a LogIS
    static bool IsInitialized();
    static void ShutDown();
    //@}

    static Universe *CreateUniverse();

    //! \name DS debug API
    //@{
    //\todo static const util::ItemStream *GetLogIS();
    static const util::ItemStream *GetVizIS();
    static const util::ItemStream *GetProfIS();
    static const util::ItemStream *GetStatIS();
    static const util::ItemStream *GetParamIS();

    static void UpdateLog();
    static void UpdateViz();
    static void UpdateProf();

    static void ClearLog();
    static void ClearViz();
    static void ClearProf();

    static void ClearStats();
    static void ClearParams();
    //@}

    static geo::ShapeLibrary &GetShapeLibrary() { return s_ShapeLib; }
    static geo::ShapeFactory &GetShapeFactory() { return s_ShapeFactory; }
    static geo::ObjectFactory &GetObjectFactory() { return s_ObjectFactory; }
    static util::LogStream &GetLogStream() { return s_LogStream; }

private:

    //! Processes control command results \todo UNUSED...
    static bool ProcessReturn();

private:
    //! \name BS static stuff
    //@{
    enum EStatus {
        eStatusNoInitialized,
        eStatusInitialized,
        eStatusShutDown,
        eStatusPanic
    };

    static EStatus s_Status;
    static ds::IDarkSide *s_pIDarkSide;
    static ds::Channel *s_pCC; //!< DS Control channel

    static geo::ShapeLibrary s_ShapeLib;
    static geo::ShapeFactory s_ShapeFactory;
    static geo::ObjectFactory s_ObjectFactory;
    //@}

    //! \name Debug Stuff
    //@{
    static util::LogStream s_LogStream;
    //@}

    //! \name Friend Classes (allowed to use BSG)
    //@{
    friend class Universe;
    //@}
};

} // namespace S2

// Some debugging macros
#define BS_ERROR( x ) BSG::GetLogStream() << util::BError() << x << util::EError()
#define BS_WARNING( x ) BSG::GetLogStream() << util::BWarning() << x << util::EWarning()
#define BS_INFO( x ) BSG::GetLogStream() << util::BInfo() << x << util::EInfo()
#define BS_BCMD( x ) BSG::GetLogStream() << util::BCmd(x)
#define BS_ECMD( b_result ) BSG::GetLogStream() << ((b_result)?"OK":"Error") << util::ECmd()

// Some useful macros
#ifdef PROFILE_FINAL
#  define BS_ASSERT( x )
#  define BS_LOG( x, ... )
#  define BS_LOG_WARNING( x, ... )
#  define BS_LOG_ERROR( x, ... )
#else
#  include <assert.h>
#  include <stdio.h>
#  define BS_ASSERT( x ) assert(x)
#  define BS_LOG( x, ... ) { printf("<BS LOG> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define BS_LOG_WARNING( x, ... ) { printf("<BS WARNING> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#  define BS_LOG_ERROR( x, ... ) { printf("<BS ERROR> "); printf( x, ##__VA_ARGS__ ); printf("\n"); }
#endif

#endif // S2_BS_BSG_H
