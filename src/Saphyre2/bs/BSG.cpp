#include <Saphyre2/bs/BSG.h>

#include <Saphyre2/ds/IDarkSide.h>
#include <Saphyre2/ds/DarkSideST.h>
#include <Saphyre2/ds/Commands.h>

#include <Saphyre2/bs/Universe.h>

#include <Geo/geo.h>

namespace S2 {

// API Static attributes
BSG::EStatus BSG::s_Status;
ds::IDarkSide *BSG::s_pIDarkSide = 0;
ds::Channel *BSG::s_pCC = 0;

geo::ShapeLibrary BSG::s_ShapeLib;
geo::ShapeFactory BSG::s_ShapeFactory(s_ShapeLib);
geo::ObjectFactory BSG::s_ObjectFactory(&s_ShapeFactory);
util::LogStream BSG::s_LogStream;

//-------- API Init/Shutdown

/*! API intialization method... just creates some internal and very
  secret things, so call it before using S2 and don't ask anything ok?
*/
bool BSG::Init( const std::string &log_file )
{
    // Check S2API status
    if( s_Status == eStatusInitialized
        || s_Status == eStatusPanic )
        return false;

    // Init libraries
    bool bOk = geo::Init();

    // Create Shape Libraries
    s_ShapeLib.Reserve( 1<<26, util::ItemStream::eRealloc_Identifiers ); //IMPORTANT: ShapeLib can realloc ids BUT NOT data or shape pointers will be invalidated, must be large enough, otherwise we'll assert

    // Create LogStream
    if( log_file != "" )
        s_LogStream.Open( log_file, "---------- S2::BS Log File ------------" );
    else
        s_LogStream.Open( "S2_BS.log.xml", "---------- S2::BS Log File ------------" );

    // Create S2DS and Connect
    s_LogStream.BeginSection("Connecting...");
      s_pIDarkSide = new S2::ds::DarkSideST;
      s_pIDarkSide->Init();
      s_pCC = s_pIDarkSide->Connect();
      if( 0 == s_pCC )
          s_LogStream << util::BError() << "Connect() failed" << util::EError();
      else if( !ProcessReturn() )
          s_LogStream << util::BError() << "Processing ReturnIS failed" << util::EError();
    s_LogStream.EndSection();

    s_Status = eStatusInitialized;

    return true;
}

bool BSG::IsInitialized()
{
    return s_Status == eStatusInitialized;
}

/*! API finishing method. Call it before exiting your app or you'll be
  cursed till the end of all days.
*/
void BSG::ShutDown()
{
    // Check S2API status
    if( s_Status != eStatusInitialized )
        return;

    // Send disconnect command
    s_LogStream.BeginSection("Disconnecting...");
    s_pIDarkSide->Disconnect();
    s_LogStream.EndSection();
    delete s_pIDarkSide;
    s_pIDarkSide = 0;
    s_LogStream.Close();

    // Shutdown libraries
    geo::ShutDown();
    s_ShapeLib.Clear();

    s_Status = eStatusShutDown;
}


Universe *BSG::CreateUniverse()
{
    Universe *pUni = new Universe();
    ds::Channel *pCommandChannel = s_pIDarkSide->OpenChannel();
    if( pCommandChannel )
    {
        pUni->SetChannel(pCommandChannel);
        return pUni;
    }
    else
        return 0;
}

//---- DS debug API
const util::ItemStream *BSG::GetVizIS()
{
    return s_pIDarkSide->GetVizIS();
}

const util::ItemStream *BSG::GetProfIS()
{
    return s_pIDarkSide->GetProfIS();
}

const util::ItemStream *BSG::GetStatIS()
{
    return s_pIDarkSide->GetStatIS();
}

const util::ItemStream *BSG::GetParamIS()
{
    return s_pIDarkSide->GetParamIS();
}

void BSG::UpdateLog() { s_pIDarkSide->UpdateLog(); }
void BSG::UpdateViz() { s_pIDarkSide->UpdateViz(); }
void BSG::UpdateProf() { s_pIDarkSide->UpdateProf(); }

void BSG::ClearLog() { s_pIDarkSide->ClearLog(); }
void BSG::ClearViz() { s_pIDarkSide->ClearViz(); }
void BSG::ClearProf() { s_pIDarkSide->ClearProf(); }

void BSG::ClearStats() { s_pIDarkSide->ClearStats(); }
void BSG::ClearParams() { s_pIDarkSide->ClearParams(); }

bool BSG::ProcessReturn()
{
    bool bResult(true);
    /*
    for( ds::ReturnIt it=s_pRETS->Begin();
         !bError && it.IsValid();
         ++it )
    {
        if( it.IsSimple() )
        {
            switch( it.GetType() )
            {
            case ds::eRet_Error: BS_ERROR("Result = Error"); break;
            case ds::eRet_Ok: BS_INFO("Result = OK"); break;
            default: BS_INFO("Result = Unknown"); break;
            }
        }
        else if( it.IsComplex() )
        {
            BS_INFO("Result = Complex");
            switch( it.GetType() )
            {
            case ds::eRet_Error_Data: BS_ERROR("Result = Error"); break;
            case ds::eRet_Ok_Data: BS_INFO("Result = OK"); break;
            case ds::eRet_NewEntity:
                BS_INFO("Result = NewEntity");
                ProcessNewEntity( it );
                break;
            case ds::eRet_Update:
                BS_INFO("Result = Update");
                ProcessUpdate( it );
                break;
            default: BS_INFO("Result = Unknown"); break;
            }
        }
        else
            bError = true;
    }
    s_pIDarkSide->DiscardReturn();
    */
    return bResult;
}

} // namespace S2
