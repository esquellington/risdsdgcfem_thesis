#include <Saphyre2/ds/DarkSideST.h>
#include <Saphyre2/ds/DSG.h>

namespace S2 { namespace ds
{

bool DarkSideST::Init()
{
    DSG::Init();
    DSG::SetVizIS( &m_VizIS );
    DSG::SetProfIS( &m_ProfIS );
    DSG::SetStatIS( &m_StatIS );
    DSG::SetParamIS( &m_ParamIS );

    m_MoP.Init();
    
    return true;
}

void DarkSideST::ShutDown()
{
    DSG::ShutDown();
    
    m_MoP.ShutDown();
}

//---- Control methods
Channel *DarkSideST::Connect()
{
    return &m_ControlChannel;
}

void DarkSideST::Disconnect()
{
}

Channel *DarkSideST::OpenChannel()
{
    Channel *pChannel = new Channel(this);
    m_Channels.push_back(pChannel);
    return pChannel;
}

//---- Command channel internal methods (Assume valid p_channel)
void DarkSideST::CloseChannel( Channel *p_channel )
{
    if( p_channel == &m_ControlChannel )
        Disconnect();
    else
    {
        m_Channels.remove(p_channel);
        delete p_channel;
    }
}

bool DarkSideST::FlushCommands( Channel *p_channel )
{
    if( m_MoP.Execute(p_channel->m_CMDS,p_channel->m_RETS) )
    {
        p_channel->m_CMDS.Clear();
        return true;
    }
    else
    {
        DS_INFO("[ABORTED] due to error in MoP::Execute()");
        DS_INFO("[TODO] dump Channel cmd/ret streams on LogStream");
        DS_ASSERT(false);
        return false;
    }
}

void DarkSideST::DiscardReturn( Channel *p_channel )
{
    p_channel->m_RETS.Clear();
}

} } // namespace S2::ds
