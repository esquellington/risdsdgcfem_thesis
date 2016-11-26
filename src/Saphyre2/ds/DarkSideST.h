#ifndef S2_DS_DARK_SIDE_ST_H
#define S2_DS_DARK_SIDE_ST_H

#include <Saphyre2/ds/IDarkSide.h>
#include <Saphyre2/ds/MoP.h>

#include <util/ItemStream.h>
#include <list>

namespace S2 { namespace ds {

//! Single Threaded DarkSide Implementation
class DarkSideST: public IDarkSide
{
public:
    DarkSideST()
    : m_ControlChannel(this)
    , m_VizIS((1<<22),(1<<15)) //1<<20 = 1mb
    , m_ProfIS((1<<20),(1<<15))
    , m_StatIS((1<<20),(1<<10))
    , m_ParamIS((1<<20),(1<<10))
    {}

    bool Init();
    void ShutDown();

    //! \name Control methods
    //@{
    Channel *Connect();
    void Disconnect();
    Channel *OpenChannel();
    //@}

    //! \name Debug
    //@{
    const util::ItemStream *GetVizIS() const { return &m_VizIS; }
    const util::ItemStream *GetProfIS() const { return &m_ProfIS; }
    const util::ItemStream *GetStatIS() const { return &m_StatIS; }
    const util::ItemStream *GetParamIS() const { return &m_ParamIS; }

    void UpdateLog() { m_MoP.UpdateLog(); }
    void UpdateViz() { m_MoP.UpdateViz(); }
    void UpdateProf() { m_MoP.UpdateProf(); }

    void ClearLog() { } //Mm_LogIS.Clear(); }
    void ClearViz() { m_VizIS.Clear(); }
    void ClearProf() { m_ProfIS.Clear(); }

    void ClearStats() { m_StatIS.Clear(); }
    void ClearParams() { m_ParamIS.Clear(); }
    //@}

protected:
    //! \name Command channel internal methods
    //@{
    bool FlushCommands( Channel *p_channel );
    void CloseChannel( Channel *p_channel );
    void DiscardReturn( Channel *p_channel );
    //@}

private:
    MoP m_MoP;
    Channel m_ControlChannel;
    typedef std::list<Channel*> ChannelContainer;
    ChannelContainer m_Channels;
    util::ItemStream m_VizIS;
    util::ItemStream m_ProfIS;
    util::ItemStream m_StatIS;
    util::ItemStream m_ParamIS;
};

} } // namespace S2::ds

#endif // S2_DS_DARK_SIDE_ST_H
