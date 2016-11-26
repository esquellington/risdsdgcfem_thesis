#ifndef S2_DS_I_DARK_SIDE_H
#define S2_DS_I_DARK_SIDE_H

#include <util/ItemStream.h>
#include "Commands.h"

namespace S2 { namespace ds {

typedef util::ItemStream::ItemIt ReturnIt;
class Channel;

//! Dark Side Interface
/*! Abstract interface to Saphyre's DS

Offers interfaces for connection/disconnection to S2DS and
request/answer queues access hiding all implementation details.

\note Subclasses may implement multithreading and network access to
the DS transparently.
*/
class IDarkSide
{
public:
    virtual ~IDarkSide() {}
    virtual bool Init() = 0;
    virtual void ShutDown() = 0;

    //! \name Control methods
    //@{
    virtual Channel *Connect() = 0;     //!< Returns the Control channel
    virtual void Disconnect() = 0;
    virtual Channel *OpenChannel() = 0; //!< Opens a new Command channel
    //@}

    //! \name Debug, Profiling and Parameter Tweaking
    //@{
    //\todo virtual const util::ItemStream *GetLogIS() const = 0;
    virtual const util::ItemStream *GetVizIS() const = 0;
    virtual const util::ItemStream *GetProfIS() const = 0;
    virtual const util::ItemStream *GetStatIS() const = 0;
    virtual const util::ItemStream *GetParamIS() const = 0;

    virtual void UpdateLog() = 0;
    virtual void UpdateViz() = 0;
    virtual void UpdateProf() = 0;

    virtual void ClearLog() = 0;
    virtual void ClearViz() = 0;
    virtual void ClearProf() = 0;

    virtual void ClearStats() = 0;
    virtual void ClearParams() = 0;
    //@}

protected:
    //! \name Command channel internal methods
    //@{
    friend class Channel;
    virtual bool FlushCommands( Channel *p_channel ) = 0;
    virtual void DiscardReturn( Channel *p_channel ) = 0;
    virtual void CloseChannel( Channel *p_channel ) = 0;
    //@}
};


/*! Command/Control channel proxy

  Simplifies API code by grouping Command and Return item streams
  and offering a simpler Write/Read API to them.
*/
class Channel
{
public:
    static const unsigned int cMaxBytesCmdsData = (1<<20);
    static const unsigned int cMaxBytesRetsData = (1<<20);
    static const unsigned int cMaxBytesStringPool = (1<<20);

public:
    Channel( IDarkSide *p_ds )
    : m_pDS(p_ds)
    , m_CMDS(cMaxBytesCmdsData,cMaxBytesStringPool)
    , m_RETS(cMaxBytesRetsData,cMaxBytesStringPool)
    , m_LastCommandId(0)
    {}
    ~Channel() {}

public:

    //!\name Commands and their params
    //@{
    inline CommandID BeginCommand( uint16 type )
    {
        m_LastCommandId = mal::Max(1,m_LastCommandId+1); //Valid CID are [1,2^31-1]
        m_CMDS.BeginComplex(m_LastCommandId,type);
        return m_LastCommandId;
    }
    inline void EndCommand() { m_CMDS.EndComplex(); }

    // Simple types (req item_type_id<T> trait and T=POD)
    template <typename T> inline void Write( int32 id, const T& data ) { m_CMDS.Write(id,data); }
    template <typename T> inline void Write( const char *name, const T& data ) { m_CMDS.Write(name,data); }

    // Pointer types (req item_type_id<T> trait and T=POD)
    template <typename T> inline void WritePtr( int32 id, const T* ptr ) { m_CMDS.WritePtr(id,ptr); }
    template <typename T> inline void WritePtr( const char *name, const T* ptr ) { m_CMDS.WritePtr(name,ptr); }

    // Deep-Copy an Item, possibly from another ItemStream
    inline void WriteItem( int32 id, const util::ItemStream::ItemIt &it ) { m_CMDS.WriteItem(id,it); }
    inline void WriteItem( const char *name, const util::ItemStream::ItemIt &it ) { m_CMDS.WriteItem(name,it); }

    // Simple type arrays (req item_type_id<T> trait and T=POD)
    template <typename T> inline void WriteArray( int32 id, const T* data, unsigned int num_elems )
    { m_CMDS.WriteArray(id,data,num_elems); }
    template <typename T> inline void WriteArray( const char *name, const T* data, unsigned int num_elems )
    { m_CMDS.WriteArray(name,data,num_elems); }

    // Complex data
    inline void BeginComplex( int32 id, uint16 type ) { m_CMDS.BeginComplex(id,type); }
    inline void BeginComplex( const char *name, uint16 type ) { m_CMDS.BeginComplex(name,type); }
    inline void EndComplex() { m_CMDS.EndComplex(); }
    //@}

    //! \todo Sync-wait for a specific CommandID return
    inline void WaitReturn( CommandID cid ) const {}
    inline ReturnIt GetReturnIt() const { return m_RETS.Begin(); }

    inline bool FlushCommands() { return m_pDS->FlushCommands(this); }
    inline void DiscardReturn() { m_pDS->DiscardReturn(this); }
    inline void Close() { m_pDS->CloseChannel(this); }

public:
    IDarkSide *m_pDS;
    util::ItemStream m_CMDS;
    util::ItemStream m_RETS;
    CommandID m_LastCommandId;
};


} } // namespace S2::ds

#endif // S2_DS_I_DARK_SIDE_H
