#ifndef S2_DS_MOP_H
#define S2_DS_MOP_H

#include <Saphyre2/ds/Commands.h>
#include <util/ItemStream.h>

namespace S2 { namespace ds {

// Fwd declarations
class IDynamicSystemHierarchy;
class GeomFactory;
    
//! Saphyre2 MoP
/*! MoP of the S2 Dark Side... receives requests from the Bright Side
  and orchestrates the Dark Side to fulfill them.
*/
class MoP
{
public:
    typedef util::ItemStream CommandStream;
    typedef util::ItemStream ReturnStream;
    typedef util::ItemStream::ItemIt CommandIt;
    typedef util::ItemStream::ItemIt ParamIt;
    typedef util::ItemStream::ItemIt ReturnIt;    
    
public:
    MoP();
    ~MoP();
    
    //! \name Public API for IDarkSide
    //@{
    void Init();
    void ShutDown();
    bool Execute( const CommandStream &cmds, ReturnStream &rets );    
    //@}

    //! \name Debug API for IDarkSide
    //@{
    void UpdateLog();
    void UpdateViz();
    void UpdateProf();
    //@}
    
private:
    bool Execute( const CommandIt &cmdit, ReturnStream &rets );

    //! \name Control commands \todo useless... consider removing them
    //@{
    /*
    bool Ctrl_Log( CommandID cid, const ParamIt &pit, ReturnStream &rets );      //!< Enable/Disable/Config
    bool Ctrl_Viz( CommandID cid, const ParamIt &pit, ReturnStream &rets );      //!< Enable/Disable/Config
    bool Ctrl_Prof( CommandID cid, const ParamIt &pit, ReturnStream &rets );     //!< Enable/Disable/Config
    */
    //@}
    
    //! \name Entity commands
    //@{
    bool Cmd_Create( CommandID cid, const ParamIt &pit, ReturnStream &rets );
    bool Cmd_Destroy( CommandID cid, const ParamIt &pit, ReturnStream &rets );
    bool Cmd_Edit( CommandID cid, const ParamIt &pit, ReturnStream &rets );
    bool Cmd_Update( CommandID cid, const ParamIt &pit, ReturnStream &rets );   
    bool Cmd_Internal( CommandID cid, const ParamIt &pit, ReturnStream &rets );   
    //@}

    //! \name Entity queries
    //@{    
    bool Query_RayCast( CommandID cid, const ParamIt &pit, ReturnStream &rets );
    //bool Query_Intersection( CommandID cid, const ParamIt &pit, ReturnStream &rets );
    //@}
    
    //! \name Entity debug commands
    //@{
    //! cmd_type = Log/Viz/Profile
    bool Dbg_Conf( CommandID cid, const ParamIt &pit, ReturnStream &rets );
    //! Stats
    bool Dbg_QueryStats( CommandID cid, const ParamIt &pit );
    //! Params
    bool Dbg_QueryParams( CommandID cid, const ParamIt &pit );
    bool Dbg_SyncParams( CommandID cid, const ParamIt &pit );
    //@}
    
private:
    //! \name Multiverse Structure related attributes
    //@{    
    IDynamicSystemHierarchy *m_pRootDSH; //!< Multiverse root DSH node
    GeomFactory *m_pGeomFactory;
    //@}
};

} } // namespace S2::ds

#endif // S2_DS_MOP_H
