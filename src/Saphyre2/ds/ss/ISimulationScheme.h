#ifndef S2_DS_SS_ISS_H
#define S2_DS_SS_ISS_H

#include "../Config.h"
#include <util/VizStream.h>

// Fwd declarations from Geo
namespace geo { namespace bp { class IBroadPhase; struct Proxy; class Pair; }}

namespace S2 { namespace ds {

// Fwd declarations
class IDynamicSystemHierarchy;
class IConnector;
class IConstraint;
class IInteraction;
class IGeom;

class ISimulationScheme
{
public:
    enum ESimulationSchemeProxyUserFlags
    {
        eSSPUF_None   = 0,
        eSSPUF_Entity = (1<<0)
    };

public:
    ISimulationScheme() {}
    virtual ~ISimulationScheme() {}

    virtual bool BindDSH( IDynamicSystemHierarchy* p_dsh ) = 0; //!< Returns false if incompatible
    virtual IDynamicSystemHierarchy* GetDSH() const = 0;

    virtual void Update( Real dt ) = 0;

    //! \name Bound DSH Notifications
    //@{
    virtual void NotifyAddChild( IDynamicSystemHierarchy* p_child ) = 0;
    virtual void NotifyAddConnector( IConnector* p_connector ) = 0;
    virtual void NotifyAddConstraint( IConstraint* p_constraint ) = 0;
    virtual void NotifyAddGeom( IGeom* p_geom ) = 0;

    virtual void NotifyRemoveChild( IDynamicSystemHierarchy* p_child ) = 0;
    virtual void NotifyRemoveConnector( IConnector* p_connector ) = 0;
    virtual void NotifyRemoveConstraint( IConstraint* p_constraint ) = 0;
    virtual void NotifyRemoveGeom( IGeom* p_geom ) = 0;

    virtual void NotifyPairBP( const geo::bp::Proxy* p_other, const geo::bp::Pair* p_pair ) = 0;
    //@}

    //! \name BP access for external queries
    //@{
    virtual const geo::bp::IBroadPhase* GetBP() const = 0;
    //@}

    //! \name Debug
    //@{
    //virtual void DoLog( util::ItemStream& ls ) const {}//= 0;
    virtual void DoViz( util::VizStream& vs ) const {}//= 0;
    //@}
};

} } // namespace S2::ds

#endif // S2_DS_SS_ISS_H
