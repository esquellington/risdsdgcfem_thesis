#ifndef S2_DS_DSH_LEAF_PARTICLE_H_
#define S2_DS_DSH_LEAF_PARTICLE_H_

#include <Saphyre2/ds/dsh/GLeafDSH_Model.h>
#include <Saphyre2/ms/Particle.h>
#include <Mal/GConversion.h>
#include <util/ItemStream.h>
#include <util/VizStream.h>
#include <vector>

#include <Geo/bv/GAABB.h>
#include <Geo/bp/Pair.h>
#include <Geo/bp/Proxy.h>
#include <Geo/util/Viz.h>
#include <Geo/mp/TestContact.h>

#include <Saphyre2/ms/Particle_Connectors.h>

namespace S2 { namespace ds {

template <unsigned N, EDSHType DSHType, ms::EModelType ModelT>
class GLeafDSH_Particle: public GLeafDSH_Model< ms::GParticle<N,ModelT> >
{
public:
    typedef GLeafDSH_Model< ms::GParticle<N,ModelT> > base_type;
    
public:
    GLeafDSH_Particle( uint32 uid, IDynamicSystemHierarchy *p_parent )
    : base_type(uid, p_parent)
    {}

    EDSHType GetDSHType() const { return DSHType; }

    void SyncGO()
    {
        if( 0 != base_type::m_pGO )
            base_type::m_pGO->SetTransform(
                typename base_type::transform_type( base_type::m_Model.GetPos(),
                                                    base_type::transform_type::mat_type::Identity() ) );
    }
    
    void RecomputeBV()
    {        
        if( 0 != base_type::m_pGO )
            base_type::m_pGO->ComputeBV( *base_type::m_pBV );
        else //TEMPORAL...
            base_type::m_pBV->template As< geo::bv::GAABB<N> >().SetPosHalfSizes(
                base_type::m_Model.GetPos(),
                typename base_type::model_type::vec_type(1) );
    }

    
    bool Create( const ParamIt &pit )
    {
        // Init entity
        base_type::m_Model.SetMass( pit.Find("mass").Get<Real>() );
        base_type::m_Model.SetPos( pit.Find("pos").Get< typename base_type::model_type::vec_type >() );
        base_type::m_Model.SetVel( pit.Find("vel").Get< typename base_type::model_type::vec_type >() );        
        return true;
    }

    bool Edit( const ParamIt &pit )
    {
        ParamIt pit2 = pit.Find("pos");
        if( pit2.IsValid() )
            base_type::m_Model.SetPos( pit2.Get< typename base_type::model_type::vec_type >() );
        pit2 = pit.Find("vel");
        if( pit2.IsValid() )
            base_type::m_Model.SetVel( pit2.Get< typename base_type::model_type::vec_type >() );
        pit2 = pit.Find("force");
        if( pit2.IsValid() )
            base_type::m_Model.ApplyForce( pit2.Get< typename base_type::model_type::vec_type >() );
        pit2 = pit.Find("impulse");
        if( pit2.IsValid() )
            base_type::m_Model.ApplyImpulse( pit2.Get< typename base_type::model_type::vec_type >() );
        return true;
    }

    void DoViz( util::VizStream &vs ) const
    {
        // Shape
        if( base_type::GetVizFlags().Test(eDbg_DSH_Shape) && 0 != base_type::GetGO() )
            geo::VizObject(base_type::GetGO(),vs);

        // SS
        if( base_type::GetVizFlags().Test(eDbg_DSH_SimScheme) && 0 != base_type::GetSS() )
            base_type::GetSS()->DoViz(vs);

        /* State
        if( base_type::GetVizFlags().Test(eDbg_State) )
        {
            // Pos
            vs.BeginComplex(1,util::eType_VizPoint);
            {
                vs.Write("pos", mal::CastDimension<3,float,N>(base_type::m_Model.GetPos()) );
                vs.BeginComplex("style",util::eType_VizStyle);
                {
                    vs.Write("color",Vec4f(1,1,1,1));
                    vs.Write("pen_size",5.0f);
                }
                vs.EndComplex();
            }
            vs.EndComplex();
            // Vel
            vs.BeginComplex(1,util::eType_VizVec);
            {
                vs.Write("pos", mal::CastDimension<3,float,N>(base_type::m_Model.GetPos()) );
                vs.Write("vec", mal::CastDimension<3,float,N>(base_type::m_Model.GetVel()) );
                vs.BeginComplex("style",util::eType_VizStyle);
                {
                    vs.Write("color",Vec4f(0,1,0,1));
                    vs.Write("pen_size",2.0f);
                }
                vs.EndComplex();
            }
            vs.EndComplex();
        }
        */
        
        // BV
        /*
        if( base_type::GetVizFlags().Test(eDbg_DSH_BV) && 0 != base_type::GetBV() )
            geo::VizBoundingVolume(base_type::GetBV(),vs);
        */
    }
};

//!\name Common explicit instantiations
//@{
//typedef GLeafDSH_Particle<2,eDSH_Leaf_Particle2D,ms::eModel_Particle2D> LeafDSH_Particle2D;

//TEMPORAL: only exists because of NotifyPairBP() specialization,
//which may disappear when contacts are handled at SS level.
class LeafDSH_Particle2D: public GLeafDSH_Particle<2,eDSH_Leaf_Particle2D,ms::eModel_Particle2D>
{
public:
    typedef GLeafDSH_Particle<2,eDSH_Leaf_Particle2D,ms::eModel_Particle2D> base_type;

public:
    LeafDSH_Particle2D( uint32 uid, IDynamicSystemHierarchy *p_parent )
    : base_type(uid, p_parent)
    {}

    /*! Should check exact contact and create or persist a
     *  ContactConstraint that would be solved in
     *  SolveAcc/SolveVel/SolvePos, NOT here
     */
    void NotifyPairBP( const geo::bp::Proxy *p_other, const geo::bp::Pair *p_pair )
    {
        IEntity *pOtherEntity( p_other->GetObj<IEntity*>() );
        DS_ASSERT( eEntity_Geom == pOtherEntity->GetEntityType() );
        IGeom *pGeom( static_cast<IGeom*>( pOtherEntity ) );
        
        switch( p_pair->m_State )
        {
        case geo::bp::Pair::eNew:
        case geo::bp::Pair::ePersistent:
            {                
                // detect contact
                //\todo Create contact constraint, store in DSH, solve...
                
                //geo::np::ContactData2 cd;
                geo::np::ContactCache2 *p_cc(0);
                // \todo an IGeom2 should return a geo::IObject2, but an IGeom should return an geo::IObject, consider 2 different GetGO() methods..
                //if( geo::mp::TestContact( base_type::GetGO(), pGeom->GetGO(), m_CD, p_cc ) )
                if( geo::mp::TestContact( static_cast<const geo::IObject2*>(base_type::GetGO()),
                                          static_cast<const geo::IObject2*>(pGeom->GetGO()),
                                          m_CD, p_cc ) ) //TEMP: casts are ugly
                {
                    // SolveVel
                    Real vel_n( mal::Dot( base_type::m_Model.GetVel(), m_CD.m_AvgNormal ) );
                    if( vel_n < Real(0) )
                    {
                        Real coeff_restitution(0.5f);
                        base_type::m_Model.ApplyImpulse( -( base_type::m_Model.GetMass()
                                                            * (1+coeff_restitution)
                                                            * vel_n )
                                                         * m_CD.m_AvgNormal );
                        //base_type::m_Model.AccImpulseOnVel();
                    }
                    // SolvePos
                    base_type::m_Model.ApplyDisplacement( m_CD.m_AvgDepth*m_CD.m_AvgNormal );
                }                
            }
            break;
        case geo::bp::Pair::eVanished: break;
        default: break;            
        }
    }

    void DoViz( util::VizStream &vs ) const
    {
        base_type::DoViz( vs );
        geo::VizContactData( m_CD, vs );
    }

#ifdef __S2_DS_ENABLE_STATS
    void QueryStats( Flags32 flags, ReturnStream &rets ) const
        {
            rets.BeginComplex("stats", eType_Property_Object );
            {
                rets.BeginComplex( "<Energy>", eType_Property_Group );
                {
                    Real T( m_Model.ComputeKineticEnergy() );
                    Real V( m_Model.ComputePotentialEnergy( Vec2(0,-1) ) ); //TEMPORAL: Hardcoded gravity!!
                    rets.Write( "kinetic", T );
                    rets.Write( "potential", V );
                    rets.Write( "total", T+V );
                }
                rets.EndComplex();
            }
            rets.EndComplex();
        }
#endif //__S2_DS_ENABLE_STATS
    
private:
    geo::np::ContactData2 m_CD;
};

typedef GLeafDSH_Particle<3,eDSH_Leaf_Particle3D,ms::eModel_Particle3D> LeafDSH_Particle3D;
//@}
    
} } // namespace S2::ds

#endif // S2_DS_DSH_LEAF_PARTICLE_H_
