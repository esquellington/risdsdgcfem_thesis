#ifndef S2_DS_DSH_LEAF_PARTICLE_SYSTEM_H_
#define S2_DS_DSH_LEAF_PARTICLE_SYSTEM_H_

#include <Saphyre2/ds/dsh/GLeafDSH_Model.h>
#include <Saphyre2/ms/ParticleSystem.h>
#include <Mal/GConversion.h>
#include <util/ItemStream.h>

#include <Geo/shape/GSphereSetShape.h>
#include <Geo/util/Viz.h>

namespace S2 { namespace ds {

template <unsigned N, EDSHType DSHType, ms::EModelType ModelT>
class GLeafDSH_ParticleSystem: public GLeafDSH_Model< ms::GParticleSystem<N,ModelT> >
{
public:
    typedef GLeafDSH_Model< ms::GParticleSystem<N,ModelT> > base_type;

public:
    GLeafDSH_ParticleSystem( uint32 uid, IDynamicSystemHierarchy *p_parent )
    : base_type(uid, p_parent)
    {}

    EDSHType GetDSHType() const { return DSHType; }

    void SyncGO()
    {
        if( 0 != base_type::m_pGO )
        {
            // Perform SyncGO using only GObjectSDOF<N,GVec<N>> base class,
            // knowing that ANY bound GO must have proper dimension
            // and SDOF_T
            //\note Assumes that transform == identity, as all
            //information is in the SDOF
            base_type::m_pGO->SetTransform( base_type::transform_type::Identity() );
            base_type::m_Model.GetDOF( static_cast< geo::GObjectSDOF< N, mal::GVec<Real,N> >* >(base_type::m_pGO)->GetVecDOF_WriteOnly() );
        }
    }

    void RecomputeBV()
    {
        if( 0 != base_type::m_pGO && 0 != base_type::m_pBV )
            base_type::m_pGO->ComputeBV( *base_type::m_pBV );
    }

    bool Create( const ParamIt &pit )
    {
        // Init entity
        //\todo pit.Find("flags").Get<uint32>() => flags could say if
        //structure must be static (efficient) or variable (dynamic
        //particle alloc)
        unsigned int num_particles = pit.Find("count").Get<uint32>();
        base_type::m_Model.SetDscr( num_particles, pit.Find("mass").Get<Real>() );

        const typename base_type::model_type::vec_type *vec_pos =
            pit.Find("pos_i").GetArrayPtr< typename base_type::model_type::vec_type >();
        for( unsigned int i=0; i<num_particles; i++ )
        {
            base_type::m_Model.SetPos(i,vec_pos[i]);
            base_type::m_Model.SetVel(i,Vec2(0,0));
        }
        return true;
    }

    bool Edit( const ParamIt &pit )
    {
        ParamIt pit_vectp = pit.Find("vec_tp");
        if( !pit_vectp.IsValid() ) return false;

        const ds::GTouchedParticle<N> *vec_tp = pit_vectp.GetArrayPtr< ds::GTouchedParticle<N> >();
        unsigned int num_tp = pit_vectp.GetArrayCount();
        for( unsigned int it_tp=0; it_tp<num_tp; it_tp++ )
        {
            switch( vec_tp[it_tp].GetType() )
            {
            case ds::GTouchedParticle<N>::eTouchedPos:
                base_type::m_Model.SetPos( vec_tp[it_tp].GetParticleId(), vec_tp[it_tp].m_Vec );
                break;
            case ds::GTouchedParticle<N>::eTouchedVel:
                base_type::m_Model.SetVel( vec_tp[it_tp].GetParticleId(), vec_tp[it_tp].m_Vec );
                break;
            case ds::GTouchedParticle<N>::eTouchedForce:
                base_type::m_Model.ApplyForce( vec_tp[it_tp].GetParticleId(), vec_tp[it_tp].m_Vec );
                break;
            case ds::GTouchedParticle<N>::eTouchedImpulse:
                base_type::m_Model.ApplyImpulse( vec_tp[it_tp].GetParticleId(), vec_tp[it_tp].m_Vec );
                break;
            default: break;
            }
        }
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

        // State
        if( base_type::GetVizFlags().Test(eDbg_State) )
        {
            for( unsigned int i=0; i<base_type::m_Model.GetNumParticles(); i++ )
            {
                // Pos
                vs.BeginComplex(1,util::eType_VizPoint);
                {
                    vs.Write("pos", mal::CastDimension<3,float,N>(base_type::m_Model.GetPos(i)) );
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
                    vs.Write("pos", mal::CastDimension<3,float,N>(base_type::m_Model.GetPos(i)) );
                    vs.Write("vec", mal::CastDimension<3,float,N>(base_type::m_Model.GetVel(i)) );
                    vs.BeginComplex("style",util::eType_VizStyle);
                    {
                        vs.Write("color",Vec4f(0,1,0,1));
                        vs.Write("pen_size",2.0f);
                    }
                    vs.EndComplex();
                }
                vs.EndComplex();
            }
        }
        // BV
        if( base_type::GetVizFlags().Test(eDbg_DSH_BV) && 0 != base_type::GetBV() )
            geo::VizBoundingVolume(base_type::GetBV(),vs);
    }
};

//!\name Common explicit instantiations
//@{
typedef GLeafDSH_ParticleSystem<2,eDSH_Leaf_ParticleSystem2D,ms::eModel_ParticleSystem2D> LeafDSH_ParticleSystem2D;
typedef GLeafDSH_ParticleSystem<3,eDSH_Leaf_ParticleSystem3D,ms::eModel_ParticleSystem3D> LeafDSH_ParticleSystem3D;
//@}

} } // namespace S2::ds

#endif // S2_DS_DSH_LEAF_PARTICLE_SYSTEM_H_
