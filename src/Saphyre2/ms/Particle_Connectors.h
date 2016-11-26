#ifndef S2_MS_PARTICLE_CONNECTORS_H
#define S2_MS_PARTICLE_CONNECTORS_H

#include "Config.h"
#include "IConnector.h"
#include "Particle.h"

namespace S2 { namespace ms {

/*! Could be made dimension-generic as GParticle<D> is...

   \note Designed to be pool-allocated: Constructed once, Initialized
   many times.
*/
class ConnectorPointOnParticle2: public GIConnectorSA<Vec2,Vec2>
{
    
public:
    inline ConnectorPointOnParticle2() : m_pModel(0) {}
    inline void Init( Particle2 *p_model, const Vec2 &pos_local )
        {
            m_pModel = p_model;
            m_PosLocal = pos_local;
        }
    
    //! \name IConnector implementation
    //@{
    void Update( Real t, Real dt )
        {
            m_X.SR<Vec2>() = m_pModel->GetPos() + m_PosLocal;
            m_dX_dT.SR<Vec2>() = m_pModel->GetVel();
            //m_dX_dT_e.SR<Vec2>() = m_pModel->GetVel() after accumulating external impulses
            //m_d2X_dT2_e.SR<Vec2>() = m_pModel->GetAcc() after accumulating external forces
        }
    void ApplyForce( const Real *f ) {}
    void ApplyImpulse( const Real *j ) {}
    void ApplyVariation( const Real *d ) {}
    //@} 
    
private:
    Vec2 m_PosLocal;
    Particle2 *m_pModel;
};

}}

#endif // S2_MS_PARTICLE_CONNECTORS_H
