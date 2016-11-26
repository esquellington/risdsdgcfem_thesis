#ifndef S2_MS_PARTICLE_H
#define S2_MS_PARTICLE_H

#include "Config.h"
#include "IModel.h"

namespace S2 { namespace ms {

/*! Mathematical model of an N-dimensional Particle
   S(t) = [ pos vel ]
*/
template <unsigned N, EModelType EModelT>
class GParticle: public IModel
{
public:
    enum EConstants { cDimension = N };
    typedef mal::GVec<Real,N> vec_type;

public:
    GParticle()
    : m_AccForce(vec_type::Zero()), m_AccImpulse(vec_type::Zero())
    {}
    ~GParticle() {}

    void SetMass( Real mass ) { m_Mass = mass; m_InvMass = Real(1.0f)/mass; }

    //!\name Model-specific methods
    //@{
    finline const vec_type& GetPos() const { return m_Pos; }
    finline const vec_type& GetVel() const { return m_Vel; }
    finline void SetPos( const vec_type &pos ) { m_Pos = pos; }
    finline void SetVel( const vec_type &vel ) { m_Vel = vel; }
    finline void ApplyForce( const vec_type &force ) { m_AccForce += force; }
    finline void ApplyImpulse( const vec_type &impulse ) { m_AccImpulse += impulse; }
    finline void ApplyDisplacement( const vec_type &displacement ) { m_Pos += displacement; } //\todo could be made async with AccDispl, or subsumed in ApplyVariation()
    //@}

    //!\name Bookkeeping methods
    //@{
    EModelType GetModelType() const { return EModelT; }
    unsigned int GetDimension() const { return cDimension; }
    //@}

    //!\name Configuration and State management methods
    //@{
    unsigned int GetNumDOF() const { return cDimension; }
    void GetDOF( Real *p_dof ) const { m_Pos.ToArray(p_dof); }

    unsigned int GetStateSize() const { return 2*cDimension; }
    void GetState( Real *p_state ) const { m_Pos.ToArray( &p_state[0] );
                                           m_Vel.ToArray( &p_state[cDimension] ); }
    void SetState( const Real *p_state ) { m_Pos.FromArray( &p_state[0] );
                                           m_Vel.FromArray( &p_state[cDimension] ); }
    void Compute_dot_S( Real *p_dot_S ) const { m_Vel.ToArray( &p_dot_S[0] );
                                                (m_InvMass*m_AccForce).ToArray( &p_dot_S[cDimension] ); }
    //@}

    //!\name Dynamics methods
    //@{
    void ApplyForce( const Real *p_force ) { m_AccForce += vec_type(p_force); }
    void GetForce( Real *p_force ) const { m_AccForce.ToArray(p_force); }
    void ApplyImpulse( const Real *p_impulse ) { m_AccImpulse += vec_type(p_impulse); }
    void GetImpulse( Real *p_impulse ) const { m_AccImpulse.ToArray(p_impulse); }
    void ResetAccumulators() { m_AccForce = vec_type::Zero(); m_AccImpulse = vec_type::Zero(); }
    void Update( Real dt )
    {
        m_AccImpulse += dt*m_AccForce;
        m_Vel += m_InvMass*m_AccImpulse;
        m_Pos += dt*m_Vel;
        ResetAccumulators();
    }
    //@}

    //!\name Misc
    //@{
    Real GetMass() const { return m_Mass; }
    Real ComputeKineticEnergy() const { return Real(0.5f)*m_Mass*m_Vel.NormSq(); }
    Real ComputePotentialEnergy( const vec_type &gravity ) const { return -m_Mass*m_Pos*gravity; }
    //@}

private:
    Real m_Mass;
    Real m_InvMass;

    vec_type m_Pos;
    vec_type m_Vel;

    vec_type m_AccForce;
    vec_type m_AccImpulse;
};

typedef GParticle<2,eModel_Particle2D> Particle2;
typedef GParticle<3,eModel_Particle3D> Particle3;

}} // namespace S2::ms

#endif // S2_MS_PARTICLE_H_
