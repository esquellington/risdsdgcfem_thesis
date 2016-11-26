#ifndef S2_MS_PARTICLE_SYSTEM_H
#define S2_MS_PARTICLE_SYSTEM_H

#include "Config.h"
#include "IModel.h"
#include <memory.h> //req by memset

namespace S2 {
namespace ms {

//! Dimension-agnostic particle system base
class IParticleSystem: public IModel
{
public:
    IParticleSystem() {}
    inline unsigned int GetNumParticles() const { return m_NumParticles; }
    //!\todo ALL mass-exclusive stuff could be here
protected:
    unsigned int m_NumParticles;
};

/*! Mathematical model of an N-dimensional Particle System
   S(t) = { [pos_0 vel_0] ... [pos_i vel_i]... }

   Uses Simplectic-Euler integration, which is the same as velocity-less Verlet for velocity-independent forces.

   \note Particles cannot be created or removed.
*/
template <unsigned N, EModelType EModelT>
class GParticleSystem: public IParticleSystem
{
public:
    enum EConstants { cDimension = N };
    typedef mal::GVec<Real,N> vec_type;

private:
    Real m_TotalMass;
    Real *m_vecInvMass;
    vec_type *m_vecPos;
    vec_type *m_vecVel;
    vec_type *m_vecAccForce;
    vec_type *m_vecAccImpulse;

public:
    GParticleSystem()
    : m_vecInvMass(0), m_vecPos(0), m_vecVel(0), m_vecAccForce(0), m_vecAccImpulse(0)
    {}
    bool SetDscr( unsigned int num_particles, Real mass )
    {
        MS_ASSERT(num_particles>0);
        m_NumParticles = num_particles;
        m_vecInvMass = new Real[num_particles];
        m_vecPos = new vec_type[num_particles];
        m_vecVel = new vec_type[num_particles];
        m_vecAccForce = new vec_type[num_particles];
        m_vecAccImpulse = new vec_type[num_particles];

        SetTotalMass(mass);
        ResetAccumulators();
        return true;
    }
    ~GParticleSystem()
    {
        if(m_vecInvMass) delete [] m_vecInvMass;
        if(m_vecPos) delete [] m_vecPos;
        if(m_vecVel) delete [] m_vecVel;
        if(m_vecAccForce) delete [] m_vecAccForce;
        if(m_vecAccImpulse) delete [] m_vecAccImpulse;
    }

    //!\name Model-specific methods
    //@{
    finline unsigned int GetNumParticles() const { return m_NumParticles; }
    finline Real GetMass( int pid ) const { return Real(1)/m_vecInvMass[pid]; }
    finline Real GetInvMass( int pid ) const { return m_vecInvMass[pid]; }
    finline const vec_type& GetPos( int pid ) const { return m_vecPos[pid]; }
    finline const vec_type& GetVel( int pid ) const { return m_vecVel[pid]; }
    finline const vec_type& GetAccForce( int pid ) const { return m_vecAccForce[pid]; }
    finline const vec_type& GetAccImpulse( int pid ) const { return m_vecAccImpulse[pid]; }
    finline void SetMass( int pid, Real mass ) { m_TotalMass += mass - GetMass(pid); m_vecInvMass[pid] = Real(1.0f)/mass; }
    finline void SetPos( int pid, const vec_type &pos ) { m_vecPos[pid] = pos; }
    finline void SetVel( int pid, const vec_type &vel ) { m_vecVel[pid] = vel; }
    finline void ApplyForce( int pid, const vec_type &force ) { m_vecAccForce[pid] += force; }
    finline void ApplyImpulse( int pid, const vec_type &impulse ) { m_vecAccImpulse[pid] += impulse; }
    //
    finline const vec_type* GetVecPos() const { return m_vecPos; }
    finline const vec_type* GetVecVel() const { return m_vecVel; }
    finline const Real* GetVecInvMass() const { return m_vecInvMass; }

    finline vec_type* GetVecPos_RW() { return m_vecPos; }
    finline vec_type* GetVecVel_RW() { return m_vecPos; }

    finline void ApplyForce( const vec_type *p_force ) { for( unsigned int i=0; i<m_NumParticles; i++ )
                                                             m_vecAccForce[i] += p_force[i]; }
    //@}

    //!\name Bookkeeping methods
    //@{
    unsigned int GetDimension() const { return cDimension; }
    EModelType GetModelType() const { return EModelT; }
    //@}

    //!\name Configuration and State management methods
    //@{
    unsigned int GetNumDOF() const { return cDimension*m_NumParticles; }
    void GetDOF( Real *p_dof ) const { memcpy( (char*)p_dof, &m_vecPos[0], m_NumParticles*sizeof(vec_type) ); }

    unsigned int GetStateSize() const { return 2*cDimension*m_NumParticles; }
    void GetState( Real *p_state ) const { memcpy( (char*)p_state, &m_vecPos[0], m_NumParticles*sizeof(vec_type) );
                                           memcpy( (char*)p_state + m_NumParticles*sizeof(vec_type), &m_vecVel[0],
                                                   m_NumParticles*sizeof(vec_type) ); }
    void SetState( const Real *p_state ) { memcpy( &m_vecPos[0], (char*)p_state, m_NumParticles*sizeof(Point2) );
                                           memcpy( &m_vecVel[0], (char*)p_state + m_NumParticles*sizeof(vec_type),
                                                   m_NumParticles*sizeof(vec_type) ); }
    void Compute_dot_S( Real *p_dot_S ) const { MS_ASSERT(false); }
    //@}

    //!\name Dynamics methods
    //@{
    void ApplyForce( const Real *p_force ) { for( unsigned int i=0; i<cDimension*m_NumParticles; i++ )
                                                 reinterpret_cast<Real*>(m_vecAccForce)[i] += p_force[i]; }
    void GetForce( Real *p_force ) const { memcpy( (char*)p_force, &m_vecAccForce[0], m_NumParticles*sizeof(vec_type) ); }
    void ApplyImpulse( const Real *p_impulse ) { /*m_AccImpulse += vec_type(p_impulse);*/ }
    void GetImpulse( Real *p_impulse ) const { /*m_AccImpulse.ToArray(p_impulse);*/ }
    void ResetAccumulators() { memset( (char*)m_vecAccForce, 0, m_NumParticles*sizeof(vec_type) );
                               memset( (char*)m_vecAccImpulse, 0, m_NumParticles*sizeof(vec_type) ); }
    void Update( Real dt )
    {
        //!\todo Handle frequent cases no-forces and/or no-impulses specifically (FASTER!!)
        for( unsigned int i=0; i<m_NumParticles; i++ )
        {
            m_vecVel[i] += m_vecInvMass[i]*(m_vecAccImpulse[i] + dt*m_vecAccForce[i]);
            m_vecPos[i] += dt*m_vecVel[i];
        }
        ResetAccumulators();
    }
    //@}

    //!\name Misc
    //@{
    Real GetMass() const { return m_TotalMass; }
    Real ComputeKineticEnergy() const { Real ke = 0;
                                        for( unsigned int i=0; i<m_NumParticles; i++ )
                                            ke += Real(0.5)*GetMass(i)*GetVel(i).NormSq();
                                        return ke; }
    Real ComputePotentialEnergy( const vec_type &gravity ) const { Real pe = 0;
                                                                   for( unsigned int i=0; i<m_NumParticles; i++ )
                                                                       pe -= GetMass(i)*GetPos(i)*gravity;
                                                                   return pe; }
    //@}

private:

    void SetTotalMass( Real mass )
    {
        m_TotalMass = mass;
        Real inv_particle_mass = Real(m_NumParticles)/mass;
        for( unsigned int i=0; i<m_NumParticles; i++ ) m_vecInvMass[i] = inv_particle_mass;
    }
};

typedef GParticleSystem<2,eModel_ParticleSystem2D> ParticleSystem2D;
typedef GParticleSystem<3,eModel_ParticleSystem3D> ParticleSystem3D;

}} // namespace S2::ms

#endif // S2_MS_PARTICLE_SYSTEM_H
