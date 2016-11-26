#ifndef S2_MS_IMODEL_H
#define S2_MS_IMODEL_H

#include "Config.h"

namespace S2 { namespace ms {

enum EModelType {
    eModel_Generic = 0,
    eModel_Particle3D,
    eModel_Particle2D,
    eModel_ParticleSystem3D,
    eModel_ParticleSystem2D,
    eModel_ParticleSystem_Verlet3D,
    eModel_ParticleSystem_Verlet2D,
    eModel_SPH_Fluid2D,
    eNumModelTypes
};

/*! Generic physics Model interface:
  - DOF q(t)
  - Mass
  - Optional:
    - dq_dt(t)
    - q(t-dt)
    - accumulate_f(t)
    - accumulate_j(t)
  - Kinetic Energy T(dq_dt(t))
*/
class IModel
{

public:
    IModel() {}
    virtual ~IModel() {}

    //!\name Bookkeeping methods
    //@{
    virtual EModelType GetModelType() const = 0;
    virtual unsigned int GetDimension() const = 0;
    //@}

    //!\name Configuration, DOF and State vector management methods
    //@{
    virtual unsigned int GetNumDOF() const = 0;
    virtual void GetDOF( Real *p_dof ) const = 0;

    virtual unsigned int GetStateSize() const = 0;
    virtual void GetState( Real *p_state ) const = 0;
    virtual void SetState( const Real *p_state ) = 0;

    virtual void Compute_dot_S( Real *p_dot_S ) const = 0;
    //@}

    //!\name Dynamics methods
    //@{
    virtual void ApplyForce( const Real *p_force ) = 0;
    virtual void GetForce( Real *p_force ) const = 0;
    virtual void ApplyImpulse( const Real *p_impulse ) = 0;
    virtual void GetImpulse( Real *p_impulse ) const = 0;
    virtual void ResetAccumulators() = 0;
    virtual void Update( Real dt ) = 0; //!< Simplest Integration
    //@}

    //!\name Misc
    //@{
    virtual Real GetMass() const = 0;
    virtual Real ComputeKineticEnergy() const { return 0; }
    //@}
};

}} // namespace S2::ms

#endif // S2_MS_IMODEL_H
