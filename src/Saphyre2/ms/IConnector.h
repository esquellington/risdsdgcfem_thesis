#ifndef S2_MS_ICONNECTOR_H
#define S2_MS_ICONNECTOR_H

#include "Config.h"
#include <Mal/GSRV.h>

namespace S2 { namespace ms {

/*! Low-level Connector Interface based on \ref INTDYN90

  Connectors can be defined on ANYTHING that has DOF
  - Kinematic Geometry
  - DSH
  - Models

  Computes:
  - x(q(t)) -> X()
  - \dot x -> dX_dT()        Current Vel
  - \dot x_e -> dX_dT_e()    Vel induced by accumulated external impulses
  - \ddot x^e -> d2X_dT2_e() Acc induced by accumulated external forces
  - \partialderiv{x}{q} -> pX_pQ()
  - (Internally) \partialderiv{x}{t}, \partialderivtwo{x}{q},
    \partialderivtwo{x}{t}, \partialderivmix{x}{t}{q},

  - All methods returning const SRV/SRM & vectors/matrices.
    - SRV/M may encode fixed-size constants (ZeroN, IdentityNxM,....)

  See \ref S2CONSTRAINTSDERIVATION for a complete derivation.
*/
class IConnector
{
public:
    /*! Connector type, according to it's host object
    enum EType {
        eStatic,          //!< Only X() returns non-null
        eKinematic,       //!< ApplyF/I/V does nothing
        eDynamic
    };

    enum EHostType {
        eHost_None,
        eHost_IGeom,
        eHost_IDSH
    };
    */

public:
    IConnector() {}
    virtual ~IConnector() {}

    //! \name Connector Interface
    //@{
    virtual unsigned int GetNumDOF() const = 0;
    virtual void Update( Real t, Real dt ) = 0;

    virtual const SRV& X() const = 0;
    virtual const SRV& dX_dT() const = 0;
    virtual const SRV& dX_dT_e() const = 0;
    virtual const SRV& d2X_dT2_e() const = 0;
    //virtual const SRM& pX_pQ() const = 0;

    virtual void ApplyForce( const Real *f ) = 0;
    virtual void ApplyImpulse( const Real *j ) = 0;
    virtual void ApplyVariation( const Real *d ) = 0;
    //@}

    /* Host-related methods
      virtual EHostType GetHostType() = 0;
      virtual void *GetHost() const = 0;
      template< typename HostT > HostT *GetHostAs() const { if(correct_type) return reinterpret_cast<HostT*>(GetHost()); else return 0; }
    */
};

/*! Generic static allocation of Connectors with compile-time defined
  input/output SR types.

  \todo GIConnectorDA could do the same, but with dynamic allocation
  to support > 1 Count... for host objects with a variable number of
  sdof (ex: a whole particle system)
*/
template< typename OutputSRT, typename InputSRT >
class GIConnectorSA: public IConnector
{
public:
    typedef OutputSRT output_sr_type;
    typedef InputSRT input_sr_type;
    typedef mal::GMat< Real,
                       OutputSRT::size_in_reals,
                       InputSRT::size_in_reals > jacobian_mat_type;

public:
    GIConnectorSA() {}
    ~GIConnectorSA() {}

    //!\name Low-level access to connector magnitudes
    //@{
    finline const output_sr_type &A() const { return m_X.template SR<output_sr_type>(); }
    finline const output_sr_type &dA_dT() const { return m_dX_dT.template SR<output_sr_type>(); }
    finline const output_sr_type &dA_dT_e() const { return m_dX_dT_e.template SR<output_sr_type>(); }
    finline const output_sr_type &d2A_dT2_e() const { return m_d2X_dT2_e.template SR<output_sr_type>(); }
    //\todo finline const jacobian_mat_type &pA_pQ() const { return m_pX_pQ.template SR<jacobian_mat_type>(); }
    //@}

    //! \name IConnector implementation
    //@{
    unsigned int GetNumDOF() const { return output_sr_type::size_in_reals; }
    const SRV& X() const { return m_X; }
    const SRV& dX_dT() const { return m_dX_dT; }
    const SRV& dX_dT_e() const { return m_dX_dT_e; }
    const SRV& d2X_dT2_e() const { return m_d2X_dT2_e; }
    //const SRM& pX_pQ() const { return m_pX_pQ; }
    //@}

protected:
    //!\name Statically allocated magnitudes
    //@{
    mal::GSRVA< output_sr_type, 1 > m_X;
    mal::GSRVA< output_sr_type, 1 > m_dX_dT;
    mal::GSRVA< output_sr_type, 1 > m_dX_dT_e;
    mal::GSRVA< output_sr_type, 1 > m_d2X_dT2_e;
    //\todo mal::GSRMA< jacobian_mat_type, 1, 1 > m_pX_pQ;
    //@}
};

// Exemples
/*
typedef GIConnectorSA<Vec3,Vec3> ConnectorPointOnParticle3;
typedef GIConnectorSA<Vec2,PosAngle2> ConnectorPointOnRigid2;
typedef GIConnectorSA<PosQuat,PosQuat> ConnectorRefSysOnRigid3;
*/

}}

#endif // S2_MS_ICONNECTOR_H
