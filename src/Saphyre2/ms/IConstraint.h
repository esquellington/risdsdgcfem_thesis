#ifndef S2_MS_ICONSTRAINT_H
#define S2_MS_ICONSTRAINT_H

#include "Config.h"
#include <Mal/GSRV.h>

namespace S2 { namespace ms {

/*! Low-level Constraint Interface based on \ref INTDYN90

  \todo USE SRV instead of Real*!!
  \todo Must support several paradigms:
  - Solve Locally for Acc, Vel and Pos
  - Export Magnitudes required to solve externally for Acc, Vel and Pos
  - Application of externally computed Force/Impulse/Variation

  Important:
  - dX_dT current time-derivative of X (eg: velocity of X)
  - dX_dT_e current time-derivative of X taking into account external forces/impulses.

  See \ref S2CONSTRAINTSDERIVATION for a complete derivation.
*/
class IConstraint
{
public:
    IConstraint() {}
    virtual ~IConstraint() {}

    //! \name Constraint Interface
    //@{
    virtual unsigned int GetNumDOF() const = 0;
    virtual void Update( Real t, Real dt ) = 0;

    virtual const Real *C() const = 0;
    virtual const Real *dC_dT() const = 0;
    virtual const Real *d2C_dT2_e() const = 0;

    virtual unsigned int GetNumInputs() const = 0;
    virtual const Real *pC_pQi( int input_idx ) const = 0;
    //@}

    //! \name Local Constraint Solving
    //@{
    virtual void SolveAcc( Real dt ) = 0;
    virtual void SolveVel( Real dt ) = 0;
    virtual void SolvePos( Real dt ) = 0;
    //@}

    //! \name Generalized Constraint Force/Impulse/Variation Application
    //@{
    virtual void ApplyForceLambda( const Real *lambda ) = 0;
    virtual void ApplyImpulseLambda( const Real *lambda ) = 0;
    virtual void ApplyVariationLambda( const Real *lambda ) = 0;
    //@}
};

/*! Single-input constraint
*/
class IUnaryConstraint: public IConstraint
{
public:
    IUnaryConstraint() {}
    virtual ~IUnaryConstraint() {}

    //! \name UC Interface
    //@{
    virtual const Real *pC_pQ1() const = 0;
    //@}

    unsigned int GetNumInputs() const { return 1; }
    const Real *pC_pQi( int input_idx ) const { MS_ASSERT(input_idx < (int)GetNumInputs()); return pC_pQ1(); }
};

/*! Two-input constraint
*/
class IBinaryConstraint: public IConstraint
{
public:
    IBinaryConstraint() {}
    virtual ~IBinaryConstraint() {}

    //! \name BC Interface
    //@{
    virtual const Real *pC_pQ1() const = 0;
    virtual const Real *pC_pQ2() const = 0;
    //@}

    unsigned int GetNumInputs() const { return 2; }
    const Real *pC_pQi( int input_idx ) const { MS_ASSERT(input_idx < (int)GetNumInputs()); return (input_idx==0) ? pC_pQ1() : pC_pQ2(); }
};

/*! Generic static allocation of Constraints with compile-time defined
  input/output SR types.

  \todo GIConstraintDA could do the same, but with dynamic allocation
  to support > 1 Count... for host objects with a variable number of
  sdof (ex: a whole particle system)
*/
template< typename OutputSRT, typename InputSRT1 >
class GIUnaryConstraintSA: public IUnaryConstraint
{
public:
    typedef OutputSRT output_sr_type;
    typedef InputSRT1 input1_sr_type;
    typedef mal::GMat< Real,
                       OutputSRT::size_in_reals,
                       InputSRT1::size_in_reals > jacobian1_mat_type;

protected:
    //!\name Statically allocated magnitudes
    //@{
    mal::GSRVA< output_sr_type, 1 > m_C;
    mal::GSRVA< output_sr_type, 1 > m_dC_dT;
    mal::GSRVA< output_sr_type, 1 > m_d2C_dT2_e;
    //\todo mal::GSRMA< jacobian1_mat_type, 1, 1 > m_pC_pQ1;
    //@}
};

template< typename OutputSRT, typename InputSRT1, typename InputSRT2 >
class GIBinaryConstraintSA: public IBinaryConstraint
{
};

} } // namespace S2::ms

#endif // S2_MS_ICONSTRAINT_H
