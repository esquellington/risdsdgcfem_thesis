#ifndef S2_MS_SPHFLUID2D_H
#define S2_MS_SPHFLUID2D_H

#include "IModel.h"
#include "SPH_Grid2D.h"

#define __ENABLE_MODEL_VIZ
#ifdef __ENABLE_MODEL_VIZ
#  include <util/VizStream.h>
#endif

#define __USE_BUFFER_DISTANCES //MUST BE ENABLED FOR KEISER!!

#define __USE_PRESSURE_RELAXATION

namespace S2 { namespace ms {

//! Saphyre2 base SPHFluid2D model
/*! This class represents a mathematical model of a 2D SPH fluid

S(t) = [ pos_i vel_i+1/2 ]

*/
class SPH_Fluid2D: public IModel
{
public:
    enum EConstants {
        cDimension = 2,
        cMaxNeighbours = 32
    };
    
    enum EApproachType {
        eApproachDG96 = 0,
        eApproachKeiser = 1,
        eApproachClavet = 2
    };
    
public:
    SPH_Fluid2D();
    ~SPH_Fluid2D();
    
    void SetDscr( unsigned int num_particles,
                  Real density3d,
                  Real thickness,
                  Vec2 initshape_AABB_pos_min,
                  Vec2 initshape_AABB_pos_max,
                  Vec2 bounds_AABB_pos_min,
                  Vec2 bounds_AABB_pos_max,
                  EApproachType approach_type );

    //!\name Bookkeeping methods
    //@{
    EModelType GetModelType() const { return eModel_SPH_Fluid2D; }
    unsigned int GetDimension() const { return 2; }
    //@}
    
    //!\name Configuration, DOF and State vector management methods
    //@{
    unsigned int GetNumDOF() const { return 2*m_NumParticles; }
    void GetDOF( Real *p_dof ) const;

    unsigned int GetStateSize() const { return 4*m_NumParticles; }
    void GetState( Real *p_state ) const;
    void SetState( const Real *p_state );
    void Compute_dot_S( Real *p_dot_S ) const;

    Real GetMass() const { return m_Mass; }
    Real GetKineticEnergy() const;
    //@}
    
    //!\name Model-Specific State Consultors (SHOULDN'T EXIST!)
    //@{
    /*
    inline const Point2 &GetPos( int idx ) const { return m_vecPos[idx];  }
    inline Vec2 GetVel( int idx ) const { return m_vecVel[idx]; }
    inline Vec2 GetAcc( int idx ) const { DS_ASSERT(false); return Vec2(0,0); }
    inline const Vec2 GetMomentum( int idx ) const { DS_ASSERT(false); return Vec2(0,0); }

    inline void SetPos( int idx, const Point2 &pos ) { DS_ASSERT(false); }
    inline void SetVel( int idx, const Vec2 &vel ) { DS_ASSERT(false); }
    inline void SetMomentum( int idx, const Vec2 &momentum ) { DS_ASSERT(false); }
    */
    //@}

    //!\name Dynamics methods
    //@{
    void ApplyForce( const Real *p_force ) {}
    void GetForce( Real *p_force ) const {}
    void ApplyImpulse( const Real *p_impulse ) {}
    void GetImpulse( Real *p_impulse ) const {}    
    void ResetAccumulators();
    void Update( Real dt );
    //@}

    void ApplyPressure( const Point2 &fluid_pos, Real radius, Real pressure );

#ifdef __ENABLE_MODEL_VIZ
    void SetVizStream( util::VizStream *p_vs ) { m_pVS = p_vs; }
#endif
    
private:
    //! \name Internal update methods
    //@{
    void UpdateSCG( const Point2 *vec_pos );
    
    void UpdateDensity_DG96( const Point2 *vec_pos );
    void UpdateDensity_Keiser( const Point2 *vec_pos );
    void UpdateDensity_Clavet( const Point2 *vec_pos );
    
    void UpdateForces_DG96( const Point2 *vec_pos, const Vec2 *vec_vel );
    void UpdateForces_Keiser( const Point2 *vec_pos, const Vec2 *vec_vel );
    
    void UpdateViscosity_Clavet( const Point2 *vec_pos, Vec2 *vec_vel );
    void UpdatePressure_Clavet( Point2 *vec_pos );
    
    void UpdateCollisions( Point2 *vec_pos, Vec2 *vec_vel );
    void UpdateBounds( Point2 *vec_pos, Vec2 *vec_vel, Real Radius );

    void UpdateStep(); //Update fixed dt
    //@}

    void UpdateDebug();
    void ComputeStatistics( unsigned int &max_neighbours, float &avg_neighbours );
    
private:
    //! \name Params (Approach-independent)
    //@{
    unsigned int m_NumParticles;
    Real m_FixedDT;
    Real m_ResidualTime;
    Real m_Density2D;        //!< density2d = thickness * density3d
    Real m_Mass;             //!< mass = surface * density2D
    Real m_ParticleMass;     //!< constant by definition
    Real m_InvParticleMass;  //!< constant by definition
    Real m_RigidRadius;      //!< no other particle can get closer
    Real m_KernelRadius;     //!< kernel support radius h
    Real m_PowKernelRadius[10]; //!< powers from h^0 to h^9
        
    Real m_RestitutionPvsB;  //!< Particle Vs Bounds restitution
    Real m_RestitutionPvsP;  //!< Particle Vs Particle restitution

    Real m_ReferenceDensity; //!< density shift: pressure = k * (density-density0)    
    //@}
    
    //! \name Approach-dependent params
    //@{
    EApproachType m_Approach;
    Real m_CoeffGas;              //!< k in ideal gas equation pressure = k * density
    Real m_CoeffCohesion_Keiser;
    Real m_CoeffViscosity_Keiser;
    Real m_SoundSpeed_DG96;        //!< Fastest wave front \sa DG96
    Real m_CoeffGasNear_Clavet;    //!< Near-pressure gas coeff 
    Real m_CoeffViscosity1_Clavet; //!< Linear term \sigma
    Real m_CoeffViscosity2_Clavet; //!< Quadratic term \beta
    //@}
       
    //! \name Per-particle state/params
    //@{    
    Point2 *m_vecPos;
    Vec2 *m_vecVel;
    Vec2 *m_vecTmp; //!tmp vel or pos, depending on approach
    Vec2 *m_vecForce;
    Real *m_vecDensity;
    //@}        
    
    typedef SPH_Grid2D scg_type;
    //! SCG contains particle idx
    scg_type m_SCG;
    float m_ReferenceNonEmptyCellFraction;
    
    //! Per-particle cached SCG data
    struct ParticleDataSCG
    {
        unsigned int m_NumNeighbours;
        uint16 *m_vecNeighbours;
#ifdef __USE_BUFFER_DISTANCES
        float *m_vecDistances;
#endif
        scg_type::CellID m_CellId;
    };
    ParticleDataSCG *m_vecPD;
    uint16 *m_bufferNeighbours;
#ifdef __USE_BUFFER_DISTANCES
    float *m_bufferDistances;
#endif

#ifdef __ENABLE_MODEL_VIZ
    util::VizStream *m_pVS;
#endif
};

}} // namespace S2::ms

#endif // S2_MS_SPHFLUID2D_H
