#include "SPH_Fluid2D.h"
#include "SPH_Kernels.h"
#include <Saphyre2/ds/DSG.h> //!\todo THIS BREAKS ENCAPSULATION, only included for Profiler and ParamStream
#include <Mal/GConversion.h>

//---- Default param values
//DG96
#define DG96_SPH_FIXED_DT 0.001  //def = ?
#define DG96_SPH_K_GAS 10.0      //def = 10.0
#define DG96_SPH_SOUND_SPEED 2.0 //def = 2.0
//Keiser
#define KEISER_SPH_FIXED_DT 0.003  //def = 0.001
#define KEISER_SPH_K_GAS 650.0     //def = 400
#define KEISER_SPH_K_COHESION 0.0  //def = ?
#define KEISER_SPH_K_VISCOSITY 35  //def = ?
//Clavet
#define CLAVET_SPH_FIXED_DT 0.003     //def = [1/30 .. 1/300] al paper
#define CLAVET_SPH_K_GAS 110.0        //def = 0.004 (non-normalized)
#define CLAVET_SPH_K_GAS_NEAR 10.0   //def = 0.01  (non-normalized)
#define CLAVET_SPH_K_VISCOSITY_1 80.0 //def = ?
#define CLAVET_SPH_K_VISCOSITY_2 80.0 //def = ?

//---- Simulation Flags
#define __DISABLE_COLLISIONS

//---- Debug Flags
#ifdef __ENABLE_MODEL_VIZ
#  define __DEBUG_INIT
#  define __DEBUG_SCG
//#  define __DEBUG_NEIGHBOURS
#  define __DEBUG_SPH_FIELD
#  define __DEBUG_ENERGY
//#  define __DEBUG_COLLISIONS
//#  define __DEBUG_KERNELS
#endif

//---- Param tweaking flags
//#define __ENABLE_MODEL_PARAMTWEAK

namespace S2 { namespace ms {

SPH_Fluid2D::SPH_Fluid2D()
: m_NumParticles(0)
, m_FixedDT(-1)
, m_ResidualTime(0)
, m_Mass(1)
, m_RestitutionPvsB(0.25f)
, m_RestitutionPvsP(0.0f)  //! \todo Crashes if > 0 !!!!!???
, m_ReferenceDensity(0)
  // Approach-dependent stuff
, m_Approach(eApproachKeiser)
, m_CoeffGas(-1)
, m_CoeffCohesion_Keiser(-1)
, m_CoeffViscosity_Keiser(-1)
, m_SoundSpeed_DG96(-1)
, m_CoeffGasNear_Clavet(-1)
, m_CoeffViscosity1_Clavet(-1)
, m_CoeffViscosity2_Clavet(-1)
  // Misc...
, m_vecPos(0)
, m_vecVel(0)
, m_vecTmp(0)
, m_vecForce(0)
, m_vecDensity(0)
, m_ReferenceNonEmptyCellFraction(0)
, m_vecPD(0)
, m_bufferNeighbours(0)
#ifdef __USE_BUFFER_DISTANCES
, m_bufferDistances(0)
#endif
#ifdef __ENABLE_MODEL_VIZ
, m_pVS(0)
#endif  
{
}


SPH_Fluid2D::~SPH_Fluid2D()
{
    if( m_vecPos ) delete m_vecPos;
    if( m_vecVel ) delete m_vecVel;
    if( m_vecTmp ) delete m_vecTmp;
    if( m_vecForce ) delete m_vecForce;
    if( m_vecDensity ) delete m_vecDensity;
    if( m_vecPD ) delete m_vecPD;
    if( m_bufferNeighbours ) delete m_bufferNeighbours;
#ifdef __USE_BUFFER_DISTANCES
    if( m_bufferDistances ) delete m_bufferDistances;
#endif
}    

void SPH_Fluid2D::SetDscr( unsigned int num_particles,
                           Real density3d,
                           Real thickness,
                           Vec2 initshape_AABB_pos_min,
                           Vec2 initshape_AABB_pos_max,
                           Vec2 bounds_AABB_pos_min,
                           Vec2 bounds_AABB_pos_max,
                           EApproachType approach_type )
{
    MS_ASSERT( m_NumParticles == 0 );

    m_Approach = approach_type;
    m_NumParticles = num_particles;

    Real shape_width = initshape_AABB_pos_max.x()-initshape_AABB_pos_min.x();
    Real shape_height = initshape_AABB_pos_max.y()-initshape_AABB_pos_min.y();
   
    m_Density2D = density3d*thickness;
    m_Mass = m_Density2D*(shape_width*shape_height);
    m_ParticleMass = m_Mass/num_particles;
    m_InvParticleMass = 1.0f/m_ParticleMass;
    
    m_vecPos = new Point2[num_particles];
    m_vecVel = new Vec2[num_particles];
    m_vecTmp = new Vec2[num_particles];
    m_vecForce = new Vec2[num_particles];
    m_vecDensity = new Real[num_particles];
    m_vecPD = new ParticleDataSCG[num_particles];
    m_bufferNeighbours = new uint16[num_particles * cMaxNeighbours];
#ifdef __USE_BUFFER_DISTANCES
    m_bufferDistances = new float[num_particles * cMaxNeighbours];
#endif

    //---- Init param values according to approach
    switch( m_Approach )
    {
    case eApproachDG96:
        m_FixedDT = DG96_SPH_FIXED_DT;
        m_CoeffGas = DG96_SPH_K_GAS;
        m_SoundSpeed_DG96 = DG96_SPH_SOUND_SPEED;
        break;
    case eApproachKeiser:
        m_FixedDT = KEISER_SPH_FIXED_DT;
        m_CoeffGas = KEISER_SPH_K_GAS;
        m_CoeffCohesion_Keiser = KEISER_SPH_K_COHESION;
        m_CoeffViscosity_Keiser = KEISER_SPH_K_VISCOSITY;
        break;
    case eApproachClavet:
        m_FixedDT = CLAVET_SPH_FIXED_DT;
        m_CoeffGas = CLAVET_SPH_K_GAS;
        m_CoeffGasNear_Clavet = CLAVET_SPH_K_GAS_NEAR;
        m_CoeffViscosity1_Clavet = CLAVET_SPH_K_VISCOSITY_1;
        m_CoeffViscosity2_Clavet = CLAVET_SPH_K_VISCOSITY_2;
        break;
    default: break;
    }
    
    //---- Gen particles inside initshape_AABB
        
    // we compute the dimensions of a uniform "particle matrix"...  
    unsigned int num_width_slots = (int)mal::Ceil( mal::Sqrt( (shape_width/shape_height) * m_NumParticles ) );
    unsigned int num_height_slots = (int)mal::Ceil( Real(m_NumParticles)/num_width_slots );

    //...compute the particle rigid radius...
    m_RigidRadius = 0.5f * Real(shape_width) / num_width_slots;
    
    //...and "raster" the particles by rows from bottom to top
    //placing them in a closest-sphere-packing hexagonal pattern where
    //even rows have num_width_slots particles and even rows have
    //num_width_slots-1 particles    
    Real particle_height_spacing = 2.0f*m_RigidRadius * mal::Sqrt(3.0f/4.0f);
    unsigned int num_particles_added = 0;
    for( unsigned int i=0; i<num_width_slots+1 && num_particles_added < m_NumParticles; i++ )
    {
        //Even-Odd rows are horizontally shifted by 1 radius
        Point2 start_pos = initshape_AABB_pos_min + Vec2( m_RigidRadius*(1+i%2),
                                                          m_RigidRadius + particle_height_spacing*i );
        for( unsigned int j=0;
             j<num_width_slots-i%2 && num_particles_added < m_NumParticles;
             j++ )
            m_vecPos[num_particles_added++] = start_pos + Vec2(j*2.0f*m_RigidRadius,0);        
    }

    // Compute kernel support radius and its powers
    m_KernelRadius = 3*m_RigidRadius;
    m_PowKernelRadius[0] = 1.0f;
    for( int i=1; i<10; i++ ) m_PowKernelRadius[i] = m_KernelRadius*m_PowKernelRadius[i-1];
    
    // Reset vel
    for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
        m_vecVel[it_particle] = Vec2(0,0);

    //---- Init SCG (ensure that cells size > kernel diameter)
    num_width_slots = (int)mal::Floor( Real(bounds_AABB_pos_max.x()-bounds_AABB_pos_min.x()) / (2.01*m_KernelRadius) );
    num_height_slots = (int)mal::Floor( Real(bounds_AABB_pos_max.y()-bounds_AABB_pos_min.y()) / (2.01*m_KernelRadius) );
    m_SCG.Init( bounds_AABB_pos_min, bounds_AABB_pos_max, num_width_slots, num_height_slots );
    
    // Add spheres to scg
    UpdateSCG( m_vecPos );
    unsigned int  max_pic;
    float avg_pic;
    m_SCG.ComputeStatistics( max_pic, avg_pic, m_ReferenceNonEmptyCellFraction );
    
    // Compute mean density and use to it as reference
    switch( m_Approach )
    {
    case eApproachDG96: UpdateDensity_DG96( m_vecPos ); break;
    case eApproachKeiser: UpdateDensity_Keiser( m_vecPos ); break;
    case eApproachClavet: UpdateDensity_Clavet( m_vecPos ); break;
    default: break;
    }
    m_ReferenceDensity = 0;
    for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
        m_ReferenceDensity += m_vecDensity[it_particle];
    m_ReferenceDensity /= m_NumParticles;
    
    // Init forces!! (required by Leapfrog)
    ResetAccumulators();
    switch( m_Approach )
    {
    case eApproachDG96: UpdateForces_DG96( m_vecPos, m_vecVel ); break;
    case eApproachKeiser: UpdateForces_Keiser( m_vecPos, m_vecVel ); break;
    case eApproachClavet: break;
    default: break;
    }
    
#ifdef __DEBUG_INIT
    UpdateDebug();
#endif

#ifdef __ENABLE_MODEL_PARAMTWEAK
    //Create and export Params
    ds::DSG::GetParameterIS()->BeginComplex(1/*"SPH_Fluid2D::Params"*/,eType_ParamGroup);
    {
        GSimpleParameter<float> param_float;

        param_float.m_Ptr = &m_FixedDT;
        param_float.m_Min = 1.0f/1000.0f;
        param_float.m_Max = 1.0f/30.0f;
        param_float.m_Step = 0.00025f;
        ds::DSG::GetParameterIS()->Write("dt",param_float);
      
        param_float.m_Ptr = &m_CoeffGas;
        param_float.m_Min = 0.0f;
        param_float.m_Max = 1000.0f;
        param_float.m_Step = 0.5f;
        ds::DSG::GetParameterIS()->Write("k_gas",param_float);
      
        switch( m_Approach )
        {
        case eApproachDG96:
            param_float.m_Ptr = &m_SoundSpeed_DG96;
            param_float.m_Min = 0.0f;
            param_float.m_Max = 20.0f;
            param_float.m_Step = 0.1f;
            ds::DSG::GetParameterIS()->Write("k_sound",param_float);
            break;
        case eApproachKeiser:
            param_float.m_Ptr = &m_CoeffCohesion_Keiser;
            param_float.m_Min = 0.0f;
            param_float.m_Max = 10.0f;
            param_float.m_Step = 0.1f;
            ds::DSG::GetParameterIS()->Write("k_cohesion",param_float);
            param_float.m_Ptr = &m_CoeffViscosity_Keiser;
            param_float.m_Min = 0.0f;
            param_float.m_Max = 100.0f;
            param_float.m_Step = 0.1f;
            ds::DSG::GetParameterIS()->Write("k_vis",param_float);
            break;
        case eApproachClavet:
            param_float.m_Ptr = &m_CoeffGasNear_Clavet;
            param_float.m_Min = 0.0f;
            param_float.m_Max = 1000.0f;
            param_float.m_Step = 0.5f;
            ds::DSG::GetParameterIS()->Write("k_gas_near",param_float);
            param_float.m_Ptr = &m_CoeffViscosity1_Clavet;
            param_float.m_Min = 0.0f;
            param_float.m_Max = 100.0f;
            param_float.m_Step = 0.1f;
            ds::DSG::GetParameterIS()->Write("k_vis1",param_float);
            param_float.m_Ptr = &m_CoeffViscosity2_Clavet;
            param_float.m_Min = 0.0f;
            param_float.m_Max = 100.0f;
            param_float.m_Step = 0.1f;
            ds::DSG::GetParameterIS()->Write("k_vis2",param_float);
            break;
        default: break;
        }
    }
    ds::DSG::GetParameterIS()->EndComplex();
#endif //__ENABLE_MODEL_PARAMTWEAK
}

//---- Configuration and State management methods
void SPH_Fluid2D::GetDOF( Real *p_dof ) const
{
    memcpy( (char*)p_dof,
            &m_vecPos[0],
            m_NumParticles*sizeof(Point2) );    
}

void SPH_Fluid2D::GetState( Real *p_state ) const
{
    memcpy( (char*)p_state,
            &m_vecPos[0],
            m_NumParticles*sizeof(Point2) );
    memcpy( (char*)p_state + m_NumParticles*sizeof(Point2),
            &m_vecVel[0],
            m_NumParticles*sizeof(Vec2) );
}

void SPH_Fluid2D::SetState( const Real *p_state )
{
    memcpy( &m_vecPos[0],
            (char*)p_state,
            m_NumParticles*sizeof(Point2) );
    memcpy( &m_vecVel[0],
            (char*)p_state + m_NumParticles*sizeof(Point2),
            m_NumParticles*sizeof(Vec2) );
}

Real SPH_Fluid2D::GetKineticEnergy() const
{
    Real ke = 0;
    for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
        ke += 0.5f*m_ParticleMass*m_vecVel[it_particle].NormSq();
    return ke;
}

void SPH_Fluid2D::Compute_dot_S( Real *p_dot_S ) const
{
    MS_ASSERT(false);
}

//---- Dynamics methods
void SPH_Fluid2D::Update( Real dt )
{
    m_ResidualTime += dt;
    m_ResidualTime = m_FixedDT;
    while( m_ResidualTime >= m_FixedDT )
    {
        UpdateStep();
        m_ResidualTime -= m_FixedDT;
    }
}    


//---- Internal update methods
void SPH_Fluid2D::UpdateSCG( const Point2 *vec_pos )
{
    DS_BPROF( "UpdateSCG" );
    
    // Clear
    m_SCG.Clear();
    
    DS_BPROF( "Classify" );
    // Classify
    for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
        m_vecPD[it_particle].m_CellId = m_SCG.AddParticle( it_particle, vec_pos[it_particle] );
    DS_EPROF();
    
    DS_BPROF( "Neighbours" );
    // Precompute neighbour lists (including the particle itself!)
    unsigned int total_neighbours = 0;
    for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
    {
        m_vecPD[it_particle].m_vecNeighbours = &m_bufferNeighbours[total_neighbours];
#ifdef __USE_BUFFER_DISTANCES
        m_vecPD[it_particle].m_vecDistances = &m_bufferDistances[total_neighbours];
        m_vecPD[it_particle].m_NumNeighbours = m_SCG.GetNeighbours2( vec_pos[it_particle], m_KernelRadius,
                                                                     vec_pos,
                                                                     m_vecPD[it_particle].m_vecNeighbours,
                                                                     m_vecPD[it_particle].m_vecDistances );
#else
        m_vecPD[it_particle].m_NumNeighbours = m_SCG.GetNeighbours( vec_pos[it_particle], m_KernelRadius,
                                                                    vec_pos,
                                                                    m_vecPD[it_particle].m_vecNeighbours );        
#endif
        total_neighbours += m_vecPD[it_particle].m_NumNeighbours;
        total_neighbours += total_neighbours & 1; //ensure parity to ensure 32-bit aligned sub-vectors        
    }
    //Detect overflow!!
    MS_ASSERT( total_neighbours < m_NumParticles*cMaxNeighbours );
    
    DS_EPROF();

    DS_EPROF();
}

// Compute density field (includes self-density)
void SPH_Fluid2D::UpdateDensity_DG96( const Point2 *vec_pos )
{
    DS_BPROF( "UpdateDensity_DG96" );
    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        Point2 pos1 = vec_pos[it_particle1];
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        Real density = 0.0f;
        
#ifdef __USE_BUFFER_DISTANCES
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
            density += m_ParticleMass
                       * Kernel_DCEW96_2D( pdscg.m_vecDistances[it_neighbour], m_PowKernelRadius );
#else
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
            density += m_ParticleMass
                       * Kernel_DCEW96_2D( (vec_pos[pdscg.m_vecNeighbours[it_neighbour]]-pos1).Norm(),
                                           m_PowKernelRadius );
#endif        
        m_vecDensity[it_particle1] = density;
    }
    
    DS_EPROF();
}

// Compute density field (includes self-density)
void SPH_Fluid2D::UpdateDensity_Clavet( const Point2 *vec_pos )
{
    DS_BPROF( "UpdateDensity_Clavet" );
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        Real density = 0.0f;        
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
            density += m_ParticleMass
                       * Kernel_Clavet_Poly2(pdscg.m_vecDistances[it_neighbour], m_PowKernelRadius);
        m_vecDensity[it_particle1] = density;
    }
    
    DS_EPROF();
}

void SPH_Fluid2D::UpdateDensity_Keiser( const Point2 *vec_pos )
{
    DS_BPROF( "UpdateDensity_Keiser" );
    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        Point2 pos1 = vec_pos[it_particle1];
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        Real density = 0.0f;        
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
            density += m_ParticleMass
                       * Kernel_Poly6_2D( pdscg.m_vecDistances[it_neighbour], m_PowKernelRadius );
        m_vecDensity[it_particle1] = density;
    }
    
    DS_EPROF();
}

void SPH_Fluid2D::ApplyPressure( const Point2 &fluid_pos, Real radius, Real pressure )
{
    // Assume SCG updated, find "all" (upper bouded to fit vec_pid[])
    // particles in "radius" around fluid_pos
    uint16 vec_pid[128];
    int num_affected_particles = m_SCG.GetParticlesInSphere( fluid_pos,
                                                             radius,
                                                             m_vecPos,
                                                             vec_pid, 128 );
    for( int it_particle=0; it_particle < num_affected_particles; it_particle++ )
    {
        Point2 r = m_vecPos[vec_pid[it_particle]] - fluid_pos;
        Real dist = r.Norm();
        m_vecForce[vec_pid[it_particle]] += (r / dist)
                                            * (1-dist/radius)
                                            * pressure;
    }
}

void SPH_Fluid2D::UpdateForces_DG96( const Point2 *vec_pos, const Vec2 *vec_vel )
{
    DS_BPROF( "UpdateForces_DG96" );
    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        Point2 pos1 = vec_pos[it_particle1];
        Real pressure1 = m_CoeffGas * (m_vecDensity[it_particle1] - m_ReferenceDensity);
        
        // External forces
        m_vecForce[it_particle1] += Vec2(0,-9.8f*m_ParticleMass);
        
        // Internal forces
        Vec2 force(0,0);
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
        {
            unsigned int it_particle2 = pdscg.m_vecNeighbours[it_neighbour];
            if( it_particle1 != it_particle2 )
            {
                Vec2 diff_pos = pos1 - vec_pos[it_particle2];
#ifdef __USE_BUFFER_DISTANCES
                Real dist = m_vecPD[it_particle1].m_vecDistances[it_neighbour];
#else
                Real dist = diff_pos.Norm();
#endif
                Vec2 u;
                u = (dist > 0.00001f) ? (diff_pos/dist) : Vec2(0,0);

                // Pressure Attraction/Repulsion
                Real pressure2 = m_CoeffGas * (m_vecDensity[it_particle2] - m_ReferenceDensity);
                force += -m_ParticleMass
                         * ( pressure1/mal::Sq(m_vecDensity[it_particle1])
                             + pressure2/mal::Sq(m_vecDensity[it_particle2]) )
                         * dWdr_Kernel_DCEW96_2D(dist,m_PowKernelRadius)
                         * u;
                // Viscosity (only when approaching)
                Vec2 vel12 = vec_vel[it_particle1]-vec_vel[it_particle2];
                Real mu = (0.5f*m_KernelRadius)*(vel12*diff_pos) / (mal::Sq(dist) + m_PowKernelRadius[2]/400.0f);
                if( mu < 0 )
                {
                    Real density12 = 0.5f*(m_vecDensity[it_particle1] + m_vecDensity[it_particle2]);
                    force += -m_ParticleMass
                             * ( (-m_SoundSpeed_DG96*mu + 2.0f*mal::Sq(mu)) / density12 )
                             * dWdr_Kernel_DCEW96_2D(dist,m_PowKernelRadius)
                             * u;
                }
            }
        }

        // Mass factor for particle1
        force *= m_ParticleMass;
        
        // Cohesion
        force += Vec2(0,0);
        
        // Total
        m_vecForce[it_particle1] += force;
    }

    DS_EPROF();
}


void SPH_Fluid2D::UpdateForces_Keiser( const Point2 *vec_pos, const Vec2 *vec_vel )
{
    DS_BPROF( "UpdateForces_Keiser" );
    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        Point2 pos1 = vec_pos[it_particle1];
        Real pressure1 = m_CoeffGas * (m_vecDensity[it_particle1]/m_ReferenceDensity - 1)
                         - m_CoeffCohesion_Keiser * mal::Sq(m_vecDensity[it_particle1]);
        Real volume1 = m_ParticleMass / m_vecDensity[it_particle1];
        
        // External forces
        m_vecForce[it_particle1] += Vec2(0,-9.8f*m_ParticleMass);
        
        // Internal forces
        Vec2 force(0,0);
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
        {
            unsigned int it_particle2 = pdscg.m_vecNeighbours[it_neighbour];
            if( it_particle1 != it_particle2 )
            {
                Vec2 diff_pos = pos1 - vec_pos[it_particle2];
                Real dist = m_vecPD[it_particle1].m_vecDistances[it_neighbour];
                Vec2 u;
                u = (dist > 0.00001f) ? (diff_pos/dist) : Vec2(0,0);

                // Pressure Attraction/Repulsion
                Real pressure2 = m_CoeffGas * (m_vecDensity[it_particle2]/m_ReferenceDensity - 1)
                                 - m_CoeffCohesion_Keiser * mal::Sq(m_vecDensity[it_particle2]);
                Real volume2 = m_ParticleMass / m_vecDensity[it_particle2];
                force += -volume2
                         * 0.5f*(pressure1+pressure2)
                         //* GradientKernel_Poly6_2D(dist,u,m_PowKernelRadius);
                         * dWdr_Kernel_DCEW96_2D(dist,m_PowKernelRadius) * u; //(Visualment millor que Poly6)
                
                // Viscosity (RADIAL version)
                Vec2 vel12 = vec_vel[it_particle1]-vec_vel[it_particle2];
                Real vel_r = vel12*u;
                force += -m_CoeffViscosity_Keiser
                         * volume2
                         * (-vel_r)
                         * dWdr_Kernel_DCEW96_2D(dist,m_PowKernelRadius) //(Estable!)
                         //* LaplacianKernel_Viscosity_Keiser_2D(dist,m_PowKernelRadius) //(Inestable!)
                         * u;
                
                /* Viscosity (Non-radial!) (ALWAYS? approaching AND separating!??)
                Vec2 vel12 = vec_vel[it_particle1]-vec_vel[it_particle2];                 
                force += -m_CoeffViscosity_Keiser
                         * volume2
                         * vel12
                         * LaplacianKernel_Viscosity_Keiser_2D(dist,m_PowKernelRadius);
                         */
            }
        }

        // Volume factor for particle1
        force *= volume1;

        // Cohesion
        force += Vec2(0,0);
        
        // Total
        m_vecForce[it_particle1] += force;
    }
        
    DS_EPROF();
}


/*! Implemented according to applyViscosity in Clavet with minor
    changes in notation and vector signs.
*/
void SPH_Fluid2D::UpdateViscosity_Clavet( const Point2 *vec_pos, Vec2 *vec_vel )
{
    DS_BPROF( "UpdateViscosity_Clavet" );
    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {        
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
        {
            unsigned int it_particle2 = pdscg.m_vecNeighbours[it_neighbour];
            Real dist = m_vecPD[it_particle1].m_vecDistances[it_neighbour];
            if( it_particle1 < it_particle2
                && dist < m_KernelRadius
                && dist > 0.0001f )
            {
                Vec2 diff_pos = vec_pos[it_particle1] - vec_pos[it_particle2];
                Vec2 u =  diff_pos / dist;
                Vec2 vel12 = vec_vel[it_particle1] - vec_vel[it_particle2];
                Real vel_r = vel12*u;
                if( vel_r < 0 )
                {
                    Vec2 impulse = m_FixedDT
                                   * (1.0f-dist/m_KernelRadius)
                                   * ( m_CoeffViscosity1_Clavet*(-vel_r)
                                       + m_CoeffViscosity2_Clavet*mal::Sq(vel_r) )
                                   * u;
                    vec_vel[it_particle1] += 0.5f*impulse;
                    vec_vel[it_particle2] -= 0.5f*impulse;
                }
            }
        }
    }

    DS_EPROF();
}

/*! Compute and apply pressure impulses using density/near-density.
  
  \note We assume that the neighbour lists do NOT change, but the
  neighbour positions do, and therefore we must recompute distances.
*/
void SPH_Fluid2D::UpdatePressure_Clavet( Point2 *vec_pos )
{
    DS_BPROF( "UpdatePressure_Clavet" );
    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        // Compute density
        Real density = 0;
        Real near_density = 0;        
        Point2 pos1 = vec_pos[it_particle1];
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];        
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
        {
            unsigned int it_particle2 = pdscg.m_vecNeighbours[it_neighbour];
            Vec2 diff_pos = pos1 - vec_pos[it_particle2];
            Real dist = diff_pos.Norm();
            if( dist < m_KernelRadius )
            {
                /* Unnormalized kernel
                Real tmp = 1.0f - dist/m_KernelRadius;
                Real tmp_sq = mal::Sq(tmp);
                density += m_ParticleMass*tmp_sq;
                near_density += m_ParticleMass*tmp_sq*tmp;
                */

                // Using normalized kernel
                density += m_ParticleMass * Kernel_Clavet_Poly2(dist,m_PowKernelRadius);
                near_density += m_ParticleMass * Kernel_DCEW96_2D(dist,m_PowKernelRadius);
            }
        }
        // Compute pressure
        Real pressure = m_CoeffGas * (density - m_ReferenceDensity);
        Real near_pressure = m_CoeffGasNear_Clavet * near_density;
        // Apply pressure displacements
        Vec2 acc_disp(0,0);
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
        {
            unsigned int it_particle2 = pdscg.m_vecNeighbours[it_neighbour];
            Vec2 diff_pos = pos1 - vec_pos[it_particle2];
            Real dist = diff_pos.Norm();
            if( it_particle1 < it_particle2 //AL PAPER NO FILTRA CASOS SIMETRICS!!?? Pero si no ho faig explota al principi!!
                && dist < m_KernelRadius )
            {
                Real tmp = 1.0f - dist/m_KernelRadius;
                Vec2 disp = mal::Sq(m_FixedDT)
                            * ( pressure*tmp + near_pressure*mal::Sq(tmp) )
                            * (diff_pos/dist);
                acc_disp += 0.5f*disp;
                vec_pos[it_particle2] -= 0.5f*disp;
            }
        }
        vec_pos[it_particle1] += acc_disp;
    }

    DS_EPROF();
}


void SPH_Fluid2D::UpdateCollisions( Point2 *vec_pos, Vec2 *vec_vel )
{
    DS_BPROF( "UpdateCollisions" );
    // Collide particles    
    for( unsigned int it_particle1=0; it_particle1<m_NumParticles; it_particle1++ )
    {
        Point2 pos1 = vec_pos[it_particle1];
        const ParticleDataSCG &pdscg = m_vecPD[it_particle1];
        for( unsigned int it_neighbour=0; it_neighbour < pdscg.m_NumNeighbours; it_neighbour++ )
        {
            unsigned int it_particle2 = pdscg.m_vecNeighbours[it_neighbour];
            if( it_particle1 < it_particle2 )
            {                
                Vec2 diff_pos = pos1-vec_pos[it_particle2];
                Real dist = diff_pos.Norm();
                Vec2 u = diff_pos/dist;
                if( dist < 2.0f*m_RigidRadius )
                {
                    // Fix positions
                    Real relaxation_factor = 0.75f;
                    vec_pos[it_particle1] += relaxation_factor*0.5*(2.0*m_RigidRadius-dist)*u;
                    vec_pos[it_particle2] -= relaxation_factor*0.5*(2.0*m_RigidRadius-dist)*u;
                    
                    // Bounce
                    Vec2 vel12 = vec_vel[it_particle1]-vec_vel[it_particle2];
                    Real vel12_rel = vel12*u;
                    Real e_factor = (1.0f-m_RestitutionPvsP);
                    if( vel12_rel < 0 )
                    {
                        vec_vel[it_particle1] -= e_factor*vel12_rel*u;
                        vec_vel[it_particle2] += e_factor*vel12_rel*u;
                    }                    
#ifdef __DEBUG_COLLISIONS
                    m_pVS->BeginComplex(1,util::eType_VizSegment);
                    {
                        m_pVS->Write("pos1", mal::CastDimension<3,float,2>(vec_pos[it_particle1]) );
                        m_pVS->Write("pos2", mal::CastDimension<3,float,2>(vec_pos[it_particle2]) );
                        m_pVS->BeginComplex("style",util::eType_VizStyle);
                        {
                            m_pVS->Write("color",Vec4f(1,0,0,1));
                            m_pVS->Write("pen_size",2.0f);
                        }
                        m_pVS->EndComplex();
                    }
                    m_pVS->EndComplex();
#endif
                }
            }
        }
    }
    DS_EPROF();
}

void SPH_Fluid2D::UpdateBounds( Point2 *vec_pos, Vec2 *vec_vel, Real radius )
{    
    DS_BPROF( "UpdateBounds" );
    
    Point2 aabb_min, aabb_max;
    m_SCG.GetAABB(aabb_min,aabb_max);    
    for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
    {
        if( vec_pos[it_particle].x() < aabb_min.x() + radius )
        {
            vec_pos[it_particle].x() = aabb_min.x() + radius;
            vec_vel[it_particle].x() = -m_RestitutionPvsB*vec_vel[it_particle].x();
        }
        else if( vec_pos[it_particle].x() > aabb_max.x() - radius )
        {
            vec_pos[it_particle].x() = aabb_max.x() - radius;
            vec_vel[it_particle].x() = -m_RestitutionPvsB*vec_vel[it_particle].x();
        }
        
        if( vec_pos[it_particle].y() < aabb_min.y() + radius )
        {
            vec_pos[it_particle].y() = aabb_min.y() + radius;
            vec_vel[it_particle].y() = -m_RestitutionPvsB*vec_vel[it_particle].y();
        }
        else if( vec_pos[it_particle].y() > aabb_max.y() - radius )
        {
            vec_pos[it_particle].y() = aabb_max.y() - radius;
            vec_vel[it_particle].y() = -m_RestitutionPvsB*vec_vel[it_particle].y();
        }            
    }
    DS_EPROF();
}
    
void SPH_Fluid2D::UpdateStep()
{    
    //---- Clear DIS
    /*
    ds::DSG::GetVizIS()->Clear();
    ds::DSG::GetProfilerIS()->Clear();
    ds::DSG::GetProfiler().Clear();
    */

    DS_BPROF( "UpdateStep" );

    switch( m_Approach )
    {
    case eApproachDG96:
        {    
            //---- Integrate Pos (Leapfrog)
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecPos[it_particle] += m_FixedDT*m_vecVel[it_particle];
            
            //---- Interact
            // Update SCG
            UpdateSCG( m_vecPos );
            
            // Precompute density field
            UpdateDensity_DG96( m_vecPos );
            
            // Predict Vel (Leapfrog)
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecTmp[it_particle] = m_vecVel[it_particle]
                                        + 0.5f*m_FixedDT*m_InvParticleMass*m_vecForce[it_particle];
            
            // Accumulate forces
            UpdateForces_DG96( m_vecPos, m_vecTmp );
            
            //---- Integrate Vel (Leapfrog)    
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecVel[it_particle] += m_FixedDT*m_InvParticleMass*m_vecForce[it_particle];
            ResetAccumulators();
            
            //---- Collide with bounds
            UpdateBounds( m_vecPos, m_vecVel, m_RigidRadius );
        }
        break;
    case eApproachKeiser:
        {
            //---- Integrate Pos (Leapfrog)
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecPos[it_particle] += m_FixedDT*m_vecVel[it_particle];
            
            //---- Interact
            // Update SCG
            UpdateSCG( m_vecPos );
            
            // Precompute density field
            UpdateDensity_Keiser( m_vecPos );
            
            // Predict Vel (Leapfrog)
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecTmp[it_particle] = m_vecVel[it_particle]
                                        + 0.5f*m_FixedDT*m_InvParticleMass*m_vecForce[it_particle];
            
            // Accumulate forces
            UpdateForces_Keiser( m_vecPos, m_vecTmp );
            
            //---- Integrate Vel (Leapfrog)    
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecVel[it_particle] += m_FixedDT*m_InvParticleMass*m_vecForce[it_particle];
            ResetAccumulators();
            
            //---- Collide with bounds
            UpdateBounds( m_vecPos, m_vecVel, m_RigidRadius );            
        }
        break;
    case eApproachClavet:
        {
            // Update SCG
            UpdateSCG( m_vecPos );

            //UpdateDensity_Clavet(m_vecPos); // FOR DEBUG ONLY!!
            
            // Apply ext forces
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
            {
                m_vecForce[it_particle] += Vec2(0,-9.8f*m_ParticleMass);
                m_vecVel[it_particle] += (m_FixedDT*m_InvParticleMass)*m_vecForce[it_particle];
            }
            ResetAccumulators();
            
            // Apply Viscosity
            UpdateViscosity_Clavet( m_vecPos, m_vecVel );
            
            // Integrate pos saving old values
            Vec2 *tmp = m_vecPos;
            m_vecPos = m_vecTmp;
            m_vecTmp = tmp;
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecPos[it_particle] = m_vecTmp[it_particle] + m_FixedDT*m_vecVel[it_particle];           
            
            // Apply pressure (density relaxation)            
            UpdatePressure_Clavet( m_vecPos );
            
            // Collide with bounds (EFFECT ON VEL IS IGNORED!!)
            UpdateBounds( m_vecPos, m_vecVel, m_RigidRadius );

            // Compute final vel from pos difference
            for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
                m_vecVel[it_particle] = (1.0f/m_FixedDT) * (m_vecPos[it_particle]-m_vecTmp[it_particle]);
        }
        break;        
    default: break;
    }

#ifndef __DISABLE_COLLISIONS
    UpdateCollisions( m_vecPos, m_vecVel );
#endif //__DISABLE_COLLISIONS
    
    //---- Debug Viz
    UpdateDebug();

    DS_EPROF(); //UpdateStep

    /*
    ds::DSG::GetProfilerIS()->WriteArray( "ProfilerData",
                                      &ds::DSG::GetProfiler().GetEntries()[0],
                                      ds::DSG::GetProfiler().GetEntries().size() );
    */
}

void SPH_Fluid2D::ResetAccumulators()
{
    memset( m_vecForce, 0, m_NumParticles*sizeof(Vec2) );
}

void SPH_Fluid2D::ComputeStatistics( unsigned int &max_neighbours, float &avg_neighbours )
{
    unsigned int acc_neighbours = 0;
    max_neighbours = 0;
    for( unsigned int it_particle=0; it_particle < m_NumParticles; it_particle++ )
    {
        unsigned int num_neighbours = m_vecPD[it_particle].m_NumNeighbours;
        acc_neighbours += num_neighbours;
        max_neighbours = mal::Max( max_neighbours, num_neighbours );
    }
    avg_neighbours = float(acc_neighbours)/m_NumParticles;
}

//! Debug Viz
void SPH_Fluid2D::UpdateDebug()
{
    DS_BPROF("UpdateDebug");
#ifdef __DEBUG_SCG
    // SCG
    {
        // Extract SCG data and stats
        Point2 aabb_min, aabb_max;
        m_SCG.GetAABB(aabb_min,aabb_max);        
        unsigned int dim1,dim2;
        m_SCG.GetDimensions( dim1, dim2 );
        unsigned int  max_neighbours, max_pic;
        float avg_neighbours, avg_pic, non_emtpy_cell_fraction;
        ComputeStatistics( max_neighbours, avg_neighbours );
        m_SCG.ComputeStatistics( max_pic, avg_pic, non_emtpy_cell_fraction );
        char str[256];
        sprintf(str,"Grid2D[%d][%d] : Neigh(%f,%d), NECF(%3.0f%%,%3.0f%%), PiC(%f,%d)",
                dim1, dim2, avg_neighbours, max_neighbours,
                100.0 * non_emtpy_cell_fraction, 100.0*non_emtpy_cell_fraction/m_ReferenceNonEmptyCellFraction,
                avg_pic, max_pic );
        // And draw it
        m_pVS->BeginComplex(str,util::eType_VizGrid2D);
        {
            m_pVS->Write("pos", mal::CastDimension<3,float,2>(aabb_min) );
            m_pVS->Write("dir1", Vec3f(aabb_max[0]-aabb_min[0],0,0) );
            m_pVS->Write("dir2", Vec3f(0,aabb_max[1]-aabb_min[1],0) );

            m_pVS->Write("size1", dim1 );
            m_pVS->Write("size2", dim2 );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color",Vec4f(0,0.5f,0.5f,1));
                m_pVS->Write("pen_size",1.0f);
            }
            m_pVS->EndComplex();
        }
        m_pVS->EndComplex();
        
        /*Draw full cells in red
          std::vector<int> vec_cells;
          m_SCG.GetNonEmptyCells( vec_cells );
          for( int it_cell=0; it_cell < vec_cells.size(); it_cell++ )
          {
          Point2 p2d = m_SCG.GetPos( vec_cells[it_cell] );
          debug::VizPoint viz_point;
          viz_point.m_Color = Vec3f(1,0,0);
          viz_point.m_Size = 10.0f;
          viz_point.m_Pos = Vec3f(p2d.x(),p2d.y(),0);
          ds::DSG::GetVizIS()->Write( 1, viz_point );
          }
        */

        /*
        uint16 vec_pid[1024];
        int num_particles_in_sphere = m_SCG.GetParticlesInSphere( 0.5*(aabb_max+aabb_min),
                                                                  0.25, m_vecPos, vec_pid );
        for( int it_particle=0; it_particle < num_particles_in_sphere; it_particle++ )
        {
            Point2 p2d = m_vecPos[vec_pid[it_particle]];
            debug::VizPoint viz_point;
            viz_point.m_Color = Vec3f(0,1,0);
            viz_point.m_Size = 10.0f;
            viz_point.m_Pos = Vec3f(p2d.x(),p2d.y(),0);
            ds::DSG::GetVizIS()->Write( 1, viz_point );
        }
        */
    }
#endif

#ifdef __DEBUG_NEIGHBOURS
    for( int it_particle=0; it_particle<m_NumParticles; it_particle++ )
    {
        for( int it_neighbour=0; it_neighbour<m_vecPD[it_particle].m_NumNeighbours; it_neighbour++ )
        {
            m_pVS->BeginComplex(1,util::eType_VizSegment);
            {
                m_pVS->Write("pos1", mal::CastDimension<3,float,2>(m_vecPos[it_particle]) );
                m_pVS->Write("pos2",
                             mal::CastDimension<3,float,2>(m_vecPos[m_vecPD[it_particle].m_vecNeighbours[it_neighbour]]) );
                m_pVS->BeginComplex("style",util::eType_VizStyle);
                {
                    m_pVS->Write("color",Vec4f(1,0,0,1));
                    m_pVS->Write("pen_size",2.0f);
                }
                m_pVS->EndComplex();
            }
            m_pVS->EndComplex();
        }
    }
#endif
    
#ifdef __DEBUG_SPH_FIELD
    // SPH Particles
    for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
    {
        // Kernels as Circles        
        m_pVS->BeginComplex(1,util::eType_VizDisk);
        {
            m_pVS->Write("pos", mal::CastDimension<3,float,2>(m_vecPos[it_particle]) );
            m_pVS->Write("rot", Quatf::Identity() );
            m_pVS->Write("radius", m_KernelRadius );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color", Vec4f(0.25,0.25,0,1) );
                m_pVS->Write("flags", Flags32(util::eVizStyle_Wire) );
            }
            m_pVS->EndComplex();
        }
        m_pVS->EndComplex();

        // Rigid inner radius
        m_pVS->BeginComplex(1,util::eType_VizDisk);
        {
            m_pVS->Write("pos", mal::CastDimension<3,float,2>(m_vecPos[it_particle]) );
            m_pVS->Write("rot", Quatf::Identity() );
            m_pVS->Write("radius", m_RigidRadius );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color", Vec4f(1,0.5,0,1) );
                m_pVS->Write("flags", Flags32(util::eVizStyle_Wire) );
            }
            m_pVS->EndComplex();
        }
        m_pVS->EndComplex();        
        
        // (normalized) Density field        
        //Real normalized_density = log(m_vecDensity[it_particle]/m_Density2D)/log(10.0);
        Real normalized_density = m_vecDensity[it_particle]/m_Density2D;
        m_pVS->BeginComplex(1,util::eType_VizVec);
        {
            m_pVS->Write("pos", mal::CastDimension<3,float,2>(m_vecPos[it_particle]) );
            m_pVS->Write("vec", Vec3f(0,0,-0.1*normalized_density) );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color",
                             Vec4f( (normalized_density>1.0f) ? 1.0f : 0,
                                    0.0f,
                                    (normalized_density<=1.0f) ? normalized_density : 0.0f,
                                    1.0f ) );
            }
            m_pVS->EndComplex();
        }
        m_pVS->EndComplex();
        
        /*CellIds
        char str[128];
        sprintf(str,"(%d,%d)",m_vecPD[it_particle].m_CellId.i, m_vecPD[it_particle].m_CellId.j );          
        m_pVS->BeginComplex(str,util::eType_VizPoint);
        {
            m_pVS->Write("pos", mal::CastDimension<3,float,2>(m_vecPos[it_particle]) );
            m_pVS->Write("vec", Vec3f(0,0,-0.1*normalized_density) );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color",Vec4f(0,0,1,1));
            }
            m_pVS->EndComplex();
        }          
        m_pVS->EndComplex();
        */
        
        // Warning if num neighbours is MAX, because we may have missed some others!        
        if( m_vecPD[it_particle].m_NumNeighbours > cMaxNeighbours )
        {
            m_pVS->BeginComplex(1,util::eType_VizDisk);
            {
                m_pVS->Write("pos", mal::CastDimension<3,float,2>(m_vecPos[it_particle]) );
                m_pVS->Write("rot", Quatf::Identity() );
                m_pVS->Write("radius", m_KernelRadius );
                m_pVS->BeginComplex("style",util::eType_VizStyle);
                {
                    m_pVS->Write("color", Vec4f(1,0,0,1) );
                    m_pVS->Write("flags", Flags32(util::eVizStyle_Solid) );
                }
                m_pVS->EndComplex();
            }
            m_pVS->EndComplex();
        }
    }    
#endif

#ifdef __DEBUG_ENERGY
    {
        // Mechanical Energy
        Real ke = GetKineticEnergy();
        Real pe = 0;
        for( unsigned int it_particle=0; it_particle<m_NumParticles; it_particle++ )
            pe += m_ParticleMass*m_vecPos[it_particle].y();
        Real me = ke + pe;        
        char str[256];
        sprintf(str,"me/ke/pe = %f/%f/%f",me,ke,pe);
        m_pVS->BeginComplex(str,util::eType_VizVec);
        {
            m_pVS->Write("pos", Vec3f(-1,0,0) );
            m_pVS->Write("vec", Vec3f(0,me,0) );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color",Vec4f(0.7f,0.7f,0.7f,1));
            }
            m_pVS->EndComplex();
        }
        m_pVS->EndComplex();
    }
#endif

#ifdef __DEBUG_KERNELS
    {
        // DG96 kernel
        float dr = 0.0001f;
        float kernel_integral_dg96 = 0.0f;
        for( float r=0.0f; r<m_KernelRadius; r += dr )
            kernel_integral_dg96 += Kernel_DCEW96_2D(r,m_PowKernelRadius)
                                    * MAL_CONSTANT_PI * (mal::Sq(r+dr) - mal::Sq(r));
        float kernel_integral_keiser_poly6 = 0.0f;
        for( float r=0.0f; r<m_KernelRadius; r += dr )
            kernel_integral_keiser_poly6 += Kernel_Poly6_2D(r,m_PowKernelRadius)
                                            * MAL_CONSTANT_PI * (mal::Sq(r+dr) - mal::Sq(r));
        float kernel_integral_keiser_spiky = 0.0f;
        for( float r=0.0f; r<m_KernelRadius; r += dr )
            kernel_integral_keiser_spiky += Kernel_Spiky_2D(r,m_PowKernelRadius)
                                            * MAL_CONSTANT_PI * (mal::Sq(r+dr) - mal::Sq(r));

        char str[256];
        sprintf(str,"Kernels (DG96,KeiserPoly6,KeiserSpiky) = %f %f %f",
                kernel_integral_dg96, kernel_integral_keiser_poly6, kernel_integral_keiser_spiky);        
        m_pVS->BeginComplex(str,util::eType_VizVec);
        {
            m_pVS->Write("pos", Vec3f(-2,0,0) );
            m_pVS->Write("vec", Vec3f(0,kernel_integral_dg96,0) );
            m_pVS->BeginComplex("style",util::eType_VizStyle);
            {
                m_pVS->Write("color",Vec4f(0.7f,0.7f,0.7f,1));
            }
            m_pVS->EndComplex();
        }
        m_pVS->EndComplex();
    }    
#endif
    
    DS_EPROF();
}

} } // namespace S2::ms
