#ifndef S2_BS_PARTICLESYS2D_H
#define S2_BS_PARTICLESYS2D_H

#include <Saphyre2/bs/ISyncEntity.h>

namespace S2
{

//! ParticleSys2D class
/*! App-level ParticleSys2D interface.
*/
class ParticleSys2D: public ISyncEntity
{
public:
    static const unsigned int cDefaultMaxParticles = 64;
    static const unsigned int cMaxForceApplicationsPerDT = 16;

    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        uint32 m_NumParticles;
        Real m_TotalMass;
        Real m_ParticleRadius;
        Real m_ParticleCoeffRestitution;
        
        Params()
        : m_Flags(0)
        , m_NumParticles(0)
        , m_TotalMass(1.0f)
        , m_ParticleRadius(0.1f)
        , m_ParticleCoeffRestitution(0.5f)
        {}
    };
    
public:
    
    EEntityType GetType() const { return eEntity_ParticleSys2D; }
    bool Update( Real dt );
    
    //! \name Initialization/Edition
    //@{
    bool SetParams( const Params &params ); //!< Returns false if error
    inline const Params &GetParams() const { return m_Params; }
        
    void DefineParticle( unsigned int pid, const Point2 &pos ); //!< \todo Could define mass/radius/coeffs if per-particle
    inline unsigned int GetNumParticles() const { return m_Params.m_NumParticles; }
    //@}

    //! \name Optimized for-each-particle edit methods
    //@{
    void SetVel( const Vec2 &vel );
    void ApplyForce( const Vec2 &f_global );
    void ApplyImpulse( const Vec2 &j_global );
    //@}
    
    //! \name Per-Particle edit methods
    //@{
    const Point2 &GetPos( int pid ) const;
    Vec2 GetVel( int pid ) const;
    
    void SetPos( int pid, const Point2 &pos );
    void SetVel( int pid, const Vec2 &vel );

    void ApplyForce( int pid, const Vec2 &f_global );
    void ApplyImpulse( int pid, const Vec2 &j_global );    
    //@}
    
private:

    //! \name ISyncEntity internal protocol
    //@{
    void BeginDef_Internal();
    bool EndDef_Internal();
    void Lock_Internal() {}
    void Unlock_Internal();
    //@}
    
    //! \name Constructor/Destructor and Methods with controlled scope
    //@{
    ParticleSys2D( ISyncEntity *p_parent );
    ~ParticleSys2D();
    //@}
    
    bool ProcessUpdate( const ds::ReturnIt &rit );
    
    friend class Universe; //!< Allows only Universe to create/destroy instances.
    
private:
    enum ETouchedAttribute {
        eTouchedNothing = 0,
        eTouchedPos     = (1<<0),
        eTouchedVel     = (1<<1),
        eTouchedForce   = (1<<2),
        eTouchedImpulse = (1<<3)
    };
    
private:
    
    //! \name Fixed Parameters
    //@{
    unsigned int m_NumParticles;
    Params m_Params;    
    //@}
    
    //! \name User State
    //@{
    Point2 *m_vecPos;
    //@}
    
    //! \name Per-particle Touch Events (pos/vel/force/impulse)
    /*!\todo specific TouchedParticle2D class sucks... Use a
      std::triad<uint16 touch_type,uint16 pid, vec2 data> array for
      touch events and, if maximum capacity is hit, END the current
      Edit(), Sync it, clear touch mask, OPEN a new Edit() and write
      the exceeding touch events there!
    */
    //@{
    unsigned int m_MaxTP;
    unsigned int m_NumTP;
    ds::TouchedParticle2D *m_vecTP;
    //@}
};

} // namespace S2

#endif // S2_PARTICLESYS2D_H
