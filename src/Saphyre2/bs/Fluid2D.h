#ifndef S2_BS_FLUID2D_H
#define S2_BS_FLUID2D_H

#include <Saphyre2/bs/ISyncEntity.h>

namespace S2
{

//!  Fluid2D class
/*! App-level Fluid2D object interface.

\note This interface should be model-agnostic (SPH, Eulerian Grid,
PiC...), maybe with extensions for specific models.

*/
class Fluid2D: public ISyncEntity
{
public:
    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        uint32 m_NumParticles;
        Real m_Density;
        Real m_Thickness;
        Vec2 m_InitShape_AABB_PosMin;
        Vec2 m_InitShape_AABB_PosMax;
        Vec2 m_Bounds_AABB_PosMin;
        Vec2 m_Bounds_AABB_PosMax;
        
        Params()
        : m_Flags(0)
        , m_NumParticles(0)
        , m_Density(1000.0f)
        , m_Thickness(0.01f)
        , m_InitShape_AABB_PosMin(0,0)
        , m_InitShape_AABB_PosMax(1,1)
        , m_Bounds_AABB_PosMin(0,0)
        , m_Bounds_AABB_PosMax(1,1)
        {}
    };
    
public:

    EEntityType GetType() const { return eEntity_Fluid2D; }
    bool Update( Real dt );
    
    //! \name Initialization
    //@{
    bool SetParams( const Params &params ); //!< Returns false if error
    inline const Params &GetParams() const { return m_Params; }
    //@}    

    //! \name Query
    //@{
    Real GetDensity() const { return m_Params.m_Density; }
    unsigned int GetNumParticles() const { return m_Params.m_NumParticles; }
    //void GetBounds_AABB( Point2 &pos_min, Point2 &pos_max );    
    //@}
    
    //! \name Per-Particle methods (optional)
    //@{
    const Point2 &GetPos( int idx ) const;
    /*
    void SetPos( int idx, const Point2 &pos );
    void SetVel( int idx, const Vec2 &vel );
    void SetAcc( int idx, const Vec2 &acc );
    
    void ApplyForce( int idx, const Vec2 &f_global );
    void ApplyImpulse( int idx, const Vec2 &j_global );
    
    Real GetMass( int idx ) const;        
    Vec2 GetVel( int idx ) const;
    Vec2 GetAcc( int idx ) const;   
    */
    //@}

    //! \name Continuous-fluid methods (mandatory)
    //@{
    /*
    void SetPos( const Point2 &fluid_pos, const Point2 &pos );
    void SetVel( const Point2 &fluid_pos, const Vec2 &vel );
    void SetAcc( const Point2 &fluid_pos, const Vec2 &acc );
    
    void ApplyForce( const Point2 &fluid_pos, const Vec2 &f_global );
    void ApplyImpulse( const Point2 &fluid_pos, const Vec2 &j_global );
    */
    bool ApplyPressure( const Point2 &fluid_pos, Real radius, Real pressure );
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
    Fluid2D( ISyncEntity *p_parent );
    ~Fluid2D();
    //@}
    
    bool ProcessUpdate( const ds::ReturnIt &rit );
    
    friend class Universe; //!< Allows only Universe to create/destroy instances.
    
private:
    enum ETouchedAttribute {
        eTouchedNothing  = 0,
        eTouchedPos      = (1<<0),
        eTouchedVel      = (1<<1),
        eTouchedForce    = (1<<2),
        eTouchedImpulse  = (1<<3),
        eTouchedPressure = (1<<4),
        eTouchedAnything = 0xFFFFFFFF
    };
    
private:
    
    //! \name Fixed Parameters
    //@{
    Params m_Params;
    //@}
    
    //! \name User State
    //@{
    Point2 *m_vecPos;
    //@}

    //! \name Single radial pressure application per frame, by now
    //@{
    ds::RadialPressureAtPoint2D m_RadialPressure;
    //@}
};

} // namespace S2

#endif // S2_FLUID2D_H
