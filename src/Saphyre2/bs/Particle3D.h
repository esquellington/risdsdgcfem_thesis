#ifndef S2_BS_PARTICLE3D_H
#define S2_BS_PARTICLE3D_H

#include <Saphyre2/bs/ISyncEntity.h>
#include <Geo/geo.h>

namespace S2
{

//! Particle3D class
/*! App-level Particle3D System interface.
*/
class Particle3D: public ISyncEntity
{
public:

    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        Real m_Mass;
        Real m_Radius;
        Real m_CoeffRestitution;
        
        Params()
        : m_Flags(0)
        , m_Mass(0.1f)
        , m_Radius(0.01f)
        , m_CoeffRestitution(0.5f)
        {}
    };
    
public:
    
    EEntityType GetType() const { return eEntity_Particle3D; }
    bool Update( Real dt );
    
    //! \name Initialization
    //@{
    bool SetParams( const Params &params ); //!< Returns false if error
    inline const Params &GetParams() const { return m_Params; }

    bool AttachShape( geo::ShapeID shape_id );
    bool AttachShape( const geo::IShape &r_shape );
    //bool AttachShape( const geo::ShapeDef &shape_def );    
    //@}

    //! \name Query
    //@{
    Real GetMass() const { return m_Params.m_Mass; }
    const Point3 &GetPos() const { return m_Pos0; }
    const Vec3 &GetVel() const { return m_Vel0; }
    //}@
    
    //! \name Edit (requires Lock()-Unlock()-Sync() to be effective)
    //@{
    //void SetMass( Real mass );
    bool SetPos( const Point3 &pos );
    bool SetVel( const Vec3 &vel );
    
    bool ApplyForce( const Vec3 &f );
    bool ApplyImpulse( const Vec3 &j );
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
    Particle3D( ISyncEntity *p_parent );
    ~Particle3D();
    //@}
    
    bool ProcessUpdate( const ds::ReturnIt &rit );
    
    friend class Universe; //!< Allows only Universe to create/destroy instances.
    
private:
    enum ETouchedAttribute {
        eTouchedNothing = 0,
        eTouchedPos     = (1<<0),
        eTouchedVel     = (1<<1),
        eTouchedForce   = (1<<2),
        eTouchedImpulse = (1<<3),
        eTouchedAnything= 0xFFFFFFFF
    };
    
private:
    
    //! \name Fixed Parameters
    //@{
    Params m_Params;    
    geo::ShapeID m_ShapeId;
    //@}
    
    //! \name User State
    //@{
    Point3 m_Pos0;
    Point3 m_Pos1;
    Vec3 m_Vel0;
    Vec3 m_Vel1;
    Vec3 m_AccForce1;
    Vec3 m_AccImpulse1;
    //@}
};

} // namespace S2

#endif // S2_PARTICLE3D_H
