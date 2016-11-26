#ifndef S2_BS_KINE2D_H
#define S2_BS_KINE2D_H

#include <Saphyre2/bs/ISyncEntity.h>
#include <Geo/geo.h>

namespace S2
{

//! Kinematic stuff common to all sub-types
class Kinematic
{
public:
    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        Params() : m_Flags(0) {}
    };
};

//! Kine2D API
/*! Generic Kinematic 2D interface:
  - Transform:
    - Position
    - Rotation
    - Scale
  - Shape-dependent DOF
    - Real array

  \note Kine does NOT store shape, only uses it for
  initialization. Users can store it in their app entities if required
  for rendering, shape edition, etc...
*/
class Kine2D: public ISyncEntity
{
public:

    //!< Static on-creation params, inited to defaults
    struct Params: public Kinematic::Params {};

public:

    EEntityType GetType() const { return eEntity_Kine2D; }
    bool Update( Real dt );

    //! \name Initialization
    //@{
    bool SetParams( const Params &params ); //!< Returns false if error
    inline const Params &GetParams() const { return m_Params; }

    bool AttachShape( geo::ShapeID shape_id );
    bool AttachShape( const geo::IShape2 &r_shape );
    //bool AttachShape( const geo::ShapeDef &shape_def );
    const geo::IShape2* GetShape() const;
    geo::IShape2* GetShape(); //\todo DANGEROUS, only enabled to test BVH in testSaphyre

    /*No reason to use a "full" geo::IObject, transform and dof
      are already in Kine2D, and shape needs not be instantiated
      unless requested by GetShape()
    const geo_object_type *GetGO();
    const geo_object_type *GetGO() const;
    */
    //@}

    //! \name Query
    //@{
    inline const Point2 &GetPos() const { return m_Transform.m_Pos; }
    inline const Mat2x2 &GetRot() const { return m_Transform.m_Rot; }
    inline const Transform2 &GetTransform() const { return m_Transform; }
    inline uint32 GetNumDOF() const { return m_NumDOF; }
    inline const Real *GetVecDOF() const { return m_vecDOF; }

    //geo::IShape *GetShapeAPI() which allows changing the DOF and params?
    //}@

    //! \name Edit (require Lock()-Unlock()-Sync() to be effective)
    //@{
    bool SetPos( const Point2 &point );
    bool SetRot( const Mat2x2 &rot );
    bool SetAngle( Real angle );
    bool SetTransform( const Transform2 &transform );
    bool SetVecDOF( const Real *vec_dof );
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
    Kine2D( ISyncEntity *p_parent );
    ~Kine2D();
    //@}

    bool ProcessUpdate( const ds::ReturnIt &rit );

    friend class Universe; //!< Allows only Universe to create/destroy instances.

private:
    enum ETouchedAttribute {
        eTouchedNothing   = 0,
        eTouchedTransform = (1<<0),
        eTouchedDOF       = (1<<1),
        eTouchedAnything  = 0xFFFFFFFF
    };

private:

    //! \name Fixed Parameters
    //@{
    Params m_Params;
    geo::ShapeID m_ShapeId;
    //@}

    //! \name User State
    //@{
    Transform2 m_Transform;
    uint32 m_NumDOF;
    Real *m_vecDOF;
    mutable geo::IShape2 *m_pShape;
    //@}
};

} // namespace S2

#endif // S2_KINE2D_H
