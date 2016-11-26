#ifndef S2_BS_KINE3D_H
#define S2_BS_KINE3D_H

#include <Saphyre2/bs/ISyncEntity.h>
#include <Geo/geo.h>

namespace S2
{

//! Kine3D API
/*! Generic Kinematic 3D interface:
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
class Kine3D: public ISyncEntity
{
public:

    //!< Static on-creation params, inited to defaults
    struct Params
    {
        Flags32 m_Flags;
        Params() : m_Flags(0) {}
    };

public:

    EEntityType GetType() const { return eEntity_Kine3D; }
    bool Update( Real dt );

    //! \name Initialization
    //@{
    bool SetParams( const Params &params ); //!< Returns false if error
    inline const Params &GetParams() const { return m_Params; }

    bool AttachShape( geo::ShapeID shape_id );
    bool AttachShape( const geo::IShape3 &r_shape );
    //bool AttachShape( const geo::ShapeDef &shape_def );
    const geo::IShape3* GetShape() const;
    geo::IShape3* GetShape(); //\todo DANGEROUS, only enabled to test BVH in testSaphyre

    /*No reason to use a "full" geo::IObject, transform and dof
      are already in Kine3D, and shape needs not be instantiated
      unless requested by GetShape()
    const geo_object_type *GetGO();
    const geo_object_type *GetGO() const;
    */
    //@}

    //! \name Query
    //@{
    inline const Point3 &GetPos() const { return m_Transform.m_Pos; }
    inline const Mat3x3 &GetRot() const { return m_Transform.m_Rot; }
    inline const Transform3 &GetTransform() const { return m_Transform; }
    inline uint32 GetNumDOF() const { return m_NumDOF; }
    inline const Real *GetVecDOF() const { return m_vecDOF; }

    //geo::IShape *GetShapeAPI() which allows changing the DOF and params?
    //}@

    //! \name Edit (require Lock()-Unlock()-Sync() to be effective)
    //@{
    bool SetPos( const Point3 &point );
    bool SetRot( const Mat3x3 &rot );
    bool SetTransform( const Transform3 &transform );
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
    Kine3D( ISyncEntity *p_parent );
    ~Kine3D();
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
    Transform3 m_Transform;
    uint32 m_NumDOF;
    Real *m_vecDOF;
    mutable geo::IShape3 *m_pShape;
    //@}
};

} // namespace S2

#endif // S2_KINE3D_H
