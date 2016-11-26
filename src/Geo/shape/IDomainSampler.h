#ifndef GEO_IDOMAIN_SAMPLER_H
#define GEO_IDOMAIN_SAMPLER_H

#include <Geo/Config.h>

namespace geo {

// fwd decl
template<unsigned D> class GIShapeD;
class IDomainSampler;
template<unsigned D> class GIDomainSamplerD;
typedef GIDomainSamplerD<2> IDomainSampler2;
typedef GIDomainSamplerD<3> IDomainSampler3;

//typedef uint32 point_on_boundary_id_type;
typedef uint32 point_on_boundary_id_type;
const point_on_boundary_id_type cInvalidPOB( ~point_on_boundary_id_type(0) ); //0xFFFFFF...

class IDomainSampler
{
public:
    IDomainSampler() {}
    virtual ~IDomainSampler() {}

    //!\name Bookkeeping
    //@{
    virtual unsigned int GetDimension() const = 0;
    virtual unsigned int GetNumPOB() const = 0;
    virtual unsigned int GetNumAllocatedPOB() const = 0; // >= max POB id/index in use
    //@}
    
    //virtual void UpdatePositions( const Real *vec_dof ) = 0; \todo Updates all POB positions, later retrieved inline with GetPosition(pob)
    //virtual void UpdateNormals( const Real *vec_dof ) = 0; \todo Maybe

    virtual void UpdateDOF( const Real *vec_dof ) = 0; //Sets new DOF, may rebuild internal stuff
    virtual void ClearPOB() = 0;
    
    //\name array-POB lifetime methods
    //@{
    virtual void APOB_Destroy( const point_on_boundary_id_type *vec_id, unsigned int count ) = 0;
    virtual void APOB_IncRef( const point_on_boundary_id_type *vec_id, unsigned int count ) = 0;
    virtual void APOB_DecRef( const point_on_boundary_id_type *vec_id, unsigned int count ) = 0;
    //@}

    //\name single-POB helper methods
    //@{
    finline void POB_Destroy( point_on_boundary_id_type id ) { APOB_Destroy( &id, 1 ); }
    finline void POB_IncRef( point_on_boundary_id_type id ) { APOB_IncRef( &id, 1 ); }
    finline void POB_DecRef( point_on_boundary_id_type id ) { APOB_DecRef( &id, 1 ); }
    //@}
    
    //TEMP
    virtual void Trace() const = 0;
};

//! Boundary element data... WiP
template <typename T, unsigned D>
struct GBoundariel
{
    enum EConstants { cDimension = D };
    typedef mal::GVec<T,D> vec_type;
    vec_type m_Pos;
    vec_type m_Normal;
    T m_Radius;
};

/* Domain:
   - Supports Interior and Boundary POB
   - Allocs POB
*/
template <unsigned D>
class GIDomainSamplerD: public IDomainSampler
{
public:
    enum EConstants { cDimension = D };
    typedef mal::GTransform<Real,D> transform_type;
    typedef mal::GVec<Real,D> vec_type;
    typedef GBoundariel<Real,D> boundariel_type;
    
public:
    GIDomainSamplerD() {}
    ~GIDomainSamplerD() {}
    
    finline unsigned int GetDimension() const { return cDimension; }

    //\todo MAYBE: virtual void Crop( const bv::BoundingVolume &bv_local ) = 0; substract bv from domain, permanently?

    //\name PointOnBoundary creation
    //@{
    virtual point_on_boundary_id_type CreatePOB( const PointOnFeature &pof ) = 0; //Create POB from valid POF
    
    virtual point_on_boundary_id_type StepPOB( point_on_boundary_id_type id,
                                               const vec_type &dir_local, Real step_length, Real eps_length ) = 0;
    /*\todo the same, but returns full boundariel info at stepped POB
    virtual point_on_boundary_id_type StepBoundariel( point_on_boundary_id_type id,
                                                      const vec_type &dir_local, Real step_length, Real eps_length,
                                                      boundariel_type &boundariel ) = 0;
    */    
    /*\todo Step count POB at once would be faster due to amortization of the virtual call boundary stepping-initialization
      \todo Maybe adaptative stepping with max_steps and max_error so that the actual num_steps <= max_steps is variable (and therefore returned)
    virtual void StepPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                          point_on_boundary_id_type id, const vec_type &dir, Real total_length, Real eps_length ) = 0;
    */
    /*! Returns POB closest to pos_local and closest to POB[id] in a [-radius,+radius] neigbourhood around POB[id]
      \note MAY return POB[id] if it's already the closest one within eps_radius.
    */
    virtual point_on_boundary_id_type ClosestPOB( point_on_boundary_id_type id,
                                                  const vec_type &pos_local,
                                                  Real neighbourhood_radius, Real neighbourhood_epsilon ) = 0;
    
    virtual point_on_boundary_id_type CreateRandomPOB() = 0;
    virtual void CreateRandomPOB( point_on_boundary_id_type *vec_id, unsigned int count ) = 0;
    //! Create structured neighours around a POB (2 in 2D, N in 3D)
    virtual void CreateNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                                     point_on_boundary_id_type id,
                                     Real neighbourhood_radius, Real neighbourhood_epsilon ) = 0;
    //\note Neighbours equire DOF as radius is assumed in global metric, and DOF may deform the domain locally (ex: FEM)
    virtual void CreateRandomNeighbourPOB( point_on_boundary_id_type *vec_id, unsigned int count,
                                           Real neighbourhood_radius, point_on_boundary_id_type id ) = 0;
    //\todo virtual void CreateRandomPID( unsigned int num_pid, GIPointInDomainD<D> **vec_ptr_pid ) = 0;

    //virtual void CreateRandomNeighbourPID( GIPointInDomainD<D> *p_pid, Real radius, GIPointInDomainD<D> **vec_ptr_neighbours ) const = 0;
    //@}

    //\name \todo array-POB queries
    //@{
    /*\todo These would avoid A LOT of virtual POB_XXXX() calls during sample matching and contact determination...
    virtual void APOB_Position( vec_type *vec_position, const point_on_boundary_id_type *vec_id, unsigned int count, const Real *vec_dof ) const = 0;
    virtual void APOB_Normal( vec_type *vec_normal, const point_on_boundary_id_type *vec_id, unsigned int count, const Real *vec_dof ) const = 0;
    */
    //@}

    //!\name all-POB queries Should be faster than array-POB and single-POB (require iterating over all VALID pob, ignoring FREE ones)
    //@{
    //virtual void AllPositions( vec_type *vec_pos ) const = 0; //positions in local coords
    //virtual void AllNormals( vec_type *vec_normal ) const = 0; //normals in local coords
    virtual void AllPositionsAndNormals( vec_type *vec_pos, vec_type *vec_normal ) const = 0;
    //@}
    
    //\name single-POB queries
    //@{
    virtual vec_type POB_Position( point_on_boundary_id_type id ) const = 0;
    virtual vec_type POB_Normal( point_on_boundary_id_type id ) const = 0;
    //\todo Not sure yet virtual Real POB_Depth( point_on_boundary_id_type id ) const = 0;
    //\todo Not sure yet virtual Real POD_RayCast( point_on_boundary_id_type id, const vec_type &ray_dir_local ) const = 0;

    /*\todo Add point-on-feature specific params      
      It's not enough with the feature_id, we need the POB barycentric
      coords relative to such a feature.
      
      POB_PointOnFeature( ) -> point_on_feature = {feature_id = (type + index) + barycentric_coords = (b0,b1,b2,b3) }
      Vtx: b0 = 1
      Edge: b0 + b1 = 1
      Triangle: b0 + b1 + b2 = 1
      Tetrahedron: b0 + b1 + b2 + b3 = 1
    */
    virtual feature_id POB_FeatureId( point_on_boundary_id_type id ) const = 0;
    virtual PointOnFeature POB_PointOnFeature( point_on_boundary_id_type id ) const = 0;
    //virtual boundariel_type POB_Boundariel( point_on_boundary_id_type id ) = 0;
    //@}
    
    /*
    virtual bool IsInDomain( const vec_type &pos_local ) const = 0; //\todo potentially expensive
    virtual bool IsInterior( const vec_type &pos_local ) const = 0; //\todo potentially expensive
    virtual bool IsOnBoundary( const vec_type &pos_local ) const = 0; //\todo potentially expensive
    */

    /*\name Set queries
      virtual bool IsFinite() const = 0; //\todo Plane is NOT! CANNOT be sampled unless finite subdomain is specified!!
      virtual bool IsClosed() const = 0;
      virtual bool IsSimplyConnected() const = 0; //!< false if there are holes
      virtual bool IsConnected() const = 0;
      virtual bool IsConvex() const = 0;
    */
    
    virtual const GIShapeD<D> *GetShape() const = 0;
};


/*!
  TEMPORAL: This is kept here as a reference, but GIPointOnBoundaryD
  does NOT NEED to exist, POB are just indexed entries in a
  GIDomainSamplerD, no need for an actual class, just POB_XXXXXX()
  methods in GIDomainSamplerD. This simplifies allocation and makes POB
  ownership explicit

  PointOnBoundary connector, parametrized on point dimension D.

  \note ALL METHODS return data in SHAPE-LOCAL COORDS, according to
  the provided vec_dof.
  
  IShape subclasses are expected to derive a specific
  GIPointOnBoundary<D> and implement its operation AND all allocation
  related functionality.

  \todo GIPointOnBoundaryD is NOT an internal class of GIDomainSamplerD
  because it MAY be useful in other contexts, but MAY become so, we'll
  see...
  
template<unsigned D>
class GIPointOnBoundaryD
{
public:
    static const unsigned int cDimension = D;
    typedef mal::GVec<Real,D> vec_type;
public:
    GIPointOnBoundaryD() {}
    virtual ~GIPointOnBoundaryD() {}
    
    virtual vec_type Position( const Real *vec_dof ) const = 0;
    virtual vec_type Normal( const Real *vec_dof ) const = 0;
    //\todo Not sure yet virtual Real Depth( const Real *vec_dof ) const = 0;
    //\todo Not sure yet virtual Real RayCast( const vec_type &ray_dir_local, const Real *vec_dof ) const = 0;
    //\todo NOT SURE if necessary, probably yes to delete it... virtual const GIDomainSamplerD<D> *GetDomainSampler() const = 0;
};
*/

} //namespace geo

#endif // GEO_IDOMAIN_SAMPLER_H
