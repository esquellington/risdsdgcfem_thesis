#ifndef S2_MS_FEM_TRIANGLE_ELEMENT_SET_2_H
#define S2_MS_FEM_TRIANGLE_ELEMENT_SET_2_H

#include "../Config.h"
#include "Params.h"
#include "TriangleElement2.h"
#include <util/GPoolDA.h>

#define __ENABLE_PLASTICITY

//#define __ENABLE_TRACE_FEM
#ifdef __ENABLE_TRACE_FEM
#  include <iostream>
#  include <Mal/GSerialization.h>
#endif

namespace S2 {

typedef mal::GVec<Real,6> Vec6;
typedef mal::GMat<Real,2,3> Mat2x3;

}

namespace S2 {
namespace ms {

/*! Generic Force-On-Particle-System interface

  \note ALL methods ACCUMULATE on input arrays to avoid temporaries,
  they should be EXTERNALLY RESET to 0 if no accumulation desired.

  \note Notation:
  dx = differential of x
  df_x = differential of function f(x,v,t) due to differential in variable x
  PD_f_x = partial derivative of f(x,v,t) wrt x

  \todo To allow both instantiation of IForce definitions and
  application to iterative solvers, IForce API is strictly const.

  Any external state and params required is stored in a specific
  ForceCache instance that will be passed to all IForce methods.
  Also, consider passing vec_x, vec_v to all methods instead setting
  them and storing them in BeginEvaluation().

  ALTERNATIVELY, save them in the cache at BeginEvaluation(). Also,
  the specific cache would be created by each specific IForce subclass
  and can contain specific data, such as the vec_node_inv_mass for
  rayleigh damping, for example. \todo All per-particle data (gravity
  => mass, electric => charge, etc...) required by a specific force
  could be handled this way.
*/
class IForce2
{
public:
    class ICache
    {
    public:
        ICache() : m_bIsValid(false) {}
        virtual ~ICache() {}
        inline bool IsValid() const { return m_bIsValid; }
    protected:
        bool m_bIsValid;
    };

public:
    IForce2() {}
    virtual ~IForce2() {}

    //\name IForce virtual API
    //@{
    virtual void BeginEvaluation( ICache* p_cache, const Vec2 *vec_prev_x, const Vec2 *vec_x, const Vec2 *vec_v, Real t ) const = 0;
    virtual void f( ICache* p_cache, Vec2 *vec_f ) const = 0; //! += Force f(x,v,t)
    virtual void df_x( ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const = 0; //! += Force differential df_x = \frac{\partial f}{partial x}(x,t) * dx
    virtual void df_v( ICache* p_cache, const Vec2 *vec_dv, Vec2 *vec_df ) const = 0; //! += Force differential df_v = \frac{\partial f}{partial v}(x,v,t) * dv
    virtual Real V( ICache* p_cache ) const = 0; //! Potential energy V(x,v,t)
    //virtual void PD_f_x( const Vec2 *vec_x, Real t ) const = 0; //! Matrix D_x(f) = \frac{\partial f}{partial x}(x,t)
    //virtual void PD_f_v( const Vec2 *vec_x, const Vec2 *vec_v, Real t ) const = 0; //! Matrix D_x(f) = \frac{\partial f}{partial v}(x,v,t)
    virtual void EndEvaluation( ICache* p_cache ) const = 0;
    //@}
};

}} //namespace S2::ms

namespace S2 {
namespace ms {
namespace fem {

/*! Constant Triangle Element Set 2D

  Supports several FEM Params::EMethod that can be changed on the fly.

  \todo Current API is const, except for per-element cached data
        ElementCache, which is mutable. To allow instancing with
        shared access to the TriangleElementSet2 mesh, ElementCache
        could be externally allocated per-instance and passed in
        BeginEvaluation, for example.
           Also, this would allow finer control of which magnitudes
        are recomputed and which ones are not when used inside
        Quasi/Modified/Exact Newton-Raphson methods.

        IMPORTANT: During NR-iterations, the f() and df_x() will be
        evaluated for several intermediate states X_k that may change
        the initial degeneration status at X_0. Non-chronological
        methods (SVD1, QR) would recompute R_k at these intermediate
        states. For our chronologically-correct methods, we *SHOULD*
        recompute the DoC_k and its associated R_k *BUT* considering
        the displacement X_k - X_0, NOT X_k - X_k-1.
*/
class TriangleElementSet2: public IForce2
{
public:
    typedef uint16 element_index_type; //\todo Could parametrize or upgrade to uint32 for large meshes
    typedef uint16 node_index_type; //\todo Could parametrize or upgrade to uint32 for large meshes
    enum EConstants { cInvalidNoC = -1 };

public:
    TriangleElementSet2();
    virtual ~TriangleElementSet2();

    inline unsigned int GetNumElements() const { return m_NumElements; }
    inline bool IsValid() const { return 0 != m_vecRefPos; } //\note otherwise, it has NOT been initialized with SetBakedData()

    //! \name Specific cache creation/Destruction (eg: needs node inv masses for rayleigh damping M term)
    //@{
    IForce2::ICache* CreateCache( const Real* vec_inv_mass );
    void DestroyCache( IForce2::ICache* p_cache );
    //@}

    //\name IForce virtual API
    //@{
    void BeginEvaluation( IForce2::ICache* p_cache, const Vec2* vec_prev_x, const Vec2* vec_x, const Vec2* vec_v, Real t ) const;
    void f( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void df_v( IForce2::ICache* p_cache, const Vec2 *vec_dv, Vec2 *vec_df ) const;
    Real V( IForce2::ICache* p_cache ) const;
    void EndEvaluation( IForce2::ICache* p_cache ) const;
    //@}

    //! Elastic force \todo Consider adding to f() implicitly OR add f_Elastic() to IForce virtual API and add vec_node_mass to BeginEvaluation() params
    void f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    //! Damping force \todo Consider adding to f() implicitly OR add f_Damping() to IForce virtual API and add vec_node_mass to BeginEvaluation() params
    void f_d( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    //! Plastic force \todo Consider adding to f() implicitly OR add f_Plastic() to IForce virtual API and add vec_node_mass to BeginEvaluation() params
    void f_p( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    // Resets plastic deformation
    void ResetPlasticity( IForce2::ICache* p_cache );

    //\name Consultors
    //@{
    inline Real TotalArea() const { return m_TotalArea; }
    inline Real Area( unsigned int eid ) const { return m_vecED[eid].m_Area; }
    void GetNID( unsigned int eid, unsigned int &nid0, unsigned int &nid1, unsigned int &nid2 ) const;

    Transform2 GetTransformE2W( const IForce2::ICache* p_cache, unsigned int eid, const Vec2* vec_pos ) const;
    unsigned int GetNoC( const IForce2::ICache* p_cache, unsigned int eid ) const;
    //@}

    //TEMP: Apply reaction constraints on inverted elements, this may be generalized into IForce (as a strain-limit) or IConstraintSet (as a general constraint eq)
    //@{
    void UpdateConstraints_ActiveSet( IForce2::ICache* p_cache, const Vec2* vec_vel, const Vec2* vec_acc );
    void ApplyConstraints_NoC_Reaction( const IForce2::ICache* p_cache, Vec2* vec_v ) const;
    //@}

public:
    struct ElementData;
    struct LinearED;

protected:
    void SetBakedData( bool b_shared,
                       uint32 num_nodes, uint32 num_elements,
                       const Vec2 *vec_ref_pos, const ElementData *vec_ed, const LinearED *vec_led );
    void ClearBakedData();

private:
    void UpdateDegeneration( IForce2::ICache* p_cache ) const;
    void UpdatePlasticity( IForce2::ICache* p_cache ) const;

    //\name Corotational R computation
    //@{
    void ComputeRotations_Id( IForce2::ICache* p_cache ) const;
    void ComputeRotations_QR( IForce2::ICache* p_cache ) const;
    void ComputeRotations_MSVD( IForce2::ICache* p_cache ) const;
    void ComputeRotations_PD_Project( IForce2::ICache* p_cache ) const;
    void ComputeRotations_PD_Reflect( IForce2::ICache* p_cache ) const;
    void ComputeRotations_DAPD( IForce2::ICache* p_cache ) const;
    //@}

    //\name f_e()
    //@{
    void L_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void C_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void C_LCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void C_CCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void H_LCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void H_CCM_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void H_NHC0_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void H_NHC1_f_e( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    //@}

    //\name V()
    //@{
    Real L_V( IForce2::ICache* p_cache ) const;
    Real C_V( IForce2::ICache* p_cache ) const;
    Real C_LCM_V( IForce2::ICache* p_cache ) const;
    Real C_CCM_V( IForce2::ICache* p_cache ) const;
    Real H_LCM_V( IForce2::ICache* p_cache ) const;
    Real H_CCM_V( IForce2::ICache* p_cache ) const;
    Real H_NHC0_V( IForce2::ICache* p_cache ) const;
    Real H_NHC1_V( IForce2::ICache* p_cache ) const;
    //@}

    //\name f_d()
    //@{
    void L_f_d( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void C_f_d( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    //@}

    //\name f_p()
    //@{
    void L_f_p( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    void C_f_p( IForce2::ICache* p_cache, Vec2 *vec_f ) const;
    //@}

    //\name df_x()
    //@{
    void L_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void C_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void C_LCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void C_LCMH_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void C_CCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void H_LCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    void H_CCM_df_x( IForce2::ICache* p_cache, const Vec2 *vec_dx, Vec2 *vec_df ) const;
    //@}

public:
    //! ElementData, constant during entire element lifetime
    struct ElementData
    {
        node_index_type m_vecNID[3];
        Mat2x2 m_InvDm; //D_m^-1
        Real m_Area;
    };
    //! LFEM specific data
    struct LinearED
    {
        Mat2x2 m_K00, m_K01, m_K02, m_K11, m_K12, m_K22; //Essential
        Mat2x2 m_K10, m_K20, m_K21; //Transposed
        /*\todo Alternatively: Mat2x2 m_matK[3][3], so that K_ij are
          directly accessible, or even explicit Mat2x2 m_K00, m_K01,
          m_K02, m_K11, m_K12, m_K22 and external transposition to
          access symmetric blocks m_K10, m_K20, m_K21.
        */
#ifdef __ENABLE_PLASTICITY
        Mat2x3 m_P0, m_P1, m_P2;
#endif
    };

protected:
    //constant representation
    node_index_type m_NumNodes;
    element_index_type m_NumElements;
    Real m_TotalArea;

    //\name params
    //@{
    Params m_Params;
    // Derived
    Real m_LameMu, m_LameLambda;
    //@}

    const Vec2 *m_vecRefPos;
    const ElementData *m_vecED;
    const LinearED *m_vecLED;

    //uint32 *m_pBuffer; //!< If exists, the baked data is NOT SHARED
};

/*! Editable Triangle Element Set 2D
  - Construction: Add/Remove elements
  - Modification: Change parameters

  \todo Despite selected EMethod, we precompute all constant data
  required by any method, as we may change it on the fly later.

  \note We follow the MeshSolidShape2 Baked/Editable data paradigm
*/
class EditableTriangleElementSet2: public TriangleElementSet2
{
public:
    EditableTriangleElementSet2();
    ~EditableTriangleElementSet2();

    void ClearEditData();
    //\todo if method is given per-element, params.m_Method == eMethod_PerElement and alloc a per-element method descriptor
    void BeginEdition( unsigned int num_nodes, const Vec2 *vec_r, unsigned int num_elements, const Params &params );
    void AddElement( unsigned int i, unsigned int j, unsigned int k );
    //\todo void AddElement_ConstitutiveModel( unsigned int nid1, unsigned int nid2, unsigned int nid3, constitutive_model_specific_params );
    bool EndEdition();
    void SetParams( const Params &params );

private:
    unsigned int m_LastElement;
    Vec2 *m_pAllocRefPos;
    ElementData *m_pAllocED;
    LinearED *m_pAllocLED;
};

}}} //namespace S2::ms::fem

#endif //S2_MS_FEM_TRIANGLE_ELEMENT_SET_2_H
