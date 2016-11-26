#ifndef S2_MS_FEM_TETRAHEDRON_ELEMENT_SET_3_H
#define S2_MS_FEM_TETRAHEDRON_ELEMENT_SET_3_H

#include "../Config.h"
#include "Params.h"
#include "TetrahedronElement3.h"
#include <util/GPoolDA.h>

//#define __ENABLE_TRACE_FEM
#ifdef __ENABLE_TRACE_FEM
#  include <iostream>
#  include <Mal/GSerialization.h>
#endif

namespace S2 {
namespace ms {

/*! Generic Force-On-Particle-System interface

  \note ALL methods ACCUMULATE on input arrays to avoid temporaries,
  they should be EXTERNALLY RESET to 0 if no accumulation desired.

  \note Notation:
  dx = differential of x
  df_x = differential of function f(x,v,t) due to differential in variable x
  PD_f_x = partial derivative of f(x,v,t) wrt x
  \todo Some forces need extra per-particle data (gravity => mass, electric => charge, etc...)
*/
class IForce3
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
    IForce3() {}
    virtual ~IForce3() {}

    //\name IForce virtual API
    //@{
    virtual void BeginEvaluation( ICache* p_cache, const Vec3 *vec_prev_x, const Vec3 *vec_x, const Vec3 *vec_v, Real t ) const = 0;
    virtual void f( ICache* p_cache, Vec3 *vec_f ) const = 0; //! += Force f(x,v,t)
    virtual void df_x( ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const = 0; //! += Force differential df_x = \frac{\partial f}{partial x}(x,t) * dx
    virtual void df_v( ICache* p_cache, const Vec3 *vec_dv, Vec3 *vec_df ) const = 0; //! += Force differential df_v = \frac{\partial f}{partial v}(x,v,t) * dv
    virtual Real V( ICache* p_cache ) const = 0; //! Potential energy V(x,v,t)
    //virtual void PD_f_x( const Vec3 *vec_x, Real t ) const = 0; //! Matrix D_x(f) = \frac{\partial f}{partial x}(x,t)
    //virtual void PD_f_v( const Vec3 *vec_x, const Vec3 *vec_v, Real t ) const = 0; //! Matrix D_x(f) = \frac{\partial f}{partial v}(x,v,t)
    virtual void EndEvaluation( ICache* p_cache ) const = 0;
    //@}
};

}} //namespace S2::ms

namespace S2 {
namespace ms {
namespace fem {

/*! Constant Tetrahedron Element Set 2D

  Supports several FEM Params::EMethod that can be changed on the fly.

*/
class TetrahedronElementSet3: public IForce3
{
public:
    typedef uint16 element_index_type; //\todo Could parametrize or upgrade to uint32 for large meshes
    typedef uint16 node_index_type; //\todo Could parametrize or upgrade to uint32 for large meshes
    //enum EConstants { cInvalidNoC = -1 };

public:
    TetrahedronElementSet3();
    virtual ~TetrahedronElementSet3();

    inline unsigned int GetNumElements() const { return m_NumElements; }
    inline bool IsValid() const { return 0 != m_vecRefPos; } //\note otherwise, it has NOT been initialized with SetBakedData()

    //! \name Specific cache creation/Destruction (eg: needs node inv masses for rayleigh damping M term)
    //@{
    IForce3::ICache* CreateCache( const Real* vec_inv_mass );
    void DestroyCache( IForce3::ICache* p_cache );
    //@}

    //\name IForce virtual API
    //@{
    void BeginEvaluation( IForce3::ICache* p_cache, const Vec3 *vec_prev_x, const Vec3 *vec_x, const Vec3 *vec_v, Real t ) const;
    void f( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void df_v( IForce3::ICache* p_cache, const Vec3 *vec_dv, Vec3 *vec_df ) const;
    Real V( IForce3::ICache* p_cache ) const;
    void EndEvaluation( IForce3::ICache* p_cache ) const;
    //@}

    //! Elastic force \todo Consider adding to f() implicitly OR add f_Elastic() to IForce virtual API and add vec_node_mass to BeginEvaluation() params
    void f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    //! Damping force \todo Consider adding to f() implicitly OR add f_Damping() to IForce virtual API and add vec_node_mass to BeginEvaluation() params
    void f_d( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    //! \todo Plastic force \todo Consider adding to f() implicitly OR add f_Plastic() to IForce virtual API and add vec_node_mass to BeginEvaluation() params
    void f_p( IForce3::ICache* p_cache, Vec3 *vec_f ) const;

    //\name Consultors
    //@{
    inline Real TotalVolume() const { return m_TotalVolume; }
    inline Real Volume( unsigned int eid ) const { return m_vecED[eid].m_Volume; }
    void GetNID( unsigned int eid, unsigned int &nid0, unsigned int &nid1, unsigned int &nid2, unsigned int &nid3 ) const;
    Transform3 GetTransformE2W( const IForce3::ICache* p_cache, unsigned int eid, const Vec3* vec_pos ) const;
    //inline unsigned int GetNoC( unsigned int eid ) const { return m_vecEC[eid].m_NoC; } //TEMP
    TetrahedronElement3::DoC GetDoC( const IForce3::ICache* p_cache, unsigned int eid ) const;
    //@}

public:
    struct ElementData;
    struct LinearED;

protected:
    void SetBakedData( bool b_shared,
                       uint32 num_nodes, uint32 num_elements,
                       const Vec3 *vec_ref_pos, const ElementData *vec_ed, const LinearED *vec_led );
    void ClearBakedData();

private:
    void UpdateDegeneration( IForce3::ICache* p_cache ) const;

    //\name Corotational R computation
    //@{
    void ComputeRotations_Id( IForce3::ICache* p_cache ) const;
    void ComputeRotations_QR( IForce3::ICache* p_cache ) const;
    void ComputeRotations_MSVD( IForce3::ICache* p_cache ) const;
    void ComputeRotations_PD_Project( IForce3::ICache* p_cache ) const;
    void ComputeRotations_PD_Reflect( IForce3::ICache* p_cache ) const;
    void ComputeRotations_DAPD( IForce3::ICache* p_cache ) const;
    //@}

    //\name f_e()
    //@{
    void L_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void C_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void C_LCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void C_CCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void H_LCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void H_CCM_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void H_NHC0_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void H_NHC1_f_e( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    //@}

    //\name V()
    //@{
    Real L_V( IForce3::ICache* p_cache ) const;
    Real C_V( IForce3::ICache* p_cache ) const;
    Real C_LCM_V( IForce3::ICache* p_cache ) const;
    Real C_CCM_V( IForce3::ICache* p_cache ) const;
    Real H_LCM_V( IForce3::ICache* p_cache ) const;
    Real H_CCM_V( IForce3::ICache* p_cache ) const;
    Real H_NHC0_V( IForce3::ICache* p_cache ) const;
    Real H_NHC1_V( IForce3::ICache* p_cache ) const;
    //@}

    //\name f_d()
    //@{
    void L_f_d( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    void C_f_d( IForce3::ICache* p_cache, Vec3 *vec_f ) const;
    //@}

    //\name df_x()
    //@{
    void L_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void C_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void C_LCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void C_LCMH_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void C_CCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void H_LCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    void H_CCM_df_x( IForce3::ICache* p_cache, const Vec3 *vec_dx, Vec3 *vec_df ) const;
    //@}

public:
    //! ElementData, constant during entire element lifetime
    struct ElementData
    {
        node_index_type m_vecNID[4];
        Mat3x3 m_InvDm; //D_m^-1
        Real m_Volume;
    };
    //! LFEM specific data
    struct LinearED
    {
        //Essential
        Mat3x3 m_K00, m_K01, m_K02, m_K03;
        Mat3x3 m_K11, m_K12, m_K13;
        Mat3x3 m_K22, m_K23;
        Mat3x3 m_K33;
        //Transposed
        Mat3x3 m_K10, m_K20, m_K21, m_K30, m_K31, m_K32;
        /*\todo Alternatively: Mat3x3 m_matK[3][3], so that K_ij are
          directly accessible.
        */
    };

protected:
    //constant representation
    node_index_type m_NumNodes;
    element_index_type m_NumElements;
    Real m_TotalVolume;

    //\name params
    //@{
    Params m_Params;
    // Derived
    Real m_LameMu, m_LameLambda;
    //@}

    const Vec3 *m_vecRefPos;
    const ElementData *m_vecED;
    const LinearED *m_vecLED;

    //uint32 *m_pBuffer; //!< If exists, the baked data is NOT SHARED
#ifdef __S2_MS_ENABLE_STATS
public:
    struct Stats
    {
        uint32 m_Num_Degenerate;
        float32 m_Sum_Degenerate_DetF;
    };
    void ComputeStats( const IForce3::ICache* p_cache, Stats& stats ) const;
#endif
};

/*! Editable Tetrahedron Element Set 3D
  - Construction: Add/Remove elements
  - Modification: Change parameters

  \todo Despite selected EMethod, we precompute all constant data
  required by any method, as we may change it on the fly later.

  \note We follow the MeshSolidShape2 Baked/Editable data paradigm
*/
class EditableTetrahedronElementSet3: public TetrahedronElementSet3
{
public:
    EditableTetrahedronElementSet3();
    ~EditableTetrahedronElementSet3();

    void ClearEditData();
    //\todo if method is given per-element, params.m_Method == eMethod_PerElement and alloc a per-element method descriptor
    void BeginEdition( unsigned int num_nodes, const Vec3 *vec_r, unsigned int num_elements, const Params &params );
    void AddElement( unsigned int i, unsigned int j, unsigned int k, unsigned int l );
    //\todo void AddElement_ConstitutiveModel( unsigned int nid1, unsigned int nid2, unsigned int nid3, constitutive_model_specific_params );
    bool EndEdition();
    void SetParams( const Params &params );

private:
    unsigned int m_LastElement;
    Vec3 *m_pAllocRefPos;
    ElementData *m_pAllocED;
    LinearED *m_pAllocLED;
};

}}} //namespace S2::ms::fem

#endif //S2_MS_FEM_TETRAHEDRON_ELEMENT_SET_3_H
