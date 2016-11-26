#ifndef GEO_MP_S_DOUBLE_DISPATCHER
#define GEO_MP_S_DOUBLE_DISPATCHER

#include <Geo/ShapeTypes.h>
#include <Geo/shape/IShape.h>
#include <Geo/np/Context.h>
#include <Geo/np/Contact.h>

namespace geo {
namespace mp {

typedef bool (* test_contact_pair_function_2_type ) ( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                                                      const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                                                      np::ContactData2& cd, np::ContactCache2* p_cc,
                                                      const np::Context* p_context );

typedef bool (* test_contact_pair_function_3_type ) ( const IShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                                                      const IShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                                                      np::ContactData3& cd, np::ContactCache3* p_cc,
                                                      const np::Context* p_context );

/*\todo
typedef bool (* test_overlap_pair_function_2_type ) ( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                                                                   const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                                                                   const np::Context* p_context,
                                                                   np::OverlapCache2* p_oc );
typedef Real (* test_proximity_pair_function_2_type ) ( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                                                                     const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                                                                     const np::Context* p_context,
                                                                     np::ProximityData& pd, np::ProximityCache2* p_pc );
*/

//\todo Merge FS2 and FS3?? OR CONSIDER independent 2D and 3D DoubleDispatchers??
class DoubleDispatcher
{
public:
    /*! Pairwise Function set that holds all query function pointers
        and flags to determine if paramaters and results must be
        flipped.
    */
    struct FunctionSet2
    {
        test_contact_pair_function_2_type m_Contact;
        /*\todo
        test_overlap_pair_function_2_type m_Overlap;
        test_proximity_pair_function_2_type m_Proximity;
        */
        Flags32 m_Flags;
        enum EFlags { eFlipContact = 1, eFlipOverlap = 2, eFlipProximity = 4 };

        void SetFunction_Contact( test_contact_pair_function_2_type tcpf, bool b_flip )
            { m_Contact = tcpf; if(b_flip) m_Flags.Enable(eFlipContact); else m_Flags.Disable(eFlipContact); }
        /*\todo
        void SetFunction_Overlap( test_overlap_pair_function_2_type tcpf, bool b_flip )...
        void SetFunction_Proximity( test_proximity_pair_function_2_type tcpf, bool b_flip )...
        */
    };
    /*! Pairwise Function set that holds all query function pointers
        and flags to determine if paramaters and results must be
        flipped.
    */
    struct FunctionSet3
    {
        test_contact_pair_function_3_type m_Contact;
        /*\todo
        test_overlap_pair_function_2_type m_Overlap;
        test_proximity_pair_function_2_type m_Proximity;
        */
        Flags32 m_Flags;
        enum EFlags { eFlipContact = 1, eFlipOverlap = 2, eFlipProximity = 4 };

        void SetFunction_Contact( test_contact_pair_function_3_type tcpf, bool b_flip )
            { m_Contact = tcpf; if(b_flip) m_Flags.Enable(eFlipContact); else m_Flags.Disable(eFlipContact); }
        /*\todo
        void SetFunction_Overlap( test_overlap_pair_function_2_type tcpf, bool b_flip )...
        void SetFunction_Proximity( test_proximity_pair_function_2_type tcpf, bool b_flip )...
        */
    };

public:
    DoubleDispatcher();
    ~DoubleDispatcher();

    //\name DD definition
    //@{
    void MakeDefault();
    void MakeDefault_Contact();
    void SetFunction_Contact( EShapeType st1, EShapeType st2, test_contact_pair_function_2_type tcpf );
    void SetFunction_Contact( EShapeType st1, EShapeType st2, test_contact_pair_function_3_type tcpf );
    /*\todo
    void MakeDefault_Overlap();
    void MakeDefault_Proximity();

      void SetFunction_Overlap( EShapeType st1, EShapeType st2, test_overlap_pair_function_2_type topf );
      void SetFunction_Proximity( EShapeType st1, EShapeType st2, test_proximity_pair_function_2_type tppf );

    */
    //@}

    bool TestContact( const IShape2* p_shape1, const Transform2& tr1, const Real* p_dof1,
                      const IShape2* p_shape2, const Transform2& tr2, const Real* p_dof2,
                      np::ContactData2& cd, np::ContactCache2* p_cc );
    bool TestContact( const IShape3* p_shape1, const Transform3& tr1, const Real* p_dof1,
                      const IShape3* p_shape2, const Transform3& tr2, const Real* p_dof2,
                      np::ContactData3& cd, np::ContactCache3* p_cc );

    /*\todo
    bool TestOverlap();
    Real TestProximity();
    */

    np::Context* GetContext() { return &m_Context; }

protected:
    FunctionSet2 m_DDFST2D[cNumShapeTypes][cNumShapeTypes];
    FunctionSet3 m_DDFST3D[cNumShapeTypes][cNumShapeTypes];
    np::Context m_Context;
};

//! Global default DD
extern DoubleDispatcher* g_pDefaultDoubleDispatcher;

}} //namespace geo::mp

#endif //GEO_MP_S_DOUBLE_DISPATCHER
