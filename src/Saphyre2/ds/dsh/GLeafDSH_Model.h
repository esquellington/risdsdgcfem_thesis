#ifndef S2_DS_DSH_GLEAFDSHMODEL_H
#define S2_DS_DSH_GLEAFDSHMODEL_H

#include <Saphyre2/ds/dsh/IDynamicSystemHierarchy.h>

namespace S2 { namespace ds {

/*! Generic Model Leaf DSH
  
  Basic implementation of ILeafDSH methods building on common Model
  functionality.

  If there is no bound ISimulationScheme, Update() is automatically
  forwarded to ModelT. Notice, however, that the ModelT may ignore
  DSH-level defined Constraints and Geoms, and most times will just
  perform an integration step.
*/
template <class ModelT>
class GLeafDSH_Model: public GIDynamicSystemHierarchyD<ModelT::cDimension,ILeafDSH>
{
public:
    enum EConstants { cDimension = ModelT::cDimension };
    typedef ModelT model_type;
    typedef GIDynamicSystemHierarchyD<ModelT::cDimension,ILeafDSH> base_type;
    
public:
    GLeafDSH_Model( uint32 uid, IDynamicSystemHierarchy *p_parent ) : base_type(uid,p_parent) {}

    void Step( Real dt ) { m_Model.Update(dt); }
    
    inline ModelT &GetModel() { return m_Model; }
    inline const ModelT &GetModel() const { return m_Model; }

    void QueryState( util::ItemStream &rets ) const
    {
        rets.BeginComplex( 666, (uint32)eRet_Update ); //\todo cid = 666??!!
        {
            rets.Write("uid",base_type::GetUID());
            Real *p_state = rets.AllocArray<Real>("state",m_Model.GetStateSize());
            m_Model.GetState(p_state);
        }
        rets.EndComplex();
    }
        
protected:
    ModelT m_Model;
};

}} // namespace S2::ds

#endif // S2_DS_DSH_GLEAFDSHMODEL_H

