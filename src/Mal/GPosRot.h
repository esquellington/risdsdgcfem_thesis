/*! \file GPosRot.h
\brief Rigid Body Position,Orientation {Trans,Rot}

\note Some UNTESTED functions fire an assert(false) to FORCE checking
their results on the first use }-D
*/

#ifndef MAL_GPOSROT_H
#define MAL_GPOSROT_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>
#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GQuatUtils.h> //!< Quat-Vec composition

namespace mal
{

template<typename T>
class GPosRot
{
public:
    typedef T real_type;
    typedef GVec<T,3> vector3_type;
    typedef GQuat<T> quat_type;
    static const unsigned int size_in_reals = vector3_type::size_in_reals + quat_type::size_in_reals ;
    
    finline static GPosRot Identity() { return GPosRot( vector3_type::Zero(), quat_type::Identity() ); }

public:
    vector3_type m_Pos;
    quat_type m_Rot;
    
public:
    //!\name Construction
    //@{
    finline GPosRot() {}
    finline GPosRot( const GPosRot &transform ) : m_Pos(transform.m_Pos), m_Rot(transform.m_Rot) {}
    finline GPosRot( const vector3_type &p, const quat_type &q ) : m_Pos(p), m_Rot(q) {}
    finline GPosRot( const T* array ) : m_Pos(&array[0]), m_Rot(&array[3]) {}
    //@}

    //!\name Access
    //@{
    //@}

    //!\name Add/Sub (NONE BY NOW, maybe needed for skinning later)
    //@{
    //@}
    
    //!\name Products
    //@{
    finline GPosRot operator*( const GPosRot &transform ) const { MAL_ASSERT(false);
                                                                  return GPosRot( (*this) * transform.m_Pos,
                                                                                  m_Rot * transform.m_Rot ); }
    //! *= NOT AVAILABLE because it is ambiguous: this = transform * this Vs this = this * transform
    //finline void operator*=( const GPosRot &transform );

    finline vector3_type operator*( const vector3_type &v ) const { MAL_ASSERT(false); return vector3_type::Zero(); }
    //@}

    //!\name Unary
    //@{
    finline GPosRot Inverse() const { MAL_ASSERT(false); return Identity(); }
    finline void Invert() { MAL_ASSERT(false); }
    //@}    

    //!\name Norms
    //@{
    //@}
    
    //!\name Utility
    //@{
    finline void ToArray( T *array ) const { m_Pos.ToArray(&array[0]); m_Rot.ToArray(&array[3]); }
    //@}
};

//---- Misc
template <typename T>
finline bool IsNaN( const GPosRot<T> &pq )
{
    return IsNaN(pq.m_Pos) || IsNaN(pq.m_Rot);
}

} // namespace mal

#endif // MAL_GPOSROT_H
