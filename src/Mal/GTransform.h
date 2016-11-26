/*! \file GTransform.h
\brief Geometric Transform {Trans,Rot}, Scale not possible by now
*/

#ifndef MAL_GTRANSFORM_H
#define MAL_GTRANSFORM_H

#include <Mal/Config.h>
#include <Mal/RealUtils.h>
#include <Mal/GVec.h>
#include <Mal/GMat.h>

namespace mal
{

/*! Rigid-body transform (pos,rot) */
template<typename T, unsigned D>
class GTransform
{
public:
    typedef T real_type;
    typedef GVec<T,D> vec_type;
    typedef GMat<T,D,D> mat_type;

    static const unsigned dimension = D;
    static const unsigned int size_in_reals = vec_type::size_in_reals + mat_type::size_in_reals;

    finline static GTransform Identity() { return GTransform( vec_type::Zero(), mat_type::Identity() ); }

public:
    vec_type m_Pos;
    mat_type m_Rot;

public:
    //!\name Construction
    //@{
    finline GTransform() {}
    finline GTransform( const GTransform &transform ) : m_Pos(transform.m_Pos), m_Rot(transform.m_Rot) {}
    finline GTransform( const vec_type &p, const mat_type &m ) : m_Pos(p), m_Rot(m) {}
    finline explicit GTransform( const T* array ) : m_Pos(&array[0]), m_Rot(&array[vec_type::size_in_reals]) {}
    //@}

    //!\name Access
    //@{
    finline const vec_type &Pos() const { return m_Pos; }
    finline const mat_type &Rot() const { return m_Rot; }
    //@}

    //!\name Add/Sub (NONE BY NOW, maybe needed for skinning later)
    //@{
    //@}

    //!\name Products
    //@{
    finline GTransform operator*( const GTransform &transform ) const { return GTransform( (*this) * transform.m_Pos,
                                                                                           m_Rot * transform.m_Rot ); }
    //! *= NOT AVAILABLE because it is ambiguous: this = transform * this Vs this = this * transform
    //finline void operator*=( const GTransform &transform );

    finline vec_type operator*( const vec_type &v ) const { return m_Rot*v + m_Pos; }
    //@}

    //!\name Unary
    //@{
    finline GTransform Inverse() const { GTransform inv;
                                         inv.m_Rot = m_Rot.Transposed();
                                         inv.m_Pos = - inv.m_Rot * m_Pos;
                                         return inv; }
    finline void Invert() { m_Rot.Transpose();
                            m_Pos = - m_Rot * m_Pos; }
    //@}

    //!\name Norms
    //@{
    //@}

    //!\name Utility
    //@{
    finline void ToArray( T *real_array ) const { m_Pos.ToArray(&real_array[0]);
                                                  m_Rot.ToArray(&real_array[vec_type::size_in_reals]); }
    finline void FromArray( const T *real_array ) { m_Pos.FromArray(&real_array[0]);
                                                    m_Rot.FromArray(&real_array[vec_type::size_in_reals]); }
    //@}
};

//---- Misc
template <typename T, unsigned D>
finline bool IsNaN( const GTransform<T,D> &tr )
{
    return IsNaN(tr.m_Pos) || IsNaN(tr.m_Rot);
}

template <typename T,unsigned D>
finline GTransform<T,D> Inverse( const GTransform<T,D> &tr )
{
    return tr.Inverse(); //\todo Consider explicit implementation here, no method call
}

} // namespace mal

#endif // MAL_GTRANSFORM_H
