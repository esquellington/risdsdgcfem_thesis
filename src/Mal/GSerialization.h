#ifndef MAL_GSERIALIZATION_H
#define MAL_GSERIALIZATION_H

#include <Mal/GVec.h>
#include <Mal/GQuat.h>
#include <Mal/GMat.h>
#include <Mal/GTransform.h>

namespace mal
{

template<typename OStreamT, typename T, int N>
inline OStreamT &operator<<( OStreamT &o_stream, const GVec<T,N> &v )
{
    o_stream << "( ";
    for(int i=0;i<N;i++) o_stream << v[i] << " ";
    o_stream << ")";
    return o_stream;
}

template<typename OStreamT, typename T>
inline OStreamT &operator<<( OStreamT &o_stream, const GQuat<T> &q )
{
    o_stream << "[ " << q.x << " " << q.y << " " << q.z << ", " << q.w << " ]";
    return o_stream;
}

template <typename OStreamT, typename T, int NR, int NC>
inline OStreamT &operator<<( OStreamT &o_stream, const GMat<T,NR,NC> &m )
{
    /*TEMP: Old format
    for( int i=0; i<NR; i++ )
    {
        o_stream << "[ ";
        for( int j=0; j<NC; j++ ) o_stream << m(i,j) << " ";
        o_stream << "]\n";
    }
    return o_stream;
    */
    // Ouput in Octave format for easy copy-paste
    o_stream << "[ ";
    for( int i=0; i<NR; i++ )
    {
        for( int j=0; j<NC; j++ )
            if( j==NC-1 )
                o_stream << m(i,j);
            else
                o_stream << m(i,j) << ", ";
        if( i==NR-1 )
            o_stream << " ]\n";
        else
            o_stream << " ;\n";
    }
    return o_stream;
}

template <typename OStreamT, typename T, unsigned N>
inline OStreamT &operator<<( OStreamT &o_stream, const GTransform<T,N> &tr )
{
    o_stream << "Pos\n" << tr.m_Pos << "\n";
    o_stream << "Rot\n" << tr.m_Rot << "\n";
    return o_stream;
}

} //namespace mal

#endif //MAL_GSERIALIZATION_H
