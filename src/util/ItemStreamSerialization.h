#ifndef UTIL_ITEM_STREAM_SERIALIZATION_H
#define UTIL_ITEM_STREAM_SERIALIZATION_H

#include <util/ItemStream.h>
#include <Mal/GSerialization.h>
#include <iostream>

namespace util {

inline std::ostream &operator<<( std::ostream &o_stream, const util::ItemStream::ItemIt &it_its )
/*
template <typename OstreamT>
inline OstreamT &operator<<( OstreamT &o_stream, const util::ItemStream::ItemIt &it_its )
*/
{
    if( !it_its.IsValid() )
    {
        o_stream << "Invalid Item!" << std::endl;
        return o_stream;
    }

    if( it_its.IsSimple() )
    {
        if( it_its.GetId() < 0 ) o_stream << "< Simple[\"" << it_its.GetName() << "\"] T" << it_its.GetType() << " S" << it_its.GetSize() << " >" << std::endl;
        else o_stream << "< Simple[#" << it_its.GetId() << "] T" << it_its.GetType() << " S" << it_its.GetSize() << " >" << std::endl;

        o_stream << "\tValue = ";
        switch( it_its.GetType() )
        {
        case eType_Int64: o_stream << it_its.Get<int64>(); break;
        case eType_UInt64: o_stream << it_its.Get<uint64>(); break;
        case eType_Int32: o_stream << it_its.Get<int32>(); break;
        case eType_UInt32: o_stream << it_its.Get<uint32>(); break;
        case eType_Int16: o_stream << it_its.Get<int16>(); break;
        case eType_UInt16: o_stream << it_its.Get<uint16>(); break;
        case eType_Int8: o_stream << (int)it_its.Get<int8>(); break;
        case eType_UInt8: o_stream << (int)it_its.Get<uint8>(); break;

        case eType_Bool32: o_stream << ( it_its.Get<bool>() ? "_true" : "_false" ); break;

        case eType_Float32: o_stream << it_its.Get<float32>(); break;
        case eType_Float64: o_stream << it_its.Get<float64>(); break;

        case eType_Flags8: o_stream << std::hex << it_its.Get<Flags8>().m_Bits << std::dec; break;
        case eType_Flags16: o_stream << std::hex << it_its.Get<Flags16>().m_Bits << std::dec; break;
        case eType_Flags32: o_stream << std::hex << it_its.Get<Flags32>().m_Bits << std::dec; break;

        case eType_String16: o_stream << "'" << it_its.Get<String16>().GetStr() << "'16"; break;
        case eType_String32: o_stream << "'" << it_its.Get<String32>().GetStr() << "'32"; break;
        case eType_String64: o_stream << "'" << it_its.Get<String64>().GetStr() << "'64"; break;

        case eType_Vec2f: o_stream << it_its.Get<Vec2f>(); break;
        case eType_Vec3f: o_stream << it_its.Get<Vec3f>(); break;
        case eType_Quatf: o_stream << it_its.Get<Quatf>(); break;

        default: o_stream << "Non-basic type!"; break;
        }
        o_stream << std::endl;
    }
    else if( it_its.IsArray() )
    {
        if( it_its.GetId() < 0 ) o_stream << "< Array[\"" << it_its.GetName() << "\"] T" << it_its.GetType() << " S" << it_its.GetSize() << " >" << std::endl;
        else o_stream << "< Array[#" << it_its.GetId() << "] T" << it_its.GetType() << " S" << it_its.GetSize() << " >" << std::endl;

        if( eType_Int8 == it_its.GetType() )
            o_stream << "\tStr[" << it_its.GetArrayCount() << "] = '"
                     << it_its.GetArrayPtr<char>() << "'" << std::endl;
        else
        {
            o_stream << "\tCount = " << it_its.GetArrayCount() << std::endl;
            o_stream << "\tValue = [ ";
            for( unsigned i=0; i<it_its.GetArrayCount(); i++ )
            {
                switch( it_its.GetType() )
                {
                case eType_Int64: o_stream << it_its.GetArrayPtr<int64>()[i] << " "; break;
                case eType_UInt64: o_stream << it_its.GetArrayPtr<uint64>()[i] << " "; break;
                case eType_Int32: o_stream << it_its.GetArrayPtr<int32>()[i] << " "; break;
                case eType_UInt32: o_stream << it_its.GetArrayPtr<uint32>()[i] << " "; break;
                case eType_Int16: o_stream << it_its.GetArrayPtr<int16>()[i] << " "; break;
                case eType_UInt16: o_stream << it_its.GetArrayPtr<uint16>()[i] << " "; break;
                case eType_Int8: o_stream << it_its.GetArrayPtr<int8>()[i] << " "; break;
                case eType_UInt8: o_stream << it_its.GetArrayPtr<uint8>()[i] << " "; break;

                case eType_Float32: o_stream << it_its.GetArrayPtr<float32>()[i] << " "; break;
                case eType_Float64: o_stream << it_its.GetArrayPtr<float64>()[i] << " "; break;

                case eType_Vec2f: o_stream << it_its.GetArrayPtr<Vec2f>()[i] << " "; break;
                case eType_Vec3f: o_stream << it_its.GetArrayPtr<Vec3f>()[i] << " "; break;
                case eType_Quatf: o_stream << it_its.GetArrayPtr<Quatf>()[i] << " "; break;

                default: o_stream << "Non-basic type!"; break;
                }
            }
            o_stream << "]" << std::endl;
        }
    }
    else if( it_its.IsComplex() )
    {
        if( it_its.GetId() < 0 ) o_stream << "< Complex[\"" << it_its.GetName() << "\"] T" << it_its.GetType() << " S" << it_its.GetSize() << " >" << std::endl;
        else o_stream << "< Complex[#" << it_its.GetId() << "] T" << it_its.GetType() << " S" << it_its.GetSize() << " >" << std::endl;
        for( util::ItemStream::ItemIt it_subitem=it_its.GetSubItem();
             it_subitem.IsValid();
             ++it_subitem )
            o_stream << it_subitem;
        o_stream << "</ Complex >" << std::endl;
    }
    return o_stream;
}

//! ItemStream dump
inline std::ostream &operator<<( std::ostream &o_stream, const util::ItemStream &its )
/*
template <typename OstreamT>
inline OstreamT &operator<<( OstreamT &o_stream, const util::ItemStream &its )
*/
{
    for( util::ItemStream::ItemIt it_its=its.Begin();
         it_its.IsValid();
         ++it_its )
        o_stream << it_its;
    return o_stream;
}

} //namespace util

namespace util
{

bool ItemStreamToFileTxt( const ItemStream &is, FILE *p_file );
bool ItemStreamFromFileTxt( ItemStream &is, FILE *p_file );

} //namespace util

#endif // UTIL_ITEM_STREAM_SERIALIZATION_H
