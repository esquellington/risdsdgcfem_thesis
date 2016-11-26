#include <util/LogStream.h>
//#include <util/ItemStreamSerialization.h>  //\todo See comments in ItemStream-related operators <<

namespace util {

//----- LogStream Implementation 
LogStream::LogStream()
: m_TabLevel(0)
, m_bEnabledTimeStamp(true)
{
}

LogStream::~LogStream()
{
    Close();
}

bool LogStream::Open( const std::string &file, const std::string &title )
{    
    if( IsOpen() )
        Close();

    m_Chronos.ResetTime();
    
    if( file != "" )
    {
        m_Ofstream.open(file.c_str());
        m_Ofstream << "<file name=\"" << file << "\" title= \"" << title << "\" >"
                   << std::endl;
        return true;
    }
    else
        return false;
}

void LogStream::Close()
{    
    while( m_stackTags.size() )
        EndTag();
    m_Ofstream << "</file>" << std::endl;
    
    m_Ofstream.close();
}

bool LogStream::IsOpen()
{
    return m_Ofstream.is_open();
}

void LogStream::BeginTag( const std::string &tag, const std::string &name )
{
    m_stackTags.push_back( tag );
    if( m_bEnabledTimeStamp )
        m_stackTagTimeStamps.push_back( m_Chronos.GetTime() );
    
    (*this) << "<" << tag;
    if( !name.empty() )
        (*this) << " name= \"" << name << "\" ";
    if( m_bEnabledTimeStamp )
        (*this) << " time= \"" << m_stackTagTimeStamps.back() << "\"";
    (*this) << ">";
    m_TabLevel++;
    (*this) << Endl();
}

//!< std::endl replacement
void LogStream::EndTag( bool b_compute_duration )
{
    m_TabLevel--;
    (*this) << Endl();
    if( m_bEnabledTimeStamp && b_compute_duration )
        (*this) << "<duration> " << (m_Chronos.GetTime() - m_stackTagTimeStamps.back() ) << " </duration>" << Endl();
    (*this) << "</" << m_stackTags.back() << ">" << Endl();
    m_stackTags.pop_back();
    if( m_bEnabledTimeStamp )
        m_stackTagTimeStamps.pop_back();
}

LogStream& LogStream::operator<< ( const Endl & )
{    
    m_Ofstream << std::endl;
    for( unsigned int i=0; i<m_TabLevel; i++)
        m_Ofstream << "\t";
    return *this;
}

LogStream& LogStream::operator<< ( const ItemStream &its )
{
    /* For some fucking reason this FAILS...
    m_Ofstream << its;
    return *this;
    */
    //\todo BY NOW, we just copy-paste the code from "ItemStreamSerialization.h"
    for( util::ItemStream::ItemIt it=its.Begin();
         it.IsValid();
         ++it )
        *this << it;
    return *this;
}

LogStream& LogStream::operator<< ( const ItemStream::ItemIt &it )
{
    /* For some fucking reason this FAILS...
    m_Ofstream << its;
    return *this;
    */
    //\todo BY NOW, we just copy-paste the code from "ItemStreamSerialization.h"
    if( !it.IsValid() )
    {
        m_Ofstream << "Invalid Item!" << std::endl;
        return *this;
    }
    
    if( it.GetId() < 0 )
        m_Ofstream << "Item[ " << it.GetName() << " ]" << std::endl;
    else
        m_Ofstream << "Item[ " << it.GetId() << " ]" << std::endl;
    
    m_Ofstream << "\tType = " << it.GetType() << std::endl;
    m_Ofstream << "\tSize = " << it.GetSize() << std::endl;
    if( it.IsSimple() )
    {        
        m_Ofstream << "\tValue = ";
        switch( it.GetType() )
        {
        case eType_Int32: m_Ofstream << it.Get<int32>(); break;
        case eType_UInt32: m_Ofstream << it.Get<uint32>(); break;
        case eType_Int16: m_Ofstream << it.Get<int16>(); break;
        case eType_UInt16: m_Ofstream << it.Get<uint16>(); break;
        case eType_Int8: m_Ofstream << (int)it.Get<int8>(); break;
        case eType_UInt8: m_Ofstream << (int)it.Get<uint8>(); break;
                
        case eType_Float32: m_Ofstream << it.Get<float32>(); break;
        case eType_Float64: m_Ofstream << it.Get<float64>(); break;

            /*
        case eType_Vec2f: m_Ofstream << it.Get<Vec2f>(); break;
        case eType_Vec3f: m_Ofstream << it.Get<Vec3f>(); break;
        case eType_Quatf: m_Ofstream << it.Get<Quatf>(); break;
            */
                
        default: m_Ofstream << "Non-basic type!"; break;
        }
        m_Ofstream << std::endl;
    }
    else if( it.IsArray() )
    {
        if( eType_Int8 == it.GetType() )
            m_Ofstream << "\tStr[" << it.GetArrayCount() << "] = '"
                     << it.GetArrayPtr<char>() << "'" << std::endl;
        else
        {
            m_Ofstream << "\tCount = " << it.GetArrayCount() << std::endl;
            m_Ofstream << "\tValue = [ ";
            for( unsigned int i=0; i<it.GetArrayCount(); i++ )
            {
                switch( it.GetType() )
                {
                case eType_Int32: m_Ofstream << it.GetArrayPtr<int32>()[i] << " "; break;
                case eType_UInt32: m_Ofstream << it.GetArrayPtr<uint32>()[i] << " "; break;
                case eType_Int16: m_Ofstream << it.GetArrayPtr<int16>()[i] << " "; break;
                case eType_UInt16: m_Ofstream << it.GetArrayPtr<uint16>()[i] << " "; break;
                case eType_Int8: m_Ofstream << it.GetArrayPtr<int8>()[i] << " "; break;
                case eType_UInt8: m_Ofstream << it.GetArrayPtr<uint8>()[i] << " "; break;
                    
                case eType_Float32: m_Ofstream << it.GetArrayPtr<float32>()[i] << " "; break;
                case eType_Float64: m_Ofstream << it.GetArrayPtr<float64>()[i] << " "; break;

                    /*
                      case eType_Vec2f: m_Ofstream << it.GetArrayPtr<Vec2f>()[i] << " "; break;
                      case eType_Vec3f: m_Ofstream << it.GetArrayPtr<Vec3f>()[i] << " "; break;
                      case eType_Quatf: m_Ofstream << it.GetArrayPtr<Quatf>()[i] << " "; break;
                    */
                
                default: m_Ofstream << "Non-basic type!"; break;
                }
            }
            m_Ofstream << "]" << std::endl;
        }
    }
    else if( it.IsComplex() )
    {
        m_Ofstream << "<Complex>" << std::endl;
        for( util::ItemStream::ItemIt it_subitem=it.GetSubItem();
             it_subitem.IsValid();
             ++it_subitem )
            *this << it_subitem;
        m_Ofstream << "</Complex>" << std::endl;
    }
    return *this;
}

} // namespace util
