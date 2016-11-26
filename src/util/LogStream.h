#ifndef UTIL_LOG_STREAM_H
#define UTIL_LOG_STREAM_H

#include <fstream>
#include <string>
#include <vector>

#include <util/Chronos.h>
#include <util/ItemStream.h>

namespace util {

//! \name Base Tags
//@{
class BTag { public: BTag( const std::string &tag ) : m_Tag( tag ) {} const std::string m_Tag; };
class ETag { public: ETag() {} };
class Endl { public: Endl() {} }; //!< std::endl replacement
//@}

//! \name Simple Info Tags
//@{
class BInfo: public BTag { public: BInfo() : BTag( "info" ) {} };
typedef ETag EInfo;
class BWarning: public BTag { public: BWarning() : BTag( "warning" ) {} };
typedef ETag EWarning;
class BError: public BTag { public: BError() : BTag( "error" ) {} };
typedef ETag EError;
//@}

//! \name Complex Info Tags
//@{
//! Utility Tag <tag name= "NAME" >
class BTagWithName: public BTag
{
public:
    BTagWithName( const std::string &tag, const std::string &name ) : BTag( tag ), m_Name( name ) {}
    const std::string m_Name;
};

//! Utility Tag to signal that TagDuration must be output
class ETagWithDuration: public ETag { public: ETagWithDuration() {} };

//! <cmd> Command Tag with a Command Name
class BCmd: public BTagWithName { public: BCmd( const std::string &cmd_name ) : BTagWithName("cmd",cmd_name) {} };
typedef ETagWithDuration ECmd;
//! <section> Command Tag with a Section Name
class BSection: public BTagWithName { public: BSection( const std::string &section_name ) : BTagWithName("section",section_name) {} };
typedef ETagWithDuration ESection;
//@}

//! Log Stream class to output structured (Xml by now) data.
/*! The LogStream class defines most of std::ostream output
  operators (and friend functions) and new ones for structured output.

  Uses an std::ofstream internally to do the final output on a file. */
class LogStream
{    
public:
    LogStream();
    ~LogStream();

    bool Open( const std::string &file, const std::string &title );
    void Close();
    bool IsOpen();

    //!\name Nested sections/tags management
    //@{
    void BeginTag( const std::string &tag, const std::string &name = "" );
    void EndTag( bool b_compute_duration = false );

    void BeginSection( const std::string &title ) { (*this) << BSection(title); }
    void EndSection() { (*this) << ESection(); }
    //@}

    //!\name Timing stuff
    //@{
    void SetEnabledTimeStamp( bool b_enabled ) { m_bEnabledTimeStamp = b_enabled; }
    //@}
    
    //!\name Tag serialization
    //@{
    LogStream& operator<< ( const BTag &btag ) { BeginTag( btag.m_Tag ); return *this; }
    LogStream& operator<< ( const ETag &etag ) { EndTag(); return *this; }
    LogStream& operator<< ( const BTagWithName &btwn ) { BeginTag( btwn.m_Tag, btwn.m_Name ); return *this; }
    LogStream& operator<< ( const ETagWithDuration &etwd ) { EndTag( true ); return *this; }
    LogStream& operator<< ( const Endl &e );
    //@}

    LogStream& operator<< ( const ItemStream &its );
    LogStream& operator<< ( const ItemStream::ItemIt &it );
    
    //!\name Standard ostream output operators and friend methods
    //@{
    LogStream& operator<< ( const bool& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const short& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const unsigned short& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const int& val ) { m_Ofstream << val; return *this; }  
    LogStream& operator<< ( const unsigned int& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const long& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const unsigned long& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const long long& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const unsigned long long& val ) { m_Ofstream << val; return *this; }    
    LogStream& operator<< ( const float& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const double& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const long double& val ) { m_Ofstream << val; return *this; }
    LogStream& operator<< ( const void*& val ) { m_Ofstream << val; return *this; }   
    
    friend LogStream& operator<< ( LogStream& ls, char ch );
    friend LogStream& operator<< ( LogStream& ls, signed char ch );
    friend LogStream& operator<< ( LogStream& ls, unsigned char ch );
    friend LogStream& operator<< ( LogStream& ls, const char* str );
    friend LogStream& operator<< ( LogStream& ls, const signed char* str );
    friend LogStream& operator<< ( LogStream& ls, const unsigned char* str );
    friend LogStream& operator<< ( LogStream& ls, const std::string& str );
    friend LogStream& operator<< ( LogStream& ls, void* value );    
    //@}
    
private:    
    unsigned int m_TabLevel;  //!< Current level of indentation
    std::ofstream m_Ofstream; //!< Internal ofstream used for actual output
    
    std::vector< std::string > m_stackTags;

    bool m_bEnabledTimeStamp;
    Chronos m_Chronos;
    std::vector< double > m_stackTagTimeStamps;
};

//external operator<<  overloaded functions:
inline LogStream& operator<< ( LogStream& ls, char ch )
{ ls.m_Ofstream << ch; return ls; }
inline LogStream& operator<< ( LogStream& ls, signed char ch )
{ ls.m_Ofstream << ch; return ls; }
inline LogStream& operator<< ( LogStream& ls, unsigned char ch )
{ ls.m_Ofstream << ch; return ls; }
inline LogStream& operator<< ( LogStream& ls, const char* str )
{ ls.m_Ofstream << str; return ls; } 
inline LogStream& operator<< ( LogStream& ls, const signed char* str )
{ ls.m_Ofstream << str; return ls; }
inline LogStream& operator<< ( LogStream& ls, const unsigned char* str )
{ ls.m_Ofstream << str; return ls; }
inline LogStream& operator<< ( LogStream& ls, const std::string& str )
{ ls.m_Ofstream << str; return ls; }

inline LogStream& operator<< ( LogStream& ls, void* value )
{ ls.m_Ofstream << std::hex << value << std::dec; return ls; }

} // namespace util

#endif // UTIL_LOG_STREAM_H
