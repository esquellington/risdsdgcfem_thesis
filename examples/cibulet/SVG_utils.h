#ifndef CIBULET_SVG_UTILS_H
#define CIBULET_SVG_UTILS_H

#include "Config.h"
#include <string>
#include "rapidxml/rapidxml.hpp"
#include <Geo/shape/shape.h>

/*! Load an SVG path_element from a file into an EditablePathShape2. */
bool SVG_Load_EPS2( geo::EditablePathShape2& eps2,
                    std::string file_name, std::string path_element_name )
{
    // Load the SVG from disk, find shape definition d="...." and pass it to Make_PathShape2_SVG()
    FILE *pFile = fopen( file_name.c_str(), "rb" );
    APP_ASSERT( pFile );
    fseek( pFile, 0, SEEK_END ); // \todo non-portable ?
    size_t file_size = ftell( pFile) ;
    char *file_buffer = new char[file_size+1];
    rewind( pFile );
    fread( file_buffer, 1, file_size, pFile );
    file_buffer[file_size] = 0; //IMPORTANT: FORCE EOS, otherwise rapidxml sometimes crashes!

    rapidxml::xml_document<> doc; // character type defaults to char
    doc.parse<0>(file_buffer);    // 0 means default parse flags
    const rapidxml::xml_node<> *p_xml_node(0);
    // <svg
    p_xml_node = doc.first_node("svg");
    APP_LOG_ASSERT( 0 != p_xml_node, "Node 'svg' not found" );
    // <g
    p_xml_node = p_xml_node->first_node("g");
    APP_LOG_ASSERT( 0 != p_xml_node, "Node 'g' not found" );

    // <path
    // Find requested path
    p_xml_node = p_xml_node->first_node("path");
    APP_LOG_ASSERT( 0 != p_xml_node, "Node 'path' not found" );
    bool bFound(false);
    while( !bFound && 0 != p_xml_node )
    {
        if( std::string( p_xml_node->first_attribute("id")->value() ) == path_element_name )
            bFound = true;
        else
            p_xml_node = p_xml_node->next_sibling("path");
    }
    if( bFound )
    {
        // <d
        rapidxml::xml_attribute<> *p_xml_attribute(0);
        p_xml_attribute = p_xml_node->first_attribute("d");
        APP_LOG_ASSERT( 0 != p_xml_attribute, "Attribute 'd' not found" );
        const char *p_path_definition = p_xml_attribute->value();
        //APP_LOG( "Parsing SVG path definition '%s'", p_path_definition );

        geo::Make_PathShape2_SVG( eps2, p_path_definition );
        geo::Vec2 barycenter = eps2.Barycenter_0();
        eps2.BeginEdition();
        eps2.Transform( geo::Transform2( -barycenter, geo::Mat2x2::Identity() ) ) ;
        eps2.EndEdition();

        eps2.BeginEdition();
        eps2.Transform( geo::Transform2( geo::Vec2::Zero(), geo::Mat2x2(1,0,0,-1) ) ) ;
        eps2.EndEdition();
    }
    else
        APP_LOG_ERROR("SVG path element %s not found in file %s", path_element_name.c_str(), file_name.c_str() );

    delete[] file_buffer;

    return bFound;
}

#endif //CIBULET_SVG_UTILS_H
