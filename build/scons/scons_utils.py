# File listing and filtering utilities

#Collects all files from root_dir recursively and returns them with:
#include_root == True: FULL path names (from /)
#include_root == False: LOCAL (from root_dir) 
def RecursiveFileList( root_dir, include_root = True ):
    import os.path
    file_list = []
    for root, dirs, files in os.walk( root_dir ):
        if include_root:
            for name in files:            
                file_list.append( os.path.join(root, name) )
        else:
            for name in files:
                file_list.append( os.path.join( root.split(root_dir)[1], name) )
        #for name in dirs: print "paths: ", os.path.join(root, name)
    return file_list

def FilterStringList( string_list, pattern ):
    import re
    test_re = re.compile( pattern, re.IGNORECASE )
    return filter(test_re.search, string_list)

def PrependStringList( prefix, string_list ):
    result = []
    for name in string_list:
        result.append( prefix + name )
    return result
