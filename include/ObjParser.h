#ifndef ____LUX_PLYPARSER_H____
#define ____LUX_PLYPARSER_H____
 
#include <string>
using namespace std;
//#include "AsciiParser.h"
//#include "Noise.h"
#include "VolumeGeometry.h"
 
namespace lux
{
 
class TriangleGeometry;
    
class ObjParser
{
  public:
 
    ObjParser(){}
   ~ObjParser(){}
 
    
    const void split(const std::string& s, char delim, std::vector<std::string>& tokens);
    void load(const std::string& filename, Mesh m);

    const bool IsComment(const std::string& str) const  { return (str.substr(0,2) == "# ") ? true : false; }
    const bool IsTexture(const std::string& str) const { return (str.substr(0,3) == "vt ") ? true : false; }
    const bool IsVertex(const std::string& str) const { return (str.substr(0,2) == "v ") ? true : false; }
    const bool IsNormal(const std::string& str) const { return (str.substr(0,3) == "vn ") ? true : false; }
    const bool IsFace(const std::string& str) const { return (str.substr(0,2) == "f ") ? true : false; }
    const bool IsGroup(const std::string& str) const { return (str.substr(0,2) == "g ") ? true : false; }

  private:
    int nb_vertices;
    int nb_faces;

     
};
 
}
 
#endif
