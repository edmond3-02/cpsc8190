#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <string>

namespace img
{

bool write( const std::string& filename, float* img_data, int Nx, int Ny, int Nc );

}

#endif
