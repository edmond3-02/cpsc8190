#include "Image.h"

#include <iostream>

#include <OpenImageIO/imageio.h>

using namespace OIIO;

namespace img
{

bool write( const std::string& filename, float* img_data, int Nx, int Ny, int Nc ) 
{
  std::unique_ptr<ImageOutput> out = ImageOutput::create(filename.c_str());
  if( !out ){return false;}
  ImageSpec spec (Nx, Ny, Nc, TypeDesc::FLOAT);
  spec.attribute("user", "imgproc");
  spec.attribute("writer", "imgproc:write" );
  out->open (filename.c_str(), spec);

  unsigned long scanlinesize = Nx * Nc * sizeof(img_data[0]);
  out->write_image (TypeDesc::FLOAT, img_data);
  //out->write_image (TypeDesc::FLOAT, img_data + (Ny-1)*scanlinesize, AutoStride, -scanlinesize, AutoStride);
  out->close();
  return true;

}

}
