#include "ImgProc.h"
using namespace img;
using namespace std;
using namespace OIIO;

#include "Color.h"

ImgProc::ImgProc() :
  Nx (0),
  Ny (0),
  Nc (0),
  Nsize (0),
  img_data (0)
  {}

ImgProc::~ImgProc()
{
  clear();
}

void ImgProc::clear()
{
  if( img_data != 0 ){ delete[] img_data; img_data = 0; }
  Nx = 0;
  Ny = 0;
  Nc = 0;
  Nsize = 0;
}
/*
void ImgProc::clear(int nX, int nY, int nC)
{
  clear();
  Nx = nX;
  Ny = nY;
  Nc = nC;
  Nsize = (long)Nx * (long)Ny * (long)Nc;
  img_data = new float[Nsize];
  #pragma omp parallel for
  for( long i=0;i<Nsize;i++ ){ img_data[i] = 0.0; };
}
*/
ImgProc::ImgProc(const ImgProc& v) :
  Nx (v.Nx),
  Ny (v.Ny),
  Nc (v.Nc),
  Nsize (v.Nsize)
{
  img_data = new float[Nsize];
  #pragma omp parallel for
  for( long i=0;i<Nsize;i++){ img_data[i] = v.img_data[i]; }
}

ImgProc& ImgProc::operator=(const ImgProc& v)
{
  if( this == &v )
  {
    return *this;
  }

  if( Nx != v.Nx || Ny != v.Ny || Nc != v.Nc )
  {
    clear();
    Nx = v.Nx;
    Ny = v.Ny;
    Nc = v.Nc;
    Nsize = v.Nsize;
    img_data = new float[Nsize];
  }

  #pragma omp parallel for
  for( long i=0;i<Nsize;i++)
  {
    img_data[i] = v.img_data[i];
  }
  return *this;
}

void ImgProc::operator*=( const float s )
{
#pragma omp parallel for
  for( long i=0;i<Nsize;i++ )
  {
    img_data[i] *= s;
  }
}

void ImgProc::operator+=( const float s )
{
  for( int j=0;j<Ny;j++ )
  {
#pragma omp parallel for
    for( int i=0;i<Nx;i++ )
    {
      long index = (long)i + (long)Nx * (long)j;
      for(int c=0;c<Nc;c++ )
      {
        img_data[ c + (long)Nc * index ] += s;
      }
    }
  }
}


void ImgProc::value( int i, int j, lux::Color& pixel) const
{
  //pixel.clear();
  if( img_data == 0 ){ return; }
  if( i<0 || i>=Nx ){ return; }
  if( j<0 || j>=Ny ){ return; }

  //pixel.resize(Nc);
  for( int c=0;c<Nc;c++ )
  {
    pixel[c] = img_data[index(i,j,c)];
  }
  return;
}

void ImgProc::set_value( int i, int j, const lux::Color& pixel)
{
  if( img_data == 0 ){ return; }
  if( i<0 || i>=Nx ){ return; }
  if( j<0 || j>=Ny ){ return; }
  //if( Nc > (int)pixel.size() ){ return; }

  #pragma omp parallel for
  for( int c=0;c<Nc;c++ )
  {
    img_data[index(i,j,c)] = pixel[c];
  }

  return;
}


bool ImgProc::load( const std::string& filename )
{
  auto inp = ImageInput::open(filename);
  if (!inp) {return false;}
  const ImageSpec& spec = inp->spec();
  clear();

  Nx = spec.width;
  Ny = spec.height;
  Nc = spec.nchannels;

  Nsize = (long)Nx * (long)Ny * (long)Nc;
  img_data = new float[Nsize];

  unsigned long scanlinesize = Nx * Nc * sizeof(img_data[0]);

  inp->read_image(0, 0, 0, Nc, TypeDesc::FLOAT, (char *)img_data + (Ny-1)*scanlinesize, AutoStride, -scanlinesize, AutoStride);
  inp->close ();
  return true;
}

bool ImgProc::write( const std::string& filename ) const
{
  std::unique_ptr<ImageOutput> out = ImageOutput::create(filename.c_str());
  if( !out ){return false;}
  ImageSpec spec (Nx, Ny, Nc, TypeDesc::FLOAT);
  spec.attribute("user", "imgproc");
  spec.attribute("writer", "imgproc:write" );
  out->open (filename.c_str(), spec);

  unsigned long scanlinesize = Nx * Nc * sizeof(img_data[0]);
  out->write_image (TypeDesc::FLOAT, (char *)img_data + (Ny-1)*scanlinesize, AutoStride, -scanlinesize, AutoStride);
  out->close ();
  return true;
}


bool ImgProc::create( const std::string& filename, int resX, int resY)
{

   Nx = resX;
   Ny = resY;
   Nc = 4;
   Nsize = (long)Nx * (long)Ny * (long)Nc;
   img_data = new float[Nsize];

   return true;
}

