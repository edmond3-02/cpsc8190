#ifndef IMGPROC_H
#define IMGPROC_H

#include <vector>
#include <string>
#include <OpenImageIO/imageio.h>
#include "Color.h"

//using namespace starter;

namespace img{

class ImgProc
{
 public:
	ImgProc();
	~ImgProc();

	void clear();
//	void clear(int nX, int nY, int nC);

	// image data
	int nx() const { return Nx; }
	int ny() const { return Ny; }
	int depth() const { return Nc; }
	float* raw() {return img_data;};

	// pixel access
	void value( int i, int j, lux::Color& pixel ) const;
	void set_value( int i, int j, const lux::Color& pixel );
	int index(int i, int j, int c) const { return c + Nc*(i + Nx*j); }

	// copies
	ImgProc( const ImgProc& img ); // copy constructor
	ImgProc& operator=(const ImgProc& img ); // copy assignment

	// operator overloads
  	void operator*=( const float s );
  	void operator+=( const float  s );

	// simple algorithms
//	void flip();
//	void flop();

	// file operations
	bool load( const std::string& filename );
	bool write( const std::string& filename ) const;
	bool create( const std::string& filename, int resX, int resY);

 private:
	int Nx, Ny, Nc;
	long Nsize;
	float* img_data;
	std::string filename;
};

}
#endif
