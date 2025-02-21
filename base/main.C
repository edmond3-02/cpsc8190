#include <cstdlib>
#include <iostream>
#include <memory>
#include <cmath>
#include <string>
#include <iostream>
#include <OpenImageIO/imageio.h>

#include "Color.h"
#include "Fields.h"
#include "Vector.h"
#include "Camera.h"
#include "SparseGrid.h"
#include "ProgressMeter.h"
#include "Stamp.h"
#include "RayMarcher.h"
#include "ObjParser.h"
#include "VolumeGeometry.h"
#include "ComplexVolumes.h"


using namespace lux;


bool write( const std::string& filename, float* img_data, int Nx, int Ny, int Nc );

int main(int argc, char *argv[])
{
	
	ScalarGrid grid = ScalarGrid(new SGrid<float>);
	
	int gridRes = 60;

	grid->setDefVal(0.0);

	// Parse obj
	ObjParser parser;
	Mesh mesh = createEmptyMesh();
	parser.load("models/bunny/bunny.obj", mesh);
	Vector gridDims = mesh->URC() - mesh->LLC();
	grid->init(gridRes, gridRes, gridRes, gridDims.X(), gridDims.Y(), gridDims.Z(), mesh->LLC() );
/*	
	// Ray march the level set of the mesh
	RayMarchLevelSet(mesh, grid, 4);

	// Make grid into volume
	ScalarField density = volumize(grid);
	
	float scale = 1;
	density = scale(density, Vector(scale,-scale,scale));
	density = mask(density);

	//ofstream outfile("volume.txt");
	//WriteVolumeGrid(*grid, outfile);
	//outfile.close();
	*/

	ScalarField density;
	ColorField color = constant(Color(1,1,1,1));
	
	makeMcTyson(density, color);
	
	/* bunny settings
	RenderData d;
	d.ds = .0004;
	d.densityField = density;
	d.colorField = color;
	d.snear = 0;
	d.sfar = 1;
	d.kappa = 10;
	*/


	RenderData d;
	d.ds = .001;
	d.densityField = density;
	d.colorField = color;
	d.snear = 0;
	d.sfar = 1;
	d.kappa = 10;

	d.lightPosition = {Vector(-1,-1,0),  Vector(1,0,0),  Vector(-1,0,0)};
	d.lightColor    = {Color(2,0,0,0), Color(0,2,0,0), Color(0,0,2,0)};
	d.lightPosition.resize(1);
	d.lightColor.resize(1);
	for (int i=0; i<d.lightPosition.size();i++)
	{

		ScalarGrid dsm = ScalarGrid(new SGrid<float>);
		dsm->init(gridRes, gridRes, gridRes, gridDims.X()+0.05, gridDims.Y()+0.05, gridDims.Z()+0.05, mesh->LLC() - Vector(.025, .025, .025) );
		dsm->setDefVal(0.0);
	
		d.dsmField.push_back(RayMarchDSMAccumulation(&d, density, d.lightPosition[i], d.ds, dsm));
	}

	float cam_dist = 0.65;
	

	std::string prefix = "grid";
	for(int i=1; i<2; i++)
	{
		float rot = 0.5;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), -0.1, cam_dist*cos((rot)*2*3.14159265359)); 
		Vector up(0.0, 1.0, 0.0);
		Vector dir = -eye;
		dir[1] = 0;

		d.camera.setEyeViewUp(eye, dir, up);

		int Nsize = d.Nx() * d.Ny() * d.Nc();
		float* image = new float[Nsize];

		char padded_num[6];
		sprintf(padded_num, "%0*d", 5, i);
		std::string framename = prefix + "_" + padded_num;
		
		ProgressMeter progress(d.Ny()*d.Nx(), framename);

		RenderFrame(&d, progress, image, RayMarchDSM);
		//RenderFrame(&d, progress, image, RayMarchEmission);

		write("images/" + framename + ".exr", image, d.Nx(), d.Ny(), d.Nc());
	}
	
}

bool write( const std::string& filename, float* img_data, int Nx, int Ny, int Nc ) 
{
  std::unique_ptr<OpenImageIO_v2_1::ImageOutput> out = OpenImageIO_v2_1::ImageOutput::create(filename.c_str());
  if( !out ){return false;}
  OpenImageIO_v2_1::ImageSpec spec (Nx, Ny, Nc, OpenImageIO_v2_1::TypeDesc::FLOAT);
  spec.attribute("user", "imgproc");
  spec.attribute("writer", "imgproc:write" );
  out->open (filename.c_str(), spec);

  unsigned long scanlinesize = Nx * Nc * sizeof(img_data[0]);
  out->write_image (OpenImageIO_v2_1::TypeDesc::FLOAT, (char *)img_data);
  out->close();
  return true;

}

