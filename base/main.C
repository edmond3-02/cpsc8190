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
	ColorField color = constant(Color(1,1,1,1));

	// tike
	int tike_res = 200;
	ScalarField tike;
	
	//makeTikeMyson(tike, color);
	ScalarGrid tike_grid = ScalarGrid(new SGrid<float>);
	tike_grid->init(tike_res, tike_res, tike_res, 2, 2, 2, Vector(-1, -1.2, -1));

	//StampField( tike_grid, tike, 2);
	//tike = volumize(tike_grid);
	
	
	ScalarField bunny;
	Mesh bunny_mesh;
	makeObjVolume("models/bunny/bunny.obj", bunny, bunny_mesh, 60);

	// bunny settings
/*	
	RenderData d;
	d.ds = .0006;
	d.densityField = bunny*2;
	d.colorField = color;
	float cam_dist = 0.5;
	d.snear = 0;
	d.sfar = 1;
	d.kappa = 5;
*/	

	// tike settings
/*	RenderData d;
	d.ds = .005;
	d.densityField = tike;
	d.colorField = color;
	float cam_dist = 4;
	d.snear = cam_dist - 1;
	d.sfar = cam_dist + 1;
	d.kappa = 2;
	
*/	


	ScalarField ajax;
	Mesh ajax_mesh;
	makeObjVolume("models/ajax/smallajax.obj", ajax, ajax_mesh, 60);
	float as = .015;
	ajax = scale(ajax, Vector(as, as, as));

	//Ajax render settings
	RenderData d;
	d.ds = .0005;
	d.densityField = ajax;
	d.colorField = color;
	float cam_dist = 1.2;
	d.snear = 1;
	d.sfar = 1.4;
	d.kappa = 4;


	d.densityField = d.densityField * 7;

	d.lightPosition = {Vector(-1,0,-2)*2,  Vector(2,3,2),  Vector(-.2,-2,1)};
	d.lightColor    = {Color(1,0,0,1), Color(0,1,0,1), Color(0.25,0.25,1.25,1)};
	d.lightPosition.resize(3);
	d.lightColor.resize(3);

	Mesh mesh = ajax_mesh;
	int dsmRes = 30;	
	Vector gridDims = (mesh->URC() - mesh->LLC()) * as*1.2;
	for (int i=0; i<d.lightPosition.size();i++)
	{
		ScalarGrid dsm = ScalarGrid(new SGrid<float>);
		dsm->init(dsmRes, dsmRes, dsmRes, gridDims.X(), gridDims.Y(), gridDims.Z(), mesh->LLC()*as );
		dsm->setDefVal(0.0);
	
		d.dsmField.push_back(RayMarchDSMAccumulation(&d, d.densityField, d.lightPosition[i], d.ds, dsm));
	}

	

	std::string prefix = "ajax";
	for(int i=0; i<120; i++)
	{
		float rot = i/120.0;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), 0.0, cam_dist*cos((rot)*2*3.14159265359)); 
		// Bunny cam
		// Vector eye(cam_dist*sin((rot)*2*3.14159265359), -0.1, cam_dist*cos((rot)*2*3.14159265359)); 
		// TYSON CAM
		//Vector eye(cam_dist*sin((rot)*2*3.14159265359), -0.2, cam_dist*cos((rot)*2*3.14159265359)); 
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

