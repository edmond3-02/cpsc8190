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
	int start_frame = 0;
	int end_frame   = 119;
	Vector urc=Vector(2,2,2), llc=Vector(-2,-1.2,-2);
	std::string filename = "empty";
	float modelscale = 1.0;
	float cam_dist = 1;
	Mesh mesh;
	int opt;
	RenderData d;
	d.kappa = 8;
	ColorField color = constant(Color(1,1,1,1));
	// obj:,tike,ajax,bunnys,combined,startframe,end
	while( (opt = getopt(argc, argv, "o:tabcs:e:")) != -1)
	{
		switch(opt)
		{
			case 'o':
			break;

			// TIKE 
			case 't':
			{
				filename = "tike";
				int tike_res = 200;
				ScalarField tike;
				
				makeTikeMyson(tike, color);
				ScalarGrid tike_grid = ScalarGrid(new SGrid<float>);
				tike_grid->init(tike_res, tike_res, tike_res, 2, 2, 2, Vector(-1, -1.2, -1));

				StampField( tike_grid, tike, 2);
				tike = volumize(tike_grid);	
				

				// tike settings
				d.ds = .005;
				d.densityField = tike;
				d.colorField = color;
				cam_dist = 4;
				d.snear = cam_dist - 1;
				d.sfar = cam_dist + 1;
				break;
			}
			case 'a':
			{
				std::cout << "starting ajax" << std::endl;
				ScalarField ajax;
				Mesh ajax_mesh;
				makeObjVolume("models/ajax/smallajax.obj", ajax, ajax_mesh, 80);
				urc = ajax_mesh->URC();
				llc = ajax_mesh->LLC();
				filename = "ajax";

				modelscale = .015;
				ajax = scale(ajax, Vector(modelscale, modelscale, modelscale));
				
				//Ajax render settings
				d.ds = .0005;
				d.densityField = ajax;
				d.colorField = color;
				cam_dist = 1.2;
				d.snear = 1;
				d.sfar = 1.4;

				break;
			}
			// BUNNY
			case 'b':
			{
				ScalarField bunny;
				Mesh bunny_mesh;
				makeObjVolume("models/bunny/bunny.obj", bunny, bunny_mesh, 100);
				filename = "bunny";
				//urc = bunny_mesh->URC();
				//llc = bunny_mesh->LLC();

				bunny = translate(bunny, Vector(0,.1,0));

				// bunny settings
				d.ds = .0003;
				d.densityField = bunny;
			  	d.colorField = color;
				cam_dist = 0.5;
				d.snear = 0;
	  			d.sfar = 1;
				break;
			}
			// COMBINED
			case 'c':
			// TO BE DONE
			{
				filename = "combined";
				ScalarField ajax;
				Mesh ajax_mesh;
				makeObjVolume("models/ajax/smallajax.obj", ajax, ajax_mesh, 100);
				urc = Vector(1,1,1)*.5;
				llc = Vector(-1,-1,-1)*.5;
				
				modelscale = .015;
				ajax = scale(ajax, Vector(modelscale, modelscale, modelscale));
				
				ajax = intersection(ajax, Plane(Vector(0,-0.04,0), Vector(.4,-1,-.4)));

				ScalarField bunny;
				Mesh bunny_mesh;
				makeObjVolume("models/bunny/bunny.obj", bunny, bunny_mesh, 100);
				bunny = translate(bunny, Vector(0,.04,0));
				bunny = rotate(bunny, 1.3, Vector(.2,0,1));
				bunny = translate(bunny, Vector(-0.04, 0,0));
				
				ajax = Union(ajax, bunny);
				ajax = mask(ajax);

				ajax = rotate(ajax, 1, Vector(0,0,1));
				float ascale = .15;
				ajax = scale(ajax, Vector(ascale,ascale,ascale));
				ajax = translate(ajax, Vector(-.07,-.07,0));

				int tike_res = 200;
				ScalarField tike;
				makeTikeMyson(tike, color);
				ScalarGrid tike_grid = ScalarGrid(new SGrid<float>);
				tike_grid->init(tike_res, tike_res, tike_res, .5, .6, .4, Vector(-0.25, -0.4, -0.2));

				StampField( tike_grid, tike, 2);
				tike = volumize(tike_grid);

				ColorField acolor = constant(Color(1,1,1,1));
				acolor = multiply(acolor, ajax);
				color = add(color, acolor);

				tike = Union(tike, ajax);
							
				d.ds = .0005;
				d.densityField = tike;
				d.colorField = color;
				cam_dist = 1.2;
				d.snear =.9;
				d.sfar = 1.5;
				

			
				break;
			}
			case 's':
				start_frame = std::stoi(optarg);
			break;

			case 'e':
				end_frame = std::stoi(optarg);
			break;
		}
	}
	
	std::cout << "start frame: " << start_frame << std::endl;
	std::cout << "end frame: " << end_frame << std::endl;
	
	// check that density field was set
	if(!d.densityField) return 0;
	d.densityField = d.densityField * 25;

	d.lightPosition = {Vector(-.2,-2,1),  Vector(2,3,2),  Vector(-1,0,-2)*2};
	d.lightColor    = {Color(1.2,0.2,0.2,1), Color(0,0.5,0,1), Color(0,0,1.0,1)};
	d.lightPosition.resize(3);
	d.lightColor.resize(3);

	int dsmRes = 200;	
	Vector gridDims = (urc - llc);
	for (int i=0; i<d.lightPosition.size();i++)
	{
		ScalarGrid dsm = ScalarGrid(new SGrid<float>);
		dsm->init(dsmRes, dsmRes, dsmRes, gridDims.X(), gridDims.Y(), gridDims.Z(), llc);
		dsm->setDefVal(0.0);
	
		d.dsmField.push_back(RayMarchDSMAccumulation(&d, d.densityField, d.lightPosition[i], d.ds, dsm));
	}

	

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
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

