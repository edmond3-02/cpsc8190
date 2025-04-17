#include <cstdlib>
#include <iostream>
#include <memory>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>

#include "Camera.h"
#include "Color.h"
#include "Fields.h"
#include "Image.h"
#include "ProgressMeter.h"
#include "RayMarcher.h"
#include "SparseGrid.h"
#include "Stamp.h"
#include "Vector.h"
#include "VolumeGeometry.h"

using namespace lux;


int main(int argc, char *argv[])
{
	RenderData d;

	ScalarField sphere = Sphere(Vector(), .15);
	ScalarField cube = CsgBox(Vector(0,0,0), .000000001, 8);
	ColorField red = constant(Color(1,0,0,1));
	ColorField blue = constant(Color(0,0,1,1));

	ColorField color = red * (mask(sphere) - mask(cube));
	color = color + (blue * mask(cube));

  	d.colorField = color;
	d.densityField = clamp(sphere, 0,1);

	std::string filename = "performance_test";

	float box_size = .5;
	Vector urc = Vector(box_size,box_size,box_size);
	Vector llc = Vector(-box_size,-box_size,-box_size);
	float cam_dist = 1;
	d.kappa = 5;
	d.resolution[0] = 1920;
	d.resolution[1] = 1080;
	cam_dist = .5;
	d.snear = 0;
	d.sfar = 1;
	d.ds = 0.001;
	d.Fmax = .1;
	d.Fmin = .0001;
	d. maxSteps = 10000;
	


	AABB bound(llc, urc);
	d.boundingBoxes.push_back(bound);	

	std::string prefix = filename;

	float rot = 0/120.0;
	Vector eye(cam_dist*sin((rot)*2*3.14159265359), 0.0, cam_dist*cos((rot)*2*3.14159265359)); 
	Vector up(0.0, 1.0, 0.0);
	Vector dir = -eye;
	dir[1] = 0;
	d.camera.setEyeViewUp(eye, dir, up);

	// Image allocation
	int Nsize = d.Nx() * d.Ny() * d.Nc();
	float* image = new float[Nsize];

	// Pad number for filename
	std::string framename = prefix;

	d.lightPosition = {Vector(0,-1,1)};
	d.lightColor    = {Color(1,1,1,1)};
	// Render
	d.dsmField.resize(0);
	int dsmRes = 200;	
	for (int j=0; j<d.lightPosition.size();j++)
	{
		ScalarGrid dsm = ScalarGrid(new SGrid<float>);
		dsm->init(dsmRes, dsmRes, dsmRes, 1, 1, 1, llc);
		dsm->setDefVal(0.0);
	
		d.dsmField.push_back(RayMarchDSMAccumulation(&d, d.densityField, d.lightPosition[j], d.ds, dsm));
	}

	ProgressMeter progress(d.Ny(), framename);
	
	RenderFrame(&d, progress, image);
	//RenderFrame(&d, progress, image, RayMarchDSM);

	// Write to file
	std::string name = "images/" + framename + ".exr";
	img::write(name, image, d.Nx(), d.Ny(), d.Nc());

}
