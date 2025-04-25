
#include <cstdlib>
#include <iostream>
#include <memory>
#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <random>

#include "Camera.h"
#include "Color.h"
#include "Fields.h"
#include "Image.h"
#include "LinearAlgebra.h"
#include "ObjParser.h"
#include "ProgressMeter.h"
#include "RayMarcher.h"
#include "SparseGrid.h"
#include "Stamp.h"
#include "Vector.h"
#include "VolumeGeometry.h"
#include "ComplexVolumes.h"

using namespace lux;

int main(int argc, char *argv[])
{

	int start_frame = 0;
	int end_frame	= 0;
	std::string filename = "pyrobunny";
	int opt;
	while( (opt = getopt(argc, argv, "s:e:")) != -1)
	{
		switch(opt)
		{
			case 's':
			  start_frame = std::stoi(optarg);
			  break;

			case 'e':
			  end_frame = std::stoi(optarg);
			  break;
		}
	}

	RenderData render_data;
	render_data.kappa = 10;
	render_data.resolution[0] = 1920;
	//render_data.resolution[0] = 960;
	render_data.resolution[1] = 1080;
	//render_data.resolution[1] = 540;
	//working ish
	render_data.ds = 0.002;
	render_data.ds = 0.001;

	float cam_dist = 0.4;
	float width = 0.1;

	render_data.snear = 0;
	render_data.sfar = 1;


	ScalarField bunny;
	Mesh bunny_mesh;

	/* START bunny load */
	ScalarField density;
	Mesh mesh;
	int gridRes = 100;
	int bandwidth = 12;
	ScalarGrid grid = ScalarGrid(new SGrid<float>);

	grid->setDefVal(0.0);

	// Parse obj
	ObjParser parser;
	mesh = createEmptyMesh();
	parser.load("models/bunny/bunny_fixed.obj", mesh);
	Vector murc = mesh->URC() + Vector(1,1,1)*.05;
	Vector mllc = mesh->LLC() - Vector(1,1,1)*.05;
	Vector gridDims =  murc - mllc;
	grid->init(gridRes, gridRes, gridRes, gridDims.X(), gridDims.Y(), gridDims.Z(), mllc );
	AABB bound(Vector(-width,-width*.8,-width), Vector(width,width*.8,width));
	render_data.boundingBoxes.push_back(bound);
	
	// Ray march the level set of the mesh
	RayMarchLevelSet(mesh, grid, bandwidth);

	// Make grid into volume
	density = volumize(grid);
	bunny = density;
	/* END bunny load */

	// POSITION BUNNY
	bunny = scale(bunny, Vector(1,-1,1));
	bunny = translate(bunny, Vector(0,.1,0));
	bunny = Union(bunny, Sphere(Vector(), .03));
	
	// TEST SPHERE
	//bunny = Sphere(Vector(), .07);

	// pyroclastic displacement noise
	Noise_t perlin_noise;
	perlin_noise.octaves = 2;
	perlin_noise.frequency = 80;
	NoiseMachine noise = perlin(perlin_noise);
	
	// advection noise field
	Noise_t noise_field_parms;
	noise_field_parms.octaves = 3;
	noise_field_parms.frequency = 70;
	noise_field_parms.amplitude = .13;
	NoiseMachine noise_field_machine = perlin(noise_field_parms);
	VectorField noise_field = VFNoise(noise_field_machine);

	// TESTING ONLY
	VectorGrid noise_grid = VectorGrid(new SGrid<Vector>);
	int noiseRes = 800;
	float noise_width = 0.15;
	noise_grid->init(noiseRes, noiseRes, noiseRes, noise_width*2, noise_width*2, noise_width*2, -Vector(noise_width, noise_width, noise_width));
//	StampField(noise_grid, noise_field);
//	noise_field = gridded(noise_grid);

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{
		float pct = std::min(float(i)/59.0, 1.0);

		render_data.densityField = bunny;
		render_data.densityField = Pyroclast(render_data.densityField, noise, .006 * pct, 2, 1);
		//render_data.densityField = Pyroclast(render_data.densityField, noise, .005, 15, 1.0/3.0);

		VectorField map = identity();
		for(int adv_i=60; adv_i<i; adv_i++)
		{
			map = advect(map, noise_field, .004);
		}
		render_data.densityField = warp(render_data.densityField, map);

		// pre render processing
		render_data.densityField = clamp(render_data.densityField*60, 0, 1)*60;
		render_data.colorField = constant(Color(1,1,1,1));
//		render_data.colorField = render_data.colorField + constant(-Color(0,1,1,0)) * mask(bunny);
	
		render_data.lightPosition = {Vector(0, -2, 0)};
		render_data.lightColor    = {Color(1,1,1,1)};
		
		// Render
		int dsmRes = 150;	
//		ScalarGrid constant_grid = ScalarGrid(new SGrid<float>);
//		constant_grid->init(dsmRes, dsmRes, dsmRes, width*2, width*2, width*2, -Vector(width, width, width));
//		constant_grid->setDefVal(0.0);
//		ScalarField con = constant(.1);
//		StampField(constant_grid, con);
//		render_data.densityField = clamp(Union(bunny*40,gridded(constant_grid)), 0, 1);
		
		if(render_data.lightPosition.size() > 0)
		{
			render_data.dsmField.resize(0);
			for (int j=0; j<render_data.lightPosition.size();j++)
			{
				ScalarGrid dsm = ScalarGrid(new SGrid<float>);
				dsm->init(dsmRes, dsmRes, dsmRes, width*2, width*2, width*2, -Vector(width, width, width));
				dsm->setDefVal(0.0);
			
				render_data.dsmField.push_back(RayMarchDSMAccumulation(&render_data, render_data.densityField, render_data.lightPosition[j], render_data.ds, dsm));
			}
		} 

		// Camera settings
		float rot = 0.0;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), 0, cam_dist*cos((rot)*2*3.14159265359)); 
		Vector up(0.0, 1.0, 0.0);
		Vector dir = -eye;
		render_data.camera.setEyeViewUp(eye, dir, up);

		// Image allocation
		int Nsize = render_data.Nx() * render_data.Ny() * render_data.Nc();
		float* image = new float[Nsize];

		// Pad number for filename
		char padded_num[6];
		sprintf(padded_num, "%0*d", 5, i);
		std::string framename = prefix + "_" + padded_num;
		
		ProgressMeter progress(render_data.Ny(), framename);
	
		RenderFrame(&render_data, progress, image);

		// Write to file
		std::string name = "images/" + framename + ".exr";
		img::write(name, image, render_data.Nx(), render_data.Ny(), render_data.Nc());
		delete image;
	}
}
