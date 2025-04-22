
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

	float box_size = .5;
	RenderData render_data;
	render_data.kappa = 20;
	render_data.resolution[0] = 1920;
	render_data.resolution[0] = 960;
	render_data.resolution[1] = 1080;
	render_data.resolution[1] = 540;
	render_data.ds = 0.002;
	float cam_dist = 0.3;
	float height = 0.1;
	float width = 0.5;
	Vector llc = Vector(-width, -width, -width);
	Vector urc = Vector(width, width, width);
	AABB bound(llc, urc);
	render_data.boundingBoxes.push_back(bound);

	render_data.snear = 0;
	render_data.sfar = 1;


	ScalarField bunny;
	Mesh bunny_mesh;
	//makeObjVolume("models/bunny/bunny_fixed.obj", bunny, bunny_mesh, 250, 10);
	//bunny = translate(bunny, Vector(0,.1,0));

	/* START bunny load */
	ScalarField density;
	Mesh mesh;
	int gridRes = 200;
	int bandwidth = 80;
	ScalarGrid grid = ScalarGrid(new SGrid<float>);

	grid->setDefVal(0.0);

	// Parse obj
	ObjParser parser;
	mesh = createEmptyMesh();
	parser.load("models/bunny/bunny_fixed.obj", mesh);
	Vector murc = mesh->URC() + Vector(1,1,1);
	Vector mllc = mesh->LLC() - Vector(1,1,1);
	Vector gridDims =  mesh->URC() - mesh->LLC();
	grid->init(gridRes, gridRes, gridRes, gridDims.X(), gridDims.Y(), gridDims.Z(), mesh->LLC() );
	
	// Ray march the level set of the mesh
	RayMarchLevelSet(mesh, grid, bandwidth);
	grid->setDefVal(-10000.0);

	// Make grid into volume
	density = volumize(grid);
	bunny = density;
	/* END bunny load */

	bunny = scale(bunny, Vector(1,-1,1));
	bunny = translate(bunny, Vector(0,.1,0));
	//bunny = Sphere(Vector(), .07);

		
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
	noise_grid->init(noiseRes, noiseRes, noiseRes, noise_width*2, noise_width*2, noise_width*2, -Vector(width, width, width));
	//StampField(noise_grid, noise_field);
	noise_field = gridded(noise_grid);

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{
		Noise_t perlin_noise;
		perlin_noise.octaves = 2;
		perlin_noise.frequency = 80;
		NoiseMachine noise = perlin(perlin_noise);
		float pct = float(i)/120.0;

		render_data.densityField = bunny;
		render_data.densityField = Pyroclast(render_data.densityField, noise, .002, 9, 1.0);
		//bunny = Pyroclast(bunny, noise, .015, 15, 1.0/3.0);
		//bunny = Pyroclast(bunny, noise, .015, 9, 1);

		VectorField map = identity();
		VectorGrid map_grid = VectorGrid(new SGrid<Vector>);
		int mapRes = 700;
		float map_width = 0.15;
		map_grid->init(mapRes, mapRes, mapRes, map_width*2, map_width*2, map_width*2, -Vector(width, width, width));
		for(int adv_i=60; adv_i<i; adv_i++)
		{
			map = advect(map, noise_field, .004);
		}
		//render_data.densityField = warp(render_data.densityField, map);

		// pre render processing
		render_data.densityField = clamp(render_data.densityField*20, 0, 1)*20;
	//	render_data.densityField = scale(render_data.densityField, Vector(1,-1,1));
		render_data.colorField = constant(Color(1,1,1,1));
	
		render_data.lightPosition = {Vector(0, -2, 0)};
		render_data.lightColor    = {Color(1,1,1,1)};
		
		// Render
		int dsmRes = 100;	
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
		float rot = 0;
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
	}
}
