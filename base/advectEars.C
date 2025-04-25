
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
	bool vis_bounds = false;
	std::string filename = "bunnyears";
	int opt;
	while( (opt = getopt(argc, argv, "s:e:v")) != -1)
	{
		switch(opt)
		{
			case 's':
			  start_frame = std::stoi(optarg);
			  break;

			case 'e':
			  end_frame = std::stoi(optarg);
			  break;
			case 'v':
			  vis_bounds = true;
			  break;
		}
	}

	RenderData render_data;
	render_data.kappa = 20;
	render_data.resolution[0] = 1920;
	//render_data.resolution[0] = 960;
	render_data.resolution[1] = 1080;
	//render_data.resolution[1] = 540;
	//working ish
	render_data.ds = 0.0008;
	//render_data.ds = 0.0005;

	float cam_dist = 0.4;
	float width = 0.1;

	render_data.snear = 0;
	render_data.sfar = 1;


	ScalarField bunny;
	Mesh bunny_mesh;

	/* START bunny load */
	ScalarField density;
	Mesh mesh;
	int gridRes = 200;
	int bandwidth = 10;
	ScalarGrid grid = ScalarGrid(new SGrid<float>);

	grid->setDefVal(0.0);

	// Parse obj
	ObjParser parser;
	mesh = createEmptyMesh();
	parser.load("models/bunny/bunny_fixed.obj", mesh);
	Vector murc = mesh->URC() + Vector(1,1,1)*.05;
	Vector mllc = mesh->LLC() - Vector(1,1,1)*.05;
	murc = mesh->URC();
	mllc = mesh->LLC();
	Vector gridDims =  murc - mllc;
	grid->init(gridRes, gridRes, gridRes, gridDims.X(), gridDims.Y(), gridDims.Z(), mllc );
	
	// Ray march the level set of the mesh
	RayMarchLevelSet(mesh, grid, bandwidth);

	// Make grid into volume
	density = volumize(grid);
	bunny = density;
	/* END bunny load */

	// POSITION BUNNY
	bunny = scale(bunny, Vector(1,-1,1));
	bunny = translate(bunny, Vector(0,.1,0));


	Vector body_llc = Vector(-.09, -.03, -.04);
	Vector body_urc = Vector(.065, .07, .06);
	Vector head_llc = Vector(-.1, -.1, -.06);
	Vector head_urc = Vector(.0, .0, .06);

	if(vis_bounds)
	{
		int res = 50;
		ScalarGrid body_con_grid = ScalarGrid(new SGrid<float>);
		body_con_grid->init(res, res, res, width*1.5, width*1.2, width*1.2, -Vector(width+.07, width-.03, width-.01)/2);
		body_con_grid->setDefVal(0.0);
		// Body bounding box
		Vector body_dim = body_urc - body_llc;
		body_con_grid->init(res, res, res, body_dim.X(), body_dim.Y(), body_dim.Z(), body_llc);
		
		// Head bounding box
		ScalarGrid head_con_grid = ScalarGrid(new SGrid<float>);
		head_con_grid->init(res, res, res, width*1.5, width*1.2, width*1.2, -Vector(width+.07, width-.03, width-.01)/2);
		head_con_grid->setDefVal(0.0);
		Vector head_dim = head_urc - head_llc;
		head_con_grid->init(res, res, res, head_dim.X(), head_dim.Y(), head_dim.Z(), head_llc);
	
		// display bounds	
		ScalarField con = constant(.05);
		StampField(body_con_grid, con);
		StampField(head_con_grid, con);
		con = Union(gridded(body_con_grid), gridded(head_con_grid));
		render_data.densityField = clamp(Union(bunny*40,con), 0, 1);
		render_data.colorField = constant(Color(1,0,0,1))+ constant(Color(0,1,1,0)) * mask(con);
	} else 
	{
	
		// advection noise field
		Noise_t noise_field_parms;
		noise_field_parms.octaves = 3;
		noise_field_parms.frequency = 80;
		noise_field_parms.amplitude = .15;
		noise_field_parms.roughness = .8;
		NoiseMachine noise_field_machine = perlin(noise_field_parms);
		VectorField noise_field = VFNoise(noise_field_machine);
	
		// TESTING ONLY
		VectorGrid noise_grid = VectorGrid(new SGrid<Vector>);
		int noiseRes = 800;
		float noise_width = 0.15;
		Vector ear_urc = Vector(.0, -.03, .06);
		Vector noiseDim = ear_urc - head_llc;
		noise_grid->init(noiseRes, noiseRes, noiseRes, noise_width*2, noise_width*2, noise_width*2, -Vector(noise_width, noise_width, noise_width));
		noise_grid->init(noiseRes, noiseRes, noiseRes, noiseDim.X(), noiseDim.Y(), noiseDim.Z(), head_llc);
		
		noise_field = noise_field - identity();
		StampField(noise_grid, noise_field);
		noise_field = gridded(noise_grid) + identity();
			
		VectorField map = identity();
		map = advect(map, noise_field, .08);
		map = advectBFECC(map, noise_field, .08, 2);
		ScalarField blend_mask = clamp(Plane(Vector(0,-.03,0), Vector(-.5,3,1))* 40, 0, 1);
		//bunny = bunny * clamp(blend_mask - constant(.4), 0, 1);
		bunny = which(bunny, warp(bunny, map), blend_mask);
//		bunny = warp(bunny, map);

	
		render_data.densityField = bunny;
		render_data.densityField = clamp(render_data.densityField*60, 0, 1)*60;
		render_data.colorField = constant(Color(1,1,1,1));
		render_data.colorField = which(constant(Color(0,0,1,1)), constant(Color(1,0,0,1)), blend_mask);
	}

	render_data.lightPosition = {Vector(0, -2, 0)};
	render_data.lightColor    = {Color(1,1,1,1)};
	AABB body_bound(body_llc, body_urc);
	AABB head_bound(head_llc, head_urc);
	render_data.boundingBoxes.push_back(body_bound);
	render_data.boundingBoxes.push_back(head_bound);


	int dsmRes = 300;
	if(render_data.lightPosition.size() > 0)
	{
		render_data.dsmField.resize(0);
		for (int j=0; j<render_data.lightPosition.size();j++)
		{
			ScalarGrid dsm = ScalarGrid(new SGrid<float>);
			Vector dsmDim = body_urc - head_llc;
			dsm->init(dsmRes, dsmRes, dsmRes, dsmDim.X(), dsmDim.Y(), dsmDim.Z(), head_llc);
			dsm->setDefVal(0.0);
		
			render_data.dsmField.push_back(RayMarchDSMAccumulation(&render_data, render_data.densityField, render_data.lightPosition[j], render_data.ds, dsm));
		}
	} 

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{
		// Camera settings
		float rot = i/119.0;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), -0.08, cam_dist*cos((rot)*2*3.14159265359)); 
		Vector up(0.0, 1.0, 0.0);
		Vector dir = -eye;
		render_data.camera.setEyeViewUp(eye, dir, up);


		// Pad number for filename
		char padded_num[6];
		sprintf(padded_num, "%0*d", 5, i);
		std::string framename = prefix + "_" + padded_num;
		std::string name = "images/" + framename + ".exr";
		ProgressMeter progress(render_data.Ny(), framename);
	
		// Image allocation
		int Nsize = render_data.Nx() * render_data.Ny() * render_data.Nc();
		float* image = new float[Nsize];
		// forces clearing of image
		#pragma omp parallel for
		for(int j=0; j<Nsize; j++) image[j] = 0.0;
		
		// Render
		RenderFrame(&render_data, progress, image);

		// Write to file
		img::write(name, image, render_data.Nx(), render_data.Ny(), render_data.Nc());
		delete image;
	}
}
