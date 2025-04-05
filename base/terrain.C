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

	int total_frames = 0;
	int start_frame = 0;
	int end_frame	= 0;
	bool renderdsm	= false;
	std::string filename = "terrain";
	int opt;
	while( (opt = getopt(argc, argv, "s:e:d")) != -1)
	{
		switch(opt)
		{
			case 's':
			  start_frame = std::stoi(optarg);
			  break;

			case 'e':
			  end_frame = std::stoi(optarg);
			  break;
			
			case 'd':
			  renderdsm = true;
			  break;
		}
	}

	float box_size = .5;
	float modelscale = 1.0;
	RenderData render_data;
	render_data.kappa = 10;
	render_data.resolution[0] = 1920;
	render_data.resolution[1] = 1080;
	// test at .004
	render_data.ds = 0.002;
	float cam_dist = 4;
	
	render_data.colorField = constant(Color(1,1,1,1));

	Noise_t terr_noise;
	terr_noise.translate = Vector(-.4,0,-.5);
	terr_noise.octaves = 4;
	terr_noise.frequency = 1;
	terr_noise.fjump = 2;
	terr_noise.terrain_up_gamma = 2;
	terr_noise.terrain_up_amplitude = 15;
	terr_noise.terrain_down_gamma = 1.0/2.0;
	terr_noise.terrain_down_amplitude = .1;
	terr_noise.terrain_gamma_transition = .0;
	NoiseMachine noise = terrain(terr_noise);

	float height = 1.5;
	float width = 1.5;
	render_data.snear = cam_dist - width;
	render_data.sfar = cam_dist + width;
	ScalarField land = Plane(Vector(0,0,0), Vector(0,-1,0));
	land = Pyroclast(land, noise, .4, 0, 1);


	// CALDERA
	Vector caldera_cen = Vector(0,-.7,.73);
	ScalarField caldera = RadialPyroclast(Vector(),.25,.3,  2,3,.2,  0,0,1.5);
	ScalarField hole = Cylinder(Vector(0,1,0), 0.1);
	caldera = Union(caldera, hole);
	ScalarGrid cal_grid = ScalarGrid(new SGrid<float>);
	cal_grid->setDefVal(-1.0);
	float cal_res = 100;
	float cal_width=.4, cal_height=1.5;
	cal_grid->init(cal_res, cal_res, cal_res, cal_width*2, cal_height, cal_width*2, -Vector(cal_width, cal_height/2, cal_width));
	StampField(cal_grid, caldera);
	caldera = gridded(cal_grid);
	caldera = translate(caldera, caldera_cen);

	// LAND BRIDGE
	Noise_t bridge_noise;
	bridge_noise.translate = Vector(-.4,0,-.5);
	bridge_noise.octaves = 2;
	bridge_noise.frequency = 8;
	bridge_noise.fjump = 2;
	bridge_noise.terrain_up_gamma = 1.;
	bridge_noise.terrain_up_amplitude = 1;
	bridge_noise.terrain_down_gamma = .5;
	bridge_noise.terrain_down_amplitude = .2;
	bridge_noise.terrain_gamma_transition = .1;
	NoiseMachine bridge_noise_m = terrain(bridge_noise);

	ScalarField bridge = Torus(Vector(0,0,0), Vector(0,0,1), 0.3, 0.05);
	bridge = Pyroclast(bridge, bridge_noise_m, .002, 4, .5);
	ScalarGrid bridge_grid = ScalarGrid(new SGrid<float>);
	bridge_grid->setDefVal(0.0);
	float bridge_res = 200;
	float bridge_width=.5, bridge_height=.5;
	bridge_grid->init(bridge_res, bridge_res, bridge_res, bridge_width*2, bridge_height+bridge_height/4.0, bridge_width*2, -Vector(bridge_width, bridge_height, bridge_width));
	StampField(bridge_grid, bridge);
	bridge = gridded(bridge_grid)*10;
	bridge = translate(bridge, Vector(.4, 0.1, .8));


	// CAVE 
	Noise_t cave_noise;
	cave_noise.translate = Vector(-.4,0,-.5);
	cave_noise.octaves = 2;
	cave_noise.frequency = 4;
	cave_noise.fjump = 2;
	NoiseMachine cave_noise_m = perlin(cave_noise);

	float cave_len = 1.5;
	ScalarField cave = Ellipse(Vector(0,0,0), Vector(0,1,0), cave_len/2, .15);
	cave = Pyroclast(cave, cave_noise_m, -1., 2, 1);
	ScalarGrid cave_grid = ScalarGrid(new SGrid<float>);
	cave_grid->setDefVal(-0.1);
	float cave_res = 200;
	float cave_width=.25, cave_height=cave_len/2.0;
	cave_grid->init(cave_res, cave_res, cave_res, cave_width*2, cave_height*2, cave_width*2, -Vector(cave_width, cave_height, cave_width));
	StampField(cave_grid, cave);
	cave = gridded(cave_grid)*10;
	cave = rotate(cave, -3.1415/3.0, Vector(1, 0, 1));
	cave = translate(cave, Vector(.8, -0.3, -.4)) * 10;
	
	// AJAX
	ScalarField ajax;
	Mesh ajax_mesh;
	makeObjVolume("models/ajax/smallajax.obj", ajax, ajax_mesh, 250);
	modelscale = .01;
	ajax = scale(ajax, Vector(modelscale, modelscale, modelscale));
	ajax = rotate(ajax, 0.5, Vector(1, 0, 1).unitvector());
	ajax = translate(ajax, Vector(-.8, -0.1, .4));


	// LAND
	ScalarGrid land_grid = ScalarGrid(new SGrid<float>);
	land_grid->setDefVal(0.0);
	float land_res = 400;
	land_grid->init(land_res, land_res, land_res, width*2, height+height/4.0, width*2, -Vector(width, height, width));
	StampField(land_grid, land);
	land = gridded(land_grid);
	ScalarField cutter = Union(cave, caldera);
	land = cutout(land, cutter);
	land = Union(land, bridge);

	// SMOKE
	
	ScalarField smoke = constant(0.0);
	AnchorChain particles;
	float smoke_w=1, smoke_h=1, smoke_d=1;
	int num_p = 50;
	std::uniform_real_distribution<> ran(-1.0, 1.0);
	std::mt19937 gen(0);
	for(int i=0; i<num_p; i++)
	{ 
		float pct = i/float(num_p);
		float angle = pct * .5 * 3.1415;
		float y = sin(angle)*.3;
		float x = (1-cos(angle*2.0))*.4;
		
	        Noise_t noise;
		noise.octaves = 3;
		noise.translate = Vector(0,0,i);
		noise.frequency = 20 - 17*pct;
		noise.fjump = 2;
		noise.roughness = pct;
		noise.amplitude = 1.0 - pct;
		noise.falloff = 1.5*pct + .5;
		float min_r=.05, max_r=.25;
		noise.radius = min_r + (max_r - min_r)*pct;
		Vector offset = Vector(ran(gen), ran(gen), ran(gen)) * .08 * pct;
		noise.P = Vector(x, y, 0) + offset;
		particles.push_back(noise);

		smoke = Union(smoke, clamp(Sphere(noise.P, noise.radius), 0, 1)*((num_p-i)*2));
	}
	int smoke_res = 300;
	ScalarGrid smoke_grid = ScalarGrid(new SGrid<float>);
	smoke_grid->init(smoke_res, smoke_res, smoke_res, smoke_w+.2, smoke_h+.2, smoke_d+.2, Vector(-.1, -.1, -.1));
	StampNoise(smoke_grid, particles);
	smoke = clamp(gridded(smoke_grid), 0, 1) * 30;
	//smoke = Union(smoke, mask(Sphere(Vector(0,0,0), .05))*20);
	//smoke = Union(smoke, mask(Sphere(Vector(1,0,0), .05))*20);
	smoke = scale(smoke, Vector(1,-1,1));
	smoke = rotate(smoke, 3.14 * .75, Vector(0,1,0));
	smoke = translate(smoke, Vector(0, -.45, .75));

	ScalarField all = land;
	all = Union(land, ajax);

	// X cut
	all = intersection(all, Plane(Vector(width,0,0), Vector(1,0,0)));
	all = intersection(all, Plane(Vector(-width,0,0), Vector(-1,0,0)));
	// Z cut
	all = intersection(all, Plane(Vector(0,0,width), Vector(0,0,1)));
	all = intersection(all, Plane(Vector(0,0,-width), Vector(0,0,-1)));	

	all = clamp(all*5, 0., 1.) * 50;
	all = Union(all, smoke);
	//all = smoke;
	render_data.densityField = all;

	render_data.lightPosition = {Vector(-4,-4,4),  Vector(3,-2,3),  Vector(-2,2,-1)*2};
	render_data.lightColor    = {Color(.7,.7,1,1), Color(1.0,.75,.85,1), Color(0,.3,0.0,1)};
	render_data.lightPosition.resize(2);
	render_data.lightColor.resize(2);
	
	// Render
	int dsmRes = 500;	
	if(renderdsm)
	{
		render_data.dsmField.resize(0);
		for (int j=0; j<render_data.lightPosition.size();j++)
		{
			ScalarGrid dsm = ScalarGrid(new SGrid<float>);
			dsm->init(dsmRes, dsmRes, dsmRes, width*2, height+height/4.0, width*2, -Vector(width, height, width));
			dsm->setDefVal(0.0);
		
			render_data.dsmField.push_back(RayMarchDSMAccumulation(&render_data, render_data.densityField, render_data.lightPosition[j], render_data.ds, dsm));
		}
	} 

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{

		// Camera settings
		float rot = i/4.0;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), -2, cam_dist*cos((rot)*2*3.14159265359)); 
//		Vector eye(cam_dist*sin((rot)*2*3.14159265359), -0, cam_dist*cos((rot)*2*3.14159265359)); 
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
		
		ProgressMeter progress(render_data.Ny()*render_data.Nx(), framename);
	
		if(renderdsm)
		  RenderFrame(&render_data, progress, image, RayMarchDSM);
		else
		  RenderFrame(&render_data, progress, image, RayMarchEmission);

		// Write to file
		std::string name = "images/" + framename + ".exr";
		img::write(name, image, render_data.Nx(), render_data.Ny(), render_data.Nc());
	}
}
