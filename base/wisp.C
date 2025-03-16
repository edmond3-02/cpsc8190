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
#include "LinearAlgebra.h"
#include "ProgressMeter.h"
#include "RayMarcher.h"
#include "SparseGrid.h"
#include "Stamp.h"
#include "Vector.h"

using namespace lux;
Noise_t gen_wisp(int i);

int main(int argc, char *argv[])
{

	int total_frames = 500;
	int start_frame = 0;
	int end_frame   = 0;
	bool renderdsm  = false;
	std::string filename = "wisp";
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
	Vector urc = Vector(box_size,box_size,box_size);
	Vector llc = Vector(-box_size,-box_size,-box_size);
	float modelscale = 1.0;
	float cam_dist = 1;
	RenderData render_data;
	render_data.kappa = 10;
	render_data.resolution[0] = 960;
	render_data.resolution[1] = 540;
	cam_dist = 2;
	render_data.snear = 1.4;
	render_data.sfar = 2.6;
	render_data.ds = 0.003;
	
  	render_data.colorField = constant(Color(1,1,1,1));

	render_data.lightPosition = {Vector(0,-1,1),  Vector(.3,1,.3),  Vector(-2,2,-1)*2};
	render_data.lightColor    = {Color(1,1,1,1), Color(0,0,.5,1), Color(0,.3,0.0,1)};
	render_data.lightPosition.resize(2);
	render_data.lightColor.resize(2);

	int dsmRes = 250;	
	Vector gridDims = (urc - llc);
		

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{

		// Camera settings
		float rot = 0/120.0;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), 0.0, cam_dist*cos((rot)*2*3.14159265359)); 
		Vector up(0.0, 1.0, 0.0);
		Vector dir = -eye;
		dir[1] = 0;
		render_data.camera.setEyeViewUp(eye, dir, up);

		// Image allocation
		int Nsize = render_data.Nx() * render_data.Ny() * render_data.Nc();
		float* image = new float[Nsize];

		// Pad number for filename
		char padded_num[6];
		sprintf(padded_num, "%0*d", 5, i);
		std::string framename = prefix + "_" + padded_num;
		

		AnchorChain particles;
                Noise_t p = gen_wisp(i);
		particles.push_back(p);


		int res = 500;
		ScalarGrid grid = ScalarGrid(new SGrid<float>);
		grid->init(res, res, res, box_size*2, box_size*2, box_size*2, llc);

                StampPointWisps(grid, particles);
//                Blur(*grid);
		ScalarField wisp = volumize(grid);


                
		render_data.densityField = clamp(wisp, 0.0, 0.5) * 5.0;
                
		// Render
		if(renderdsm)
		{
			render_data.dsmField.resize(0);
			for (int j=0; j<render_data.lightPosition.size();j++)
			{
				ScalarGrid dsm = ScalarGrid(new SGrid<float>);
				dsm->init(dsmRes, dsmRes, dsmRes, gridDims.X(), gridDims.Y(), gridDims.Z(), llc);
				dsm->setDefVal(0.0);
			
				render_data.dsmField.push_back(RayMarchDSMAccumulation(&render_data, render_data.densityField, render_data.lightPosition[j], render_data.ds, dsm));
			}
		} 
		
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


Noise_t gen_wisp(int i)
{
	
	float octaves = i / 125;
	octaves = remap(octaves, 0, 3, 2., 4.);

	float frequency = (i % 125) / 25;
	frequency = remap(frequency, 0, 4, 1., 4.);

	float fjump = (i % 25) / 5;
	fjump = remap(fjump, 0, 4, 1., 3.);

	float clump = i % 5;
	clump = remap(clump, 4, 0, 0.25, 2.0);


	Noise_t noise;
	noise.octaves = 2 * octaves;
	noise.frequency = frequency;
	noise.fjump = fjump;
        noise.pscale = (5 - (i%5))*0.1;
	noise.wispOctaves = 3;
        noise.wispFreq = 2;
	noise.nbWisps = 5000000;
	noise.wispRadialGroup = clump;

	return noise;
}
