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
#include "Noise.h"
#include "NoiseMachine.h"
#include "PerlinNoise.h"
#include "ProgressMeter.h"
#include "RayMarcher.h"
#include "SparseGrid.h"
#include "Stamp.h"
#include "Vector.h"
#include "VolumeGeometry.h"

using namespace lux;

Noise_t gen_pnoise(int i);

int main(int argc, char *argv[])
{

	int start_frame = 0;
	int end_frame   = 0;
	bool renderdsm = false;
	std::string filename = "noise_sphere";
	int opt;
	while( (opt = getopt(argc, argv, "ds:e:")) != -1)
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
	float cam_dist = 1;
	RenderData d;
	d.kappa = 10;
	d.resolution[0] = 960;
	d.resolution[1] = 540;
	cam_dist = 2;
	d.snear = 0;
	d.sfar = 2.6;
	d.ds = 0.01;
	
	ColorField color = constant(Color(1,1,1,1));
  	d.colorField = color;

	Noise_t	param;
	param.octaves = 2.0;
	param.frequency = 1.5;
		
	NoiseMachine noise = perlin(param);
					
	d.lightPosition = {Vector(0,-1,1),  Vector(.3,1,.3),  Vector(-2,2,-1)*2};
	d.lightColor    = {Color(1,1,1,1), Color(0,0,.5,1), Color(0,.3,0.0,1)};
	d.lightPosition.resize(1);
	d.lightColor.resize(1);

	int dsmRes = 100;	
	Vector gridDims = (urc - llc);
	AABB bound(llc, urc);
	d.boundingBoxes.push_back(bound);	

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{

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
		char padded_num[6];
		sprintf(padded_num, "%0*d", 5, i);
		std::string framename = prefix + "_" + padded_num;
		

		ScalarField noise;


		AnchorChain particles;
		particles.push_back(gen_pnoise(i));


		int res = 500;
		ScalarGrid grid = ScalarGrid(new SGrid<float>);
		grid->init(res, res, res, box_size*2, box_size*2, box_size*2, llc);

		StampNoise(grid, particles);

		noise = volumize(grid);

		d.densityField = clamp(noise, 0.0, 1.0);
		d.densityField = d.densityField * 40;

		// Render
		d.dsmField.resize(0);
		for (int j=0; j<d.lightPosition.size();j++)
		{
			ScalarGrid dsm = ScalarGrid(new SGrid<float>);
			dsm->init(dsmRes, dsmRes, dsmRes, gridDims.X(), gridDims.Y(), gridDims.Z(), llc);
			dsm->setDefVal(0.0);
		
			d.dsmField.push_back(RayMarchDSMAccumulation(&d, d.densityField, d.lightPosition[j], d.ds, dsm));
		}
	
		std::cout << "dsms" << d.dsmField.size() << std::endl;
		ProgressMeter progress(d.Ny()*d.Nx(), framename);
		
		//RenderFrame(&d, progress, image);
		RenderFrame(&d, progress, image, RayMarchDSM);

		// Write to file
		std::string name = "images/" + framename + ".exr";
		img::write(name, image, d.Nx(), d.Ny(), d.Nc());

	}
}



Noise_t gen_pnoise(int i)
{

	float fade = i / 125;
	fade = remap(fade, 0, 3, 1, 16.0);

	float frequency = (i % 125) / 25;
	frequency = remap(frequency, 0, 4, 3., 5.);

	float fjump = (i % 25) / 5;
	fjump = remap(fjump, 0, 4, 1., 3.);

	float octaves = i % 5;
	octaves = remap(octaves, 0, 4, 2., 4.);

        //fade = float(i)/2.0;
        //frequency = 4;
        //fjump = 2;
        //octaves = 4;

        Noise_t noise;
	noise.octaves = octaves;
	noise.frequency = frequency;
	noise.fjump = fjump;
	noise.falloff = fade;
	noise.radius = 0.3;

	return noise;
}
