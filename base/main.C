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

	// I wish this wasn't needed. it is though.
  	std::unique_ptr<OpenImageIO_v2_1::ImageOutput> out = OpenImageIO_v2_1::ImageOutput::create("blah");
	
	int start_frame = 0;
	int end_frame   = 0;
	std::string filename = "noise_sphere";
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

	Vector urc = Vector(.6,.6,.6);
	Vector llc = Vector(-.6,-.6,-.6);
	float modelscale = 1.0;
	float cam_dist = 1;
	Mesh mesh;
	RenderData d;
	d.kappa = 10;
	d.resolution[0] = 960;
	d.resolution[1] = 540;
	
	ColorField color = constant(Color(1,1,1,1));

	Noise_t	param;
	param.octaves = 2.0;
	param.frequency = 1.5;
		
	NoiseMachine noise = perlin(param);
					
	ScalarField sphere = clamp(Sphere(Vector(0,0,0), .6),0,1);
	sphere = pow(sphere, 2);

	d.densityField = sphere * scale(SFNoise(noise), Vector(.2,.2,.2));
  	d.colorField = color;
	cam_dist = 2;
	d.snear = 1.4;
			d.sfar = 2.6;
	d.ds = 0.002;


	
	std::cout << "start frame: " << start_frame << std::endl;
	std::cout << "end frame: " << end_frame << std::endl;
	
	// check that density field was set
	d.densityField = d.densityField * 20;

	d.lightPosition = {Vector(0,1,1),  Vector(.3,-1,.3),  Vector(-2,2,-1)*2};
	d.lightColor    = {Color(1,1,1,1), Color(0,0,.5,1), Color(0,.3,0.0,1)};
	d.lightPosition.resize(2);
	d.lightColor.resize(2);

	int dsmRes = 100;	
	Vector gridDims = (urc - llc);
		

	std::string prefix = filename;
	for(int i=start_frame; i<=end_frame; i++)
	{
		Noise_t anim = gen_pnoise(i);

		//anim.time = float(i)/10;

		noise->setParameters(anim);
		//std::cout << "octaves = " << anim.octaves << std::endl;
		/*d.dsmField.resize(0);
		for (int i=0; i<d.lightPosition.size();i++)
		{
			ScalarGrid dsm = ScalarGrid(new SGrid<float>);
			dsm->init(dsmRes, dsmRes, dsmRes, gridDims.X(), gridDims.Y(), gridDims.Z(), llc);
			dsm->setDefVal(0.0);
		
			d.dsmField.push_back(RayMarchDSMAccumulation(&d, d.densityField, d.lightPosition[i], d.ds, dsm));
		}*/

		float rot = 0/120.0;
		Vector eye(cam_dist*sin((rot)*2*3.14159265359), 0.0, cam_dist*cos((rot)*2*3.14159265359)); 
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

		//RenderFrame(&d, progress, image, RayMarchDSM);
		RenderFrame(&d, progress, image, RayMarchEmission);

		img::write("images/" + framename + ".exr", image, d.Nx(), d.Ny(), d.Nc());
	}
}



Noise_t gen_pnoise(int i)
{
	Noise_t noise;
	
	// range size
	int rs = 96;
	noise.octaves = (i < rs*.5) ? remap(i, 0, rs*.5, 1, 4) : remap(i, rs, rs*.5, 3, 4);
	noise.frequency = (i < rs *1.5) ? remap(i, rs, rs*1.5, 1, 4) : remap(i, rs*2, rs*1.5, 2, 4);
	noise.fjump = (i < rs*2.5) ? remap(i, rs*2, rs*2.5, 2, 3) : remap(i, rs*3, rs*2.5, 1, 3);

	return noise;
}
