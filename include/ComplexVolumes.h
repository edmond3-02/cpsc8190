#ifndef __COMPLEX_VOLUMES_H__
#define __COMPLEX_VOLUMES_H__


#include <string>

#include "Fields.h"
#include "VolumeGeometry.h"

using namespace lux;

void makeObjVolume(std::string filename, ScalarField& density, Mesh& mesh, int gridRes=50, int bandwidth=0);

void makeTikeMyson(ScalarField& density, ColorField& color);


#endif
