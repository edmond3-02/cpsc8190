#ifndef __COMPLEX_VOLUMES_H__
#define __COMPLEX_VOLUMES_H__


#include <string>

#include "Fields.h"

using namespace lux;


void makeMcTyson(ScalarField& density, ColorField& color);

void makeObjVolume(std::string filename);

#endif
