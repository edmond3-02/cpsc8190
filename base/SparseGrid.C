#include "SparseGrid.h"

namespace lux
{

ScalarGrid::ScalarGrid() :  boost::shared_ptr<SGrid<float> >() {}
ScalarGrid::ScalarGrid( SGrid<float>* f ) :  boost::shared_ptr<SGrid<float> >( f ) {}
ScalarGrid::~ScalarGrid(){}

}
