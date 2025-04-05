#include "SparseGrid.h"

namespace lux
{

ScalarGrid::ScalarGrid() :  boost::shared_ptr<SGrid<float> >() {}
ScalarGrid::ScalarGrid( SGrid<float>* f ) :  boost::shared_ptr<SGrid<float> >( f ) {}
ScalarGrid::~ScalarGrid(){}

VectorGrid::VectorGrid() :  boost::shared_ptr<SGrid<Vector> >() {}
VectorGrid::VectorGrid( SGrid<Vector>* f ) :  boost::shared_ptr<SGrid<Vector> >( f ) {}
VectorGrid::~VectorGrid(){}
}
