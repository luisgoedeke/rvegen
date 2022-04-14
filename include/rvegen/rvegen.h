#ifndef RVEGEN_H
#define RVEGEN_H


#ifndef NDEBUG
#define RVE_DEBUG
#endif

//shape basis
#include "include/shape_base.h"

//shapes
#include "include/circle.h"
#include "include/cylinder.h"
#include "include/ellipse.h"
#include "include/ellipsoid.h"

//rve generator
#include "include/rve_generator.h"

//shape input
#include "include/rve_shape_input.h"

//gmsh writer
#include "include/write_gmsh_geo.h"


#endif // RVEGEN_H
