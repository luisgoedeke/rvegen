#ifndef WRITE_GMSH_GEO_H
#define WRITE_GMSH_GEO_H

#include <iostream>

#include "circle.h"
#include "ellipse.h"
#include "cylinder.h"

namespace rvegen {

template<typename _T>
class write_gmsh_geo
{
public:
    using value_typ = _T;
    using size_type = std::size_t;

    write_gmsh_geo() {}

    template<typename _RVE>
    inline void write_file(std::fstream & __file, _RVE const& __rve){
        const auto dimension{__rve.dimension()};

        if(dimension == 2){
            write_file_2D(__file, __rve);
        }else if (dimension == 3) {
            write_file_3D(__file, __rve);
        }

    }

private:
    template<typename _RVE>
    inline void write_file_2D(std::fstream & __file, _RVE const& __rve){
        const auto [x_box, y_box, z_box]{__rve.box()};
        const auto& shapes{__rve.shapes()};

        __file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
        __file<<"Rectangle(1) = {0, 0, 0,"<<x_box<<","<<y_box<<", 0};"<<std::endl;

        //Rectangle + 4 lines ???? idk...
        size_type start = 5;

        for(size_t i{0}; i<shapes.size(); ++i){
            if(dynamic_cast<circle<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<circle<value_typ>*>(shapes[i].get())};
                __file<<"Circle("<<start+i<<") = {"<<data(0)<<","<<data(1)<<", 0,"<<data.radius()<<", 0, 2*Pi};"<<std::endl;
            }
        }


        for(size_t i{0};i<shapes.size();++i){
            __file<<"Curve Loop("<<i+2<<") = {"<<5+i<<"};"<<std::endl;
        }
        for(size_t i{0};i<shapes.size();++i){
            __file<<"Plane Surface("<<i+2<<") = {"<<i+2<<"};"<<std::endl;
        }

        for(size_t i{0};i<shapes.size();++i){
            __file<<"BooleanIntersection{ Surface{1}; }{ Surface{"<<i+2<<"}; Delete; }"<<std::endl;
            __file<<"BooleanDifference{ Surface{1}; Delete; }{ Surface{"<<i+2<<"}; }"<<std::endl;
        }
    }

    template<typename _RVE>
    inline void write_file_3D(std::fstream & __file, _RVE const& __rve){
        const auto [x_box, y_box, z_box]{__rve.box()};
        const auto& shapes{__rve.shapes()};

        __file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
        __file<<"Mesh.CharacteristicLengthMin = 0.05;"<<std::endl;
        __file<<"Mesh.CharacteristicLengthMax = 0.1;"<<std::endl;
        __file<<"Box(1) = {0, 0, 0,"<<x_box<<","<<y_box<<", "<<z_box<<"};"<<std::endl;

        //Rectangle + 4 lines ???? idk...
        size_type start = 2;

        for(size_t i{0}; i<shapes.size(); ++i){
            if(dynamic_cast<cylinder<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<cylinder<value_typ>*>(shapes[i].get())};
                __file<<"Cylinder("<<start+i<<") = {"<<data(0)<<","<<data(1)<<","<<data(2)<<", 0, 0, "<<data.height()<<", "<<data.radius()<<", 2*Pi};"<<std::endl;
            }
        }
/*
       for(size_t i{0};i<shapes.size();++i){
            __file<<"Curve Loop("<<i+2<<") = {"<<5+i<<"};"<<std::endl;
        }
        for(size_t i{0};i<shapes.size();++i){
            __file<<"Plane Surface("<<i+2<<") = {"<<i+2<<"};"<<std::endl;
        }
*/
        for(size_t i{0};i<shapes.size();++i){
            __file<<"BooleanIntersection{ Volume{1}; }{ Volume{"<<i+2<<"}; Delete; }"<<std::endl;
        }

        for(size_t i{0};i<shapes.size();++i){
            __file<<"BooleanDifference{ Volume{1}; Delete; }{ Volume{"<<i+2<<"}; }"<<std::endl;
        }


    }
};

}

#endif // WRITE_GMSH_GEO_H
