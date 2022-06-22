#ifndef WRITE_GMSH_GEO_H
#define WRITE_GMSH_GEO_H

#include <iostream>
#include <iomanip>

#include "circle.h"
#include "ellipse.h"
#include "cylinder.h"
#include "sphere.h"
#include "ellipsoid.h"
#include "rectangle.h"
#include "box_bounding.h"

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


    template<typename _RVE>
    inline void write_bounding_boxes(std::fstream & __file, _RVE const& __rve){
        const auto dimension{__rve.dimension()};
        if(dimension == 2){
            write_bounding_boxes_2D(__file, __rve);
        }else if (dimension == 3) {
            write_bounding_boxes_3D(__file, __rve);
        }
    }

private:
    template<typename _RVE>
    inline void write_bounding_boxes_2D(std::fstream & __file, _RVE const& __rve){
        const auto& shapes{__rve.shapes()};
        const size_type start = 5 + shapes.size();
        for(size_t i{0}; i<shapes.size(); ++i){
            const auto& box_ptr{*static_cast<rectangle_bounding<value_typ>*>(shapes[i].get()->bounding_box())};
            const auto& top = box_ptr.top_point();
            const auto& bottom = box_ptr.bottom_point();
            __file<<"Rectangle("<<start+i<<") = {"<<bottom[0]<<", "<<bottom[1]<<", 0, "<<top[0]-bottom[0]<<", "<<top[1]-bottom[1]<<"};"<<std::endl;
        }
    }

    template<typename _RVE>
    inline void write_bounding_boxes_3D(std::fstream & __file, _RVE const& __rve){
        const auto& shapes{__rve.shapes()};
        const size_type start = 2 + shapes.size();
        for(size_t i{0}; i<shapes.size(); ++i){
            const auto& box_ptr{*static_cast<box_bounding<value_typ>*>(shapes[i].get()->bounding_box())};
            const auto& top = box_ptr.top_point();
            const auto& bottom = box_ptr.bottom_point();
            __file<<"Box("<<start+i<<") = {"<<bottom[0]<<", "<<bottom[1]<<", "<<bottom[2]<<", "<<top[0]-bottom[0]<<", "<<top[1]-bottom[1]<<", "<<top[2]-bottom[2]<<"};"<<std::endl;
        }
    }

    template<typename _RVE>
    inline void write_file_2D(std::fstream & __file, _RVE const& __rve){
        const auto [x_box, y_box, z_box]{__rve.box()};
        const auto& shapes{__rve.shapes()};

        __file<<"//Number of shapes: "<<__rve.get_number_of_shapes()<<std::endl;
        __file<<"//Volume fraction: "<<__rve.get_vol_frac_inclusion()<<std::endl;

        __file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
        __file<<"Rectangle(1) = {0, 0, 0, "<<x_box<<", "<<y_box<<", 0};"<<std::endl;

        //Rectangle + 4 lines ???? idk...
        size_type curve_loop{1}, start=0;
        size_type basic_entities{4};//lines, circles, ...?
        size_type surfaces{1};

        for(size_t i{0}; i<shapes.size(); ++i){
            if(dynamic_cast<rectangle<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<rectangle<value_typ>*>(shapes[i].get())};
                const auto& w = data.width();
                const auto& bottom = data.height();
                ++surfaces;
                __file<<"Rectangle("<<surfaces<<") = {"<<-data.width()*0.5<<", "<<-data.height()*0.5<<", 0, "<<data.width()<<", "<<data.height()<<"};"<<std::endl;
                __file<<"Rotate {{0, 0, 1}, {0, 0, 0}, "<<data.rotation()<<"*2*Pi} {Surface{"<<surfaces<<"};}"<<std::endl;
                __file<<"Translate{"<<data(0)<<", "<<data(1)<<", 0} {Surface{"<<surfaces<<"};}"<<std::endl;
                basic_entities += 4;
                ++curve_loop;
            }
        }

        for(size_t i{0}; i<shapes.size(); ++i){
             if(dynamic_cast<rectangle<value_typ>*>(shapes[i].get())){
                 //do nothing
             }else if(dynamic_cast<circle<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<circle<value_typ>*>(shapes[i].get())};
                ++basic_entities;
                ++surfaces;
                __file<<"Circle("<<basic_entities<<") = {"<<data(0)<<","<<data(1)<<", 0, "<<data.radius()<<", 0, 2*Pi};"<<std::endl;
                __file<<"Curve Loop("<<curve_loop+1<<") = {"<<basic_entities<<"};"<<std::endl;
                __file<<"Surface("<<surfaces<<") = {"<<curve_loop+1<<"};"<<std::endl;
                curve_loop+=2;
            }else if(dynamic_cast<ellipse<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<ellipse<value_typ>*>(shapes[i].get())};
                ++basic_entities;
                ++surfaces;
                __file<<"Ellipse("<<basic_entities<<") = {"<<data(0)<<","<<data(1)<<", 0, "<<data.radius_a()<<","<<data.radius_b()<<", 0, 2*Pi};"<<std::endl;
                __file<<"Curve Loop("<<curve_loop+1<<") = {"<<basic_entities<<"};"<<std::endl;
                __file<<"Surface("<<surfaces<<") = {"<<curve_loop+1<<"};"<<std::endl;
                curve_loop+=2;
                __file<<"Rotate {{0, 0, 1}, {"<<data(0)<<", "<<data(1)<<", 0}, "<<data.rotation()<<"*2*Pi} {Surface{"<<surfaces<<"};}"<<std::endl;
            }else{
                throw std::runtime_error("write_gmsh_geo::write_file_2D(): no matching shape type");
            }
        }


//        for(size_t i{0};i<shapes.size();++i){
//            if(!dynamic_cast<rectangle<value_typ>*>(shapes[i].get())){
//            __file<<"Curve Loop("<<curve_loop*4+i+start<<") = {"<<curve_loop*4+start+i<<"};"<<std::endl;
//            }
//        }

//        for(size_t i{0};i<shapes.size();++i){
//            if(!dynamic_cast<rectangle<value_typ>*>(shapes[i].get())){
//            __file<<"Plane Surface("<<curve_loop*4+i+start<<") = {"<<curve_loop*4+start+i<<"};"<<std::endl;
//            }
//        }

        for(size_t i{0};i<shapes.size();++i){
            __file<<"BooleanIntersection{ Surface{1}; }{ Surface{"<<i+2<<"}; Delete; }"<<std::endl;
        }

        for(size_t i{0};i<shapes.size();++i){
            __file<<"BooleanDifference{ Surface{1}; Delete; }{ Surface{"<<i+2<<"}; }"<<std::endl;
        }

        for(size_t i{0};i<shapes.size();++i){
            __file<<"BooleanDifference{ Surface{1}; Delete; }{ Surface{"<<i+2<<"}; }"<<std::endl;
        }
    }

    template<typename _RVE>
    inline void write_file_3D(std::fstream & __file, _RVE const& __rve){
        const auto [x_box, y_box, z_box]{__rve.box()};
        const auto& shapes{__rve.shapes()};

        __file<<"//Number of shapes: "<<__rve.get_number_of_shapes()<<std::endl;
        __file<<"//Volume fraction: "<<__rve.get_vol_frac_inclusion()<<std::endl;

        __file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
        __file<<"Mesh.CharacteristicLengthMin = 0;"<<std::endl;
        __file<<"Mesh.CharacteristicLengthMax = 0.05;"<<std::endl;
        __file<<"Box(1) = {0, 0, 0,"<<x_box<<","<<y_box<<", "<<z_box<<"};"<<std::endl;

        //Rectangle + 4 lines ???? idk...
        const size_type start = 2;

        for(size_t i{0}; i<shapes.size(); ++i){
            if(dynamic_cast<cylinder<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<cylinder<value_typ>*>(shapes[i].get())};
                __file<<"Cylinder("<<start+i<<") = {"<<data(0)<<","<<data(1)<<","<<data(2)<<", 0," <<data.height()<<", "<<"0"<<", "<<data.radius()<<", 2*Pi};"<<std::endl;

            }else if(dynamic_cast<ellipsoid<value_typ>*>(shapes[i].get())){
                const auto& data{*static_cast<ellipsoid<value_typ>*>(shapes[i].get())};
                __file<<"Sphere("<<start+i<<") = {0, 0, 0, 1};"<<std::endl;
                __file<<"Dilate{{0, 0, 0}, {"<<data.radius_a()<<", "<<data.radius_b()<<", "<<data.radius_c()<<"}} {Volume{"<<start+i<<"};}"<<std::endl;
                __file<<"Rotate {{1, 0, 0}, {0, 0, 0}, -"<<data.rotation_x()<<"*2*Pi} {Volume{"<<start+i<<"};}"<<std::endl;
                __file<<"Rotate {{0, 1, 0}, {0, 0, 0}, -"<<data.rotation_y()<<"*2*Pi} {Volume{"<<start+i<<"};}"<<std::endl;
                __file<<"Rotate {{0, 0, 1}, {0, 0, 0}, -"<<data.rotation_z()<<"*2*Pi} {Volume{"<<start+i<<"};}"<<std::endl;
                __file<<"Translate{"<<data(0)<<", "<<data(1)<<", "<<data(2)<<"} {Volume{"<<start+i<<"};}"<<std::endl;
            }else if (dynamic_cast<sphere<value_typ>*>(shapes[i].get())){
                 const auto& data{*static_cast<sphere<value_typ>*>(shapes[i].get())};
                __file<<"Sphere("<<start+i<<") = {"<<data(0)<<", "<<data(1)<<", "<<data(2)<<", "<<data.radius()<<"};"<<std::endl;
            }else{
                throw std::runtime_error("write_gmsh_geo::write_file_3D(): no matching shape type");
            }
        }

        for(size_t i{0};i<shapes.size();++i){
//            __file<<"BooleanIntersection{ Volume{1}; }{ Volume{"<<i+2<<"}; Delete; }"<<std::endl;
        }

        for(size_t i{0};i<shapes.size();++i){
//            __file<<"BooleanDifference{ Volume{1}; Delete; }{ Volume{"<<i+2<<"}; }"<<std::endl;
        }
    }
};

}

#endif // WRITE_GMSH_GEO_H
