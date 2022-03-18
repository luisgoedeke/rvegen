#ifndef CHECK_DISTANCE_H
#define CHECK_DISTANCE_H

#include <vector>

#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"


namespace rvegen {

template<typename T>
bool ellipsecontrol (ellipse<T> const& lhs, ellipse<T> const& rhs){
    using value_type = T;

    int Anzahl{100};
    value_type phi;

    for (int i=1; i <= 100; i++){
        phi = 2*M_PI*(i/100);
        value_type r = (rhs.radius_b()*rhs.radius_b()/(rhs.radius_a()-rhs.focus()*cos(phi)));

        //Bestimmung der Koordinaten in Koordiantensystem mit Ursprung in dem linken Brennpunkt von rhs
        value_type pointxrhs = cos(phi)*r;
        value_type pointyrhs = sin(phi)*r;

        //Umrechnung der Koordinaten in allgemeines Koordinatensystem
        value_type pointx = cos(2*M_PI-rhs.rotation()*M_PI)*pointxrhs+sin(2*M_PI-rhs.rotation()*M_PI)*pointyrhs-rhs.focusp_l_x();
        value_type pointy = -sin(2*M_PI-rhs.rotation()*M_PI)*pointxrhs+cos(2*M_PI-rhs.rotation()*M_PI)*pointyrhs-rhs.focusp_l_y();

        //Umrechnung der Koordinaten in Koordinatensystem mit Ursprung im Schnittpunkt der Hauptachsen von lhs
        value_type pointxlhs = cos(lhs.rotation()*M_PI)*pointx-sin(lhs.rotation()*M_PI)*pointy-lhs(0);
        value_type pointylhs = sin(lhs.rotation()*M_PI)*pointx-cos(lhs.rotation()*M_PI)*pointy-lhs(1);

        if((pointxlhs*pointxlhs)/(lhs.radius_a()*lhs.radius_a())+(pointylhs*pointylhs)/(lhs.radius_b()*lhs.radius_b()) <= 1) {
            return true;
        }
    }
    return false;
}

template<typename T>
auto check_distance(circle<T> const& lhs, circle<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const auto dx = lhs(0) - rhs(0);
    const auto dy = lhs(1) - rhs(1);
    return a * a > (dx * dx + dy * dy);
}

template<typename T>
auto check_distance(cylinder<T> const& lhs, cylinder<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const auto dx = lhs(0) - rhs(0);
    const auto dy = lhs(1) - rhs(1);
    const auto start_lhs = lhs(2);
    const auto end_lhs = lhs(2)+lhs.height();
    const auto start_rhs = rhs(2);
    const auto end_rhs = rhs(2)+rhs.height();
    if (a * a > (dx * dx + dy * dy)){
        return !((end_lhs < start_rhs)||start_lhs > end_rhs);
    }
    return false;
}

template<typename T>
auto check_distance(ellipse<T> const& lhs, ellipse<T> const& rhs){
    const auto a{lhs.radius_a() + rhs.radius_a()};
    const auto b{lhs.radius_b() + rhs.radius_b()};
    const double dx = (lhs(0) - rhs(0));
    const double dy = (lhs(1) - rhs(1));
    const double lhsx{lhs(0)};
    const double lhsy{lhs(1)};
    const double rhsx{rhs(0)};
    const double rhsy{rhs(1)};
    const double lhsvx{lhs.radius_a()*cos(lhs.rotation()*M_PI)};
    const double lhsvy{lhs.radius_a()*sin(lhs.rotation()*M_PI)};
    const double rhsvx{rhs.radius_a()*cos(rhs.rotation()*M_PI)};
    const double rhsvy{rhs.radius_a()*sin(rhs.rotation()*M_PI)};
    const double lhshw{lhs.radius_b()/lhs.radius_a()};
    const double rhshw{rhs.radius_b()/rhs.radius_a()};

    //if (a * a > (dx * dx + dy * dy)) {
    //    return true;
    //}
    //else {
        return ellipsecontrol(lhs, rhs);
    //};
}


template<template<class> class _PTR, typename _T>
bool check_distance(std::vector<_PTR<shape_base<_T>>> const& __shapes, circle<_T> const& __circle){
    for(auto& shape : __shapes){
        if(dynamic_cast<circle<_T>*>(shape.get())){
            //check distance
            if(check_distance(__circle, *static_cast<circle<_T>*>(shape.get()))){
                //collision
                return true;
            }
        }
    }
    return false;
}

template<template<class> class _PTR, typename _T>
bool check_distance(std::vector<_PTR<shape_base<_T>>> const& __shapes, cylinder<_T> const& __cylinder){
    for(auto& shape : __shapes){
        if(dynamic_cast<cylinder<_T>*>(shape.get())){
            //check distance
            if(check_distance(__cylinder, *static_cast<cylinder<_T>*>(shape.get()))){
                //collision
                return true;
            }
        }
    }
    return false;
}

template<template<class> class _PTR, typename _T>
bool check_distance(std::vector<_PTR<shape_base<_T>>> const& __shapes, ellipse<_T> const& __ellipse){
    for(auto& shape : __shapes){
        if(dynamic_cast<ellipse<_T>*>(shape.get())){
            //check distance
            if(check_distance(__ellipse, *static_cast<ellipse<_T>*>(shape.get()))){
                //collision
                return true;
            }
        }
    }
    return false;
}


//template<typename T>
//bool check_distance(std::vector<circle<T>> const& vec_circle, circle<T> const& circle){
//    for(size_t i{0}; i<vec_circle.size(); ++i){
//        //check distance
//        if(check_distance(circle, vec_circle[i])){
//            //collision
//            return true;
//        }
//    }
//    return false;
//}


//template<typename T>
//bool check_distance(std::vector<ellipse<T>> const& vec_ellipse, ellipse<T> const& ellipse){
//    for(size_t i{0}; i<vec_ellipse.size(); ++i){
//        //check distance
//        if(check_distance(ellipse, vec_ellipse[i])){
//            //collision
//            return true;
//        }
//    }
//    return false;
//}

}
#endif // CHECK_DISTANCE_H
