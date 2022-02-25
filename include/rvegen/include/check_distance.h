#ifndef CHECK_DISTANCE_H
#define CHECK_DISTANCE_H

#include <vector>

#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"


namespace rvegen {

template<typename T>
auto check_distance(circle<T> const& lhs, circle<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const double dx = lhs(0) - rhs(0);
    const double dy = lhs(1) - rhs(1);
    return a * a > (dx * dx + dy * dy);
}

template<typename T>
auto check_distance_flat(cylinder<T> const& lhs, cylinder<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const double dx = lhs(0) - rhs(0);
    const double dy = lhs(1) - rhs(1);
    const double start_lhs = lhs(2);
    const double end_lhs = lhs(2)+lhs.height();
    const double start_rhs = rhs(2);
    const double end_rhs = rhs(2)+rhs.height();
    if (a * a > (dx * dx + dy * dy)){
        return !((end_lhs < start_rhs)||start_lhs > end_rhs);
    }
    return false;
}

template<typename T>
auto check_distance_space(cylinder<T> const& lhs, cylinder<T> const& rhs){
    const double start_lhs = lhs(2);
    const double end_lhs = lhs(2)+lhs.height();
    const double start_rhs = rhs(2);
    const double end_rhs = rhs(2)+rhs.height();
    return ((end_lhs < start_rhs)||start_lhs > end_rhs);
}


template<typename T>
auto check_distance(ellipse<T> const& lhs, ellipse<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const double dx = abs(lhs(0) - rhs(0));
    const double dy = abs(lhs(1) - rhs(1));
    if((lhs.radius_a() + rhs.radius_a())<dx && (lhs.radius_b() + rhs.radius_b())<dy){
        return true;
    }else{
        return false;
    }
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
            if(check_distance_flat(__cylinder, *static_cast<cylinder<_T>*>(shape.get()))){
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
