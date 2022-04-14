#ifndef CHECK_DISTANCE_H
#define CHECK_DISTANCE_H

#include <vector>
#include <iostream>
#include <array>

#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"
#include "ellipsoid.h"

namespace rvegen {



template<typename T>
bool ellipsecontrol (ellipse<T> const& lhs, ellipse<T> const& rhs){
    using value_type = T;

    //Mittelpunkt
    //Mittelpunkt von lhs im Koordinatensystem von rhs
    value_type centerx = lhs(0)-rhs(0);
    value_type centery = lhs(1)-rhs(1);

    //Rotation mit Winkel von rhs
    value_type centerxrhs = cos(rhs.rotation()*M_PI)*centerx+sin(rhs.rotation()*M_PI)*centery;
    value_type centeryrhs = -sin(rhs.rotation()*M_PI)*centerx+cos(rhs.rotation()*M_PI)*centery;

    value_type x = centerxrhs*centerxrhs;
    value_type a = rhs.radius_a()*rhs.radius_a();
    value_type y = centeryrhs*centeryrhs;
    value_type b = rhs.radius_b()*rhs.radius_b();


    if ((x/a)+(y/b) <= 1){
        return true;
    }

    //Mittelpunkt
    //Mittelpunkt von rhs im Koordinatensystem von lhs
    centerx = rhs(0)-lhs(0);
    centery = rhs(1)-lhs(1);

    //Rotation mit Winkel von lhs
    centerxrhs = cos(lhs.rotation()*M_PI)*centerx+sin(lhs.rotation()*M_PI)*centery;
    centeryrhs = -sin(lhs.rotation()*M_PI)*centerx+cos(lhs.rotation()*M_PI)*centery;

    x = centerxrhs*centerxrhs;
    a = lhs.radius_a()*lhs.radius_a();
    y = centeryrhs*centeryrhs;
    b = lhs.radius_b()*lhs.radius_b();


    if ((x/a)+(y/b) <= 1){
        return true;
    }

    for (int i=1; i <= 360; i++){
        value_type const phi = M_PI*i*2/360;

        //Bestimmung der Koordinaten im Koordiantensystem von lhs
        value_type pointxrhs = (rhs(0)+rhs.radius_a()*cos(phi)*cos(rhs.rotation()*M_PI)-rhs.radius_b()*sin(phi)*sin(rhs.rotation()*M_PI))-lhs(0);
        value_type pointyrhs = (rhs(1)+rhs.radius_a()*cos(phi)*sin(rhs.rotation()*M_PI)+rhs.radius_b()*sin(phi)*cos(rhs.rotation()*M_PI))-lhs(1);

        //Umrechnung der Koordinaten um Drehwinkel von lhs
        value_type pointxlhs = cos(lhs.rotation()*M_PI)*pointxrhs+sin(lhs.rotation()*M_PI)*pointyrhs;
        value_type pointylhs = -sin(lhs.rotation()*M_PI)*pointxrhs+cos(lhs.rotation()*M_PI)*pointyrhs;

        value_type x = pointxlhs*pointxlhs;
        value_type a = lhs.radius_a()*lhs.radius_a();
        value_type y = pointylhs*pointylhs;
        value_type b = lhs.radius_b()*lhs.radius_b();


        if (((x/a)+(y/b)) <= 1){
            return true;
        }
    }
    return false;
}

template<typename T>
bool rectanglecontrol (ellipse<T> const& lhs, ellipse<T> const& rhs){
    using value_type = T;

    //Mittelpunkt
    //Mittelpunkt von lhs im Koordinatensystem von rhs
    value_type centerx = lhs(0)-rhs(0);
    value_type centery = lhs(1)-rhs(1);

    //Rotation mit Winkel von rhs
    value_type centerxrhs = cos(rhs.rotation()*M_PI)*centerx+sin(rhs.rotation()*M_PI)*centery;
    value_type centeryrhs = -sin(rhs.rotation()*M_PI)*centerx+cos(rhs.rotation()*M_PI)*centery;

    value_type x = centerxrhs*centerxrhs;
    value_type a = rhs.radius_a()*rhs.radius_a();
    value_type y = centeryrhs*centeryrhs;
    value_type b = rhs.radius_b()*rhs.radius_b();


    if ((x/a)+(y/b) <= 1){
        return true;
    }

    //Mittelpunkt
    //Mittelpunkt von rhs im Koordinatensystem von lhs
    centerx = rhs(0)-lhs(0);
    centery = rhs(1)-lhs(1);

    //Rotation mit Winkel von lhs
    centerxrhs = cos(lhs.rotation()*M_PI)*centerx+sin(lhs.rotation()*M_PI)*centery;
    centeryrhs = -sin(lhs.rotation()*M_PI)*centerx+cos(lhs.rotation()*M_PI)*centery;

    x = centerxrhs*centerxrhs;
    a = lhs.radius_a()*lhs.radius_a();
    y = centeryrhs*centeryrhs;
    b = lhs.radius_b()*lhs.radius_b();


    if ((x/a)+(y/b) <= 1){
        return true;
    }

    value_type pointslhs[4][2];
    value_type pointsrhs[4][2];

    //lhs without rotation
    pointslhs[0][0] = +lhs.radius_a();
    pointslhs[0][1] = +lhs.radius_b();
    pointslhs[1][0] = -lhs.radius_a();
    pointslhs[1][1] = +lhs.radius_b();
    pointslhs[2][0] = -lhs.radius_a();
    pointslhs[2][1] = -lhs.radius_b();
    pointslhs[3][0] = +lhs.radius_a();
    pointslhs[3][1] = -lhs.radius_b();

    //rhs without rotation
    pointsrhs[0][0] = +rhs.radius_a();
    pointsrhs[0][1] = +rhs.radius_b();
    pointsrhs[1][0] = -rhs.radius_a();
    pointsrhs[1][1] = +rhs.radius_b();
    pointsrhs[2][0] = -rhs.radius_a();
    pointsrhs[2][1] = -rhs.radius_b();
    pointsrhs[3][0] = +rhs.radius_a();
    pointsrhs[3][1] = -rhs.radius_b();

    //rotation lhs

    for (int i = 0; i <= 3; i++) {

        value_type alpha = lhs.rotation()*M_PI;

        value_type x = pointslhs[i][0];
        value_type y = pointslhs[i][1];

        pointslhs[i][0] = (x*cos(alpha)-y*sin(alpha))+lhs(0);
        pointslhs[i][1] = (x*sin(alpha)+y*cos(alpha))+lhs(1);


    }

    for (int i = 0; i <= 3; i++) {

        value_type alpha = rhs.rotation()*M_PI;

        value_type x = pointsrhs[i][0];
        value_type y = pointsrhs[i][1];

        pointsrhs[i][0] = (x*cos(alpha)-y*sin(alpha))+rhs(0);
        pointsrhs[i][1] = (x*sin(alpha)+y*cos(alpha))+rhs(1);

    }

    for (int i = 0; i <= 3; i++){

        value_type x1 = pointslhs[i][0];
        value_type y1 = pointslhs[i][1];
        value_type x2 = pointslhs[i+1][0];
        value_type y2 = pointslhs[i+1][1];

        if (i ==3){
            x2 = pointslhs[0][0];
            y2 = pointslhs[0][1];
        }

        for (int j = 0; j <= 3; j++ ){


            value_type x3 = pointsrhs[j][0];
            value_type y3 = pointsrhs[j][1];
            value_type x4 = pointsrhs[j+1][0];
            value_type y4 = pointsrhs[j+1][1];

            if (j ==3){
                x4 = pointsrhs[0][0];
                y4 = pointsrhs[0][1];
            }

            //y=mx+b


            value_type m1 = ((y2-y1)/(x2-x1));
            value_type m2 = ((y4-y3)/(x4-x3));

            value_type b1 = y1-m1*x1;
            value_type b2 = y3-m2*x3;

            value_type xschnitt = (b2-b1)/(m1-m2);

            if ((((xschnitt >= x1) && (xschnitt <= x2)) || ((xschnitt >= x2) && (xschnitt <= x1)))
                    && (((xschnitt >= x3) && (xschnitt <= x4)) || ((xschnitt >= x4) && (xschnitt <= x3)))){
                     return true;
                }
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

    if (rectanglecontrol(lhs, rhs)) {
        return ellipsecontrol(lhs, rhs);
    }
    else {
        return false;
    };
}

template<typename T>
auto check_distance(ellipsoid<T> const& lhs, ellipsoid<T> const& rhs){

    return false;
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

template<template<class> class _PTR, typename _T>
bool check_distance(std::vector<_PTR<shape_base<_T>>> const& __shapes, ellipsoid<_T> const& __ellipsoid){
    for(auto& shape : __shapes){
        if(dynamic_cast<ellipsoid<_T>*>(shape.get())){
            //check distance
            if(check_distance(__ellipsoid, *static_cast<ellipsoid<_T>*>(shape.get()))){
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
