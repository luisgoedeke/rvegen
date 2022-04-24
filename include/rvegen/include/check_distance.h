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
void vectorrotation3d (T v[], T rotation, int axis) {
    using value_type = T;

    if (axis == 0){
        value_type y = v[1];
        value_type z = v[2];

        v[1] = y*cos(2*M_PI*rotation)-z*sin(2*M_PI*rotation);
        v[2] = y*sin(2*M_PI*rotation)+z*cos(2*M_PI*rotation);
    }

    if (axis == 1){
        value_type x = v[0];
        value_type z = v[2];

        v[0] = x*cos(2*M_PI*rotation)+z*sin(2*M_PI*rotation);
        v[2] = -x*sin(2*M_PI*rotation)+z*cos(2*M_PI*rotation);
    }

    if (axis == 2){
        value_type x = v[0];
        value_type y = v[1];

        v[0] = x*cos(2*M_PI*rotation)-y*sin(2*M_PI*rotation);
        v[1] = x*sin(2*M_PI*rotation)+y*cos(2*M_PI*rotation);
    }
}

/*template<typename T>
void vectortranslation3d (T v[], T x, T y, T z) {

    v[0] = v[0]+x;
    v[1] = v[1]+y;
    v[2] = v[2]+z;
}
*/

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
bool ellipsoidcontrol (ellipsoid<T> const& lhs, ellipsoid<T> const& rhs){
    using value_type = T;

    //Mittelpunkt
    //Mittelpunkt von lhs im Koordinatensystem von rhs
    value_type center_x = lhs(0)-rhs(0);
    value_type center_y = lhs(1)-rhs(1);
    value_type center_z = lhs(2)-rhs(2);

    //Rotation um die x-Achse mit Winkel von rhs
    value_type rot_x_center_x_rhs = center_x;
    value_type rot_x_center_y_rhs = center_y*cos(rhs.rotation_x()*2*M_PI)-center_z*sin(rhs.rotation_x()*2*M_PI);
    value_type rot_x_center_z_rhs = center_y*sin(rhs.rotation_x()*2*M_PI)+center_z*cos(rhs.rotation_x()*2*M_PI);

    //Rotation um die y-Achse mit Winkel von rhs
    value_type rot_y_center_x_rhs = rot_x_center_x_rhs*cos(rhs.rotation_y()*2*M_PI)+rot_x_center_z_rhs*sin(rhs.rotation_y()*2*M_PI);
    value_type rot_y_center_y_rhs = rot_x_center_y_rhs;
    value_type rot_y_center_z_rhs = rot_x_center_x_rhs*sin(rhs.rotation_y()*2*M_PI)+rot_x_center_z_rhs*cos(rhs.rotation_y()*2*M_PI);

    //Rotation um die z-Achse mit Winkel von rhs
    value_type rot_z_center_x_rhs = rot_y_center_x_rhs*cos(rhs.rotation_z()*2*M_PI)-rot_y_center_y_rhs*sin(rhs.rotation_z()*2*M_PI);
    value_type rot_z_center_y_rhs = rot_y_center_x_rhs*sin(rhs.rotation_z()*2*M_PI)+rot_y_center_y_rhs*cos(rhs.rotation_z()*2*M_PI);
    value_type rot_z_center_z_rhs = rot_y_center_z_rhs;

    value_type x = rot_z_center_x_rhs*rot_z_center_x_rhs;
    value_type a = rhs.radius_a()*rhs.radius_a();
    value_type y = rot_z_center_y_rhs*rot_z_center_y_rhs;
    value_type b = rhs.radius_b()*rhs.radius_b();
    value_type z = rot_z_center_z_rhs*rot_z_center_z_rhs;
    value_type c = rhs.radius_c()*rhs.radius_c();


    if (((x/a)+(y/b)+(z/c)) <= 1){
        std::cout<<"Mittelpunkt"<<std::endl;
        return true;
    }

    //Mittelpunkt von rhs im Koordinatensystem von lhs
    center_x = rhs(0)-lhs(0);
    center_y = rhs(1)-lhs(1);
    center_z = rhs(2)-lhs(2);

    //Rotation um die x-Achse mit Winkel von lhs
    value_type rot_x_center_x_lhs = center_x;
    value_type rot_x_center_y_lhs = center_y*cos(lhs.rotation_x()*2*M_PI)-center_z*sin(lhs.rotation_x()*2*M_PI);
    value_type rot_x_center_z_lhs = center_y*sin(lhs.rotation_x()*2*M_PI)+center_z*cos(lhs.rotation_x()*2*M_PI);

    //Rotation um die y-Achse mit Winkel von lhs
    value_type rot_y_center_x_lhs = rot_x_center_x_lhs*cos(lhs.rotation_y()*2*M_PI)+rot_x_center_z_lhs*sin(lhs.rotation_y()*2*M_PI);
    value_type rot_y_center_y_lhs = rot_x_center_y_lhs;
    value_type rot_y_center_z_lhs = rot_x_center_x_lhs*sin(lhs.rotation_y()*2*M_PI)+rot_x_center_z_lhs*cos(lhs.rotation_y()*2*M_PI);

    //Rotation um die z-Achse mit Winkel von lhs
    value_type rot_z_center_x_lhs = rot_y_center_x_lhs*cos(lhs.rotation_z()*2*M_PI)-rot_y_center_y_lhs*sin(lhs.rotation_z()*2*M_PI);
    value_type rot_z_center_y_lhs = rot_y_center_x_lhs*sin(lhs.rotation_z()*2*M_PI)+rot_y_center_y_lhs*cos(lhs.rotation_z()*2*M_PI);
    value_type rot_z_center_z_lhs = rot_y_center_z_lhs;

    x = rot_z_center_x_lhs*rot_z_center_x_lhs;
    a = lhs.radius_a()*lhs.radius_a();
    y = rot_z_center_y_lhs*rot_z_center_y_lhs;
    b = lhs.radius_b()*lhs.radius_b();
    z = rot_z_center_z_lhs*rot_z_center_z_lhs;
    c = lhs.radius_c()*lhs.radius_c();


    if (((x/a)+(y/b)+(z/c)) <= 1){
        std::cout<<"Mittelpunkt"<<std::endl;
        return true;
    }

    value_type f0[3] = {rhs(0)-lhs(0), rhs(1)-lhs(1), rhs(2)-lhs(2)};
    value_type f1[3] = {rhs.radius_a(), 0, 0};
    value_type f2[3] = {0, rhs.radius_b(), 0};
    value_type f3[3] = {0, 0, rhs.radius_c()};

    //Vektoren rotieren in Koordinatensystem von lhs
    //x-Achse
    vectorrotation3d(f0, rhs.rotation_x()-lhs.rotation_x(), 0);
    vectorrotation3d(f1, rhs.rotation_x()-lhs.rotation_x(), 0);
    vectorrotation3d(f2, rhs.rotation_x()-lhs.rotation_x(), 0);
    vectorrotation3d(f3, rhs.rotation_x()-lhs.rotation_x(), 0);

    //y-Achse
    vectorrotation3d(f0, rhs.rotation_y()-lhs.rotation_y(), 1);
    vectorrotation3d(f1, rhs.rotation_y()-lhs.rotation_y(), 1);
    vectorrotation3d(f2, rhs.rotation_y()-lhs.rotation_y(), 1);
    vectorrotation3d(f3, rhs.rotation_y()-lhs.rotation_y(), 1);

    //z-Achse
    vectorrotation3d(f0, rhs.rotation_z()-lhs.rotation_z(), 2);
    vectorrotation3d(f1, rhs.rotation_z()-lhs.rotation_z(), 2);
    vectorrotation3d(f2, rhs.rotation_z()-lhs.rotation_z(), 2);
    vectorrotation3d(f3, rhs.rotation_z()-lhs.rotation_z(), 2);

    value_type phi;
    value_type theta;

    for (int i=1; i<= 360; i++){
        for (int j=-90; j<= 90; j++) {

            phi = i*2*M_PI/360;
            theta = j*M_PI_4/90;

            x = ((f0[0]+f1[0]*cos(theta)*cos(phi)+f2[0]*cos(theta)*sin(phi)+f3[0]*sin(theta))*(f0[0]+f1[0]*cos(theta)*cos(phi)+f2[0]*cos(theta)*sin(phi)+f3[0]*sin(theta)));
            a = lhs.radius_a()*lhs.radius_a();
            y = ((f0[1]+f1[1]*cos(theta)*cos(phi)+f2[1]*cos(theta)*sin(phi)+f3[1]*sin(theta))*(f0[1]+f1[1]*cos(theta)*cos(phi)+f2[1]*cos(theta)*sin(phi)+f3[1]*sin(theta)));
            b = lhs.radius_b()*lhs.radius_b();
            z = ((f0[2]+f1[2]*cos(theta)*cos(phi)+f2[2]*cos(theta)*sin(phi)+f3[2]*sin(theta))*(f0[2]+f1[2]*cos(theta)*cos(phi)+f2[2]*cos(theta)*sin(phi)+f3[2]*sin(theta)));
            c = lhs.radius_c()*lhs.radius_c();

            if (((x/a)+(y/b)+(z/c)) <= 1){
                std::cout<<"Kollision"<<std::endl;
                return true;
            };

        }
    }
    return false;
}

template<typename T>
bool boxcontrol (ellipsoid<T> const& lhs, ellipsoid<T> const& rhs){
/*    using value_type = T;

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

    return false; */
    return true;

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

    if (boxcontrol(lhs, rhs)) {
        return ellipsoidcontrol(lhs, rhs);
    }
    else {
    return false;
    };
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
