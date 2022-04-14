#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <cstdlib>
#include <math.h>
#include "shape_base.h"

namespace rvegen {

template<typename T = double>
class ellipsoid : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    ellipsoid():
        _point({0,0,0}),
        _radius_a(0),
        _radius_b(0),
        _radius_c(0),
        _rotation_x(0),
        _rotation_y(0),
        _rotation_z(0)



    {}

    ellipsoid(const value_type x, const value_type y, const value_type z, const value_type radius_a, const value_type radius_b, const value_type radius_c, const value_type rotation_x, const value_type rotation_y, const value_type rotation_z):
        _point({x,y,z}),
        _radius_a(radius_a),
        _radius_b(radius_b),
        _radius_c(radius_c),
        _rotation_x(rotation_x),
        _rotation_y(rotation_y),
        _rotation_z(rotation_z)
    {}

    //copy constructor

    value_type operator()(const size_type idx)const{
        return _point[idx];
    }

    value_type& operator()(const size_type idx){
        return _point[idx];
    }

    value_type radius_a()const{
        return _radius_a;
    }

    value_type& radius_a(){
        return _radius_a;
    }


    value_type radius_b()const{
        return _radius_b;
    }

    value_type& radius_b(){
        return _radius_b;
    }

    value_type radius_c()const{
        return _radius_c;
    }

    value_type& radius_c(){
        return _radius_c;
    }

    value_type rotation_x()const{
        return _rotation_x;
    }

    value_type& rotation_x(){
        return _rotation_x;
    }

    value_type rotation_y()const{
        return _rotation_y;
    }

    value_type& rotation_y(){
        return _rotation_y;
    }

    value_type rotation_z()const{
        return _rotation_z;
    }

    value_type& rotation_z(){
        return _rotation_z;
    }

    value_type area()const{
        return _radius_a*_radius_b*M_PI;
    }

    value_type volume()const{
        return (4/3*_radius_a*_radius_b*_radius_c*M_PI);
    }

    //bsp function
    virtual void print() const override {
        std::cout<<"Hallo bin ein Ellipsoid"<<std::endl;
    }

private:
    std::array<value_type, 3> _point;
    value_type _radius_a; // in x direction
    value_type _radius_b; // in y direction
    value_type _radius_c; // in z direction
    value_type _rotation_x;
    value_type _rotation_y;
    value_type _rotation_z;
};
}
#endif // ELLIPSOID_H
