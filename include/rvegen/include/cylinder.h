#ifndef CYLINDER_H
#define CYLINDER_H

#include <cstdlib>
#include <math.h>
#include "shape_base.h"

namespace rvegen {

template<typename T = double>
class cylinder : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    cylinder():
        _point({0}),
        _radius(0),
        _height(0)
    {}

    cylinder(const value_type x, const value_type y, const value_type z, const value_type radius, const value_type height):
        _point({x,y,z}),
        _radius(radius),
        _height(height)
    {}

    cylinder(cylinder const& __data):
        _point(__data._point),
        _radius(__data._radius),
        _height(__data._height)
    {}

    virtual ~cylinder(){}

    value_type operator()(const size_type idx)const{
        return _point[idx];
    }

    value_type& operator()(const size_type idx){
        return _point[idx];
    }

    value_type radius()const{
        return _radius;
    }

    value_type& radius(){
        return _radius;
    }

    value_type height()const{
        return _height;
    }

    value_type& height(){
        return _height;
    }

    value_type area()const{
        return _radius*_radius*M_PI;
    }

    value_type volume() const {
        return _radius*_radius*M_PI*_height;
    }

    value_type max_expansion()const{
        if (_height/2 > _radius){
            return _height;
        }
        else{
            return _radius;
        }
    }

    value_type& max_expansion(){
        if (_height/2 > _radius){
            return _height;
        }
        else{
            return _radius;
        }
    }

    std::array<T,3> get_middle_point()const override{
        return {_point[0], _point[1], _point[2]};
    }

    void set_middle_point(std::array<T,3> middle_point){
        _point[0] = middle_point[0];
        _point[1] = middle_point[1];
        _point[2] = middle_point[2];
    }

    void move(value_type x, value_type y, value_type z)const{
        _point[0] = x;
        _point[1] = y;
        _point[2] = z;
    }

    virtual void make_bounding_box() override {

    }
private:
    std::array<value_type, 3> _point;
    value_type _radius;
    value_type _height;
};

}
#endif // CYLINDER_H
