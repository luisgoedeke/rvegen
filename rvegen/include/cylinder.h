#ifndef CYLINDER_H
#define CYLINDER_H

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
        _data({0}),
        _radius(0),
        _height(0)
    {}

    cylinder(const value_type x, const value_type y, const value_type radius, const value_type height):
        _data({x,y}),
        _radius(radius),
        _height(height)
    {}

    value_type operator()(const size_type idx)const{
        return _data[idx];
    }

    value_type& operator()(const size_type idx){
        return _data[idx];
    }

    value_type radius()const{
        return _radius;
    }

    value_type& radius(){
        return _radius;
    }

    value_type area()const{
        return _radius*_radius*M_PI;
    }

    //bsp function
    virtual void print() const override {
        std::cout<<"Hallo bin ein Cylinder"<<std::endl;
    }
private:
    std::array<value_type, 2> _data;
    value_type _radius;
    value_type _height;
};

}
#endif // CYLINDER_H
