#ifndef CIRCLE_H
#define CIRCLE_H

#include <cstdlib>
#include <math.h>
#include "shape_base.h"

namespace rvegen {

template<typename T = double>
class circle : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    circle():
        _point({0,0}),
        _radius(0)
    {}

    circle(const value_type x, const value_type y, const value_type radius):
        _point({x,y}),
        _radius(radius)
    {}

    circle(circle const& __circle):
        _point(__circle._point),
        _radius(__circle._radius)
    {}

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

    virtual value_type area()const override{
        return _radius*_radius*M_PI;
    }

    //bsp function
    virtual void print() const override {
        std::cout<<"Hallo bin ein Kreis"<<std::endl;
    }


private:
    std::array<value_type, 2> _point;
    value_type _radius;
};

}

#endif // CIRCLE_H
