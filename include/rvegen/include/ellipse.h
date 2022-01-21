#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <cstdlib>
#include <math.h>

#include "shape_base.h"

namespace rvegen {

template<typename T = double>
class ellipse : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    ellipse():
        _data({0,0}),
        _radius_a(0),
        _radius_b(0)
    {}

    ellipse(const value_type x, const value_type y, const value_type radius_a, const value_type radius_b):
        _data({x,y}),
        _radius_a(radius_a),
        _radius_b(radius_b)
    {}

    //copy constructor

    value_type operator()(const size_type idx)const{
        return _data[idx];
    }

    value_type& operator()(const size_type idx){
        return _data[idx];
    }

    value_type radius_b()const{
        return _radius_b;
    }

    value_type& radius_b(){
        return _radius_b;
    }

    value_type radius_a()const{
        return _radius_a;
    }

    value_type& radius_a(){
        return _radius_a;
    }

    virtual value_type area()const override{
        return _radius_a*_radius_b*M_PI;
    }

    //bsp function
    virtual void print() const override {
        std::cout<<"Hallo bin ein Ellipse"<<std::endl;
    }

private:
    std::array<value_type, 2> _data;
    value_type _radius_a; // in x direction
    value_type _radius_b; // in y direction
};

}
#endif // ELLIPSE_H
