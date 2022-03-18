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
        _point({0,0}),
        _radius_a(0),
        _radius_b(0),
        _rotation(0),
        _focus{0},
        _focusp_l({0,0}),
        _focusp_r({0,0})



    {}

    ellipse(const value_type x, const value_type y, const value_type radius_a, const value_type radius_b, const value_type rotation):
        _point({x,y}),
        _radius_a(radius_a),
        _radius_b(radius_b),
        _rotation(rotation),
        _focus{sqrt((radius_a*radius_a)-(radius_b*radius_b))},
        _focusp_l({(_point[0]-cos(_rotation*M_PI)*_focus),(_point[1]-sin(_rotation*M_PI)*_focus)}),
        _focusp_r({(_point[0]+cos(_rotation*M_PI)*_focus),(_point[1]+sin(_rotation*M_PI)*_focus)})
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

    value_type rotation()const{
        return _rotation;
    }

    value_type& rotation(){
        return _rotation;
    }

    value_type focus()const{
        return _focus;
    }

    value_type& focus(){
        return _focus;
    }

    value_type focusp_l_x()const{
        return _focusp_l[0];
    }

    value_type& focusp_l_x(){
        return _focusp_l[0];
    }

    value_type focusp_l_y()const{
        return _focusp_l[1];
    }

    value_type& focusp_l_y(){
        return _focusp_l[1];
    }

    value_type focusp_r_x()const{
        return _focusp_r[0];
    }

    value_type& focusp_r_x(){
        return _focusp_r[0];
    }

    value_type focusp_r_y()const{
        return _focusp_r[1];
    }

    value_type& focusp_r_y(){
        return _focusp_r[1];
    }



    virtual value_type area()const override{
        return _radius_a*_radius_b*M_PI;
    }

    //bsp function
    virtual void print() const override {
        std::cout<<"Hallo bin ein Ellipse"<<std::endl;
    }

private:
    std::array<value_type, 2> _point;
    value_type _radius_a; // in x direction
    value_type _radius_b; // in y direction
    value_type _rotation;
    value_type _focus;
    std::array<value_type, 2> _focusp_l; //left focus if rotation is 0
    std::array<value_type, 2> _focusp_r; //right focus if rotation is 0

};

}
#endif // ELLIPSE_H
