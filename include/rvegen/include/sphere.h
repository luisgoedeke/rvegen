#ifndef SPHERE_H
#define SPHERE_H

#include <cstdlib>
#include <math.h>
#include "shape_base.h"
#include "box_bounding.h"

namespace rvegen {

template <typename T>
class sphere : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    sphere():
        _point({0,0,0}),
        _radius(0)
    {}

    sphere(const value_type x, const value_type y, const value_type z, const value_type radius):
        _point({x,y,z}),
        _radius(radius)
    {}

    sphere(sphere const& __sphere):
        _point(__sphere._point),
        _radius(__sphere._radius)
    {}

    virtual ~sphere(){}

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
        return 0;
    }

    virtual value_type volume()const override{
        return _radius*_radius*_radius*M_PI*4/3;
    }

    virtual std::array<T,3> max_expansion()const{
        return {_radius, _radius, _radius};
    }

    std::array<T,3> get_middle_point()const override{
        return {_point[0], _point[1], _point[2]};
    }

    void set_middle_point(std::array<T,3> middle_point){
        _point[0] = middle_point[0];
        _point[1] = middle_point[1];
        _point[2] = middle_point[2];
    }

    virtual void make_bounding_box() override{
        auto box_ptr = std::make_unique<box_bounding<value_type>>();
        box_ptr.get()->top_point() = {_point[0] + _radius,  _point[1] + _radius, _point[2] + _radius};
        box_ptr.get()->bottom_point() = {_point[0] - _radius,  _point[1] - _radius, _point[2] - _radius};
        this->_bounding_box = std::move(box_ptr);
    }

    constexpr inline auto const& point()const{
        return _point;
    }

    constexpr inline auto& point(){
        return _point;
    }

private:
    std::array<value_type, 3> _point;
    value_type _radius;
};
}
#endif // SPHERE_H
