#ifndef CYLINDER_H
#define CYLINDER_H

#include <cstdlib>
#include <math.h>
#include "shape_base.h"
#include "box_bounding.h"

namespace rvegen {

template<typename T = double>
class cylinder_lf : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    cylinder_lf():
        _point({0}),
        _radius(0),
        _height(0)
    {}

    cylinder_lf(const value_type x, const value_type y, const value_type z, const value_type radius, const value_type height):
        _point({x, y, z+0.5*height}),
        _radius(radius),
        _height(height)
    {}

    cylinder_lf(cylinder_lf const& __data):
        _point(__data._point),
        _radius(__data._radius),
        _height(__data._height)
    {}

    virtual ~cylinder_lf(){}

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

    virtual value_type area()const{
        return _radius*_radius*M_PI;
    }

    virtual value_type volume() const {
        return _radius*_radius*M_PI*_height;
    }

    virtual std::array<T,3> max_expansion()const{
        return {_radius, _radius, 0.5*_height};
    }

    std::array<T,3> get_middle_point()const override{
        return {_point[0], _point[1], _point[2]};
    }

    void set_middle_point(std::array<T,3> middle_point){
        _point[0] = middle_point[0];
        _point[1] = middle_point[1];
        _point[2] = 0.5*_height;
    }

    constexpr inline auto const& point()const{
        return _point;
    }

    constexpr inline auto& point(){
        return _point;
    }

    virtual void make_bounding_box() override{
        auto box_ptr = std::make_unique<box_bounding<value_type>>();
        box_ptr.get()->top_point() = {_point[0] + _radius,  _point[1] + _radius, _point[2] + 0.5*_height};
        box_ptr.get()->bottom_point() = {_point[0] - _radius,  _point[1] - _radius, _point[2] - 0.5*_height};
        this->_bounding_box = std::move(box_ptr);
    }

private:
    std::array<value_type, 3> _point;
    value_type _radius;
    value_type _height;
};

template<typename T = double>
class cylinder_sf : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    cylinder_sf():
        _point({0}),
        _radius(0),
        _height(0)
    {}

    cylinder_sf(const value_type x, const value_type y, const value_type z, const value_type radius, const value_type height):
        _point({x, y, z}),
        _radius(radius),
        _height(height)
        {}

    cylinder_sf(cylinder_sf const& __data):
        _point(__data._point),
        _radius(__data._radius),
        _height(__data._height)
    {}

    virtual ~cylinder_sf(){}

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

    virtual value_type area()const{
        return _radius*_radius*M_PI;
    }

    virtual value_type volume() const {
        return _radius*_radius*M_PI*_height;
    }

    virtual std::array<T,3> max_expansion()const{
        return {_radius, _radius, 0.5*_height};
    }

    std::array<T,3> get_middle_point()const override{
        return {_point[0], _point[1], _point[2]};
    }

    void set_middle_point(std::array<T,3> middle_point){
        _point[0] = middle_point[0];
        _point[1] = middle_point[1];
        _point[2] = middle_point[2];
    }

    constexpr inline auto const& point()const{
        return _point;
    }

    constexpr inline auto& point(){
        return _point;
    }

    virtual void make_bounding_box() override{
        auto box_ptr = std::make_unique<box_bounding<value_type>>();
        box_ptr.get()->top_point() = {_point[0] + _radius,  _point[1] + _radius, _point[2] + 0.5*_height};
        box_ptr.get()->bottom_point() = {_point[0] - _radius,  _point[1] - _radius, _point[2] - 0.5*_height};
        this->_bounding_box = std::move(box_ptr);
    }

private:
   std::array<value_type, 3> _point;
   value_type _radius;
   value_type _height;
};





}
#endif // CYLINDER_H
