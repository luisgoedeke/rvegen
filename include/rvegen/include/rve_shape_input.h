#ifndef RVE_SHAPE_INPUT_H
#define RVE_SHAPE_INPUT_H

#include <vector>
#include <random>
#include <memory>
#include <map>

#include "shape_base.h"
#include "rectangle.h"
#include "circle.h"
#include "ellipse.h"

namespace rvegen {

template <typename T>
class distribution
{
public:
    using value_type = T;

    distribution() {}

    virtual value_type operator()() = 0;
};

template <typename T>
class uniform_real_distribution : public distribution<T>
{
public:
    using value_type = T;

    uniform_real_distribution(value_type __a, value_type __b):
        _dis(__a, __b)
    {}

    value_type operator()() override {
        std::random_device _rd;
        std::mt19937_64 _gen(_rd());
        return _dis(_gen);
    }

private:
    std::uniform_real_distribution<T> _dis;
};

class rve_shape_input
{
public:
    using value_type = double;

    rve_shape_input():
        _volume_fraction(0),
        _random_position(true),
        _random_geometry(false)
    {}

    rve_shape_input(value_type const __volume_fraction, bool const __random_position, bool const __random_geometry = false):
        _volume_fraction(__volume_fraction),
        _random_position(__random_position),
        _random_geometry(__random_geometry)
    {}

    virtual ~rve_shape_input() {}

    inline auto set_volume_fraction(bool __val){
        _volume_fraction = __val;
    }

    inline auto get_volume_fraction()const{
        return _volume_fraction;
    }

    inline auto set_random_geometry(bool __val){
        _random_geometry = __val;
    }

    inline auto is_random_geometry()const{
        return _random_geometry;
    }

    inline auto set_random_position(bool __val){
        _random_position = __val;
    }

    inline auto is_random_position()const{
        return _random_position;
    }

    virtual value_type min_area() const = 0;

    virtual std::unique_ptr<shape_base<value_type>> new_shape() const = 0;

    virtual inline void setup_position(value_type __x_max, value_type __y_max, value_type __z_max = 0){
        _distributions.insert({"pos_x", std::make_unique<rvegen::uniform_real_distribution<value_type>>(0.0, __x_max)});
        _distributions.insert({"pos_y", std::make_unique<rvegen::uniform_real_distribution<value_type>>(0.0, __y_max)});
        if(__z_max != 0){
            _distributions.insert({"pos_z", std::make_unique<rvegen::uniform_real_distribution<value_type>>(0.0, __z_max)});
        }
    }

    inline auto insert_distribution(std::string __key, std::unique_ptr<distribution<value_type>> && __distributions){
        _distributions.insert({__key, std::move(__distributions)});
    }

protected:
    std::map<std::string, std::unique_ptr<distribution<value_type>>> _distributions;
private:
    value_type _volume_fraction;
    bool _random_position;
    bool _random_geometry;
};

class rectangle_input : public rve_shape_input
{
public:
    using value_type = double;

    rectangle_input(){}

    rectangle_input(bool const __random_position, bool const __random_geometry,  value_type const __volume_fraction):
        rve_shape_input(__volume_fraction, __random_position),
        _random_geometry(__random_geometry)
    {}

    virtual ~rectangle_input(){}

    virtual std::unique_ptr<shape_base<value_type>> new_shape() const override {
        auto& pos_x = *this->_distributions.at("pos_x").get();
        auto& pos_y = *this->_distributions.at("pos_y").get();
        auto& width = *this->_distributions.at("width").get();
        auto& height = *this->_distributions.at("height").get();
        auto& rot = *this->_distributions.at("rotation").get();

        auto shape = std::make_unique<rectangle<value_type>>();
        shape.get()->width() = width();
        shape.get()->height() = height();
        shape.get()->point() = {pos_x(), pos_y()};
        shape.get()->rotation() = rot();
        return shape;
    }

    virtual value_type min_area() const override{
        return 0.1;//this->_distributions.at("width")->a()*this->_distributions.at("height")->a()
    }

private:
    bool _random_geometry;
};



class circle_input : public rve_shape_input
{
public:
    using value_type = double;

    circle_input():
        _random_radius(),
        _min_radius(),
        _max_radius()
    {}

    circle_input(bool const __random_position, bool const __random_radius, value_type const _min_radius, value_type const _max_radius, value_type const __volume_fraction):
        rve_shape_input(__volume_fraction, __random_position),
        _random_radius(__random_radius),
        _min_radius(_min_radius),
        _max_radius(_max_radius)
    {}

    virtual ~circle_input(){}

    inline auto set_random_radius(bool __val){
        _random_radius = __val;
    }

    inline auto is_random_radius()const{
        return _random_radius;
    }

    inline auto get_radius_min()const{
        return _min_radius;
    }

    inline auto get_radius_max()const{
        return _max_radius;
    }

    inline auto set_radius_min(double const __val){
        _min_radius = __val;
    }

    inline auto set_radius_max(double const __val){
        _max_radius = __val;
    }

    virtual value_type min_area() const override{
        return _min_radius*_min_radius*M_PI;
    }

    virtual std::unique_ptr<shape_base<value_type>> new_shape() const override {
        auto& pos_x = *this->_distributions.at("pos_x").get();
        auto& pos_y = *this->_distributions.at("pos_y").get();
        auto& radius = *this->_distributions.at("radius").get();

        auto shape = std::make_unique<circle<value_type>>();
        shape.get()->radius() = radius();
        shape.get()->point() = {pos_x(), pos_y()};
        return shape;
    }
private:
    bool _random_radius;
    value_type _min_radius;
    value_type _max_radius;
};

class ellipse_input : public rve_shape_input
{
public:
    using value_type = double;

    ellipse_input():
        _random_radius(),
        _min_radius_a(),
        _max_radius_a(),
        _min_radius_b(),
        _max_radius_b(),
        _min_rotation(),
        _max_rotation()
    {}

    ellipse_input(bool const __random_position, bool const __random_radius, value_type const _min_radius_a, value_type const _max_radius_a, value_type const _min_radius_b, value_type const _max_radius_b, value_type const _min_rotation, value_type const _max_rotation, value_type const __volume_fraction):
        rve_shape_input(__volume_fraction, __random_position),
        _random_radius(__random_radius),
        _min_radius_a(_min_radius_a),
        _max_radius_a(_max_radius_a),
        _min_radius_b(_min_radius_b),
        _max_radius_b(_max_radius_b),
        _min_rotation(_min_rotation),
        _max_rotation(_max_rotation)
    {}

    virtual ~ellipse_input(){}

    inline auto set_random_radius(bool __val){
        _random_radius = __val;
    }

    inline auto is_random_radius()const{
        return _random_radius;
    }

    inline auto get_radius_min_a()const{
        return _min_radius_a;
    }

    inline auto get_radius_max_a()const{
        return _max_radius_a;
    }

    inline auto get_radius_min_b()const{
        return _min_radius_b;
    }

    inline auto get_radius_max_b()const{
        return _max_radius_b;
    }

    inline auto get_min_rotation()const{
        return _min_rotation;
    }

    inline auto get_max_rotation()const{
        return _max_rotation;
    }

    inline auto set_radius_min_a(double const __val){
        _min_radius_a = __val;
    }

    inline auto set_radius_max_a(double const __val){
        _max_radius_a = __val;
    }

    inline auto set_radius_min_b(double const __val){
        _min_radius_b = __val;
    }

    inline auto set_radius_max_b(double const __val){
        _max_radius_b = __val;
    }

    inline auto set_min_rotation(double const __val){
        _min_rotation = __val;
    }

    inline auto set_max_rotation(double const __val){
        _max_rotation = __val;
    }

    virtual value_type min_area() const override {
        return 0.1;
    }

    virtual std::unique_ptr<shape_base<value_type>> new_shape() const override {
        auto& pos_x = *this->_distributions.at("pos_x").get();
        auto& pos_y = *this->_distributions.at("pos_y").get();
        auto& radius_a = *this->_distributions.at("radius_a").get();
        auto& radius_b = *this->_distributions.at("radius_b").get();
        auto& rotation = *this->_distributions.at("rotation").get();

        auto shape = std::make_unique<ellipse<value_type>>();
        shape.get()->radius_a() = radius_a();
        shape.get()->radius_b() = radius_b();
        shape.get()->rotation() = rotation();
        shape.get()->point() = {pos_x(), pos_y()};
        return shape;
    }
private:
    bool _random_radius;
    value_type _min_radius_a;
    value_type _max_radius_a;
    value_type _min_radius_b;
    value_type _max_radius_b;
    value_type _min_rotation;
    value_type _max_rotation;
};

class cylinder_input : public rve_shape_input
{
public:
    using value_type = double;

    cylinder_input():
        _random_radius(),
        _random_height(),
        _min_radius(),
        _max_radius(),
        _min_height(),
        _max_height()
    {}

    cylinder_input(bool const __random_position, bool const __random_radius, bool const __random_height, value_type const _min_radius, value_type const _max_radius, value_type const _min_height, value_type const _max_height,value_type const __volume_fraction):
        rve_shape_input(__volume_fraction, __random_position),
        _random_radius(__random_radius),
        _min_radius(_min_radius),
        _max_radius(_max_radius),
        _random_height(__random_height),
        _min_height(_min_height),
        _max_height(_max_height)
    {}

    virtual ~cylinder_input(){}

    inline auto set_random_radius(bool __val){
        _random_radius = __val;
    }

    inline auto is_random_radius()const{
        return _random_radius;
    }

    inline auto get_radius_min()const{
        return _min_radius;
    }

    inline auto get_radius_max()const{
        return _max_radius;
    }

    inline auto set_radius_min(double const __val){
        _min_radius = __val;
    }

    inline auto set_radius_max(double const __val){
        _max_radius = __val;
    }

    inline auto set_random_height(bool __val){
        _random_height = __val;
    }

    inline auto is_random_height()const{
        return _random_height;
    }

    inline auto get_height_min()const{
        return _min_height;
    }

    inline auto get_height_max()const{
        return _max_height;
    }

    inline auto set_height_min(double const __val){
        _min_height = __val;
    }

    inline auto set_height_max(double const __val){
        _max_height = __val;
    }

    virtual value_type min_area() const {
        return 0.1;
    }

    virtual std::unique_ptr<shape_base<value_type>> new_shape() const override {

        auto shape = std::make_unique<rectangle<value_type>>();

        return shape;
    }
private:
    bool _random_radius;
    bool _random_height;
    value_type _min_radius;
    value_type _max_radius;
    value_type _min_height;
    value_type _max_height;
};

class ellipsoid_input : public rve_shape_input
{
public:
    using value_type = double;

    ellipsoid_input():
        _random_radius(),
        _min_radius_a(),
        _max_radius_a(),
        _min_radius_b(),
        _max_radius_b(),
        _min_radius_c(),
        _max_radius_c(),
        _min_rotation_x(),
        _max_rotation_x(),
        _min_rotation_y(),
        _max_rotation_y(),
        _min_rotation_z(),
        _max_rotation_z()

    {}

    ellipsoid_input(bool const __random_position, bool const __random_radius, value_type const _min_radius_a, value_type const _max_radius_a, value_type const _min_radius_b, value_type const _max_radius_b, value_type const _min_radius_c, value_type const _max_radius_c, value_type const _min_rotation_x, value_type const _max_rotation_x, value_type const _min_rotation_y, value_type const _max_rotation_y, value_type const _min_rotation_z, value_type const _max_rotation_z,  value_type const __volume_fraction):
        rve_shape_input(__volume_fraction, __random_position),
        _random_radius(__random_radius),
        _min_radius_a(_min_radius_a),
        _max_radius_a(_max_radius_a),
        _min_radius_b(_min_radius_b),
        _max_radius_b(_max_radius_b),
        _min_radius_c(_min_radius_c),
        _max_radius_c(_max_radius_c),
        _min_rotation_x(_min_rotation_x),
        _max_rotation_x(_max_rotation_x),
        _min_rotation_y(_min_rotation_y),
        _max_rotation_y(_max_rotation_y),
        _min_rotation_z(_min_rotation_z),
        _max_rotation_z(_max_rotation_z)
    {}

    virtual ~ellipsoid_input(){}

    inline auto set_random_radius(bool __val){
        _random_radius = __val;
    }

    inline auto is_random_radius()const{
        return _random_radius;
    }

    inline auto get_radius_min_a()const{
        return _min_radius_a;
    }

    inline auto get_radius_max_a()const{
        return _max_radius_a;
    }

    inline auto get_radius_min_b()const{
        return _min_radius_b;
    }

    inline auto get_radius_max_b()const{
        return _max_radius_b;
    }

    inline auto get_radius_min_c()const{
        return _min_radius_c;
    }

    inline auto get_radius_max_c()const{
        return _max_radius_c;
    }

    inline auto get_min_rotation_x()const{
        return _min_rotation_x;
    }

    inline auto get_max_rotation_x()const{
        return _max_rotation_x;
    }

    inline auto get_min_rotation_y()const{
        return _min_rotation_y;
    }

    inline auto get_max_rotation_y()const{
        return _max_rotation_y;
    }

    inline auto get_min_rotation_z()const{
        return _min_rotation_z;
    }

    inline auto get_max_rotation_z()const{
        return _max_rotation_z;
    }

    inline auto set_radius_min_a(double const __val){
        _min_radius_a = __val;
    }

    inline auto set_radius_max_a(double const __val){
        _max_radius_a = __val;
    }

    inline auto set_radius_min_b(double const __val){
        _min_radius_b = __val;
    }

    inline auto set_radius_max_b(double const __val){
        _max_radius_b = __val;
    }

    inline auto set_radius_min_c(double const __val){
        _min_radius_c = __val;
    }

    inline auto set_radius_max_c(double const __val){
        _max_radius_c = __val;
    }

    inline auto set_min_rotation_x(double const __val){
        _min_rotation_x = __val;
    }

    inline auto set_max_rotation_x(double const __val){
        _max_rotation_x = __val;
    }

    inline auto set_min_rotation_y(double const __val){
        _min_rotation_y = __val;
    }

    inline auto set_max_rotation_y(double const __val){
        _max_rotation_y = __val;
    }

    inline auto set_min_rotation_z(double const __val){
        _min_rotation_z = __val;
    }

    inline auto set_max_rotation_z(double const __val){
        _max_rotation_z = __val;
    }

    virtual value_type min_area() const {
        return 0.1;
    }

    virtual std::unique_ptr<shape_base<value_type>> new_shape() const override {

        auto shape = std::make_unique<rectangle<value_type>>();

        return shape;
    }
private:
    bool _random_radius;
    value_type _min_radius_a;
    value_type _max_radius_a;
    value_type _min_radius_b;
    value_type _max_radius_b;
    value_type _min_radius_c;
    value_type _max_radius_c;
    value_type _min_rotation_x;
    value_type _max_rotation_x;
    value_type _min_rotation_y;
    value_type _max_rotation_y;
    value_type _min_rotation_z;
    value_type _max_rotation_z;

};
}
#endif // RVE_SHAPE_INPUT_H
