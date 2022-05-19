#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include <cstdlib>
#include <math.h>
#include <memory>

#include "shape_base.h"
#include "box_bounding.h"
#include "../../../external/Eigen/Eigen"

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

    ellipsoid(ellipsoid const& data):
        _point(data._point),
        _radius_a(data._radius_a),
        _radius_b(data._radius_b),
        _radius_c(data._radius_c),
        _rotation_x(data._rotation_x),
        _rotation_y(data._rotation_y),
        _rotation_z(data._rotation_z)
    {}

    virtual ~ellipsoid(){}

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

    void move(value_type x, value_type y, value_type z)const{
        _point[0] = x;
        _point[1] = y;
        _point[2] = z;
    }

    virtual void make_bounding_box() override {
        using Matrix33  = Eigen::Matrix<value_type,3,3>;
        using Vector3   = Eigen::Vector3<value_type>;
        using AngleAxis = Eigen::AngleAxis<value_type>;

        //https://math.stackexchange.com/questions/3926884/smallest-axis-aligned-bounding-box-of-hyper-ellipsoid
        Eigen::DiagonalMatrix<value_type,3> D{1./(_radius_a*_radius_a), 1./(_radius_b*_radius_b), 1./(_radius_c*_radius_c)};

        Matrix33 rot;
        rot = AngleAxis(_rotation_z, Vector3::UnitZ())
                * AngleAxis(_rotation_y, Vector3::UnitY())
                * AngleAxis(_rotation_x, Vector3::UnitX());

        Matrix33 A = (rot*D*rot.transpose()).inverse();
        std::array<value_type, 3> max, min;
        for(std::size_t i{0}; i<3; ++i){
            max[i] = _point[i] +  std::sqrt(A(i,i));
            min[i] = _point[i] -  std::sqrt(A(i,i));
        }

        auto box_ptr = std::make_unique<box_bounding<value_type>>();
        box_ptr.get()->top_point() = max;
        box_ptr.get()->bottom_point() = min;
        this->_bounding_box = std::move(box_ptr);
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
