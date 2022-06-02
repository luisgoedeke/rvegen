#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <cstdlib>
#include <math.h>
#include <memory>
#include "shape_base.h"
#include "bounding_box_base.h"
#include "rectangle_bounding.h"
#include "../../../external/Eigen/Eigen"
namespace rvegen {

template<typename T = double>
class ellipse : public shape_base<T>
{
public:
    using value_type = T;
    using size_type = std::size_t;

    ellipse():
        shape_base<T>(),
        _point({0,0}),
        _radius_a(0),
        _radius_b(0),
        _rotation(0)
    {}

    ellipse(const value_type x, const value_type y, const value_type radius_a, const value_type radius_b, const value_type rotation):
        shape_base<T>(),
        _point({x,y}),
        _radius_a(radius_a),
        _radius_b(radius_b),
        _rotation(rotation)
    {}

    ellipse(ellipse const& __data):
        shape_base<T>(),
        _point(__data._point),
        _radius_a(__data._radius_a),
        _radius_b(__data._radius_b),
        _rotation(__data._rotation)
    {}

    virtual ~ellipse(){}

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

    constexpr inline auto const& point()const{
        return _point;
    }

    constexpr inline auto& point(){
        return _point;
    }

    std::array<T,3> get_middle_point()const override{
        return {_point[0], _point[1], 0};
    }

    void set_middle_point(std::array<T,3> middle_point){
        _point[0] = middle_point[0];
        _point[1] = middle_point[1];
    }

    virtual void make_bounding_box() override{
        using Matrix22  = Eigen::Matrix<value_type,2,2>;
        using Vector2   = Eigen::Vector2<value_type>;
        using AngleAxis = Eigen::AngleAxis<value_type>;

        Eigen::DiagonalMatrix<value_type,2> D{1./(_radius_a*_radius_a), 1./(_radius_b*_radius_b)};

        Eigen::Rotation2D<value_type> rot(2*_rotation*M_PI);

        Matrix22 A = (rot.toRotationMatrix()*D*rot.toRotationMatrix().transpose()).inverse();
        std::array<value_type, 2> max, min;
        for(std::size_t i{0}; i<2; ++i){
            max[i] = _point[i] +  std::sqrt(A(i,i));
            min[i] = _point[i] -  std::sqrt(A(i,i));
        }

        auto box_ptr = std::make_unique<rectangle_bounding<value_type>>();
        box_ptr.get()->top_point() = max;
        box_ptr.get()->bottom_point() = min;
        this->_bounding_box = std::move(box_ptr);
    }

    virtual value_type area()const override{
        return _radius_a*_radius_b*M_PI;
    }

    value_type max_expansion()const{
        if (_radius_a > _radius_b){
            return _radius_a;
        }
        else{
            return _radius_b;
        }
    }

    value_type& max_expansion(){
        if (_radius_a > _radius_b){
            return _radius_a;
        }
        else{
            return _radius_b;
        }
    }


    void move(value_type x, value_type y)const{
        _point[0] = x;
        _point[1] = y;
    }

private:
    std::array<value_type, 2> _point;
    value_type _radius_a; // in x direction
    value_type _radius_b; // in y direction
    value_type _rotation;
};

}
#endif // ELLIPSE_H
