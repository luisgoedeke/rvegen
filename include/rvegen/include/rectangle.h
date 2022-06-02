#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "shape_base.h"
#include "rectangle_bounding.h"
#include "../../../external/Eigen/Eigen"

namespace rvegen {

template <typename T>
class rectangle : public shape_base<T>
{
public:
    using value_type = T;

    rectangle():
        _width(),
        _height(),
        _point(),
        _rotation()
    {}

    virtual ~rectangle(){}

    constexpr inline auto operator()(std::size_t i)const{
        return _point[i];
    }

    virtual value_type area() const {
        return _width*_height;
    }

    virtual void make_bounding_box() override {
        using Matrix  = Eigen::Matrix<value_type,2,2>;
        using Vector   = Eigen::Vector2<value_type>;
        using AngleAxis = Eigen::AngleAxis<value_type>;
        Eigen::Rotation2D<value_type> rot(2*_rotation*M_PI);
        const auto w{_width*0.5}, h{_height*0.5};
        std::array<Vector, 4> points{Vector{w,h}, Vector{-w,h}, Vector{w,-h}, Vector{-w,-h}};
        value_type maxx{std::numeric_limits<value_type>::min()}, maxy{std::numeric_limits<value_type>::min()};
        for(auto& vec : points){
            vec = rot*vec;
            if(vec[0]>maxx){
                maxx = vec[0];
            }
            if(vec[1]>maxy){
                maxy = vec[1];
            }
        }

        auto box_ptr = std::make_unique<rectangle_bounding<value_type>>();
        box_ptr.get()->top_point() = {_point[0]+maxx, _point[1]+maxy};
        box_ptr.get()->bottom_point() = {_point[0]-maxx, _point[1]-maxy};
        this->_bounding_box = std::move(box_ptr);
    }

    constexpr inline auto const& height()const{
        return _height;
    }

    constexpr inline auto& height(){
        return _height;
    }

    constexpr inline auto const& width()const{
        return _width;
    }

    constexpr inline auto& width(){
        return _width;
    }

    constexpr inline auto const& point()const{
        return _point;
    }

    constexpr inline auto& point(){
        return _point;
    }

    constexpr inline auto const& rotation()const{
        return _rotation;
    }

    constexpr inline auto& rotation(){
        return _rotation;
    }

    std::array<T,3> get_middle_point()const override{
    return {_point[0], _point[1], _point[2]};
    }

    void set_middle_point(std::array<T,3> middle_point){
        _point[0] = middle_point[0];
        _point[1] = middle_point[1];
    }

    value_type max_expansion()const{
        if (_width > _height){
            return _width;
        }
        else{
            return _height;
        }
    }

private:
    value_type _width;
    value_type _height;
    std::array<value_type, 2> _point;
    value_type _rotation;
};

}

#endif // RECTANGLE_H
