#ifndef RECTANGLE_H
#define RECTANGLE_H

#include <array>
#include "bounding_box_base.h"

namespace rvegen {

template <typename T>
class rectangle_bounding : public bounding_box_base<T>
{
public:
    using value_type = T;

    rectangle_bounding():
        bounding_box_base<T>(),
        _top(),
        _bottom()
    {}

    rectangle_bounding(std::array<value_type, 2> __top, std::array<value_type, 2> __bottom):
        bounding_box_base<T>(),
        _top(__top),
        _bottom(__bottom)
    {}

    rectangle_bounding(value_type __top_x, value_type __top_y, value_type __bottom_x, value_type __bottom_y):
        bounding_box_base<T>(),
        _top({__top_x,__top_y}),
        _bottom({__bottom_x,__bottom_y})
    {}

    virtual ~rectangle_bounding(){}

    constexpr inline auto const& top_point()const{
        return _top;
    }

    constexpr inline auto const& bottom_point()const{
        return _bottom;
    }

    constexpr inline auto& top_point(){
        return _top;
    }

    constexpr inline auto& bottom_point(){
        return _bottom;
    }

private:
    std::array<value_type, 2> _top;
    std::array<value_type, 2> _bottom;
};

}


#endif // RECTANGLE_H
