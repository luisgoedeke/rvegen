#ifndef RECTANGLE_BOUNDING_H
#define RECTANGLE_BOUNDING_H

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
        _top({__top[0], __top[1], 0}),
        _bottom({__bottom[0], __bottom[1], 0})
    {}

    rectangle_bounding(value_type __top_x, value_type __top_y, value_type __bottom_x, value_type __bottom_y):
        bounding_box_base<T>(),
        _top({__top_x, __top_y, 0}),
        _bottom({__bottom_x, __bottom_y, 0})
    {}

    virtual ~rectangle_bounding(){}

    virtual inline std::array<value_type, 3> const& top_point()const{
        return _top;
    }

    virtual inline std::array<value_type, 3> const& bottom_point()const{
        return _bottom;
    }

    virtual inline std::array<value_type, 3>& top_point(){
        return _top;
    }

    virtual inline std::array<value_type, 3>& bottom_point(){
        return _bottom;
    }

private:
    std::array<value_type, 3> _top;
    std::array<value_type, 3> _bottom;
};

}


#endif // RECTANGLE_BOUNDING_H
