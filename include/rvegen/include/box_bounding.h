#ifndef BOX_BOUNDING_H
#define BOX_BOUNDING_H

#include <array>
#include "bounding_box_base.h"

namespace rvegen {

template <typename T>
class box_bounding : public bounding_box_base<T>
{
public:
    using value_type = T;

    box_bounding():
        bounding_box_base<T>(),
        _top(),
        _bottom()
    {}

    box_bounding(std::array<value_type, 3> __top, std::array<value_type, 3> __bottom):
        bounding_box_base<T>(),
        _top(__top),
        _bottom(__bottom)
    {}

    box_bounding(value_type __top_x, value_type __top_y, value_type __top_z, value_type __bottom_x, value_type __bottom_y, value_type __bottom_z):
        bounding_box_base<T>(),
        _top({__top_x,__top_y,__top_z}),
        _bottom({__bottom_x,__bottom_y,__bottom_z})
    {}

    virtual ~box_bounding(){}

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

#endif // BOX_BOUNDING_H
