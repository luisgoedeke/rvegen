#ifndef BOUNDING_BOX_BASE_H
#define BOUNDING_BOX_BASE_H

#include <array>

namespace rvegen {

template <typename T>
class bounding_box_base
{

    using value_type = T;
public:
    bounding_box_base() {}

    //virtual bool collision(bounding_box_base * __data) const = 0;

    virtual inline std::array<value_type, 3> const& top_point()const = 0;

    virtual inline std::array<value_type, 3> const& bottom_point()const = 0;

};



}
#endif // BOUNDING_BOX_BASE_H
