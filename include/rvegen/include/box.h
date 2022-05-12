#ifndef BOX_H
#define BOX_H

#include "shape_base.h"

namespace rvegen {

template <typename T>
class box : public shape_base<T>
{
public:
    using value_type = T;

    box():
        _top(),
        _bottom(),
        _rotation()
    {}

    virtual ~box(){}

private:
    std::array<value_type, 3> _top;
    std::array<value_type, 3> _bottom;
    std::array<value_type, 3> _rotation;
};

}

#endif // BOX_H
