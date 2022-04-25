#ifndef SHAPE_BASE_H
#define SHAPE_BASE_H

#include <memory>
#include "bounding_box_base.h"

namespace rvegen {

template<typename T>
class shape_base
{
public:
    using value_type = T;

    shape_base():
        _bounding_box()
    {}

    virtual ~shape_base() {}

    //functionen wie area
    virtual value_type area() const = 0;

    //bsp function
//    virtual void print() const = 0;

    bounding_box_base<T>* bounding_box()const{
        return _bounding_box.get();
    }

    virtual void make_bounding_box() = 0;

protected:
    std::unique_ptr<bounding_box_base<value_type>> _bounding_box;
};

}
#endif // SHAPE_BASE_H
