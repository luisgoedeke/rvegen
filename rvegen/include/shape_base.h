#ifndef SHAPE_BASE_H
#define SHAPE_BASE_H

namespace rvegen {

template<typename T>
class shape_base
{
public:
    using value_type = T;

    shape_base() {}

    virtual ~shape_base() {}

    //functionen wie area
    virtual value_type area() const = 0;

    //bsp function
    virtual void print() const = 0;

private:
};

}
#endif // SHAPE_BASE_H
