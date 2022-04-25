#ifndef SPHERE_H
#define SPHERE_H

namespace rvegen {

template <typename T>
class sphere : public shape_base<T>
{
public:
    using value_type = T;

    sphere() {}
    virtual ~sphere() {}

private:
    value_type _radius;
    std::array<value_type, 3> _point;
};
}
#endif // SPHERE_H
