#ifndef POLYGON_H
#define POLYGON_H

#include <cstdlib>
#include <vector>
#include <array>
#include "shape_base.h"

namespace rvegen {

template <typename T, std::size_t Dim>
class polygon : public shape_base<T>
{
public:
    using value_type = T;

    polygon():_data()
    {}

    virtual ~polygon(){}

    virtual value_type area() const {
        //https://de.wikipedia.org/wiki/Polygon
        const auto size{_data.size()};
        value_type area{0};
        for(std::size_t i{0}; i<size; ++i){
            std::size_t j{(i+1)%size};
            area += (_data[i][0]*_data[j][1] -  _data[j][0]*_data[i][1]);
        }
        return 0.5*area;
    }

    virtual void make_bounding_box() {

    }

private:
    std::vector<std::array<value_type, Dim>> _data;
};


}
#endif // POLYGON_H
