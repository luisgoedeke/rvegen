#ifndef RVE_SHAPE_INPUT_H
#define RVE_SHAPE_INPUT_H

namespace rvegen {

class rve_shape_input
{
public:
    using value_type = double;

    rve_shape_input():
        _volume_fraction(0),
        _random_position(true),
        _random_geometry(false)
    {}

    rve_shape_input(value_type const __volume_fraction, bool const __random_position, bool const __random_geometry = false):
        _volume_fraction(__volume_fraction),
        _random_position(__random_position),
        _random_geometry(__random_geometry)
    {}

    virtual ~rve_shape_input() {}

    inline auto set_volume_fraction(bool __val){
        _volume_fraction = __val;
    }

    inline auto get_volume_fraction()const{
        return _volume_fraction;
    }

    inline auto set_random_geometry(bool __val){
        _random_geometry = __val;
    }

    inline auto is_random_geometry()const{
        return _random_geometry;
    }

    inline auto set_random_position(bool __val){
        _random_position = __val;
    }

    inline auto is_random_position()const{
        return _random_position;
    }

    virtual void print() const = 0;

private:
    value_type _volume_fraction;
    bool _random_position;
    bool _random_geometry;
};



class circle_input : public rve_shape_input
{
public:
    using value_type = double;

    circle_input():
        _random_radius(),
        _min_radius(),
        _max_radius()
    {}

    circle_input(bool const __random_position, bool const __random_radius, value_type const _min_radius, value_type const _max_radius, value_type const __volume_fraction):
        rve_shape_input(__volume_fraction, __random_position),
        _random_radius(__random_radius),
        _min_radius(_min_radius),
        _max_radius(_max_radius)
    {}

    ~circle_input(){}

    inline auto set_random_radius(bool __val){
        _random_radius = __val;
    }

    inline auto is_random_radius()const{
        return _random_radius;
    }

    inline auto get_radius_min()const{
        return _min_radius;
    }

    inline auto get_radius_max()const{
        return _max_radius;
    }

    inline auto set_radius_min(double const __val){
        _min_radius = __val;
    }

    inline auto set_radius_max(double const __val){
        _max_radius = __val;
    }

    virtual void print() const override{
        std::cout<<"circle"<<std::endl;
    }

private:
    bool _random_radius;
    value_type _min_radius;
    value_type _max_radius;
};
}
#endif // RVE_SHAPE_INPUT_H
