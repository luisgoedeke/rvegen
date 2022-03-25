#ifndef RVE_GENERATOR_H
#define RVE_GENERATOR_H

#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <memory>
#include <map>

#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"
#include "check_distance.h"
#include "rve_shape_input.h"

namespace rvegen {

enum rveType{
    Periodic,
    Random,
    OnlyInside,
};

//base point is at (0,0,0)
//          ________(x,y,z) --> box (breite, länge, höhe)
//         /      / |
//        /______/  |
//        |      |  |
//        |      | /
// (0,0,0)|______|/


template <typename _Distribution>
class rve_generator
{
public:
    using value_type = double;
    using size_type = std::size_t;


    rve_generator(rveType const __rve_type = rveType::OnlyInside, std::size_t const __dimension = 3, value_type const _x = 1, value_type const _y = 1, value_type const _z = 1):
        _rve_type(__rve_type),
        _dim(__dimension),
        _max_iter(500000),
        _vol_frac_inclusion(0),
        _box{_x,_y,_z}
    {}

    template<typename _Generator>
    constexpr inline auto compute(std::vector<std::unique_ptr<rvegen::rve_shape_input>> const& __input, _Generator& __random_generator){
        //check sum volume fraction <= 1

        if(__input.size() == 1){
            //single inclusion type
            //check which type of inclusion and do some single inclusion stuff
            if(dynamic_cast<circle_input*>(__input[0].get())){
                //check dimension
                if(_dim != 2){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_circle(*static_cast<circle_input*>(__input[0].get()), __random_generator);
            }

            else if(dynamic_cast<ellipse_input*>(__input[0].get())){
                //check dimension
                if(_dim != 2){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_ellipse(*static_cast<ellipse_input*>(__input[0].get()), __random_generator);
            }

            else if (dynamic_cast<cylinder_input*>(__input[0].get())){
                if (_dim !=3){
                   throw std::runtime_error("no matching dimensions");
                }
                compute_single_cylinder(*static_cast<cylinder_input*>(__input[0].get()), __random_generator);
            }
        }else{
            //do some more complex stuff here with different inclusions
        }
    }

    constexpr inline auto& dimension(){
        return _dim;
    }

    constexpr inline auto dimension()const{
        return _dim;
    }


    constexpr inline auto box()const{
        return _box;
    }

    constexpr inline auto const& shapes()const{
        return _shapes;
    }

    constexpr inline auto rve_type()const{
        return _rve_type;
    }

    constexpr inline auto& rve_type(){
        return _rve_type;
    }

private:
    template<typename _Generator>
    constexpr inline auto compute_single_circle(circle_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_circle_only_inside(circle_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_circle_periodic(circle_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_circle_random(circle_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipse(ellipse_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipse_only_inside(ellipse_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipse_periodic(ellipse_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipse_random(ellipse_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder(cylinder_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_only_inside(cylinder_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_periodic(cylinder_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_random(cylinder_input const& __input, _Generator& __random_generator);

private:
    rveType _rve_type;
    size_type _dim;
    size_type _max_iter;
    value_type _vol_frac_inclusion;
    std::array<value_type, 3> _box;
    std::vector<std::unique_ptr<shape_base<value_type>>> _shapes;
};


template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_circle(circle_input const& __input, _Generator& __random_generator){

#ifdef RVE_DEBUG
    //print only in debug mode
    //do some stuff
    std::cout<<"Do some random circle stuff"<<std::endl;
    std::cout<<"radius min      "<<__input.get_radius_min()<<std::endl;
    std::cout<<"radius max      "<<__input.get_radius_max()<<std::endl;
    std::cout<<"volume fraction "<<__input.get_volume_fraction()<<std::endl;
#endif

    if(_rve_type == rveType::OnlyInside){
        compute_single_circle_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_circle_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_circle_random(__input, __random_generator);
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipse(ellipse_input const& __input, _Generator& __random_generator){

#ifdef RVE_DEBUG
    //print only in debug mode
    //do some stuff
    std::cout<<"Do some random ellipse stuff"<<std::endl;
    std::cout<<"radius min a     "<<__input.get_radius_min_a()<<std::endl;
    std::cout<<"radius max a     "<<__input.get_radius_max_a()<<std::endl;
    std::cout<<"radius min b     "<<__input.get_radius_min_b()<<std::endl;
    std::cout<<"radius max b     "<<__input.get_radius_max_b()<<std::endl;
    std::cout<<"volume fraction "<<__input.get_volume_fraction()<<std::endl;
#endif
    if(_rve_type == rveType::OnlyInside){
        compute_single_ellipse_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_ellipse_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_ellipse_random(__input, __random_generator);
    }
}

template<typename T>
bool intersection (ellipse<T> const& e, T x1, T y1, T x2, T y2){

    using value_type = T;

    //Translation ins Koordinatensystem der Ellipse
    value_type x1t = x1-e(0);
    value_type y1t = y1-e(1);
    value_type x2t = x2-e(0);
    value_type y2t = y2-e(1);

    //Rotation ins Koordinatensystem der Ellipse
    value_type x1r = cos(e.rotation()*M_PI)*x1t+sin(e.rotation()*M_PI)*y1t;
    value_type y1r = -sin(e.rotation()*M_PI)*x1t+cos(e.rotation()*M_PI)*y1t;
    value_type x2r = cos(e.rotation()*M_PI)*x2t+sin(e.rotation()*M_PI)*y2t;
    value_type y2r = -sin(e.rotation()*M_PI)*x2t+cos(e.rotation()*M_PI)*y2t;

    value_type m = ((y2r-y1r)/(x2r-x1r));
    value_type d = y1r-m*x1r;

    value_type a = e.radius_a();
    value_type b = e.radius_b();

    //Prüfung ob Term unter Wurzel der p-q-Formel >= 0 ist => Schnittpunkte sind vorhanden

    value_type p = (2*m*d*a*a)/(b*b+m*m*a*a);
    value_type q = (a*a*d*d-a*a*b*b)/(b*b+m*m*a*a);

    if (((p/2)*(p/2)-q) >=0){
        return true;
    }

    return false;
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder(cylinder_input const& __input, _Generator& __random_generator){

#ifdef RVE_DEBUG
    //print only in debug mode
    //do some stuff
    std::cout<<"Do some random cylinder stuff"<<std::endl;
    std::cout<<"radius min      "<<__input.get_radius_min()<<std::endl;
    std::cout<<"radius max      "<<__input.get_radius_max()<<std::endl;
    std::cout<<"height min      "<<__input.get_height_min()<<std::endl;
    std::cout<<"height max      "<<__input.get_height_max()<<std::endl;
    std::cout<<"volume fraction "<<__input.get_volume_fraction()<<std::endl;
#endif

    if(_rve_type == rveType::OnlyInside){
        compute_single_cylinder_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_cylinder_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_cylinder_random(__input, __random_generator);
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_circle_random(circle_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);

    size_type iter{0};
    const value_type area{minR*minR*M_PI};
    const size_type max_circles{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_circles);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const double radius{dis_radius(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const circle<value_type> circle_(x, y, radius);

        if(!check_distance(_shapes, circle_)){
            _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_));
            _vol_frac_inclusion += circle_.area();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipse_random(ellipse_input const& __input, _Generator& __random_generator){
    const value_type minRa{__input.get_radius_min_a()}, maxRa{__input.get_radius_max_a()}, minRb{__input.get_radius_min_b()}, maxRb{__input.get_radius_max_b()}, minRot{__input.get_min_rotation()}, maxRot{__input.get_max_rotation()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius_a(minRa, maxRa);
    _Distribution dis_radius_b(minRb, maxRb);
    _Distribution dis_rotation(minRot, maxRot);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);

    size_type iter{0};
    const value_type area{minRa*minRb*M_PI};
    const size_type max_ellipses{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_ellipses);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const double radius_a{dis_radius_a(__random_generator)};
        const double radius_b{dis_radius_b(__random_generator)};
        const double rotation{dis_rotation(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const ellipse<value_type> ellipse_(x, y, radius_a, radius_b, rotation);

        if(!check_distance(_shapes, ellipse_)){
            _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_));
            _vol_frac_inclusion += ellipse_.area();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }

}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_random(cylinder_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()}, minH{__input.get_height_min()}, maxH{__input.get_height_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_height(minH, maxH);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);
    _Distribution dis_z(-_box[2],_box[2]);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type area{minR*minR*M_PI};
    const value_type volume {minR*minR*M_PI*minH};
    const size_type max_cylinders{static_cast<size_type>((dx*dy*dz*volume_faction)/volume)};

    size_type iter{0};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_cylinders);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new cylinder
        const double radius{dis_radius(__random_generator)};
        const double height{dis_height(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        const cylinder<value_type> cylinder_(x, y, z, radius, height);

        if(!check_distance(_shapes, cylinder_)){
            _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_));
            _vol_frac_inclusion += cylinder_.volume();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_circle_periodic(circle_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);

    const value_type dx{_box[0]}, dy{_box[1]};

    const value_type area{minR*minR*M_PI};
    const size_type max_circles{static_cast<size_type>((dx*dy*volume_faction)/area)};

    size_type iter{0};

    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_circles);

    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const value_type radius{dis_radius(__random_generator)};
        const value_type x{dis_x(__random_generator)};
        const value_type y{dis_y(__random_generator)};
        const circle<value_type> circle_(x, y, radius);

        bool check_distance_{check_distance(_shapes, circle_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            if((circle_(0)-circle_.radius()) < 0){
                is_outside = true;
                const circle<value_type> circle_periodic(x+dx,y,radius);
                check_distance_periodic = check_distance(_shapes, circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            //right side
            if(circle_(0)+circle_.radius()>dx){
                is_outside = true;
                const circle<value_type> circle_periodic(x-dx,y,radius);
                check_distance_periodic = check_distance(_shapes, circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            //bottom
            if(circle_(1)-circle_.radius()<0){
                is_outside = true;
                const circle<value_type> circle_periodic(x,y+dy,radius);
                check_distance_periodic = check_distance(_shapes, circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            //top
            if(circle_(1)+circle_.radius()>dy){
                is_outside = true;
                const circle<value_type> circle_periodic(x,y-dy,radius);
                check_distance_periodic = check_distance(_shapes, circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_));
                _vol_frac_inclusion += circle_.area();
            }
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(dx*dy))) < 0.005){
            break;
        }
        ++iter;
    }

    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipse_periodic(ellipse_input const& __input, _Generator& __random_generator){
    const value_type minRa{__input.get_radius_min_a()}, maxRa{__input.get_radius_max_a()}, minRb{__input.get_radius_min_b()}, maxRb{__input.get_radius_max_b()}, minRot{__input.get_min_rotation()}, maxRot{__input.get_max_rotation()};
    const value_type volume_faction{__input.get_volume_fraction()};


    _Distribution dis_radius_a(minRa, maxRa);
    _Distribution dis_radius_b(minRb, maxRb);
    _Distribution dis_rotation(minRot, maxRot);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);

    const value_type dx{_box[0]}, dy{_box[1]};

    const value_type area{minRa*minRb*M_PI};
    const size_type max_ellipses{static_cast<size_type>((dx*dy*volume_faction)/area)};

    size_type iter{0};

    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_ellipses);

    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const value_type radius_a{dis_radius_a(__random_generator)};
        const value_type radius_b{dis_radius_b(__random_generator)};
        const value_type rotation{dis_rotation(__random_generator)};
        const value_type x{dis_x(__random_generator)};
        const value_type y{dis_y(__random_generator)};
        const ellipse<value_type> ellipse_(x, y, radius_a, radius_b, rotation);

        bool check_distance_{check_distance(_shapes, ellipse_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            if(intersection(ellipse_, 0.0, 0.0, 0.0, 1.0)){
                is_outside = true;
                const ellipse<value_type> ellipse_periodic(x+dx,y,radius_a, radius_b, rotation);
                check_distance_periodic = check_distance(_shapes, ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            //right side
            if(intersection(ellipse_, 1.0, 0.0, 1.0, 1.0)){
                is_outside = true;
                const ellipse<value_type> ellipse_periodic(x-dx,y,radius_a, radius_b, rotation);
                check_distance_periodic = check_distance(_shapes, ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            //bottom
            if(intersection(ellipse_, 0.0, 0.0, 1.0, 0.0)){
                is_outside = true;
                const ellipse<value_type> ellipse_periodic(x,y+dy,radius_a, radius_b, rotation);
                check_distance_periodic = check_distance(_shapes, ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            //top
            if(intersection(ellipse_, 0.0, 1.0, 1.0, 1.0)){
                is_outside = true;
                const ellipse<value_type> ellipse_periodic(x,y-dy,radius_a, radius_b, rotation);
                check_distance_periodic = check_distance(_shapes, ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_));
                _vol_frac_inclusion += ellipse_.area();
            }
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(dx*dy))) < 0.005){
            break;
        }
        ++iter;
    }

    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_periodic(cylinder_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()}, minH{__input.get_height_min()}, maxH{__input.get_height_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_height(minH, maxH);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);
    _Distribution dis_z(0, _box[2]);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type area{minR*minR*M_PI};
    const value_type volume {minR*minR*M_PI*minH};
    const size_type max_cylinders{static_cast<size_type>((dx*dy*dz*volume_faction)/volume)};

    size_type iter{0};

    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_cylinders);

    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new cylinder
        const double radius{dis_radius(__random_generator)};
        const double height{dis_height(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        const cylinder<value_type> cylinder_(x, y, z, radius, height);

        bool check_distance_{check_distance(_shapes, cylinder_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            bool is_outside_left{false};
            bool is_outside_right{false};
            bool is_outside_bottom{false};
            bool is_outside_top{false};
            //left side
            if((cylinder_(0)-cylinder_.radius()) < 0){
                is_outside = true;
                is_outside_left = true;
                const cylinder<value_type> cylinder_periodic(x+dx,y,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //right side
            if(cylinder_(0)+cylinder_.radius()>dx){
                is_outside = true;
                is_outside_right = true;
                const cylinder<value_type> cylinder_periodic(x-dx,y,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //bottom
            if(cylinder_(1)-cylinder_.radius()<0){
                is_outside = true;
                is_outside_bottom = true;
                const cylinder<value_type> cylinder_periodic(x,y+dy,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //top
            if(cylinder_(1)+cylinder_.radius()>dy){
                is_outside = true;
                bool is_outside_top = true;
                const cylinder<value_type> cylinder_periodic(x,y-dy,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }

            //front
             if(cylinder_(2)+cylinder_.height()>dz){
                 is_outside = true;
                 const cylinder<value_type> cylinder_periodic(x,y,0,radius,(height-(dz-z)));
                 check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                 if(check_distance_periodic == false){
                     _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                 }
                 //front+left
                 if(is_outside_left == true){
                     const cylinder<value_type> cylinder_periodic(x+dx,y,0,radius,(height-(dz-z)));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //front+right
                 if(is_outside_right == true){
                     const cylinder<value_type> cylinder_periodic(x-dx,y,0,radius,(height-(dz-z)));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //front+bottom
                 if(is_outside_bottom == true){
                     const cylinder<value_type> cylinder_periodic(x,y+dy,0,radius,(height-(dz-z)));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //front+top
                 if(is_outside_top == true){
                     const cylinder<value_type> cylinder_periodic(x,y-dy,0,radius,(height-(dz-z)));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
             }

             //back
             if(cylinder_(2)+cylinder_.height()<0){
                 is_outside = true;
                 const cylinder<value_type> cylinder_periodic(x,y,_box[2],radius,(height+z));
                 check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                 if(check_distance_periodic == false){
                     _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                 }
                 //back+left
                 if(is_outside_left == true){
                     const cylinder<value_type> cylinder_periodic(x+dx,y,0,radius,(height+z));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //back+right
                 if(is_outside_right == true){
                     const cylinder<value_type> cylinder_periodic(x-dx,y,0,radius,(height+z));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //back+bottom
                 if(is_outside_bottom == true){
                     const cylinder<value_type> cylinder_periodic(x,y+dy,0,radius,(height+z));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //back+top
                 if(is_outside_top == true){
                     const cylinder<value_type> cylinder_periodic(x,y-dy,0,radius,(height+z));
                     check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
             }

            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_));
                _vol_frac_inclusion += cylinder_.volume();
            }
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(dx*dy))) < 0.005){
            break;
        }
        ++iter;
    }

    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_circle_only_inside(circle_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_x(maxR, (_box[0] - maxR));
    _Distribution dis_y(maxR, (_box[1] - maxR));

    size_type iter{0};
    const value_type area{minR*minR*M_PI};
    const size_type max_circles{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_circles);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const double radius{dis_radius(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const circle<value_type> circle_(x, y, radius);

        if(!check_distance(_shapes, circle_)){
            _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_));
            _vol_frac_inclusion += circle_.area();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipse_only_inside(ellipse_input const& __input, _Generator& __random_generator){
    const value_type minRa{__input.get_radius_min_a()}, maxRa{__input.get_radius_max_a()}, minRb{__input.get_radius_min_b()}, maxRb{__input.get_radius_max_b()}, minRot{__input.get_min_rotation()}, maxRot{__input.get_max_rotation()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius_a(minRa, maxRa);
    _Distribution dis_radius_b(minRb, maxRb);
    _Distribution dis_rotation(minRot, maxRot);
    _Distribution dis_x(maxRa, _box[0]-maxRa);
    _Distribution dis_y(maxRb, _box[1]-maxRa);

    size_type iter{0};
    const value_type area{minRa*minRb*M_PI};
    const size_type max_ellipses{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_ellipses);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new ellipse
        const double radius_a{dis_radius_a(__random_generator)};
        const double radius_b{dis_radius_b(__random_generator)};
        const double rotation{dis_rotation(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const ellipse<value_type> ellipse_(x, y, radius_a, radius_b, rotation);

        if(!check_distance(_shapes, ellipse_)){
            _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_));
            _vol_frac_inclusion += ellipse_.area();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_only_inside(cylinder_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()}, minH{__input.get_height_min()}, maxH{__input.get_height_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_height(minH, maxH);
    _Distribution dis_x(maxR, (_box[0] - maxR));
    _Distribution dis_y(maxR, (_box[1] - maxR));
    _Distribution dis_z(0, _box[2]);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type area{minR*minR*M_PI};
    const value_type volume {minR*minR*M_PI*minH};
    const size_type max_cylinders{static_cast<size_type>((dx*dy*dz*volume_faction)/volume)};

    size_type iter{0};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_cylinders);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new cylinder
        const double radius{dis_radius(__random_generator)};
        const double height{dis_height(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        const cylinder<value_type> cylinder_(x, y, z, radius, height);


        if((!check_distance(_shapes, cylinder_)) && (cylinder_.height() + cylinder_(2) < _box[2]) && (cylinder_(2) > 0)){
            _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_));
            _vol_frac_inclusion += cylinder_.volume();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}



//template<typename Type, typename Distribution>
//class rve_generator_3D;

//template<typename Type, typename Distribution>
//class rve_generator_2D;

//template<typename T, typename Distribution>
//class rve_generator_2D<circle<T>, Distribution>
//{
//public:
//    using value_type = T;
//    using size_type = std::size_t;
//    using Type = circle<value_type>;

//    rve_generator_2D(value_type const volume_faction_inclusion, value_type const x=1, value_type const y=1)
//        :_data(),_point(x,y),_volume_faction_inclusion(volume_faction_inclusion),_max_iter(5000000),_area(0) {}

//    rve_generator_2D(value_type const volume_faction_inclusion, point<value_type, 2> const& point)
//        :_data(),_point(point),_volume_faction_inclusion(volume_faction_inclusion),_max_iter(5000000),_area(0) {}

//    template<typename Generator>
//    void random_only_inside(Generator& gen, value_type const min_radius, value_type const max_radus){
//         Distribution dis_radius(min_radius,max_radus);
//         Distribution dis_x(max_radus, (_point(0) - max_radus));
//         Distribution dis_y(max_radus, (_point(1) - max_radus));

//         size_type iter{0};

//         const double area{min_radius*min_radius*M_PI};
//         const size_type anzahl_der_kreise{static_cast<size_type>((_point(0)*_point(1)*_volume_faction_inclusion)/area)};


//         bool finished{false};
//         for(std::size_t i{0};i<10;++i){
//             std::cout<<"Attempt: "<<i<<std::endl;
//             _data.clear();
//             _data.reserve(anzahl_der_kreise);//it is faster to push back new elements
//             iter = 0;
//             _area = 0;
//             while (true) {
//                 //neuer kreis
//                 const double radius{dis_radius(gen)};

//                 const auto x{dis_x(gen)};
//                 const auto y{dis_y(gen)};
//                 const circle<> circle_(x,y,radius);

//                 if(!check_distance(_data,circle_)){
//                     _data.push_back(circle_);
//                     _area += circle_.area();
//                     //break;
//                 }

//                 if(iter == _max_iter){
//                     break;
//                 }

//                 if((_volume_faction_inclusion - (_area/(_point(0)*_point(1)))) < 0.005){
//                     std::cout<<"Erfolgreich beendet"<<std::endl;
//                     std::cout<<"Volume faction: "<<_area/(_point(0)*_point(1))<<std::endl;
//                     finished = true;
//                     break;
//                 }
//                 ++iter;
//             }

//             if(finished){
//                 break;
//             }
//         }

//         if(!finished){
//             std::cout<<"Max iteration reached"<<std::endl;
//             std::cout<<"Volume faction: "<<_area/(_point(0)*_point(1))<<std::endl;

//         }
//    }


//    template<typename Generator>
//    void random(Generator& gen, value_type const min_radius, value_type const max_radus){
//        Distribution dis_radius(min_radius,max_radus);
//        Distribution dis_x(0, _point(0));
//        Distribution dis_y(0, _point(1));

//        const auto dx{_point(0)}, dy{_point(1)};
//        size_type iter{0};

//        const double area{min_radius*min_radius*M_PI};
//        const size_type anzahl_der_kreise{static_cast<size_type>((_point(0)*_point(1)*_volume_faction_inclusion)/area)};

//        bool finished{false};
//        for(std::size_t i{0};i<10;++i){
//            std::cout<<"Attempt: "<<i<<std::endl;
//            _data.clear();
//            _data.reserve(anzahl_der_kreise);//it is faster to push back new elements
//            iter = 0;
//            _area = 0;
//            while (true) {
//                //neuer kreis
//                const double radius{dis_radius(gen)};
//                const auto x{dis_x(gen)};
//                const auto y{dis_y(gen)};
//                const circle<> circle_(x,y,radius);

//                bool check_distance_{check_distance(_data, circle_)};
//                if(check_distance_ == false){
//                    //check if outside of rve
//                    bool check_distance_periodic{true};
//                    bool is_outside{false};
//                    //left side
//                    if(circle_(0)-circle_.radius()<0){
//                        is_outside = true;
//                        const circle<> circle_periodic(x+dx,y,radius);
//                        check_distance_periodic = check_distance(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    //right side
//                    if(circle_(0)+circle_.radius()>dx){
//                        is_outside = true;
//                        const circle<> circle_periodic(x-dx,y,radius);
//                        check_distance_periodic = check_distance(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    //bottom
//                    if(circle_(1)-circle_.radius()<0){
//                        is_outside = true;
//                        const circle<> circle_periodic(x,y+dy,radius);
//                        check_distance_periodic = check_distance(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    //top
//                    if(circle_(1)+circle_.radius()>dy){
//                        is_outside = true;
//                        const circle<> circle_periodic(x,y-dy,radius);
//                        check_distance_periodic = check_distance(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    if(check_distance_periodic == false || is_outside == false){
//                        _data.push_back(circle_);
//                        _area += circle_.area();
//                    }
//                }

//                if(iter == _max_iter){break;}

//                if((_volume_faction_inclusion - (_area/(_point(0)*_point(1)))) < 0.005){
//                    std::cout<<"Erfolgreich beendet"<<std::endl;
//                    std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;
//                    finished = true;
//                    break;
//                }
//                ++iter;
//            }

//            if(finished){
//                break;
//            }
//        }

//        if(!finished){
//            std::cout<<"Max iteration reached"<<std::endl;
//            std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;
//        }
//    }

//    auto max_iter()const{
//        return _max_iter;
//    }

//    auto& max_iter(){
//        return _max_iter;
//    }

//    template<typename File>
//    void write_gmsh_geo_file(File & file){
//        const auto dx{_point(0)}, dy{_point(1)};
//        file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
//        file<<"Rectangle(1) = {0, 0, 0,"<<dx<<","<<dy<<", 0};"<<std::endl;
//        for(size_t i{0};i<_data.size();++i){
//            file<<"Circle("<<5+i<<") = {"<<_data[i](0)<<","<<_data[i](1)<<", 0,"<< _data[i].radius()<<", 0, 2*Pi};"<<std::endl;
//        }

//        for(size_t i{0};i<_data.size();++i){
//            file<<"Curve Loop("<<i+2<<") = {"<<5+i<<"};"<<std::endl;
//        }
//        for(size_t i{0};i<_data.size();++i){
//            file<<"Plane Surface("<<i+2<<") = {"<<i+2<<"};"<<std::endl;
//        }

//        for(size_t i{0};i<_data.size();++i){
//            file<<"BooleanIntersection{ Surface{1}; }{ Surface{"<<i+2<<"}; Delete; }"<<std::endl;
//            file<<"BooleanDifference{ Surface{1}; Delete; }{ Surface{"<<i+2<<"}; }"<<std::endl;
//        }
//    }

//private:
//    std::vector<Type> _data;
//    point<value_type, 2> _point;
//    value_type _volume_faction_inclusion;
//    size_type _max_iter;
//    value_type _area;
//};


/*
////------------------------------------------------------------------------------------------------
////-------------------------------------------cylinder---------------------------------------------
////------------------------------------------------------------------------------------------------
template<typename Distribution>
class rve_generator_3d
{

public:
    using value_type = double;
    using size_type = std::size_t;

    rve_generator_3d(rveType const __rve_type = rveType::OnlyInside, std::size_t const __dimension = 3, value_type const _x = 1, value_type const _y = 1, value_type const _z = 1):
        _rve_type(__rve_type),
        _dim(__dimension),
        _max_iter(50000),
        _vol_frac_inclusion(0),
        _box{_x,_y,_z}
    {}

    template<typename _Generator>
    constexpr inline auto compute(std::vector<std::unique_ptr<rvegen::rve_shape_input>> const& __input, _Generator& __random_generator){
        //check sum volume fraction <= 1

        if(__input.size() == 1){
            //single inclusion type
            //check which type of inclusion and do some single inclusion stuff
            if(dynamic_cast<cylinder_input*>(__input[0].get())){
                //check dimension
                if(_dim != 3){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_cylinder(*static_cast<cylinder_input*>(__input[0].get()), __random_generator);
            }
        }else{
            //do some more complex stuff here with different inclusions
        }
    }

    constexpr inline auto& dimension(){
        return _dim;
    }

    constexpr inline auto dimension()const{
        return _dim;
    }


    constexpr inline auto box()const{
        return _box;
    }

    constexpr inline auto const& shapes()const{
        return _shapes;
    }

    constexpr inline auto rve_type()const{
        return _rve_type;
    }

    constexpr inline auto& rve_type(){
        return _rve_type;
    }

private:
    template<typename _Generator>
    constexpr inline auto compute_single_cylinder(cylinder_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_only_inside(cylinder_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_periodic(cylinder_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_random(cylinder_input const& __input, _Generator& __random_generator);
private:
    rveType _rve_type;
    size_type _dim;
    size_type _max_iter;
    value_type _vol_frac_inclusion;
    std::array<value_type, 3> _box;
    std::vector<std::unique_ptr<shape_base<value_type>>> _shapes;
};

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator_3d<_Distribution>::compute_single_cylinder(cylinder_input const& __input, _Generator& __random_generator){

#ifdef RVE_DEBUG
    //print only in debug mode
    //do some stuff
    std::cout<<"Do some random cylinder stuff"<<std::endl;
    std::cout<<"radius min      "<<__input.get_radius_min()<<std::endl;
    std::cout<<"radius max      "<<__input.get_radius_max()<<std::endl;
    std::cout<<"volume fraction "<<__input.get_volume_fraction()<<std::endl;
#endif

    if(_rve_type == rveType::OnlyInside){
        compute_single_cylinder_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_cylinder_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_cylinder_random(__input, __random_generator);
    }
}


template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator_3d<_Distribution>::compute_single_cylinder_random(cylinder_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()}, minH{__input.get_height_min()}, maxH{__input.get_height_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_height(minH, maxH);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);
    _Distribution dis_z(0,0);
  //  _Distribution dis_z(0, _box[2]);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type area{minR*minR*M_PI};
    const value_type volume {minR*minR*M_PI*minH};
    const size_type max_cylinders{static_cast<size_type>((dx*dy*dz*volume_faction)/volume)};

    size_type iter{0};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_cylinders);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new cylinder
        const double radius{dis_radius(__random_generator)};
        const double height{dis_height(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        const cylinder<value_type> cylinder_(x, y, z, radius, height);

        if(!check_distance(_shapes, cylinder_)){
            _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_));
            _vol_frac_inclusion += cylinder_.volume();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator_3d<_Distribution>::compute_single_cylinder_periodic(cylinder_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()}, minH{__input.get_height_min()}, maxH{__input.get_height_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_height(minH, maxH);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);
    _Distribution dis_z(0, _box[2]);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type area{minR*minR*M_PI};
    const value_type volume {minR*minR*M_PI*minH};
    const size_type max_cylinders{static_cast<size_type>((dx*dy*dz*volume_faction)/volume)};

    size_type iter{0};

    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_cylinders);

    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new cylinder
        const double radius{dis_radius(__random_generator)};
        const double height{dis_height(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        const cylinder<value_type> cylinder_(x, y, z, radius, height);

        bool check_distance_{check_distance(_shapes, cylinder_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            if((cylinder_(0)-cylinder_.radius()) < 0){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x+dx,y,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //right side
            if(cylinder_(0)+cylinder_.radius()>dx){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x-dx,y,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //bottom
            if(cylinder_(1)-cylinder_.radius()<0){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x,y+dy,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //top
            if(cylinder_(1)+cylinder_.radius()>dy){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x,y-dy,z,radius,height);
                check_distance_periodic = check_distance(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_));
                _vol_frac_inclusion += cylinder_.volume();
            }
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(dx*dy))) < 0.005){
            break;
        }
        ++iter;
    }

    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator_3d<_Distribution>::compute_single_cylinder_only_inside(cylinder_input const& __input, _Generator& __random_generator){
    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()}, minH{__input.get_height_min()}, maxH{__input.get_height_max()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius(minR, maxR);
    _Distribution dis_height(minH, maxH);
    _Distribution dis_x(maxR, (_box[0] - maxR));
    _Distribution dis_y(maxR, (_box[1] - maxR));
    _Distribution dis_z(0, 0);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type area{minR*minR*M_PI};
    const value_type volume {minR*minR*M_PI*minH};
    const size_type max_cylinders{static_cast<size_type>((dx*dy*dz*volume_faction)/volume)};

    size_type iter{0};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_cylinders);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new cylinder
        const double radius{dis_radius(__random_generator)};
        const double height{dis_height(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        const cylinder<value_type> cylinder_(x, y, z, radius, height);


        if(!check_distance(_shapes, cylinder_)){
            _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_));
            _vol_frac_inclusion += cylinder_.volume();
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
            finished = true;
            break;
        }

        ++iter;
    }


    if(iter == _max_iter){
        throw std::runtime_error("max iterations reached");
    }
}


/*    using Type = circle<value_type>;

    rve_generator_3D(value_type const volume_faction_inclusion, value_type const x=1, value_type const y=1, value_type const z=1)
        :_data(),_point(x,y,z),_volume_faction_inclusion(volume_faction_inclusion),_max_iter(5000000),_area(0) {}

    rve_generator_3D(value_type const volume_faction_inclusion, point<value_type, 3> const& point)
        :_data(),_point(point),_volume_faction_inclusion(volume_faction_inclusion),_max_iter(5000000),_area(0) {}


    template<typename Generator>
    void random_only_inside(Generator& gen, value_type const min_radius, value_type const max_radus){
         Distribution dis_radius(min_radius,max_radus);
         Distribution dis_x(max_radus, (_point(0) - max_radus));
         Distribution dis_y(max_radus, (_point(1) - max_radus));

         size_type iter{0};

         const double area{min_radius*min_radius*M_PI};
         const size_type anzahl_der_kreise{static_cast<size_type>((_point(0)*_point(1)*_volume_faction_inclusion)/area)};


         bool finished{false};
         for(std::size_t i{0};i<10;++i){
             std::cout<<"Attempt: "<<i<<std::endl;
             _data.clear();
             _data.reserve(anzahl_der_kreise);//it is faster to push back new elements
             iter = 0;
             _area = 0;
             while (true) {
                 //neuer kreis
                 const double radius{dis_radius(gen)};

                 const auto x{dis_x(gen)};
                 const auto y{dis_y(gen)};
                 const circle<> circle_(x,y,radius);

                 if(!check_distance(_data,circle_)){
                     _data.push_back(circle_);
                     _area += circle_.area();
                     break;
                 }

                 if(iter == _max_iter){

                     break;
                 }

                 if((_volume_faction_inclusion - (_area/(_point(0)*_point(1)))) < 0.005){
                     std::cout<<"Erfolgreich beendet"<<std::endl;
                     std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;
                     finished = true;
                     break;
                 }
                 ++iter;
             }

             if(finished){
                 break;
             }
         }

         if(!finished){
             std::cout<<"Max iteration reached"<<std::endl;
             std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;

         }
    }


    template<typename Generator>
    void random(Generator& gen, value_type const min_radius, value_type const max_radus){
        Distribution dis_radius(min_radius,max_radus);
        Distribution dis_x(0, _point(0));
        Distribution dis_y(0, _point(1));

        const auto dx{_point(0)}, dy{_point(1)};
        size_type iter{0};

        const double area{min_radius*min_radius*M_PI};
        const size_type anzahl_der_kreise{static_cast<size_type>((_point(0)*_point(1)*_volume_faction_inclusion)/area)};

        bool finished{false};
       for(std::size_t i{0};i<10;++i){
            std::cout<<"Attempt: "<<i<<std::endl;
            _data.clear();
            _data.reserve(anzahl_der_kreise);//it is faster to push back new elements
            iter = 0;
            _area = 0;
            while (true) {
                //neuer kreis
                const double radius{dis_radius(gen)};
                const auto x{dis_x(gen)};
                const auto y{dis_y(gen)};
                const circle<> circle_(x,y,radius);

                bool check_distance_{check_distance(_data, circle_)};
                if(check_distance_ == false){
                    //check if outside of rve
                    bool check_distance_periodic{true};
                    bool is_outside{false};
                    //left side
                    if(circle_(0)-circle_.radius()<0){
                        is_outside = true;
                        const circle<> circle_periodic(x+dx,y,radius);
                        check_distance_periodic = check_distance(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    //right side
                    if(circle_(0)+circle_.radius()>dx){
                        is_outside = true;
                        const circle<> circle_periodic(x-dx,y,radius);
                        check_distance_periodic = check_distance(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    //bottom
                    if(circle_(1)-circle_.radius()<0){
                        is_outside = true;
                        const circle<> circle_periodic(x,y+dy,radius);
                        check_distance_periodic = check_distance(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    //top
                    if(circle_(1)+circle_.radius()>dy){
                        is_outside = true;
                        const circle<> circle_periodic(x,y-dy,radius);
                        check_distance_periodic = check_distance(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    if(check_distance_periodic == false || is_outside == false){
                        _data.push_back(circle_);
                        _area += circle_.area();
                    }
                }

                if(iter == _max_iter){break;}

                if((_volume_faction_inclusion - (_area/(_point(0)*_point(1)))) < 0.005){
                    std::cout<<"Erfolgreich beendet"<<std::endl;
                    std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;
                    finished = true;
                    break;
                }
                ++iter;
            }

            if(finished){
                break;
            }
        }

        if(!finished){
            std::cout<<"Max iteration reached"<<std::endl;
            std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;
        }
    }

    auto max_iter()const{
        return max_iter;
    }

    auto& max_iter(){
        return max_iter;
    }

    template<typename File>
    void write_gmsh_geo_file(File & file){
        file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
        file<<"Box(1) = {0, 0, 0, 1, 1, 1};"<<std::endl;
        for(size_t i{0};i<_point.size();++i){
            file<<"Cylinder("<<i+2<<") = {"<<_point[i](0)<<", "<<_point[i](1)<<",0 , 0, 0, 1 ,"<<_point[i].radius()<<", 2*Pi};"<<std::endl;
        }
        for(size_t i{0};i<_point.size();++i){
            file<<"BooleanIntersection{ Volume{1}; }{ Volume{"<<i+2<<"}; Delete; }"<<std::endl;
            file<<"BooleanDifference{ Volume{1}; Delete; }{ Volume{"<<i+2<<"}; }"<<std::endl;
        }

        //+
        file<<"Physical Volume(1) = {1};"<<std::endl;
        file<<"Physical Volume(2) = {";
        for(size_t i{0};i<_point.size()-1;++i){
            file<<i+2<<", ";
        }
        file<<_point.size()+1<<"};"<<std::endl;
    }
*/
//  private:
//   std::vector<Type> _point;
//    point<value_type, 3> _point;
//    value_type _volume_faction_inclusion;
//    size_type _max_iter;
//    value_type _area;
//    std::array<value_type, 3> _box;

}
#endif // RVE_GENERATOR_H
