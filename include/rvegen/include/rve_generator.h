#ifndef RVE_GENERATOR_H
#define RVE_GENERATOR_H

#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <memory>
#include <map>
#include <utility>
#include <algorithm>

#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"
#include "ellipsoid.h"
#include "check_distance.h"
#include "rve_shape_input.h"
#include "samepack_heuristic.h"
#include "write_gmsh_geo.h"


namespace rvegen {

enum rveType{
    Periodic,
    Random,
    OnlyInside,
    OnlyInsideSamepackHeuristic,
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


    rve_generator(rveType const __rve_type = rveType::Random, std::size_t const __dimension = 3, value_type const _x = 1, value_type const _y = 1, value_type const _z = 1):
        _rve_type(__rve_type),
        _dim(__dimension),
        _max_iter(50000),
        _vol_frac_inclusion(0),
        number_of_shapes(0),
        _box{_x,_y,_z}
    {}

    template<typename _Generator>
    constexpr inline auto compute(std::vector<std::unique_ptr<rvegen::rve_shape_input>> const& __input, _Generator& __random_generator){
        //check sum volume fraction <= 1

        if(__input.size() == 1){
            compute_inclusion_stuff(__input[0].get(), __random_generator);
        }else{
            compute_multiple_inclusion_stuff(__input, __random_generator);
        }


        return ;
        if(__input.size() == 1){
            //single inclusion type
            //check which type of inclusion and do some single inclusion stuff
            if(dynamic_cast<circle_input*>(__input[0].get())){
                //check dimension
                if(_dim != 2){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_circle(*static_cast<circle_input*>(__input[0].get()), __random_generator);
            }else if(dynamic_cast<ellipse_input*>(__input[0].get())){
                //check dimension
                if(_dim != 2){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_ellipse(*static_cast<ellipse_input*>(__input[0].get()), __random_generator);
            }else if (dynamic_cast<cylinder_lf_input*>(__input[0].get())){
                if (_dim !=3){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_cylinder_lf(*static_cast<cylinder_lf_input*>(__input[0].get()), __random_generator);
            }else if (dynamic_cast<cylinder_sf_input*>(__input[0].get())){
                if (_dim !=3){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_cylinder_sf(*static_cast<cylinder_sf_input*>(__input[0].get()), __random_generator);
            }else if (dynamic_cast<ellipsoid_input*>(__input[0].get())){
                if (_dim !=3){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_ellipsoid(*static_cast<ellipsoid_input*>(__input[0].get()), __random_generator);
            }else if (dynamic_cast<rectangle_input*>(__input[0].get())){
                if (_dim !=2){
                    throw std::runtime_error("no matching dimensions");
                }
                compute_single_rectangle(*static_cast<rectangle_input*>(__input[0].get()), __random_generator);
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

    constexpr inline int get_number_of_shapes()const{
        return number_of_shapes;
    }

    constexpr inline int& get_number_of_shapes(){
        return number_of_shapes;
    }

    constexpr inline auto get_vol_frac_inclusion()const{
        return _vol_frac_inclusion;
    }

    constexpr inline auto& get_vol_frac_inclusion(){
        return _vol_frac_inclusion;
    }

private:

    template<typename _Generator>
    constexpr inline auto compute_multiple_inclusion_stuff(std::vector<std::unique_ptr<rvegen::rve_shape_input>>const& __input_shape, _Generator& __random_generator){
        std::cout<<"compute_inclusion_stuff"<<std::endl;

        if(_rve_type == rveType::OnlyInside){
            //compute_single_circle_only_inside(__input, __random_generator);
        }else if(_rve_type == rveType::Periodic){
            //compute_single_circle_periodic(__input, __random_generator);
        }else if(_rve_type == rveType::Random){
            compute_multiple_inclusion_stuff_random(__input_shape, __random_generator);
        }else if(_rve_type == rveType::OnlyInsideSamepackHeuristic){
            compute_multiple_inclusion_stuff_only_inside_samepack_heuristic(__input_shape, __random_generator);
        }
        number_of_shapes = _shapes.size();
    }


    template<typename _Generator>
    constexpr inline auto compute_inclusion_stuff(rve_shape_input* __input_shape, _Generator& __random_generator){
        std::cout<<"compute_inclusion_stuff"<<std::endl;

        if(_rve_type == rveType::OnlyInside){
            //compute_single_circle_only_inside(__input, __random_generator);
        }else if(_rve_type == rveType::Periodic){
            //compute_single_circle_periodic(__input, __random_generator);
        }else if(_rve_type == rveType::Random){
            compute_inclusion_stuff_random(__input_shape, __random_generator);
        }else if(_rve_type == rveType::OnlyInsideSamepackHeuristic){
            compute_inclusion_stuff_only_inside_samepack_heuristic(__input_shape, __random_generator);
        }
        number_of_shapes = _shapes.size();
    }

    template<typename _Generator>
    constexpr inline auto compute_inclusion_stuff_random(rve_shape_input* __input_shape, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_multiple_inclusion_stuff_random(std::vector<std::unique_ptr<rvegen::rve_shape_input>>const& __input_shape, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_inclusion_stuff_only_inside_samepack_heuristic(rve_shape_input* __input_shape, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_multiple_inclusion_stuff_only_inside_samepack_heuristic(std::vector<std::unique_ptr<rvegen::rve_shape_input>>const& __input_shape, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_rectangle(rectangle_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_rectangle_random(rectangle_input const& __input, _Generator& __random_generator);


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
    constexpr inline auto compute_single_cylinder_lf(cylinder_lf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_sf(cylinder_sf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_lf_only_inside(cylinder_lf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_sf_only_inside(cylinder_sf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_lf_periodic(cylinder_lf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_sf_periodic(cylinder_sf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_lf_random(cylinder_lf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_cylinder_sf_random(cylinder_sf_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipsoid(ellipsoid_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipsoid_only_inside(ellipsoid_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipsoid_periodic(ellipsoid_input const& __input, _Generator& __random_generator);

    template<typename _Generator>
    constexpr inline auto compute_single_ellipsoid_random(ellipsoid_input const& __input, _Generator& __random_generator);

private:
    rveType _rve_type;
    size_type _dim;
    size_type _max_iter;
    value_type _vol_frac_inclusion;
    int number_of_shapes;
    std::array<value_type, 3> _box;
    std::vector<std::unique_ptr<shape_base<value_type>>> _shapes;
    std::vector<std::pair<value_type, std::unique_ptr<shape_base<value_type>>>> _generated_shapes;
    std::vector<std::pair<int, int>> sections;
};

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

template<typename T>
auto max(T a, T b, T c){
    if ((a>=b)&&(a>=c)){
        return a;
    }
    if ((b>=a)&&(b>=c)){
        return b;
    }
    if ((c>=a)&&(c>=b)){
        return c;
    }
    return a;

}


template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_rectangle(rectangle_input const& __input, _Generator& __random_generator){
#ifdef RVE_DEBUG
    //print only in debug mode
    //do some stuff
    std::cout<<"Do some random rectangle stuff"<<std::endl;
    std::cout<<"volume fraction "<<__input.get_volume_fraction()<<std::endl;
#endif

    if(_rve_type == rveType::OnlyInside){
        //compute_single_circle_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        //compute_single_circle_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_rectangle_random(__input, __random_generator);
    }
}

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

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_lf(cylinder_lf_input const& __input, _Generator& __random_generator){

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
        compute_single_cylinder_lf_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_cylinder_lf_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_cylinder_lf_random(__input, __random_generator);
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_sf(cylinder_sf_input const& __input, _Generator& __random_generator){

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
        compute_single_cylinder_sf_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_cylinder_sf_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_cylinder_sf_random(__input, __random_generator);
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipsoid(ellipsoid_input const& __input, _Generator& __random_generator){

#ifdef RVE_DEBUG
    //print only in debug mode
    //do some stuff
    std::cout<<"Do some random ellipsoid stuff"<<std::endl;
    std::cout<<"radius min a     "<<__input.get_radius_min_a()<<std::endl;
    std::cout<<"radius max a     "<<__input.get_radius_max_a()<<std::endl;
    std::cout<<"radius min b     "<<__input.get_radius_min_b()<<std::endl;
    std::cout<<"radius max b     "<<__input.get_radius_max_b()<<std::endl;
    std::cout<<"radius min c     "<<__input.get_radius_min_c()<<std::endl;
    std::cout<<"radius max c     "<<__input.get_radius_max_c()<<std::endl;
    std::cout<<"volume fraction "<<__input.get_volume_fraction()<<std::endl;
#endif
    if(_rve_type == rveType::OnlyInside){
        compute_single_ellipsoid_only_inside(__input, __random_generator);
    }else if(_rve_type == rveType::Periodic){
        compute_single_ellipsoid_periodic(__input, __random_generator);
    }else if(_rve_type == rveType::Random){
        compute_single_ellipsoid_random(__input, __random_generator);
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_multiple_inclusion_stuff_random(std::vector<std::unique_ptr<rvegen::rve_shape_input>>const& __input_shapes, _Generator& __random_generator){
    std::cout<<"compute_multiple_inclusion_stuff_random"<<std::endl;

    value_type volume_faction{0};
    value_type area{0};
    for(auto & shape : __input_shapes){
        volume_faction += shape.get()->get_volume_fraction();
        area += shape.get()->min_area();
    }

    if(1.0 - volume_faction <= 0){
        throw std::runtime_error("rve_generator: volume_faction >= 1");
    }

    //reserve data for faster push back of new elements
    const size_type max_inclusions{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};
    _shapes.clear();
    _shapes.reserve(max_inclusions);
    for(auto & shape : __input_shapes){
        if(shape.get()->is_random_position()){
            if(dimension() == 2){
                shape.get()->setup_position(_box[0], _box[1]);
            }else{
                shape.get()->setup_position(_box[0], _box[1], _box[2]);
            }
        }

        volume_faction = shape.get()->get_volume_fraction();
        size_type iter{0};
        bool finished{false};
        _vol_frac_inclusion = 0;
        while (iter <= _max_iter) {
            //new shape
            auto new_shape = shape.get()->new_shape();
            new_shape.get()->make_bounding_box();

            if(!collision(_shapes, new_shape.get())){
                _vol_frac_inclusion += new_shape->area();
                _shapes.emplace_back(std::move(new_shape));
            }

#ifdef RVE_DEBUG
            std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
            std::fstream test_bsp;
            test_bsp.open("test_bsp_" + std::to_string(iter) + ".geo",std::ios::out);
            write_gmsh_geo<double> gmsh;
            gmsh.write_file(test_bsp, *this, 0.0);
            //gmsh.write_bounding_boxes(test_bsp, *this);
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
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_multiple_inclusion_stuff_only_inside_samepack_heuristic(std::vector<std::unique_ptr<rvegen::rve_shape_input>>const& __input_shapes, _Generator& __random_generator){
    std::cout<<"compute_multiple_inclusion_stuff_samepack_heuristic"<<std::endl;

    value_type volume_fraction{0};
    value_type area{0};
    for(auto & shape : __input_shapes){
        volume_fraction += shape.get()->get_volume_fraction();
        area += shape.get()->min_area();
    }

    if(1.0 - volume_fraction <= 0){
        throw std::runtime_error("rve_generator: volume_fraction >= 1");
    }

    //reserve data for faster push back of new elements
    const size_type max_inclusions{static_cast<size_type>((_box[0]*_box[1]*volume_fraction)/area)};
    _shapes.clear();
    _shapes.reserve(max_inclusions);
    for(auto & shape : __input_shapes){
        if(shape.get()->is_random_position()){
            if(dimension() == 2){
                shape.get()->setup_position(_box[0], _box[1]);
            }else{
                shape.get()->setup_position(_box[0], _box[1], _box[2]);
            }
        }

        volume_fraction = shape.get()->get_volume_fraction();
        size_type iter{0};
        bool finished{false};
        _vol_frac_inclusion = 0;
        int number_of_shapes = shape.get()->get_number_of_shapes();
/*
        //generate_shapes
        for (int i = 0; i <= number_of_shapes;i++){
            auto new_shape = shape.get()->new_shape();
            new_shape.get()->make_bounding_box();
            _generated_shapes.emplace_back(std::move(new_shape));
        }

        //sort_shapes
        _sorted_shapes.reserve(number_of_shapes);
        for(auto& shape : _generated_shapes){
            _sorted_shapes.push_back({shape.get()->area, *shape.get()});
        }

        for (int i = 0; i <= number_of_shapes; i++){
            _sorted_shapes.emplace_back({i, _sorted_shapes(i)});
        }
*/
        if(iter == _max_iter){
            throw std::runtime_error("max iterations reached");
        }
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_inclusion_stuff_random(rve_shape_input* __input_shape, _Generator& __random_generator){
    std::cout<<"compute_inclusion_stuff_random"<<std::endl;

    //    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()};
    //
    //    _Distribution dis_radius(minR, maxR);
    //    _Distribution dis_x(0, _box[0]);
    //    _Distribution dis_y(0, _box[1]);

    if(__input_shape->is_random_position()){
        if(dimension() == 2){
            __input_shape->setup_position(_box[0], _box[1]);
        }else{
            __input_shape->setup_position(_box[0], _box[1], _box[2]);
        }

    }

    size_type iter{0};
    const value_type volume_faction{__input_shape->get_volume_fraction()};
    const value_type area{__input_shape->min_area()};
    const size_type max_inclusions{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_inclusions);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new shape
        auto new_shape = __input_shape->new_shape();
        new_shape.get()->make_bounding_box();

        if(!collision(_shapes, new_shape.get())){
            if (dimension() == 2){
               _vol_frac_inclusion += new_shape->area();
            }else{
               _vol_frac_inclusion += new_shape->volume();
            }
            _shapes.emplace_back(std::move(new_shape));
        }

#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
        std::fstream test_bsp;
        test_bsp.open("test_bsp_" + std::to_string(iter) + ".geo",std::ios::out);
        write_gmsh_geo<double> gmsh;
        gmsh.write_file(test_bsp, *this, 0.0);
        gmsh.write_bounding_boxes(test_bsp, *this);
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
constexpr inline auto rve_generator<_Distribution>::compute_inclusion_stuff_only_inside_samepack_heuristic(rve_shape_input* __input_shape, _Generator& __random_generator){
    std::cout<<"compute_inclusion_stuff_only_inside_samepack_heuristic"<<std::endl;

    //    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()};
    //
    //    _Distribution dis_radius(minR, maxR);
    //    _Distribution dis_x(0, _box[0]);
    //    _Distribution dis_y(0, _box[1]);

    if(__input_shape->is_random_position()){
        if(dimension() == 2){
            __input_shape->setup_position(0, 0);
        }else{
            __input_shape->setup_position(0, 0, 0);
        }

    }

    value_type max_frac = __input_shape->get_volume_fraction();

    const value_type volume_fraction{__input_shape->get_volume_fraction()};
    const value_type volume{__input_shape->min_volume()};
    const value_type area{__input_shape->min_area()};

    _vol_frac_inclusion = 0.0;
    int number_of_shapes = __input_shape->get_number_of_shapes();

    _shapes.clear();
    _generated_shapes.clear();

    //reserve data for faster push back of new elements
    if (dimension()==2){
        const size_type max_inclusions{static_cast<size_type>((_box[0]*_box[1]*volume_fraction)/area)};
        _shapes.reserve(max_inclusions);
    }
    if (dimension()==3){
        const size_type max_inclusions{static_cast<size_type>((_box[0]*_box[1]*_box[2]*volume_fraction)/volume)};
        _shapes.reserve(max_inclusions);
    }

    _generated_shapes.reserve(number_of_shapes);

    generate_shapes(number_of_shapes, dimension(), __input_shape, _generated_shapes);

    sort_shapes(_generated_shapes);

    int AnzahlBereiche = 4;
    bool right_size = check_size(AnzahlBereiche, _generated_shapes);

    set_sections(AnzahlBereiche, number_of_shapes, sections);

    int AnzahlDurchläufeGesamt = 10;
    int FehlversucheFrei = 0;
    int FehlversucheSeitenMax = 10;
    int FehlversucheFreiMax = 1000;
    int DurchläufeGravitationMax = 10;

    if (right_size){
        //Durchgänge gesamt
        for (int i=0;i<AnzahlDurchläufeGesamt;i++){
           //Bereiche
            for (int j=0;j<AnzahlBereiche;j++){
               fill_sides(_vol_frac_inclusion, max_frac, dimension(), j, FehlversucheSeitenMax, sections, _shapes, _generated_shapes);
               while ((FehlversucheFrei <= FehlversucheFreiMax) && (_vol_frac_inclusion < max_frac)){
                     if (!arrange_next(_vol_frac_inclusion, dimension(), j, sections, _shapes, _generated_shapes)){
                        FehlversucheFrei++;
                    }
                     else{
                         FehlversucheFrei = 0;
                     }
                }
               FehlversucheFrei = 0;
            }
            if(_vol_frac_inclusion < max_frac){
                for (int i=0; i<DurchläufeGravitationMax; i++){
                  add_gravity(dimension(),_shapes);
                }
            }

        }

    }else{
        std::cout<< "Bitte durch "<<AnzahlBereiche<<" teilbare Zahl angeben" <<std::endl;
    }


/*
#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_fraction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
        std::fstream test_bsp;
        test_bsp.open("test_bsp_" + std::to_string(iter) + ".geo",std::ios::out);
        write_gmsh_geo<double> gmsh;
        gmsh.write_file(test_bsp, *this);
        gmsh.write_bounding_boxes(test_bsp, *this);
#endif
*/
}



template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_rectangle_random(rectangle_input const& __input, _Generator& __random_generator){

    //    const value_type minR{__input.get_radius_min()}, maxR{__input.get_radius_max()};
    //    const value_type volume_faction{__input.get_volume_fraction()};
    //    _Distribution dis_radius(minR, maxR);
    //    _Distribution dis_x(0, _box[0]);
    //    _Distribution dis_y(0, _box[1]);

    //    size_type iter{0};
    //    const value_type area{minR*minR*M_PI};
    //    const size_type max_circles{static_cast<size_type>((_box[0]*_box[1]*volume_faction)/area)};

    //    bool finished{false};


    //    _shapes.clear();
    //    //reserve data for faster push back of new elements
    //    _shapes.reserve(max_circles);
    //    iter = 0;
    //    _vol_frac_inclusion = 0;
    //    while (iter <= _max_iter) {
    //        //new circle
    //        const double radius{dis_radius(__random_generator)};
    //        const auto x{dis_x(__random_generator)};
    //        const auto y{dis_y(__random_generator)};
    //        const circle<value_type> circle_(x, y, radius);

    //        if(!collision(_shapes, circle_)){
    //            _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_));
    //            _vol_frac_inclusion += circle_.area();
    //        }

    //#ifdef RVE_DEBUG
    //        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
    //#endif

    //        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]))) < 0.005){
    //            finished = true;
    //            break;
    //        }
    //        ++iter;
    //    }


    //    if(iter == _max_iter){
    //        throw std::runtime_error("max iterations reached");
    //    }
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
        circle<value_type> circle_(x, y, radius);

        if(!collision(_shapes, &circle_)){
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
        ellipse<value_type> ellipse_(x, y, radius_a, radius_b, rotation);
        ellipse_.make_bounding_box();

        if(!collision(_shapes, &ellipse_)){
            _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_));
            _vol_frac_inclusion += ellipse_.area();
            _shapes.back().get()->make_bounding_box();

#ifdef RVE_DEBUG
            std::fstream test_bsp;
            test_bsp.open("test_bsp_" + std::to_string(iter) + ".geo",std::ios::out);
            write_gmsh_geo<double> gmsh;
            gmsh.write_file(test_bsp, *this, 0.0);
            gmsh.write_bounding_boxes(test_bsp, *this);
#endif
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
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_lf_random(cylinder_lf_input const& __input, _Generator& __random_generator){
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
        cylinder_lf<value_type> cylinder_(x, y, z, radius, height);

        if(!collision(_shapes, &cylinder_)){
            _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_));
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
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_sf_random(cylinder_sf_input const& __input, _Generator& __random_generator){
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
        cylinder_sf<value_type> cylinder_(x, y, z, radius, height);

        if(!collision(_shapes, &cylinder_)){
            _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_));
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
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipsoid_random(ellipsoid_input const& __input, _Generator& __random_generator){
    const value_type minRa{__input.get_radius_min_a()}, maxRa{__input.get_radius_max_a()}, minRb{__input.get_radius_min_b()}, maxRb{__input.get_radius_max_b()}, minRc{__input.get_radius_min_c()}, maxRc{__input.get_radius_max_c()}, minRotx{__input.get_min_rotation_x()}, maxRotx{__input.get_max_rotation_x()}, minRoty{__input.get_min_rotation_y()}, maxRoty{__input.get_max_rotation_y()}, minRotz{__input.get_min_rotation_z()}, maxRotz{__input.get_max_rotation_z()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius_a(minRa, maxRa);
    _Distribution dis_radius_b(minRb, maxRb);
    _Distribution dis_radius_c(minRc, maxRc);
    _Distribution dis_rotation_x(minRotx, maxRotx);
    _Distribution dis_rotation_y(minRoty, maxRoty);
    _Distribution dis_rotation_z(minRotz, maxRotz);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);
    _Distribution dis_z(0, _box[2]);

    size_type iter{0};
    const value_type volume{minRa*minRb*minRc*M_PI*4/3};
    const size_type max_ellipsoids{static_cast<size_type>((_box[0]*_box[1]*_box[2]*volume_faction)/volume)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_ellipsoids);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const value_type radius_a{dis_radius_a(__random_generator)};
        const value_type radius_b{dis_radius_b(__random_generator)};
        const value_type radius_c{dis_radius_c(__random_generator)};
        const value_type rotation_x{dis_rotation_x(__random_generator)};
        const value_type rotation_y{dis_rotation_y(__random_generator)};
        const value_type rotation_z{dis_rotation_z(__random_generator)};
        const value_type x{dis_x(__random_generator)};
        const value_type y{dis_y(__random_generator)};
        const value_type z{dis_z(__random_generator)};
        ellipsoid<value_type> ellipsoid_(x, y, z, radius_a, radius_b, radius_c, rotation_x, rotation_y, rotation_z);
        ellipsoid_.make_bounding_box();

        if(!collision(_shapes, &ellipsoid_)){
            _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_));
            _vol_frac_inclusion += ellipsoid_.volume();
            _shapes.back().get()->make_bounding_box();

#ifdef RVE_DEBUG
            std::fstream test_bsp;
            test_bsp.open("test_bsp_" + std::to_string(iter) + ".geo",std::ios::out);
            write_gmsh_geo<double> gmsh;
            gmsh.write_file(test_bsp, *this, 0.0);
            gmsh.write_bounding_boxes(test_bsp, *this);
#endif
        }

        //#ifdef RVE_DEBUG
        std::cout<<"iter "<<iter<<" volume fraction "<<_vol_frac_inclusion<<" error "<<(volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1])))<<std::endl;
        //#endif

        if((volume_faction - (_vol_frac_inclusion/(_box[0]*_box[1]*_box[2]))) < 0.005){
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
        circle<value_type> circle_(x, y, radius);

        bool check_distance_{collision(_shapes, &circle_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            if((circle_(0)-circle_.radius()) < 0){
                is_outside = true;
                circle<value_type> circle_periodic(x+dx,y,radius);
                check_distance_periodic = collision(_shapes, &circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            //right side
            if(circle_(0)+circle_.radius()>dx){
                is_outside = true;
                circle<value_type> circle_periodic(x-dx,y,radius);
                check_distance_periodic = collision(_shapes, &circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            //bottom
            if(circle_(1)-circle_.radius()<0){
                is_outside = true;
                circle<value_type> circle_periodic(x,y+dy,radius);
                check_distance_periodic = collision(_shapes, &circle_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<circle<value_type>>(circle_periodic));
                }
            }
            //top
            if(circle_(1)+circle_.radius()>dy){
                is_outside = true;
                circle<value_type> circle_periodic(x,y-dy,radius);
                check_distance_periodic = collision(_shapes, &circle_periodic);
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
        ellipse<value_type> ellipse_(x, y, radius_a, radius_b, rotation);

        bool check_distance_{collision(_shapes, &ellipse_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            if(intersection(ellipse_, 0.0, 0.0, 0.0, 1.0)){
                is_outside = true;
                ellipse<value_type> ellipse_periodic(x+dx,y,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, &ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            //right side
            if(intersection(ellipse_, 1.0, 0.0, 1.0, 1.0)){
                is_outside = true;
                ellipse<value_type> ellipse_periodic(x-dx,y,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, &ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            //bottom
            if(intersection(ellipse_, 0.0, 0.0, 1.0, 0.0)){
                is_outside = true;
                ellipse<value_type> ellipse_periodic(x,y+dy,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, &ellipse_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipse<value_type>>(ellipse_periodic));
                }
            }
            //top
            if(intersection(ellipse_, 0.0, 1.0, 1.0, 1.0)){
                is_outside = true;
                ellipse<value_type> ellipse_periodic(x,y-dy,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, &ellipse_periodic);
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
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_lf_periodic(cylinder_lf_input const& __input, _Generator& __random_generator){
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
        cylinder_lf<value_type> cylinder_(x, y, z, radius, height);

        bool check_distance_{collision(_shapes, &cylinder_)};
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
                cylinder_lf<value_type> cylinder_periodic(x+dx,y,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                }
            }
            //right side
            if(cylinder_(0)+cylinder_.radius()>dx){
                is_outside = true;
                is_outside_right = true;
                cylinder_lf<value_type> cylinder_periodic(x-dx,y,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                }
            }
            //bottom
            if(cylinder_(1)-cylinder_.radius()<0){
                is_outside = true;
                is_outside_bottom = true;
                cylinder_lf<value_type> cylinder_periodic(x,y+dy,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                }
            }
            //top
            if(cylinder_(1)+cylinder_.radius()>dy){
                is_outside = true;
                bool is_outside_top = true;
                cylinder_lf<value_type> cylinder_periodic(x,y-dy,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                }
            }

            //front
            if(cylinder_(2)+cylinder_.height()>dz){
                is_outside = true;
                cylinder_lf<value_type> cylinder_periodic(x,y,0,radius,(height-(dz-z)));
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                }
                //front+left
                if(is_outside_left == true){
                    cylinder_lf<value_type> cylinder_periodic(x+dx,y,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
                //front+right
                if(is_outside_right == true){
                    cylinder_lf<value_type> cylinder_periodic(x-dx,y,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
                //front+bottom
                if(is_outside_bottom == true){
                    cylinder_lf<value_type> cylinder_periodic(x,y+dy,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
                //front+top
                if(is_outside_top == true){
                    cylinder_lf<value_type> cylinder_periodic(x,y-dy,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
            }

            //back
            if(cylinder_(2)+cylinder_.height()<0){
                is_outside = true;
                cylinder_lf<value_type> cylinder_periodic(x,y,_box[2],radius,(height+z));
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                }
                //back+left
                if(is_outside_left == true){
                    cylinder_lf<value_type> cylinder_periodic(x+dx,y,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
                //back+right
                if(is_outside_right == true){
                    cylinder_lf<value_type> cylinder_periodic(x-dx,y,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
                //back+bottom
                if(is_outside_bottom == true){
                    cylinder_lf<value_type> cylinder_periodic(x,y+dy,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
                //back+top
                if(is_outside_top == true){
                    cylinder_lf<value_type> cylinder_periodic(x,y-dy,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_periodic));
                    }
                }
            }

            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_));
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
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_sf_periodic(cylinder_sf_input const& __input, _Generator& __random_generator){
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
        cylinder_sf<value_type> cylinder_(x, y, z, radius, height);

        bool check_distance_{collision(_shapes, &cylinder_)};
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
                cylinder_sf<value_type> cylinder_periodic(x+dx,y,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                }
            }
            //right side
            if(cylinder_(0)+cylinder_.radius()>dx){
                is_outside = true;
                is_outside_right = true;
                cylinder_sf<value_type> cylinder_periodic(x-dx,y,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                }
            }
            //bottom
            if(cylinder_(1)-cylinder_.radius()<0){
                is_outside = true;
                is_outside_bottom = true;
                cylinder_sf<value_type> cylinder_periodic(x,y+dy,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                }
            }
            //top
            if(cylinder_(1)+cylinder_.radius()>dy){
                is_outside = true;
                bool is_outside_top = true;
                cylinder_sf<value_type> cylinder_periodic(x,y-dy,z,radius,height);
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                }
            }

            //front
            if(cylinder_(2)+cylinder_.height()>dz){
                is_outside = true;
                cylinder_sf<value_type> cylinder_periodic(x,y,0,radius,(height-(dz-z)));
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                }
                //front+left
                if(is_outside_left == true){
                    cylinder_sf<value_type> cylinder_periodic(x+dx,y,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
                //front+right
                if(is_outside_right == true){
                    cylinder_sf<value_type> cylinder_periodic(x-dx,y,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
                //front+bottom
                if(is_outside_bottom == true){
                    cylinder_sf<value_type> cylinder_periodic(x,y+dy,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
                //front+top
                if(is_outside_top == true){
                    cylinder_sf<value_type> cylinder_periodic(x,y-dy,0,radius,(height-(dz-z)));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
            }

            //back
            if(cylinder_(2)+cylinder_.height()<0){
                is_outside = true;
                cylinder_sf<value_type> cylinder_periodic(x,y,_box[2],radius,(height+z));
                check_distance_periodic = collision(_shapes, &cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                }
                //back+left
                if(is_outside_left == true){
                    cylinder_sf<value_type> cylinder_periodic(x+dx,y,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
                //back+right
                if(is_outside_right == true){
                    cylinder_sf<value_type> cylinder_periodic(x-dx,y,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
                //back+bottom
                if(is_outside_bottom == true){
                    cylinder_sf<value_type> cylinder_periodic(x,y+dy,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
                //back+top
                if(is_outside_top == true){
                    cylinder_sf<value_type> cylinder_periodic(x,y-dy,0,radius,(height+z));
                    check_distance_periodic = collision(_shapes, &cylinder_periodic);
                    if(check_distance_periodic == false){
                        _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_periodic));
                    }
                }
            }

            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_));
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
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipsoid_periodic(ellipsoid_input const& __input, _Generator& __random_generator){
    const value_type minRa{__input.get_radius_min_a()}, maxRa{__input.get_radius_max_a()}, minRb{__input.get_radius_min_b()}, maxRb{__input.get_radius_max_b()}, minRc{__input.get_radius_min_c()}, maxRc{__input.get_radius_max_c()}, minRotx{__input.get_min_rotation_x()}, maxRotx{__input.get_max_rotation_x()}, minRoty{__input.get_min_rotation_y()}, maxRoty{__input.get_max_rotation_y()}, minRotz{__input.get_min_rotation_z()}, maxRotz{__input.get_max_rotation_z()};
    const value_type volume_faction{__input.get_volume_fraction()};


    _Distribution dis_radius_a(minRa, maxRa);
    _Distribution dis_radius_b(minRb, maxRb);
    _Distribution dis_radius_c(minRc, maxRc);
    _Distribution dis_rotation_x(minRotx, maxRotx);
    _Distribution dis_rotation_y(minRoty, maxRoty);
    _Distribution dis_rotation_z(minRotz, maxRotz);
    _Distribution dis_x(0, _box[0]);
    _Distribution dis_y(0, _box[1]);
    _Distribution dis_z(0, _box[2]);

    const value_type dx{_box[0]}, dy{_box[1]}, dz{_box[2]};

    const value_type volume{minRa*minRb*minRc*M_PI*4/3};
    const value_type max_ellipsoids{static_cast<size_type>(dx*dy*dz*volume_faction)/volume};

    size_type iter{0};

    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_ellipsoids);

    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new circle
        const value_type radius_a{dis_radius_a(__random_generator)};
        const value_type radius_b{dis_radius_b(__random_generator)};
        const value_type radius_c{dis_radius_c(__random_generator)};
        const value_type rotation_x{dis_rotation_x(__random_generator)};
        const value_type rotation_y{dis_rotation_y(__random_generator)};
        const value_type rotation_z{dis_rotation_z(__random_generator)};
        const value_type x{dis_x(__random_generator)};
        const value_type y{dis_y(__random_generator)};
        const value_type z{dis_z(__random_generator)};
        ellipsoid<value_type> ellipsoid_(x, y, z, radius_a, radius_b, radius_c, rotation_x, rotation_y, rotation_z);

        bool check_distance_{collision(_shapes, &ellipsoid_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            /*
          if(intersection(ellipsoid_, 0.0, 0.0, 0.0, 1.0)){
                is_outside = true;
                const ellipsoid<value_type> ellipsoid_periodic(x+dx,y,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, ellipsoid_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_periodic));
                }
            }
            //right side
            if(intersection(ellipsoid_, 1.0, 0.0, 1.0, 1.0)){
                is_outside = true;
                const ellipsoid<value_type> ellipsoid_periodic(x-dx,y,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, ellipsoid_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_periodic));
                }
            }
            //bottom
            if(intersection(ellipsoid_, 0.0, 0.0, 1.0, 0.0)){
                is_outside = true;
                const ellipsoid<value_type> ellipsoid_periodic(x,y+dy,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, ellipsoid_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_periodic));
                }
            }
            //top
            if(intersection(ellipsoid_, 0.0, 1.0, 1.0, 1.0)){
                is_outside = true;
                const ellipsoid<value_type> ellipsoid_periodic(x,y-dy,radius_a, radius_b, rotation);
                check_distance_periodic = collision(_shapes, ellipsoid_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_periodic));
                }
            }

            */
            /*
            //front
             if(cylinder_(2)+cylinder_.height()>dz){
                 is_outside = true;
                 const cylinder<value_type> cylinder_periodic(x,y,0,radius,(height-(dz-z)));
                 check_distance_periodic = collision(_shapes, cylinder_periodic);
                 if(check_distance_periodic == false){
                     _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                 }
                 //front+left
                 if(is_outside_left == true){
                     const cylinder<value_type> cylinder_periodic(x+dx,y,0,radius,(height-(dz-z)));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //front+right
                 if(is_outside_right == true){
                     const cylinder<value_type> cylinder_periodic(x-dx,y,0,radius,(height-(dz-z)));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //front+bottom
                 if(is_outside_bottom == true){
                     const cylinder<value_type> cylinder_periodic(x,y+dy,0,radius,(height-(dz-z)));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //front+top
                 if(is_outside_top == true){
                     const cylinder<value_type> cylinder_periodic(x,y-dy,0,radius,(height-(dz-z)));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
             }

             //back
             if(cylinder_(2)+cylinder_.height()<0){
                 is_outside = true;
                 const cylinder<value_type> cylinder_periodic(x,y,_box[2],radius,(height+z));
                 check_distance_periodic = collision(_shapes, cylinder_periodic);
                 if(check_distance_periodic == false){
                     _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                 }
                 //back+left
                 if(is_outside_left == true){
                     const cylinder<value_type> cylinder_periodic(x+dx,y,0,radius,(height+z));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //back+right
                 if(is_outside_right == true){
                     const cylinder<value_type> cylinder_periodic(x-dx,y,0,radius,(height+z));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //back+bottom
                 if(is_outside_bottom == true){
                     const cylinder<value_type> cylinder_periodic(x,y+dy,0,radius,(height+z));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
                 //back+top
                 if(is_outside_top == true){
                     const cylinder<value_type> cylinder_periodic(x,y-dy,0,radius,(height+z));
                     check_distance_periodic = collision(_shapes, cylinder_periodic);
                     if(check_distance_periodic == false){
                         _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                     }
                 }
             }
*/

            if(check_distance_periodic == false || is_outside == false){
                _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_));
                _vol_frac_inclusion += ellipsoid_.area();
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
        circle<value_type> circle_(x, y, radius);

        if(!collision(_shapes, &circle_)){
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
        ellipse<value_type> ellipse_(x, y, radius_a, radius_b, rotation);

        if(!collision(_shapes, &ellipse_)){
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
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_lf_only_inside(cylinder_lf_input const& __input, _Generator& __random_generator){
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
        cylinder_lf<value_type> cylinder_(x, y, z, radius, height);


        if((!collision(_shapes, &cylinder_)) && (cylinder_.height() + cylinder_(2) < _box[2]) && (cylinder_(2) > 0)){
            _shapes.emplace_back(std::make_unique<cylinder_lf<value_type>>(cylinder_));
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
constexpr inline auto rve_generator<_Distribution>::compute_single_cylinder_sf_only_inside(cylinder_sf_input const& __input, _Generator& __random_generator){
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
        cylinder_sf<value_type> cylinder_(x, y, z, radius, height);


        if((!collision(_shapes, &cylinder_)) && (cylinder_.height() + cylinder_(2) < _box[2]) && (cylinder_(2) > 0)){
            _shapes.emplace_back(std::make_unique<cylinder_sf<value_type>>(cylinder_));
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
constexpr inline auto rve_generator<_Distribution>::compute_single_ellipsoid_only_inside(ellipsoid_input const& __input, _Generator& __random_generator){
    const value_type minRa{__input.get_radius_min_a()}, maxRa{__input.get_radius_max_a()}, minRb{__input.get_radius_min_b()}, maxRb{__input.get_radius_max_b()}, minRc{__input.get_radius_min_c()}, maxRc{__input.get_radius_max_c()}, minRotx{__input.get_min_rotation_x()}, maxRotx{__input.get_max_rotation_x()}, minRoty{__input.get_min_rotation_y()}, maxRoty{__input.get_max_rotation_y()}, minRotz{__input.get_min_rotation_z()}, maxRotz{__input.get_max_rotation_z()};
    const value_type volume_faction{__input.get_volume_fraction()};
    _Distribution dis_radius_a(minRa, maxRa);
    _Distribution dis_radius_b(minRb, maxRb);
    _Distribution dis_radius_c(minRc, maxRc);
    _Distribution dis_rotation_x(minRotx, maxRotx);
    _Distribution dis_rotation_y(minRoty, maxRoty);
    _Distribution dis_rotation_z(minRotz, maxRotz);

    value_type maxR = max(maxRa, maxRb, maxRc);

    _Distribution dis_x(maxR, _box[0]-maxR);
    _Distribution dis_y(maxR, _box[1]-maxR);
    _Distribution dis_z(maxR, _box[1]-maxR);

    size_type iter{0};
    const value_type volume{minRa*minRb*minRc*M_PI*4/3};
    const size_type max_ellipsoids{static_cast<size_type>((_box[0]*_box[1]*_box[2]*volume_faction)/volume)};

    bool finished{false};


    _shapes.clear();
    //reserve data for faster push back of new elements
    _shapes.reserve(max_ellipsoids);
    iter = 0;
    _vol_frac_inclusion = 0;
    while (iter <= _max_iter) {
        //new ellipsoid
        const double radius_a{dis_radius_a(__random_generator)};
        const double radius_b{dis_radius_b(__random_generator)};
        const double radius_c{dis_radius_c(__random_generator)};
        const double rotation_x{dis_rotation_x(__random_generator)};
        const double rotation_y{dis_rotation_y(__random_generator)};
        const double rotation_z{dis_rotation_z(__random_generator)};
        const auto x{dis_x(__random_generator)};
        const auto y{dis_y(__random_generator)};
        const auto z{dis_z(__random_generator)};
        ellipsoid<value_type> ellipsoid_(x, y, z, radius_a, radius_b, radius_c, rotation_x, rotation_y, rotation_z);

        if(!collision(_shapes, &ellipsoid_)){
            _shapes.emplace_back(std::make_unique<ellipsoid<value_type>>(ellipsoid_));
            _vol_frac_inclusion += ellipsoid_.volume();
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

//                 if(!collision(_data,circle_)){
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

//                bool check_distance_{collision(_data, circle_)};
//                if(check_distance_ == false){
//                    //check if outside of rve
//                    bool check_distance_periodic{true};
//                    bool is_outside{false};
//                    //left side
//                    if(circle_(0)-circle_.radius()<0){
//                        is_outside = true;
//                        const circle<> circle_periodic(x+dx,y,radius);
//                        check_distance_periodic = collision(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    //right side
//                    if(circle_(0)+circle_.radius()>dx){
//                        is_outside = true;
//                        const circle<> circle_periodic(x-dx,y,radius);
//                        check_distance_periodic = collision(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    //bottom
//                    if(circle_(1)-circle_.radius()<0){
//                        is_outside = true;
//                        const circle<> circle_periodic(x,y+dy,radius);
//                        check_distance_periodic = collision(_data, circle_periodic);
//                        if(check_distance_periodic == false){
//                            _data.push_back(circle_periodic);
//                        }
//                    }
//                    //top
//                    if(circle_(1)+circle_.radius()>dy){
//                        is_outside = true;
//                        const circle<> circle_periodic(x,y-dy,radius);
//                        check_distance_periodic = collision(_data, circle_periodic);
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

        if(!collision(_shapes, cylinder_)){
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

        bool check_distance_{collision(_shapes, cylinder_)};
        if(check_distance_ == false){
            //check if outside of rve
            bool check_distance_periodic{true};
            bool is_outside{false};
            //left side
            if((cylinder_(0)-cylinder_.radius()) < 0){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x+dx,y,z,radius,height);
                check_distance_periodic = collision(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //right side
            if(cylinder_(0)+cylinder_.radius()>dx){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x-dx,y,z,radius,height);
                check_distance_periodic = collision(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //bottom
            if(cylinder_(1)-cylinder_.radius()<0){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x,y+dy,z,radius,height);
                check_distance_periodic = collision(_shapes, cylinder_periodic);
                if(check_distance_periodic == false){
                    _shapes.emplace_back(std::make_unique<cylinder<value_type>>(cylinder_periodic));
                }
            }
            //top
            if(cylinder_(1)+cylinder_.radius()>dy){
                is_outside = true;
                const cylinder<value_type> cylinder_periodic(x,y-dy,z,radius,height);
                check_distance_periodic = collision(_shapes, cylinder_periodic);
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


        if(!collision(_shapes, cylinder_)){
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

                 if(!collision(_data,circle_)){
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

                bool check_distance_{collision(_data, circle_)};
                if(check_distance_ == false){
                    //check if outside of rve
                    bool check_distance_periodic{true};
                    bool is_outside{false};
                    //left side
                    if(circle_(0)-circle_.radius()<0){
                        is_outside = true;
                        const circle<> circle_periodic(x+dx,y,radius);
                        check_distance_periodic = collision(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    //right side
                    if(circle_(0)+circle_.radius()>dx){
                        is_outside = true;
                        const circle<> circle_periodic(x-dx,y,radius);
                        check_distance_periodic = collision(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    //bottom
                    if(circle_(1)-circle_.radius()<0){
                        is_outside = true;
                        const circle<> circle_periodic(x,y+dy,radius);
                        check_distance_periodic = collision(_data, circle_periodic);
                        if(check_distance_periodic == false){
                            _data.push_back(circle_periodic);
                        }
                    }
                    //top
                    if(circle_(1)+circle_.radius()>dy){
                        is_outside = true;
                        const circle<> circle_periodic(x,y-dy,radius);
                        check_distance_periodic = collision(_data, circle_periodic);
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
