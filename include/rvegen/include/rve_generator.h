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
#include "check_distance.h"
#include "rve_shape_input.h"

namespace rvegen {

enum rveType{
    Periodic,
    NonPeriodic,
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
        _max_iter(500000000),
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
    constexpr inline auto compute_single_circle_nonperiodic(circle_input const& __input, _Generator& __random_generator);
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
    }else if(_rve_type == rveType::NonPeriodic){
        compute_single_circle_nonperiodic(__input, __random_generator);
    }
}

template <typename _Distribution>
template<typename _Generator>
constexpr inline auto rve_generator<_Distribution>::compute_single_circle_nonperiodic(circle_input const& __input, _Generator& __random_generator){

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





////------------------------------------------------------------------------------------------------
////-------------------------------------------cylinder---------------------------------------------
////------------------------------------------------------------------------------------------------
//template<typename T, typename Distribution>
//class rve_generator_3D<cylinder<T>, Distribution>
//{
//public:
//    using value_type = T;
//    using size_type = std::size_t;
//    using Type = circle<value_type>;

//    rve_generator_3D(value_type const volume_faction_inclusion, value_type const x=1, value_type const y=1, value_type const z=1)
//        :_data(),_point(x,y,z),_volume_faction_inclusion(volume_faction_inclusion),_max_iter(5000000),_area(0) {}

//    rve_generator_3D(value_type const volume_faction_inclusion, point<value_type, 3> const& point)
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
//                     std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;
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
//             std::cout<<"Volume fraction: "<<_area/(_point(0)*_point(1))<<std::endl;

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
//        return max_iter;
//    }

//    auto& max_iter(){
//        return max_iter;
//    }

//    template<typename File>
//    void write_gmsh_geo_file(File & file){
//        file<<"SetFactory(\"OpenCASCADE\");"<<std::endl;
//        file<<"Box(1) = {0, 0, 0, 1, 1, 1};"<<std::endl;
//        for(size_t i{0};i<_data.size();++i){
//            file<<"Cylinder("<<i+2<<") = {"<<_data[i](0)<<", "<<_data[i](1)<<",0 , 0, 0, 1 ,"<<_data[i].radius()<<", 2*Pi};"<<std::endl;
//        }
//        for(size_t i{0};i<_data.size();++i){
//            file<<"BooleanIntersection{ Volume{1}; }{ Volume{"<<i+2<<"}; Delete; }"<<std::endl;
//            file<<"BooleanDifference{ Volume{1}; Delete; }{ Volume{"<<i+2<<"}; }"<<std::endl;
//        }

//        //+
//        file<<"Physical Volume(1) = {1};"<<std::endl;
//        file<<"Physical Volume(2) = {";
//        for(size_t i{0};i<_data.size()-1;++i){
//            file<<i+2<<", ";
//        }
//        file<<_data.size()+1<<"};"<<std::endl;
//    }

//private:
//    std::vector<Type> _data;
//    point<value_type, 3> _point;
//    value_type _volume_faction_inclusion;
//    size_type _max_iter;
//    value_type _area;
//};

}
#endif // RVE_GENERATOR_H
