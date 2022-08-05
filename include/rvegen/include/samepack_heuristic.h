#ifndef SAMEPACK_HEURISTIC_H
#define SAMEPACK_HEURISTIC_H

#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <memory>
#include <map>
#include <utility>
#include <algorithm>
#include <ctime>
#include <random>

#include "rve_generator.h"
#include "shape_base.h"
#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"
#include "ellipsoid.h"
#include "check_distance.h"
#include "rve_shape_input.h"
#include "write_gmsh_geo.h"

namespace rvegen {

template<typename T>
void generate_shapes(int number, int dimension, rve_shape_input* __input, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& shapes){
    using value_type = T;
    for (int i=1; i <= number; i++){
        auto new_shape = __input->new_shape();
        if (dimension ==2){
            shapes.emplace_back(new_shape.get()->area(),std::move(new_shape));
        }else{
            shapes.emplace_back(new_shape.get()->volume(),std::move(new_shape));
        }
    }
}

template<typename T>
void sort_shapes(std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& shapes){
    using value_type = T;
    auto lambda = [](const auto & __a, const auto & __b){
        return __a.first > __b.first;
    };
    std::sort(shapes.begin(), shapes.end(), lambda);
}

template<typename T>
bool check_size(int a, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& shapes){
    using value_type = T;
    if((shapes.size() % a) > 0){
        std::cout << "Anzahl der Geometrien lässt sich nicht ohne Rest teilen, bitte neue Anzahl vorgeben!" <<std::endl;
        return false;
    }
    return true;
}

template<typename T>
void set_sections(int a, int b, std::vector<std::pair<T, T>>& sections){
    using value_type = T;
    int GroesseBereich = b/a;

    for (int i=0; i<a; i++){
        sections.emplace_back(std::make_pair(i*GroesseBereich, i*GroesseBereich+GroesseBereich-1));
    }
}

template<typename T>
void adjust_sections(T section, std::vector<std::pair<T, T>>& sections){
    using value_type = T;

    std::get<1>(sections[section])--;

    for (int i=section+1; i<sections.size(); i++){
        std::get<0>(sections[i])--;
        std::get<1>(sections[i])--;
    }
}

template<typename T>
void fill_sides(T& _vol_frac_inclusion, T max_frac, int dimension, int j, int FehlversucheSeitenMax, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    int FehlversucheSeiten = 0;

    while ((FehlversucheSeiten <= FehlversucheSeitenMax) && ((std::get<1>(sections[j])-std::get<0>(sections[j]))>0) && (_vol_frac_inclusion < max_frac)){
        if (!arrange_bottom(dimension, j, sections, _shapes, sorted_shapes)){
            FehlversucheSeiten++;
        }
        else{
            if (dimension == 2){
                _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
            }
            if (dimension == 3){
                _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
            }
            FehlversucheSeiten = 0;
        }
    }
    FehlversucheSeiten = 0;
    while ((FehlversucheSeiten <= FehlversucheSeitenMax) && ((std::get<1>(sections[j])-std::get<0>(sections[j]))>0) && (_vol_frac_inclusion < max_frac)){
        if (!arrange_right(dimension, j, sections, _shapes, sorted_shapes)){
            FehlversucheSeiten++;
        }
        else{
            if (dimension == 2){
                _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
            }
            if (dimension == 3){
                _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
            }
            FehlversucheSeiten = 0;
        }
    }
    FehlversucheSeiten = 0;
    while ((FehlversucheSeiten <= FehlversucheSeitenMax) && ((std::get<1>(sections[j])-std::get<0>(sections[j]))>0) && (_vol_frac_inclusion < max_frac)){
        if (!arrange_top(dimension, j, sections, _shapes, sorted_shapes)){
            FehlversucheSeiten++;
        }
        else{
            if (dimension == 2){
                _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
            }
            if (dimension == 3){
                _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
            }
            FehlversucheSeiten = 0;
        }
    }
    FehlversucheSeiten = 0;
    while ((FehlversucheSeiten <= FehlversucheSeitenMax) && ((std::get<1>(sections[j])-std::get<0>(sections[j]))>0) && (_vol_frac_inclusion < max_frac)){
        if (!arrange_left(dimension, j, sections, _shapes, sorted_shapes)){
            FehlversucheSeiten++;
        }
        else{
            if (dimension == 2){
                _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
            }
            if (dimension == 3){
                _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
            }
            FehlversucheSeiten = 0;
        }
    }

    if (dimension == 3){
        FehlversucheSeiten = 0;
        while ((FehlversucheSeiten <= FehlversucheSeitenMax) && ((std::get<1>(sections[j])-std::get<0>(sections[j]))>0) && (_vol_frac_inclusion < max_frac)){
            if (!arrange_front(dimension, j, sections, _shapes, sorted_shapes)){
                FehlversucheSeiten++;
            }
            else{
                if (dimension == 2){
                    _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
                }
                if (dimension == 3){
                    _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
                }
                FehlversucheSeiten = 0;
            }
        }
        FehlversucheSeiten = 0;
        while ((FehlversucheSeiten <= FehlversucheSeitenMax) && ((std::get<1>(sections[j])-std::get<0>(sections[j]))>0)&& (_vol_frac_inclusion < max_frac)){
            if (!arrange_back(dimension, j, sections, _shapes, sorted_shapes)){
                FehlversucheSeiten++;
            }
            else{
                if (dimension == 2){
                    _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
                }
                if (dimension == 3){
                    _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
                }
                FehlversucheSeiten = 0;
            }
        }
    }
}

template<typename T>
bool arrange_bottom(int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    uniform_real_distribution<value_type> distribution;

    srand(time(0));
    int pos = (rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1));
    int position = pos+std::get<0>(sections[section]);

    //new coordinates

    std::array<value_type,3> pos_old = (std::get<1>(sorted_shapes[position]))->get_middle_point();

    std::array<value_type, 3 > min, max;

    std::array<value_type, 3> max_exp = (std::get<1>(sorted_shapes[position]))->max_expansion();

    if (dimension == 2){
        min = {max_exp[0], max_exp[1], 0};
        max = {1-max_exp[0], 3*max_exp[1], 0};
    }

    if (dimension == 3){
        min = {max_exp[0], max_exp[1], max_exp[2]};
        max = {1-max_exp[0], 1-max_exp[1], 3*max_exp[2]};
    }

    std::array<value_type, 3> pos_new;

    for (int i=0; i<3;i++){
        distribution.set_parameter(min[i], max[i]);
        pos_new[i] = distribution();
    }

    bool is_outside = false;
    
    for (int i=0; i<3; i++){
        if(pos_new[i]-max_exp[i] < 0 || pos_new[i]+max_exp[i] > 1){
            is_outside = true;
        }
    }

    std::get<1>(sorted_shapes[position])->set_middle_point(pos_new);

    if((!collision(_shapes, std::get<1>(sorted_shapes[position]).get()))&&(is_outside == false)){
        std::get<1>(sorted_shapes[position])->make_bounding_box();
        _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position])));
        sorted_shapes.erase(sorted_shapes.begin()+position);
        adjust_sections(section, sections);
        return true;
    }
    else{
        std::get<1>(sorted_shapes[position])->set_middle_point(pos_old);
        return false;
    }
}

template<typename T>
bool arrange_right(int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    uniform_real_distribution<value_type> distribution;

    srand(time(0));
    int position = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

    //new coordinates

    auto pos_old = (std::get<1>(sorted_shapes[position]))->get_middle_point();

    std::array<value_type, 3> min, max;

    std::array<value_type, 3> max_exp = (std::get<1>(sorted_shapes[position]))->max_expansion();

    if (dimension == 2){
        min = {1-3*max_exp[0], max_exp[1], 0};
        max = {1-max_exp[0], 1-max_exp[1], 0};
    }

    if (dimension == 3){
        min = {1-3*max_exp[0], max_exp[1], max_exp[2]};
        max = {1-max_exp[0], 1-max_exp[1], 1-max_exp[2]};
    }

    std::array<value_type, 3> pos_new;

    for (int i=0; i<3;i++){
        distribution.set_parameter(min[i], max[i]);
        pos_new[i] = distribution();
    }

    std::get<1>(sorted_shapes[position])->set_middle_point(pos_new);

    if(!collision(_shapes, std::get<1>(sorted_shapes[position]).get())){
        std::get<1>(sorted_shapes[position])->make_bounding_box();
        _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position])));
        sorted_shapes.erase(sorted_shapes.begin()+position);
        adjust_sections(section, sections);
        return true;
    }else{
        std::get<1>(sorted_shapes[position])->set_middle_point(pos_old);
        return false;
    }

}

template<typename T>
bool arrange_top(int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    uniform_real_distribution<value_type> distribution;

    srand(time(0));
    int position = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

    //new coordinates

    auto pos_old = (std::get<1>(sorted_shapes[position]))->get_middle_point();

    std::array<value_type, 3> min, max;

    std::array<value_type, 3> max_exp = (std::get<1>(sorted_shapes[position]))->max_expansion();

    if (dimension == 2){
        min = {max_exp[0], 1-3*max_exp[1], 0};
        max = {1-max_exp[0], 1-max_exp[1], 0};
    }

    if (dimension == 3){
        min = {max_exp[0], max_exp[1], 1-3*max_exp[2]};
        max = {1-max_exp[0], 1-max_exp[1], 1-max_exp[2]};
    }

    std::array<value_type, 3> pos_new;

    for (int i=0; i<3;i++){
        distribution.set_parameter(min[i], max[i]);
        pos_new[i] = distribution();
    }

    bool is_outside = false;

    for (int i=0; i<3; i++){
        if(pos_new[i]-max_exp[i] < 0 || pos_new[i]+max_exp[i] > 1){
            is_outside = true;
        }
    }

    std::get<1>(sorted_shapes[position])->set_middle_point(pos_new);

    if((!collision(_shapes, std::get<1>(sorted_shapes[position]).get()))&&(is_outside == false)){
        std::get<1>(sorted_shapes[position])->make_bounding_box();
        _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position])));
        sorted_shapes.erase(sorted_shapes.begin()+position);
        adjust_sections(section, sections);
        return true;
    }else{
        std::get<1>(sorted_shapes[position])->set_middle_point(pos_old);
        return false;
    }

}

template<typename T>
bool arrange_left(int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    uniform_real_distribution<value_type> distribution;

    srand(time(0));
    int position = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

    //new coordinates

    auto pos_old = (std::get<1>(sorted_shapes[position]))->get_middle_point();

    std::array<value_type, 3> min, max;

    std::array<value_type, 3> max_exp = (std::get<1>(sorted_shapes[position]))->max_expansion();

    if (dimension == 2){
        min = {max_exp[0], max_exp[1], 0};
        max = {3*max_exp[0], 1-max_exp[1], 0};
    }

    if (dimension == 3){
        min = {max_exp[0], max_exp[1], max_exp[2]};
        max = {3*max_exp[0], 1-max_exp[1], 1-max_exp[2]};
    }

    std::array<value_type, 3> pos_new;

    for (int i=0; i<3;i++){
        distribution.set_parameter(min[i], max[i]);
        pos_new[i] = distribution();
    }

    std::get<1>(sorted_shapes[position])->set_middle_point(pos_new);

    if(!collision(_shapes, std::get<1>(sorted_shapes[position]).get())){
        std::get<1>(sorted_shapes[position])->make_bounding_box();
        _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position])));
        sorted_shapes.erase(sorted_shapes.begin()+position);
        adjust_sections(section, sections);
        return true;
    }else{
        std::get<1>(sorted_shapes[position])->set_middle_point(pos_old);
        return false;
    }

}

template<typename T>
bool arrange_front(int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    uniform_real_distribution<value_type> distribution;

    srand(time(0));
    int position = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

    //new coordinates

    auto pos_old = (std::get<1>(sorted_shapes[position]))->get_middle_point();

    std::array<value_type, 3> min, max;

    std::array<value_type, 3> max_exp = (std::get<1>(sorted_shapes[position]))->max_expansion();

    if (dimension == 2){
        min = {max_exp[0], max_exp[1], 0};
        max = {3*max_exp[0], 1-max_exp[1], 0};
    }

    if (dimension == 3){
        min = {max_exp[0], max_exp[1], max_exp[2]};
        max = {1-max_exp[0], 3*max_exp[1], 1-max_exp[2]};
    }

    std::array<value_type, 3> pos_new;

    for (int i=0; i<3;i++){
        distribution.set_parameter(min[i], max[i]);
        pos_new[i] = distribution();
    }

    std::get<1>(sorted_shapes[position])->set_middle_point(pos_new);

    if(!collision(_shapes, std::get<1>(sorted_shapes[position]).get())){
        std::get<1>(sorted_shapes[position])->make_bounding_box();
        _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position])));
        sorted_shapes.erase(sorted_shapes.begin()+position);
        adjust_sections(section, sections);
        return true;
    }else{
        std::get<1>(sorted_shapes[position])->set_middle_point(pos_old);
        return false;
    }
}

template<typename T>
bool arrange_back(int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;
    uniform_real_distribution<value_type> distribution;

    srand(time(0));
    int position = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

    //new coordinates

    auto pos_old = (std::get<1>(sorted_shapes[position]))->get_middle_point();

    std::array<value_type, 3> min, max;

    std::array<value_type, 3> max_exp = (std::get<1>(sorted_shapes[position]))->max_expansion();

    if (dimension == 2){
        min = {max_exp[0], max_exp[1], 0};
        max = {3*max_exp[0], 1-max_exp[1], 0};
    }

    if (dimension == 3){
        min = {max_exp[0], 1-3*max_exp[1], max_exp[2]};
        max = {1-max_exp[0], 1-max_exp[1], 1-max_exp[2]};
    }

    std::array<value_type, 3> pos_new;

    for (int i=0; i<3;i++){
        distribution.set_parameter(min[i], max[i]);
        pos_new[i] = distribution();
    }

    std::get<1>(sorted_shapes[position])->set_middle_point(pos_new);

    if(!collision(_shapes, std::get<1>(sorted_shapes[position]).get())){
        std::get<1>(sorted_shapes[position])->make_bounding_box();
        _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position])));
        sorted_shapes.erase(sorted_shapes.begin()+position);
        adjust_sections(section, sections);
        return true;
    }else{
        std::get<1>(sorted_shapes[position])->set_middle_point(pos_old);
        return false;
    }
}

template<typename T>
bool arrange_next(T& _vol_frac_inclusion, int dimension, int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& sorted_shapes){
    using value_type = T;

    if ((std::get<1>(sections[section])-std::get<0>(sections[section]))>0){
        uniform_real_distribution<value_type> distribution;

        srand(time(0));
        int position_new_shape = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

        //new coordinates

        auto pos_old = (std::get<1>(sorted_shapes[position_new_shape]))->get_middle_point();

        int position_exisiting_shape = rand() % _shapes.size();

        std::array<value_type, 3> max_exp1 = std::get<1>(sorted_shapes[position_new_shape]).get()->max_expansion();
        std::array<value_type, 3> max_exp2 = _shapes[position_exisiting_shape].get()->max_expansion();

        int pos_max_1 = *std::max_element(max_exp1.begin(), max_exp1.begin() + 3);
        int pos_max_2 = *std::max_element(max_exp2.begin(), max_exp2.begin() + 3);
        int pos_min_1 = *std::min_element(max_exp1.begin(), max_exp1.begin() + 3);
        int pos_min_2 = *std::min_element(max_exp2.begin(), max_exp2.begin() + 3);


        value_type max_expansion = std::max(max_exp1[pos_max_1], max_exp2[pos_max_2]);
        value_type min_expansion = std::min(max_exp1[pos_min_1], max_exp2[pos_min_2]);

        distribution.set_parameter(0.0, 1.0);
        value_type phi, theta;
        if (dimension == 2){
            phi = 2*M_PI*distribution();
        }
        if (dimension == 3){
            phi = 2*M_PI*distribution();
            theta = M_PI*distribution();
        }

        distribution.set_parameter(2*min_expansion, 2*max_expansion);
        //         distribution.set_parameter(2*0.05, 2*0.07);
        value_type r = distribution();

        std::array<value_type, 3> point_existing_shape = _shapes[position_exisiting_shape].get()->get_middle_point();
        std::array<value_type, 3> point_new_shape;

        if (dimension == 2){
            point_new_shape[0] = point_existing_shape[0]+cos(phi)*r;
            point_new_shape[1] = point_existing_shape[1]+sin(phi)*r;
            point_new_shape[2] = 0;
        }

        if (dimension == 3){
            point_new_shape[0] = point_existing_shape[0]+r*sin(theta)*cos(phi);
            point_new_shape[1] = point_existing_shape[1]+r*sin(theta)*sin(phi);
            point_new_shape[2] = point_existing_shape[2]+r*cos(theta);
        }

        std::get<1>(sorted_shapes[position_new_shape])->set_middle_point(point_new_shape);
        point_new_shape = std::get<1>(sorted_shapes[position_new_shape])->get_middle_point();

        std::array<value_type, 6>min_max; // 0 = xmin, 1 = xmax, 2 = ymin,...

        if (dimension == 2){
            min_max[0] = point_new_shape[0]-max_exp1[0];
            min_max[1] = point_new_shape[0]+max_exp1[0];
            min_max[2] = point_new_shape[1]-max_exp1[1];
            min_max[3] = point_new_shape[1]+max_exp1[1];
            min_max[4] = 0;
            min_max[5] = 0;
        }

        if (dimension == 3){
            min_max[0] = point_new_shape[0]-max_exp1[0];
            min_max[1] = point_new_shape[0]+max_exp1[0];
            min_max[2] = point_new_shape[1]-max_exp1[1];
            min_max[3] = point_new_shape[1]+max_exp1[1];
            min_max[4] = point_new_shape[2]-max_exp1[2];
            min_max[5] = point_new_shape[2]+max_exp1[2];
        }


        if( !collision(_shapes, std::get<1>(sorted_shapes[position_new_shape]).get()) && !((min_max[0] < 0)||(min_max[1] > 1)||(min_max[2] < 0) || (min_max[3] > 1) || (min_max[4] < 0 ) || (min_max[5] > 1))){
            std::get<1>(sorted_shapes[position_new_shape])->make_bounding_box();
            _shapes.emplace_back(std::move(std::get<1>(sorted_shapes[position_new_shape])));
            sorted_shapes.erase(sorted_shapes.begin()+position_new_shape);
            adjust_sections(section, sections);
            if (dimension == 2){
                _vol_frac_inclusion += (_shapes.back()).get()->area()/1;
            }
            if (dimension == 3){
                _vol_frac_inclusion += (_shapes.back()).get()->volume()/1;
            }
            return true;
        }else{
            std::get<1>(sorted_shapes[position_new_shape])->set_middle_point(pos_old);
            return false;
        }
    }
    else{
        return false;
    }
}

template<typename T>
void add_gravity(int dimension, std::vector<std::unique_ptr<shape_base<T>>>& _shapes){
    using value_type = T;

    int axis;
    std::random_device rd;
    std::mt19937_64 gen(rd());

    if (dimension == 2){
        axis = rd()%2;
    }

    if (dimension == 3){
        axis = rd()%3;
    }

    std::vector<std::pair<value_type, std::unique_ptr<shape_base<T>>>> _shapes_height;

    generate_altitude_sorted(axis, dimension, _shapes, _shapes_height);

    for (int i=0; i < _shapes_height.size(); i++){
        move_geometry(axis, dimension, i, _shapes_height);
    }

    write_shapes_back(_shapes, _shapes_height);

}

template<typename T>
void write_shapes_back(std::vector<std::unique_ptr<shape_base<T>>>& _shapes, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& _shapes_height){
    while (!_shapes_height.empty()){
        _shapes.emplace_back(std::move(std::get<1>(_shapes_height[0])));
        _shapes_height.erase(_shapes_height.begin());
    }
}

template<typename T>
void generate_altitude_sorted(int axis, int dimension, std::vector<std::unique_ptr<shape_base<T>>>& _shapes, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& _shapes_height){
    using value_type = T;

    while (!_shapes.empty()){
        std::array<T,3> coordinates = _shapes[0].get()->get_middle_point();
        value_type height = coordinates[axis];
        _shapes_height.emplace_back(height, std::move(_shapes[0]));
        _shapes.erase(_shapes.begin());
    }
    altitude_sort(_shapes_height);
}

template<typename T>
void altitude_sort( std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& _shapes_height){
    using value_type = T;
    auto lambda = [](const auto & __a, const auto & __b){
        return __a.first < __b.first;
    };
    std::sort(_shapes_height.begin(), _shapes_height.end(), lambda);
}

template<typename T>
void move_geometry(int axis, int dimension, int vector_position, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& _shapes_height){
    using value_type = T;

    std::vector<std::unique_ptr<shape_base<T>>> shapes_with_overlapping_boxes;
    //Position der Koordinate für die Höhe
    //    int position_height;

    if (dimension == 2){
        // y= Höhe
        //        position_height = 1;

        for (int i=0; i < vector_position; i++){
            //Alle Geometrien finden, die sich "unterhalb" der zu bewegenden Geometrie befinden und diese in den neuen Vektor einordnen, ist derzeit auskommentiert da es zu nicht nachvollziehbaren Kollisionen kam
            if(bounding_overlap_check(axis,*static_cast<rectangle_bounding<T>*> (std::get<1>(_shapes_height[i])->bounding_box()), *static_cast<rectangle_bounding<T>*> (std::get<1>(_shapes_height[vector_position])->bounding_box()))){
                shapes_with_overlapping_boxes.emplace_back(std::move(std::get<1>(_shapes_height[i])));
                _shapes_height.erase(_shapes_height.begin()+i);
                i--;
                vector_position--;
            }
        }
    }

    if (dimension == 3){
        //z= Höhe
        //        position_height = 2;

        for (int i=0; i < vector_position; i++){
            //Alle Geometrien finden, die sich "unterhalb" der zu bewegenden Geometrie befinden und diese in den neuen Vektor einordnen, ist derzeit auskommentiert da es zu nicht nachvollziehbaren Kollisionen kam
            if(bounding_overlap_check(axis,*static_cast<box_bounding<T>*> (std::get<1>(_shapes_height[i])->bounding_box()), *static_cast<box_bounding<T>*> (std::get<1>(_shapes_height[vector_position])->bounding_box()))){
                shapes_with_overlapping_boxes.emplace_back(std::move(std::get<1>(_shapes_height[i])));
                _shapes_height.erase(_shapes_height.begin()+i);
                i--;
                vector_position--;
            }
        }
    }

    //Maximale Versuche um das Objekt nach unten zu verschieben
    int max_attempts = 3;
    int attempts = 0;

    while(attempts < max_attempts){

        if(!shapes_with_overlapping_boxes.empty()){
            std::array<value_type,3> expansion_max = std::get<1>(_shapes_height[vector_position])->max_expansion();

            std::array<T,3> current_position = std::get<1>(_shapes_height[vector_position])->get_middle_point();
            std::array<T,3> position_highest_shape = shapes_with_overlapping_boxes[shapes_with_overlapping_boxes.size()-1]->get_middle_point();

            uniform_real_distribution<value_type> distribution;
            distribution.set_parameter(position_highest_shape[axis], current_position[axis]);

            value_type height_old = current_position[axis];
            std::array<T,3> bot_p = std::get<1>(_shapes_height[vector_position])->bounding_box()->bottom_point();
            std::array<T,3> top_p = shapes_with_overlapping_boxes[shapes_with_overlapping_boxes.size()-1]->bounding_box()->top_point();
            value_type distance_bounding = bot_p[axis]-top_p[axis];
            value_type height_new;

            if (distance_bounding > 0){
                height_new = current_position[axis]-distance_bounding;

            }else{
                height_new = distribution();
            }


                current_position[axis] = height_new;
/*
                for (int i=0; i < vector_position; i++){
                    //Alle Geometrien finden, die sich "unterhalb" der zu bewegenden Geometrie befinden und diese in den neuen Vektor einordnen, ist derzeit auskommentiert da es zu nicht nachvollziehbaren Kollisionen kam
//                    if(bounding_overlap_check(axis,*static_cast<rectangle_bounding<T>*> (std::get<1>(_shapes_height[i])->bounding_box()), *static_cast<rectangle_bounding<T>*> (std::get<1>(_shapes_height[vector_position])->bounding_box()))){
                        shapes_with_overlapping_boxes.emplace_back(std::move(std::get<1>(_shapes_height[i])));
                        _shapes_height.erase(_shapes_height.begin()+i);
                        i--;
                        vector_position--;
                    }

                for (int i=_shapes_height.size()-1; i > vector_position; i--){
                //Alle Geometrien finden, die sich "unterhalb" der zu bewegenden Geometrie befinden und diese in den neuen Vektor einordnen, ist derzeit auskommentiert da es zu nicht nachvollziehbaren Kollisionen kam
//                    if(bounding_overlap_check(axis,*static_cast<rectangle_bounding<T>*> (std::get<1>(_shapes_height[i])->bounding_box()), *static_cast<rectangle_bounding<T>*> (std::get<1>(_shapes_height[vector_position])->bounding_box()))){
                    shapes_with_overlapping_boxes.emplace_back(std::move(std::get<1>(_shapes_height[i])));
                    _shapes_height.erase(_shapes_height.begin()+i);

                }
*/
                //OnlyInside
                if((current_position[axis]-expansion_max[axis])>0 && (current_position[axis]-expansion_max[axis])<1)  {
                  std::get<1>(_shapes_height[vector_position])->set_middle_point(current_position);
                  std::get<1>(_shapes_height[vector_position])->make_bounding_box();
                  if(collision(shapes_with_overlapping_boxes, std::get<1>(_shapes_height[vector_position]).get())){
                      current_position[axis] = height_old;
                      std::get<1>(_shapes_height[vector_position])->set_middle_point(current_position);
                      std::get<1>(_shapes_height[vector_position])->make_bounding_box();
                  }
                  attempts++;
                }
                attempts++;
            }
            else{
                    attempts++;
                }
            }else{
                //Behandlung für unterste Reihe
                std::array<T,3> current_position = std::get<1>(_shapes_height[vector_position])->get_middle_point();
                value_type height_old = current_position[axis];
                std::array<value_type,3> expansion_max = std::get<1>(_shapes_height[vector_position])->max_expansion();
                uniform_real_distribution<value_type> distribution;
                distribution.set_parameter(1.1*expansion_max[axis], 1.2*expansion_max[axis]);
                value_type height_new = distribution();
                current_position[axis] = std::min(height_old, height_new);
                std::get<1>(_shapes_height[vector_position])->set_middle_point(current_position);
                std::get<1>(_shapes_height[vector_position])->make_bounding_box();
                attempts++;
            }
        }else{
            //Behandlung für unterste Reihe
            std::array<T,3> current_position = std::get<1>(_shapes_height[vector_position])->get_middle_point();
            value_type height_old = current_position[axis];
            std::array<value_type,3> expansion_max = std::get<1>(_shapes_height[vector_position])->max_expansion();
            uniform_real_distribution<value_type> distribution;
            distribution.set_parameter(1.1*expansion_max[axis], 1.2*expansion_max[axis]);
            value_type height_new = distribution();
            current_position[axis] = std::min(height_old, height_new);
            std::get<1>(_shapes_height[vector_position])->set_middle_point(current_position);
            attempts++;
        }
    }
    while(!shapes_with_overlapping_boxes.empty()){
        std::array<T,3> middle_point = (shapes_with_overlapping_boxes[0]).get()->get_middle_point();
        value_type middle_point_height = middle_point[axis];
        _shapes_height.emplace_back(middle_point_height, std::move(shapes_with_overlapping_boxes[0]));
        shapes_with_overlapping_boxes.erase(shapes_with_overlapping_boxes.begin());
    }
    altitude_sort(_shapes_height);
}

}

#endif // SAMEPACK_HEURISTIC_H
