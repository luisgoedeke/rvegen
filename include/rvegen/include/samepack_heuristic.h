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

#include "rve_generator.h"
#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"
#include "ellipsoid.h"
#include "check_distance.h"
#include "rve_shape_input.h"
#include "write_gmsh_geo.h"

namespace rvegen {

template<typename T>
void generate_shapes(int number, rve_shape_input* __input, std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& shapes){
    using value_type = T;
    for (int i=1; i <= number; i++){
       auto new_shape = __input->new_shape();
       new_shape.get()->make_bounding_box();
       shapes.emplace_back(new_shape.get()->area(),std::move(new_shape));
    }
}

template<typename T>
void sort_shapes(std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& shapes){
    using value_type = T;
    auto lambda = [](const auto & __a, const auto & __b){
        return __a.first < __b.first;
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
void arrange_bottom(int section, std::vector<std::pair<int, int>>& sections, std::vector<std::unique_ptr<shape_base<T>>>& _shapes,std::vector<std::pair<T, std::unique_ptr<shape_base<T>>>>& shapes, rve_shape_input* input_shape){
    using value_type = T;

    int position = rand() % (std::get<1>(sections[section])-std::get<0>(sections[section])+1);

    //new coordinates

    value_type x_old = (std::get<1>(shapes[position]))->(0);
    value_type y_old = (std::get<1>(shapes[position]))->(1);

    /*
    value_type max_expansion = (std::get<1>(shapes[position]))->max_expansion();
    std::mt19937 random;
    std::uniform_real_distribution<value_type> dis_x(max_expansion, 1.0-max_expansion);
    std::uniform_real_distribution<value_type> dis_y(max_expansion, 2*max_expansion);
    value_type x = dis_x(random);
    value_type y = dis_y(random);

    std::get<1>(shapes[position])->move(x, y);

    if(!collision(_shapes, std::get<1>(shapes[position]).get())){
    _shapes.emplace_back(std::get<1>(shapes[position]));
    shapes[position].erase();
    }else {
        std::get<1>(shapes[position])->move(x_old, y_old);
    }
*/
    }

template<typename T>
void arrange_right(){
    using value_type = T;
    }

template<typename T>
void arrange_top(){
    using value_type = T;
    }

template<typename T>
void arrange_left(){
    using value_type = T;
    }

template<typename T>
void arrange_front(){
    using value_type = T;
    }

template<typename T>
void arrange_back(){
    using value_type = T;
    }



}



#endif // SAMEPACK_HEURISTIC_H
