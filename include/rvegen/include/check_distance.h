#ifndef CHECK_DISTANCE_H
#define CHECK_DISTANCE_H

#include <vector>
#include <iostream>
#include <array>
#include <algorithm>
#include <numeric>
#include "../../../external/Eigen/Eigen"

#include "circle.h"
#include "cylinder.h"
#include "ellipse.h"
#include "ellipsoid.h"
#include "sphere.h"
#include "rectangle.h"
#include "box.h"
#include "rectangle_bounding.h"

namespace rvegen {

//circle -- rectangle
//https://stackoverflow.com/questions/401847/circle-rectangle-collision-detection-intersection


//bounding_box_collision
//rectangle -- rectangle    ok
//box -- box                ok


//collision
//2D
//rectangle -- rectangle
//circle -- rectangle
//ellipse -- rectangle
//circle -- circle          ok
//circle -- ellipse
//ellipse -- ellipse        ok
//polygon -- polygon        ok

//3D
//box -- box
//cylinder -- box
//ellipsoid -- box
//sphere -- box
//cylinder -- cylinder      ok
//ellipsoid -- cylinder
//sphere -- cylinder
//ellipsoid -- ellipsoid    ok
//sphere -- ellipsoid
//sphere -- sphere          ok


template<typename T>
bool collision_details(box_bounding<T> const& lhs, box_bounding<T> const& rhs);

template<typename T>
bool collision_details(rectangle_bounding<T> const& lhs, rectangle_bounding<T> const& rhs);

template<typename T>
bool collision_details(ellipse<T> const& lhs, ellipse<T> const& rhs);

template<typename T>
auto collision_details(circle<T> const& lhs, circle<T> const& rhs);

template<typename T>
auto collision_details(rectangle<T> const& lhs, rectangle<T> const& rhs);

template<typename T>
bool collision_details(ellipsoid<T> const& lhs, ellipsoid<T> const& rhs);

template<typename T>
auto collision_details(sphere<T> const& lhs, sphere<T> const& rhs);

template<typename T>
auto collision_details(cylinder<T> const& lhs, cylinder<T> const& rhs);

template<typename T>
auto collision_details(rectangle<T> const& lhs, rectangle<T> const& rhs);

template<typename T>
auto collision_details(rectangle<T> const& lhs, circle<T> const& rhs);

template<typename T>
auto collision_details(circle<T> const& lhs, rectangle<T> const& rhs);

template<typename T>
bool collision_details(ellipse<T> const& lhs, circle<T> const& rhs);

template<typename T>
bool collision_details(circle<T> const& lhs, ellipse<T> const& rhs);

template<typename T>
bool collision_details(ellipse<T> const& lhs, rectangle<T> const& rhs);

template<typename T>
bool collision_details(rectangle<T> const& lhs, ellipse<T> const& rhs);

// Returns true if the point p lies inside the polygon
template<typename _Array, typename _Point>
bool isInside(_Array const& polygon, _Point const& p){
    auto size{polygon.size()};

    // There must be at least 3 vertices in polygon[]
    if (size < 3) return false;

    // Given three collinear points p, q, r, the function checks if
    // point q lies on line segment 'pr'
    auto onSegment = [](auto const& p, auto const& q, auto const& r)
    {
        if (q[0] <= std::max(p[0], r[0]) && q[0] >= std::min(p[0], r[0]) &&
                q[1] <= std::max(p[1], r[1]) && q[1] >= std::min(p[1], r[1]))
            return true;
        return false;
    };

    // To find orientation of ordered triplet (p, q, r).
    // The function returns following values
    // 0 --> p, q and r are collinear
    // 1 --> Clockwise
    // 2 --> Counterclockwise
    auto orientation = [](auto const& p, auto const& q, auto const& r)
    {
        int val = (q[1] - p[1]) * (r[0] - q[0]) -
                (q[0] - p[0]) * (r[1] - q[1]);

        if (std::abs(val) <= 5e-8) return 0; // collinear
        return (val > 0)? 1: 2; // clock or counterclock wise
    };

    auto doIntersect = [&](auto const& p1, auto const& q1, auto const& p2, auto const& q2)
    {
        // Find the four orientations needed for general and
        // special cases
        int o1 = orientation(p1, q1, p2);
        int o2 = orientation(p1, q1, q2);
        int o3 = orientation(p2, q2, p1);
        int o4 = orientation(p2, q2, q1);

        // General case
        if (o1 != o2 && o3 != o4)
            return true;

        // Special Cases
        // p1, q1 and p2 are collinear and p2 lies on segment p1q1
        if (o1 == 0 && onSegment(p1, p2, q1)) return true;

        // p1, q1 and p2 are collinear and q2 lies on segment p1q1
        if (o2 == 0 && onSegment(p1, q2, q1)) return true;

        // p2, q2 and p1 are collinear and p1 lies on segment p2q2
        if (o3 == 0 && onSegment(p2, p1, q2)) return true;

        // p2, q2 and q1 are collinear and q1 lies on segment p2q2
        if (o4 == 0 && onSegment(p2, q1, q2)) return true;

        return false; // Doesn't fall in any of the above cases
    };

    // Create a point for line segment from p to infinite
    _Point extreme = {std::numeric_limits<double>::max(), p[1]};

    // Count intersections of the above line with sides of polygon
    //int count = 0, i = 0;
    //    std::size_t count{0};
    //    for(std::size_t i{0}/*, j{polygon.size()-1}*/; i < polygon.size(); /*j =*/ i++){
    //        int j = (i+1)%polygon.size();
    //        std::cout<<i<<" "<<j<<std::endl;
    //        if (doIntersect(polygon[i], polygon[j], p, extreme)){
    //            if (orientation(polygon[i], p, polygon[j]) == 0){
    //                return onSegment(polygon[i], p, polygon[j]);
    //            }
    //            count++;
    //        }
    //    }
    // Count intersections of the above line with sides of polygon
    int count = 0, i = 0;
    do
    {
        int next = (i+1)%polygon.size();

        // Check if the line segment from 'p' to 'extreme' intersects
        // with the line segment from 'polygon[i]' to 'polygon[next]'
        if (doIntersect(polygon[i], polygon[next], p, extreme))
        {
            // If the point 'p' is collinear with line segment 'i-next',
            // then check if it lies on segment. If it lies, return true,
            // otherwise false
            if (orientation(polygon[i], p, polygon[next]) == 0)
                return onSegment(polygon[i], p, polygon[next]);

            count++;
        }
        i = next;
    } while (i != 0);

    // Return true if count is odd, false otherwise
    return count&1; // Same as (count%2 == 1)
    // Return true if count is odd, false otherwise
    return count&1; // Same as (count%2 == 1)
}


//template<typename _Array, typename _Point>
//bool collision_details(_Array const& _points, _Point const& _point){
//    //https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon/2922778#2922778
//    bool c{false};
//    for(std::size_t i{0}, j{_points.size()-1}; i < _points.size(); j = i++){
//        if(((_points[i][1]>_point[1]) != (_points[j][1]>_point[1])) &&
//                (_point[0] < (_points[j][0] - _points[i][0]) * (_point[1] - _points[i][1]) / (_points[j][1] - _points[i][1]) + _points[i][0]))
//            c = !c;
//    }
//    return c;
//}


// Linear transform to find the orthogonal vector of the edge
template<typename Vector>
Vector calculate_normalised_projection_axis(const Vector &current_point, const Vector &next_point) {
    const double axis_x = -(next_point[1] - current_point[1]);
    const double axis_y =   next_point[0] - current_point[0];
    const double magnitude = std::hypot(axis_x, axis_y);

    Vector axis_normalised;
    axis_normalised[0] = axis_x / magnitude;
    axis_normalised[1] = axis_y / magnitude;

    return axis_normalised;
}

// Project the vertices of each polygon onto a axis
template<typename Polygon, typename Vector>
void compute_projections(const Polygon &bounds_a, const Polygon &bounds_b, const Vector &axis_normalised, std::vector<double> &projections_a, std::vector<double> &projections_b) {
    projections_a.clear();
    projections_b.clear();

    for (size_t i = 0; i < bounds_a.size(); i++) {
        const double projection_a = axis_normalised.dot(bounds_a[i]);
        const double projection_b = axis_normalised.dot(bounds_b[i]);
        projections_a.push_back(projection_a);
        projections_b.push_back(projection_b);
    }
}

// Check if the projections of two polygons overlap
bool is_overlapping(const std::vector<double> &projections_a, const std::vector<double> &projections_b) {
    const double max_projection_a = *std::max_element(projections_a.begin(), projections_a.end());
    const double min_projection_a = *std::min_element(projections_a.begin(), projections_a.end());
    const double max_projection_b = *std::max_element(projections_b.begin(), projections_b.end());
    const double min_projection_b = *std::min_element(projections_b.begin(), projections_b.end());

    // True if projection overlaps but does not necessarily mean the polygons are intersecting yet
    return !(max_projection_a < min_projection_b || max_projection_b < min_projection_a);
}

// Check if two convex polygons intersect
template<typename Polygon>
bool separating_axis_intersect(const Polygon &bounds_a, const Polygon &bounds_b){
    //https://github.com/winstxnhdw/2d-separating-axis-theorem
    //works
    using Vector = typename Polygon::value_type;

    std::vector<double> projections_a;
    std::vector<double> projections_b;
    projections_a.reserve(bounds_a.size());
    projections_b.reserve(bounds_b.size());

    for (size_t i = 0; i < bounds_a.size(); i++) {
        const Vector current_point = bounds_a[i];
        const Vector next_point = bounds_a[(i + 1) % bounds_a.size()];
        const Vector axis_normalised = calculate_normalised_projection_axis(current_point, next_point);
        compute_projections(bounds_a, bounds_b, axis_normalised, projections_a, projections_b);

        if (!is_overlapping(projections_a, projections_b)) return false;
    }

    for (size_t i = 0; i < bounds_b.size(); i++) {
        const Vector current_point = bounds_b[i];
        const Vector next_point = bounds_b[(i + 1) % bounds_b.size()];
        const Vector axis_normalised = calculate_normalised_projection_axis(current_point, next_point);
        compute_projections(bounds_a, bounds_b, axis_normalised, projections_a, projections_b);

        if (!is_overlapping(projections_a, projections_b)) return false;
    }
    // Intersects if all projections overlap
    return true;
}


template<typename T>
bool collision_details(ellipse<T> const& lhs, circle<T> const& rhs){
    return collision_details(rhs, lhs);
}

template<typename T>
bool collision_details(circle<T> const& lhs, ellipse<T> const& rhs){
    return false;
}

template<typename T>
bool collision_details(ellipse<T> const& lhs, rectangle<T> const& rhs){
    return collision_details(rhs, lhs);
}

template<typename T>
bool collision_details(rectangle<T> const& lhs, ellipse<T> const& rhs){
    return false;
}


template<typename T>
auto collision_details(circle<T> const& lhs, rectangle<T> const& rhs){
    return collision_details(rhs, lhs);
}

template<typename T>
auto collision_details(rectangle<T> const& lhs, circle<T> const& rhs){
    using value_type = T;
    //https://gist.github.com/snorpey/8134c248296649433de2

    const auto lhs_w{lhs.width()*0.5}, lhs_h{lhs.height()*0.5};
    const auto rectReferenceX = lhs(0) - lhs_w;
    const auto rectReferenceY = lhs(1) - lhs_h;

    // Rotate circle's center point back
    const auto unrotatedCircleX = std::cos(-lhs.rotation()*2*M_PI) * ( rhs(0) - lhs(0) ) - std::sin(-lhs.rotation()*2*M_PI) * ( rhs(1) - lhs(1) ) + lhs(0);
    const auto unrotatedCircleY = std::sin(-lhs.rotation()*2*M_PI) * ( rhs(0) - lhs(0) ) + std::cos(-lhs.rotation()*2*M_PI) * ( rhs(1) - lhs(1) ) + lhs(1);

    // Closest point in the rectangle to the center of circle rotated backwards(unrotated)
    value_type closestX, closestY;
    // Find the unrotated closest x point from center of unrotated circle
    if(unrotatedCircleX < rectReferenceX){
        closestX = rectReferenceX;
    }else if(unrotatedCircleX > rectReferenceX + lhs.width()){
        closestX = rectReferenceX + lhs.width();
    }else{
        closestX = unrotatedCircleX;
    }

    // Find the unrotated closest y point from center of unrotated circle
    if(unrotatedCircleY < rectReferenceY){
        closestY = rectReferenceY;
    }else if( unrotatedCircleY > rectReferenceY + lhs.height()){
        closestY = rectReferenceY + lhs.height();
    }else{
        closestY = unrotatedCircleY;
    }

    // Determine collision
    const auto dx = unrotatedCircleX - closestX;
    const auto dy = unrotatedCircleY - closestY;
    return ((dx*dx+dy*dy) < rhs.radius()*rhs.radius());
}


//
template<typename T>
auto collision_details(rectangle<T> const& lhs, rectangle<T> const& rhs){
    using value_type = T;
    using Matrix     = Eigen::Matrix<value_type,2,2>;
    using Vector     = Eigen::Vector<value_type,2>;

    const Matrix rhs_rot = Eigen::Rotation2D<value_type>(rhs.rotation()*2*M_PI).toRotationMatrix();
    const Matrix lhs_rot = Eigen::Rotation2D<value_type>(lhs.rotation()*2*M_PI).toRotationMatrix();
    const Vector x0_rhs{rhs.point().begin()};
    const Vector x0_lhs{lhs.point().begin()};
    const value_type lhs_w{lhs.width()*0.5}, lhs_h{lhs.height()*0.5};
    const value_type rhs_w{rhs.width()*0.5}, rhs_h{rhs.height()*0.5};

    std::array<Vector, 4> lhs_points, rhs_points;
    lhs_points[0] = lhs_rot*Vector{ lhs_w,  lhs_h} + x0_lhs;
    lhs_points[1] = lhs_rot*Vector{-lhs_w,  lhs_h} + x0_lhs;
    lhs_points[2] = lhs_rot*Vector{-lhs_w, -lhs_h} + x0_lhs;
    lhs_points[3] = lhs_rot*Vector{ lhs_w, -lhs_h} + x0_lhs;

    rhs_points[0] = rhs_rot*Vector{ rhs_w,  rhs_h} + x0_rhs;
    rhs_points[1] = rhs_rot*Vector{-rhs_w,  rhs_h} + x0_rhs;
    rhs_points[2] = rhs_rot*Vector{-rhs_w, -rhs_h} + x0_rhs;
    rhs_points[3] = rhs_rot*Vector{ rhs_w, -rhs_h} + x0_rhs;

    return separating_axis_intersect(lhs_points, rhs_points);
}

template<typename T>
bool collision_details(ellipse<T> const& lhs, ellipse<T> const& rhs){
    using value_type = T;
    using Matrix     = Eigen::Matrix<value_type,2,2>;
    using Vector     = Eigen::Vector<value_type,2>;

    constexpr std::size_t max_inc_theta = 360;

    //left hand side ellipe fixed
    //right hand side ellipe is discretised

    const value_type a_lhs{lhs.radius_a()}, b_lhs{lhs.radius_b()};
    const value_type x_lhs{lhs(0)}, y_lhs{lhs(1)};
    const value_type alpha_lhs{lhs.rotation()};
    const Vector x0_lhs{x_lhs, y_lhs};

    const value_type x_rhs{rhs(0)}, y_rhs{rhs(1)};
    const value_type a_rhs{rhs.radius_a()}, b_rhs{rhs.radius_b()};
    const value_type alpha_rhs{rhs.rotation()};
    const Vector x0_rhs{x_rhs, y_rhs};

    const Matrix rhs_rot = Eigen::Rotation2D<value_type>(2*alpha_rhs*M_PI).toRotationMatrix();
    const Matrix lhs_rot = Eigen::Rotation2D<value_type>(2*alpha_lhs*M_PI).toRotationMatrix();

    {
        //Center point
        const Vector pos = lhs_rot.transpose()*(x0_rhs - x0_lhs);
        if((pos(0)*pos(0)/(a_lhs*a_lhs) + pos(1)*pos(1)/(b_lhs*b_lhs)) <= 1){
            return true;
        }
    }

    //theta
    for(std::size_t j{0}; j<max_inc_theta; ++j){
        const value_type theta{j/value_type(max_inc_theta) * 360};
        const value_type s_theta{std::sin(theta*M_PI/180.)};
        const value_type c_theta{std::cos(theta*M_PI/180.)};
        //local rhs to global point
        const Vector pos_rhs = rhs_rot*Vector{a_rhs*c_theta, b_rhs*s_theta} + x0_rhs;
        //global to local lhs point
        const Vector pos = lhs_rot.transpose()*(pos_rhs - x0_lhs);
        if((pos(0)*pos(0)/(a_lhs*a_lhs) + pos(1)*pos(1)/(b_lhs*b_lhs)) <= 1.0){
            return true;
        }
    }
    return false;
}


template<typename T>
bool collision_details(box_bounding<T> const& lhs, box_bounding<T> const& rhs){
    const auto& lhs_top{lhs.top_point()};
    const auto& rhs_top{rhs.top_point()};
    const auto& lhs_bottom{lhs.bottom_point()};
    const auto& rhs_bottom{rhs.bottom_point()};

    return ((lhs_bottom[0] <= rhs_top[0] && lhs_top[0] >= rhs_bottom[0]) &&
            (lhs_bottom[1] <= rhs_top[1] && lhs_top[1] >= rhs_bottom[1]) &&
            (lhs_bottom[2] <= rhs_top[2] && lhs_top[2] >= rhs_bottom[2]));
}

template<typename T>
bool collision_details(rectangle_bounding<T> const& lhs, rectangle_bounding<T> const& rhs){
    const auto& lhs_top{lhs.top_point()};
    const auto& rhs_top{rhs.top_point()};
    const auto& lhs_bottom{lhs.bottom_point()};
    const auto& rhs_bottom{rhs.bottom_point()};

    return ((lhs_bottom[0] <= rhs_top[0] && lhs_top[0] >= rhs_bottom[0]) &&
            (lhs_bottom[1] <= rhs_top[1] && lhs_top[1] >= rhs_bottom[1]));
}

template<typename T>
bool collision_details(ellipsoid<T> const& lhs, ellipsoid<T> const& rhs){
    using value_type = T;
    using Matrix33   = Eigen::Matrix<value_type,3,3>;
    using Vector3    = Eigen::Vector3<value_type>;

    constexpr std::size_t max_inc_phi = 720;
    constexpr std::size_t max_inc_theta = 360;

    //left hand side ellipsoid is fixed
    //right hand side ellipsoid is discretised

    //a + interface := a + a*0.1
    const value_type a_lhs{lhs.radius_a() + lhs.radius_a()*0.1}, b_lhs{lhs.radius_b() + lhs.radius_b()*0.1}, c_lhs{lhs.radius_c() + lhs.radius_c()*0.1};
    const value_type x_lhs{lhs(0)}, y_lhs{lhs(1)}, z_lhs{lhs(2)};
    const value_type alpha_lhs{lhs.rotation_x()}, beta_lhs{lhs.rotation_y()}, gamma_lhs{lhs.rotation_z()};
    const Vector3 x0_lhs{x_lhs, y_lhs, z_lhs};

    const value_type x_rhs{rhs(0)}, y_rhs{rhs(1)}, z_rhs{rhs(2)};
    const value_type a_rhs{rhs.radius_a()}, b_rhs{rhs.radius_b()}, c_rhs{rhs.radius_c()};
    const value_type alpha_rhs{rhs.rotation_x()}, beta_rhs{rhs.rotation_y()}, gamma_rhs{rhs.rotation_z()};
    const Vector3 x0_rhs{x_rhs, y_rhs, z_rhs};

    Matrix33 rot;
    rot = Eigen::AngleAxisd(gamma_rhs, Vector3::UnitZ())
            * Eigen::AngleAxisd(beta_rhs, Vector3::UnitY())
            * Eigen::AngleAxisd(alpha_rhs, Vector3::UnitX());

    Matrix33 rot_lhs;
    rot_lhs = Eigen::AngleAxisd(gamma_lhs, Vector3::UnitZ())
            * Eigen::AngleAxisd(beta_lhs, Vector3::UnitY())
            * Eigen::AngleAxisd(alpha_lhs, Vector3::UnitX());

    {
        //Center point
        const Eigen::Vector3d pos = rot_lhs.transpose()*(x0_rhs - x0_lhs);
        if((pos(0)*pos(0)/(a_lhs*a_lhs) + pos(1)*pos(1)/(b_lhs*b_lhs) + pos(2)*pos(2)/(c_lhs*c_lhs)) <= 1){
            return true;
        }
    }


    //phi
    for(std::size_t i{0}; i<max_inc_phi; ++i){
        const value_type phi{i/value_type(max_inc_phi) * 360.};
        const value_type s_phi{std::sin(phi*M_PI/180.)};
        const value_type c_phi{std::cos(phi*M_PI/180.)};
        //theta
        for(std::size_t j{0}; j<max_inc_theta; ++j){
            const value_type theta{j/value_type(max_inc_theta) * 180};
            const value_type s_theta{std::sin(theta*M_PI/180.)};
            const value_type c_theta{std::cos(theta*M_PI/180.)};
            //local rhs to global point
            const Vector3 pos_rhs = rot*Vector3{a_rhs*s_theta*c_phi, c_rhs*s_theta*s_phi, b_rhs*c_theta} + x0_rhs;
            //global to local lhs point
            const Vector3 pos = rot_lhs.transpose()*(pos_rhs - x0_lhs);
            if((pos(0)*pos(0)/(a_lhs*a_lhs) + pos(1)*pos(1)/(b_lhs*b_lhs) + pos(2)*pos(2)/(c_lhs*c_lhs)) <= 1.0){
                return true;
            }
        }
    }
    return false;
}


template<typename T>
auto collision_details(sphere<T> const& lhs, sphere<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const auto dx = lhs(0) - rhs(0);
    const auto dy = lhs(1) - rhs(1);
    const auto dz = lhs(2) - rhs(2);
    return a * a >= (dx * dx + dy * dy + dz * dz);
}

template<typename T>
auto collision_details(circle<T> const& lhs, circle<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const auto dx = lhs(0) - rhs(0);
    const auto dy = lhs(1) - rhs(1);
    return a * a >= (dx * dx + dy * dy);
}

template<typename T>
auto collision_details(cylinder<T> const& lhs, cylinder<T> const& rhs){
    const auto a{lhs.radius() + rhs.radius()};
    const auto dx = lhs(0) - rhs(0);
    const auto dy = lhs(1) - rhs(1);
    const auto start_lhs = lhs(2);
    const auto end_lhs = lhs(2)+lhs.height();
    const auto start_rhs = rhs(2);
    const auto end_rhs = rhs(2)+rhs.height();
    if (a * a > (dx * dx + dy * dy)){
        return !((end_lhs < start_rhs)||start_lhs > end_rhs);
    }
    return false;
}

template<typename T>
auto collision(ellipse<T> const& lhs, ellipse<T> const& rhs){
    if (collision_details(*static_cast<rectangle_bounding<T>*>(lhs.bounding_box()), *static_cast<rectangle_bounding<T>*>(rhs.bounding_box()))){
        return collision_details(lhs, rhs);
    }else{
        return false;
    };
}

template<typename T>
auto collision(ellipsoid<T> const& lhs, ellipsoid<T> const& rhs){
    if(collision_details(*static_cast<box_bounding<T>*>(lhs.bounding_box()), *static_cast<box_bounding<T>*>(rhs.bounding_box()))){
        return collision_details(lhs, rhs);
    }else{
        return false;
    };
}

template<template<class> class _PTR, typename _T>
bool collision(std::vector<_PTR<shape_base<_T>>> const& __shapes, shape_base<_T> * __shape){
    if(dynamic_cast<rectangle<_T>*>(__shape)){
        for(const auto& shape : __shapes){
            if(dynamic_cast<rectangle<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<rectangle<_T>*>(__shape), *static_cast<rectangle<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else if(dynamic_cast<circle<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<rectangle<_T>*>(__shape), *static_cast<circle<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else if(dynamic_cast<ellipse<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<rectangle<_T>*>(__shape), *static_cast<ellipse<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else{
                throw std::runtime_error("collision(): no matching shape type");
            }
        }
    }else if(dynamic_cast<circle<_T>*>(__shape)){
        for(const auto& shape : __shapes){
            if(dynamic_cast<rectangle<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<circle<_T>*>(__shape), *static_cast<rectangle<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else if(dynamic_cast<circle<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<circle<_T>*>(__shape), *static_cast<circle<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else if(dynamic_cast<ellipse<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<circle<_T>*>(__shape), *static_cast<ellipse<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else{
                throw std::runtime_error("collision(): no matching shape type");
            }
        }
    }else if(dynamic_cast<ellipse<_T>*>(__shape)){
        for(const auto& shape : __shapes){
            if(dynamic_cast<rectangle<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<ellipse<_T>*>(__shape), *static_cast<rectangle<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else if(dynamic_cast<circle<_T>*>(shape.get())){
                //check distance
                if(collision_details(*static_cast<ellipse<_T>*>(__shape), *static_cast<circle<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else if(dynamic_cast<ellipse<_T>*>(shape.get())){
                //check distance
                if(collision(*static_cast<ellipse<_T>*>(__shape), *static_cast<ellipse<_T>*>(shape.get()))){
                    //collision
                    return true;
                }
            }else{
                throw std::runtime_error("collision(): no matching shape type");
            }
        }
    }else {
        throw std::runtime_error("collision(): no matching shape type");
    }

    return false;
}

//template<template<class> class _PTR, typename _T>
//bool collision(std::vector<_PTR<shape_base<_T>>> const& __shapes, circle<_T> const& __circle){
//    for(auto& shape : __shapes){
//        if(dynamic_cast<circle<_T>*>(shape.get())){
//            //check distance
//            if(collision_details(__circle, *static_cast<circle<_T>*>(shape.get()))){
//                //collision
//                return true;
//            }
//        }else if(dynamic_cast<ellipse<_T>*>(shape.get())){
//            //check distance
//            //            if(check_distance(__circle, *static_cast<ellipse<_T>*>(shape.get()))){
//            //                //collision
//            //                return true;
//            //            }
//        }
//    }
//    return false;
//}

//template<template<class> class _PTR, typename _T>
//bool collision(std::vector<_PTR<shape_base<_T>>> const& __shapes, cylinder<_T> const& __cylinder){
//    for(auto& shape : __shapes){
//        if(dynamic_cast<cylinder<_T>*>(shape.get())){
//            //check distance
//            if(collision_details(__cylinder, *static_cast<cylinder<_T>*>(shape.get()))){
//                //collision
//                return true;
//            }
//        }
//    }
//    return false;
//}

//template<template<class> class _PTR, typename _T>
//bool collision(std::vector<_PTR<shape_base<_T>>> const& __shapes, ellipse<_T> const& __ellipse){
//    for(auto& shape : __shapes){
//        if(dynamic_cast<ellipse<_T>*>(shape.get())){
//            //check distance
//            if(collision_details(__ellipse, *static_cast<ellipse<_T>*>(shape.get()))){
//                //collision
//                return true;
//            }
//        }
//    }
//    return false;
//}

//template<template<class> class _PTR, typename _T>
//bool collision(std::vector<_PTR<shape_base<_T>>> const& __shapes, ellipsoid<_T> const& __ellipsoid){
//    for(auto& shape : __shapes){
//        if(dynamic_cast<ellipsoid<_T>*>(shape.get())){
//            //check distance
//            if(collision_details(__ellipsoid, *static_cast<ellipsoid<_T>*>(shape.get()))){
//                //collision
//                return true;
//            }
//        }
//    }
//    return false;
//}

//template<typename T>
//bool check_distance(std::vector<circle<T>> const& vec_circle, circle<T> const& circle){
//    for(size_t i{0}; i<vec_circle.size(); ++i){
//        //check distance
//        if(check_distance(circle, vec_circle[i])){
//            //collision
//            return true;
//        }
//    }
//    return false;
//}


//template<typename T>
//bool check_distance(std::vector<ellipse<T>> const& vec_ellipse, ellipse<T> const& ellipse){
//    for(size_t i{0}; i<vec_ellipse.size(); ++i){
//        //check distance
//        if(check_distance(ellipse, vec_ellipse[i])){
//            //collision
//            return true;
//        }
//    }
//    return false;
//}


//template<typename T>
//bool rectanglecontrol(ellipse<T> const& lhs, ellipse<T> const& rhs){
////    using value_type = T;

////    //Mittelpunkt
////    //Mittelpunkt von lhs im Koordinatensystem von rhs
////    value_type centerx = lhs(0)-rhs(0);
////    value_type centery = lhs(1)-rhs(1);

////    //Rotation mit Winkel von rhs
////    value_type centerxrhs = cos(rhs.rotation()*M_PI)*centerx+sin(rhs.rotation()*M_PI)*centery;
////    value_type centeryrhs = -sin(rhs.rotation()*M_PI)*centerx+cos(rhs.rotation()*M_PI)*centery;

////    value_type x = centerxrhs*centerxrhs;
////    value_type a = rhs.radius_a()*rhs.radius_a();
////    value_type y = centeryrhs*centeryrhs;
////    value_type b = rhs.radius_b()*rhs.radius_b();


////    if ((x/a)+(y/b) <= 1){
////        return true;
////    }

////    //Mittelpunkt
////    //Mittelpunkt von rhs im Koordinatensystem von lhs
////    centerx = rhs(0)-lhs(0);
////    centery = rhs(1)-lhs(1);

////    //Rotation mit Winkel von lhs
////    centerxrhs = cos(lhs.rotation()*M_PI)*centerx+sin(lhs.rotation()*M_PI)*centery;
////    centeryrhs = -sin(lhs.rotation()*M_PI)*centerx+cos(lhs.rotation()*M_PI)*centery;

////    x = centerxrhs*centerxrhs;
////    a = lhs.radius_a()*lhs.radius_a();
////    y = centeryrhs*centeryrhs;
////    b = lhs.radius_b()*lhs.radius_b();


////    if ((x/a)+(y/b) <= 1){
////        return true;
////    }

////    value_type pointslhs[4][2];
////    value_type pointsrhs[4][2];

////    //lhs without rotation
////    pointslhs[0][0] = +lhs.radius_a();
////    pointslhs[0][1] = +lhs.radius_b();
////    pointslhs[1][0] = -lhs.radius_a();
////    pointslhs[1][1] = +lhs.radius_b();
////    pointslhs[2][0] = -lhs.radius_a();
////    pointslhs[2][1] = -lhs.radius_b();
////    pointslhs[3][0] = +lhs.radius_a();
////    pointslhs[3][1] = -lhs.radius_b();

////    //rhs without rotation
////    pointsrhs[0][0] = +rhs.radius_a();
////    pointsrhs[0][1] = +rhs.radius_b();
////    pointsrhs[1][0] = -rhs.radius_a();
////    pointsrhs[1][1] = +rhs.radius_b();
////    pointsrhs[2][0] = -rhs.radius_a();
////    pointsrhs[2][1] = -rhs.radius_b();
////    pointsrhs[3][0] = +rhs.radius_a();
////    pointsrhs[3][1] = -rhs.radius_b();

////    //rotation lhs

////    for (int i = 0; i <= 3; i++) {

////        value_type alpha = lhs.rotation()*M_PI;

////        value_type x = pointslhs[i][0];
////        value_type y = pointslhs[i][1];

////        pointslhs[i][0] = (x*cos(alpha)-y*sin(alpha))+lhs(0);
////        pointslhs[i][1] = (x*sin(alpha)+y*cos(alpha))+lhs(1);


////    }

////    for (int i = 0; i <= 3; i++) {

////        value_type alpha = rhs.rotation()*M_PI;

////        value_type x = pointsrhs[i][0];
////        value_type y = pointsrhs[i][1];

////        pointsrhs[i][0] = (x*cos(alpha)-y*sin(alpha))+rhs(0);
////        pointsrhs[i][1] = (x*sin(alpha)+y*cos(alpha))+rhs(1);

////    }

////    for (int i = 0; i <= 3; i++){

////        value_type x1 = pointslhs[i][0];
////        value_type y1 = pointslhs[i][1];
////        value_type x2 = pointslhs[i+1][0];
////        value_type y2 = pointslhs[i+1][1];

////        if (i ==3){
////            x2 = pointslhs[0][0];
////            y2 = pointslhs[0][1];
////        }

////        for (int j = 0; j <= 3; j++ ){


////            value_type x3 = pointsrhs[j][0];
////            value_type y3 = pointsrhs[j][1];
////            value_type x4 = pointsrhs[j+1][0];
////            value_type y4 = pointsrhs[j+1][1];

////            if (j ==3){
////                x4 = pointsrhs[0][0];
////                y4 = pointsrhs[0][1];
////            }

////            //y=mx+b


////            value_type m1 = ((y2-y1)/(x2-x1));
////            value_type m2 = ((y4-y3)/(x4-x3));

////            value_type b1 = y1-m1*x1;
////            value_type b2 = y3-m2*x3;

////            value_type xschnitt = (b2-b1)/(m1-m2);

////            if ((((xschnitt >= x1) && (xschnitt <= x2)) || ((xschnitt >= x2) && (xschnitt <= x1)))
////                    && (((xschnitt >= x3) && (xschnitt <= x4)) || ((xschnitt >= x4) && (xschnitt <= x3)))){
////                return true;
////            }
////        }
////    }
////    return false;
//}

//    using value_type = T;

//    //Mittelpunkt
//    //Mittelpunkt von lhs im Koordinatensystem von rhs
//    value_type centerx = lhs(0)-rhs(0);
//    value_type centery = lhs(1)-rhs(1);

//    //Rotation mit Winkel von rhs
//    value_type centerxrhs = cos(rhs.rotation()*M_PI)*centerx+sin(rhs.rotation()*M_PI)*centery;
//    value_type centeryrhs = -sin(rhs.rotation()*M_PI)*centerx+cos(rhs.rotation()*M_PI)*centery;

//    value_type x = centerxrhs*centerxrhs;
//    value_type a = rhs.radius_a()*rhs.radius_a();
//    value_type y = centeryrhs*centeryrhs;
//    value_type b = rhs.radius_b()*rhs.radius_b();


//    if ((x/a)+(y/b) <= 1){
//        return true;
//    }

//    //Mittelpunkt
//    //Mittelpunkt von rhs im Koordinatensystem von lhs
//    centerx = rhs(0)-lhs(0);
//    centery = rhs(1)-lhs(1);

//    //Rotation mit Winkel von lhs
//    centerxrhs = cos(lhs.rotation()*M_PI)*centerx+sin(lhs.rotation()*M_PI)*centery;
//    centeryrhs = -sin(lhs.rotation()*M_PI)*centerx+cos(lhs.rotation()*M_PI)*centery;

//    x = centerxrhs*centerxrhs;
//    a = lhs.radius_a()*lhs.radius_a();
//    y = centeryrhs*centeryrhs;
//    b = lhs.radius_b()*lhs.radius_b();


//    if ((x/a)+(y/b) <= 1){
//        return true;
//    }

//    for (int i=1; i <= 360; i++){
//        value_type const phi = M_PI*i*2/360;

//        //Bestimmung der Koordinaten im Koordiantensystem von lhs
//        value_type pointxrhs = (rhs(0)+rhs.radius_a()*cos(phi)*cos(rhs.rotation()*M_PI)-rhs.radius_b()*sin(phi)*sin(rhs.rotation()*M_PI))-lhs(0);
//        value_type pointyrhs = (rhs(1)+rhs.radius_a()*cos(phi)*sin(rhs.rotation()*M_PI)+rhs.radius_b()*sin(phi)*cos(rhs.rotation()*M_PI))-lhs(1);

//        //Umrechnung der Koordinaten um Drehwinkel von lhs
//        value_type pointxlhs = cos(lhs.rotation()*M_PI)*pointxrhs+sin(lhs.rotation()*M_PI)*pointyrhs;
//        value_type pointylhs = -sin(lhs.rotation()*M_PI)*pointxrhs+cos(lhs.rotation()*M_PI)*pointyrhs;

//        value_type x = pointxlhs*pointxlhs;
//        value_type a = lhs.radius_a()*lhs.radius_a();
//        value_type y = pointylhs*pointylhs;
//        value_type b = lhs.radius_b()*lhs.radius_b();


//        if (((x/a)+(y/b)) <= 1){
//            return true;
//        }
//    }
//    return false;


//template<typename T>
//bool collision(ellipsoid<T> const& lhs, ellipsoid<T> const& rhs){
//    //    //Mittelpunkt
//    //    //Mittelpunkt von lhs im Koordinatensystem von rhs
//    //    value_type center_x = lhs(0)-rhs(0);
//    //    value_type center_y = lhs(1)-rhs(1);
//    //    value_type center_z = lhs(2)-rhs(2);

//    //    //Rotation um die x-Achse mit Winkel von rhs
//    //    value_type rot_x_center_x_rhs = center_x;
//    //    value_type rot_x_center_y_rhs = center_y*cos(rhs.rotation_x()*2*M_PI)-center_z*sin(rhs.rotation_x()*2*M_PI);
//    //    value_type rot_x_center_z_rhs = center_y*sin(rhs.rotation_x()*2*M_PI)+center_z*cos(rhs.rotation_x()*2*M_PI);

//    //    //Rotation um die y-Achse mit Winkel von rhs
//    //    value_type rot_y_center_x_rhs = rot_x_center_x_rhs*cos(rhs.rotation_y()*2*M_PI)+rot_x_center_z_rhs*sin(rhs.rotation_y()*2*M_PI);
//    //    value_type rot_y_center_y_rhs = rot_x_center_y_rhs;
//    //    value_type rot_y_center_z_rhs = rot_x_center_x_rhs*sin(rhs.rotation_y()*2*M_PI)+rot_x_center_z_rhs*cos(rhs.rotation_y()*2*M_PI);

//    //    //Rotation um die z-Achse mit Winkel von rhs
//    //    value_type rot_z_center_x_rhs = rot_y_center_x_rhs*cos(rhs.rotation_z()*2*M_PI)-rot_y_center_y_rhs*sin(rhs.rotation_z()*2*M_PI);
//    //    value_type rot_z_center_y_rhs = rot_y_center_x_rhs*sin(rhs.rotation_z()*2*M_PI)+rot_y_center_y_rhs*cos(rhs.rotation_z()*2*M_PI);
//    //    value_type rot_z_center_z_rhs = rot_y_center_z_rhs;

//    //    value_type x = rot_z_center_x_rhs*rot_z_center_x_rhs;
//    //    value_type a = rhs.radius_a()*rhs.radius_a();
//    //    value_type y = rot_z_center_y_rhs*rot_z_center_y_rhs;
//    //    value_type b = rhs.radius_b()*rhs.radius_b();
//    //    value_type z = rot_z_center_z_rhs*rot_z_center_z_rhs;
//    //    value_type c = rhs.radius_c()*rhs.radius_c();


//    //    if (((x/a)+(y/b)+(z/c)) <= 1){
//    //        std::cout<<"Mittelpunkt"<<std::endl;
//    //        return true;
//    //    }

//    //    //Mittelpunkt von rhs im Koordinatensystem von lhs
//    //    center_x = rhs(0)-lhs(0);
//    //    center_y = rhs(1)-lhs(1);
//    //    center_z = rhs(2)-lhs(2);

//    //    //Rotation um die x-Achse mit Winkel von lhs
//    //    value_type rot_x_center_x_lhs = center_x;
//    //    value_type rot_x_center_y_lhs = center_y*cos(lhs.rotation_x()*2*M_PI)-center_z*sin(lhs.rotation_x()*2*M_PI);
//    //    value_type rot_x_center_z_lhs = center_y*sin(lhs.rotation_x()*2*M_PI)+center_z*cos(lhs.rotation_x()*2*M_PI);

//    //    //Rotation um die y-Achse mit Winkel von lhs
//    //    value_type rot_y_center_x_lhs = rot_x_center_x_lhs*cos(lhs.rotation_y()*2*M_PI)+rot_x_center_z_lhs*sin(lhs.rotation_y()*2*M_PI);
//    //    value_type rot_y_center_y_lhs = rot_x_center_y_lhs;
//    //    value_type rot_y_center_z_lhs = rot_x_center_x_lhs*sin(lhs.rotation_y()*2*M_PI)+rot_x_center_z_lhs*cos(lhs.rotation_y()*2*M_PI);

//    //    //Rotation um die z-Achse mit Winkel von lhs
//    //    value_type rot_z_center_x_lhs = rot_y_center_x_lhs*cos(lhs.rotation_z()*2*M_PI)-rot_y_center_y_lhs*sin(lhs.rotation_z()*2*M_PI);
//    //    value_type rot_z_center_y_lhs = rot_y_center_x_lhs*sin(lhs.rotation_z()*2*M_PI)+rot_y_center_y_lhs*cos(lhs.rotation_z()*2*M_PI);
//    //    value_type rot_z_center_z_lhs = rot_y_center_z_lhs;

//    //    x = rot_z_center_x_lhs*rot_z_center_x_lhs;
//    //    a = lhs.radius_a()*lhs.radius_a();
//    //    y = rot_z_center_y_lhs*rot_z_center_y_lhs;
//    //    b = lhs.radius_b()*lhs.radius_b();
//    //    z = rot_z_center_z_lhs*rot_z_center_z_lhs;
//    //    c = lhs.radius_c()*lhs.radius_c();


//    //    if (((x/a)+(y/b)+(z/c)) <= 1){
//    //        std::cout<<"Mittelpunkt"<<std::endl;
//    //        return true;
//    //    }

//    //    value_type f0[3] = {rhs(0)-lhs(0), rhs(1)-lhs(1), rhs(2)-lhs(2)};
//    //    value_type f1[3] = {rhs.radius_a(), 0, 0};
//    //    value_type f2[3] = {0, rhs.radius_b(), 0};
//    //    value_type f3[3] = {0, 0, rhs.radius_c()};

//    //    //Vektoren rotieren in Koordinatensystem von lhs
//    //    //x-Achse
//    //    vectorrotation3d(f0, rhs.rotation_x()-lhs.rotation_x(), 0);
//    //    vectorrotation3d(f1, rhs.rotation_x()-lhs.rotation_x(), 0);
//    //    vectorrotation3d(f2, rhs.rotation_x()-lhs.rotation_x(), 0);
//    //    vectorrotation3d(f3, rhs.rotation_x()-lhs.rotation_x(), 0);

//    //    //y-Achse
//    //    vectorrotation3d(f0, rhs.rotation_y()-lhs.rotation_y(), 1);
//    //    vectorrotation3d(f1, rhs.rotation_y()-lhs.rotation_y(), 1);
//    //    vectorrotation3d(f2, rhs.rotation_y()-lhs.rotation_y(), 1);
//    //    vectorrotation3d(f3, rhs.rotation_y()-lhs.rotation_y(), 1);

//    //    //z-Achse
//    //    vectorrotation3d(f0, rhs.rotation_z()-lhs.rotation_z(), 2);
//    //    vectorrotation3d(f1, rhs.rotation_z()-lhs.rotation_z(), 2);
//    //    vectorrotation3d(f2, rhs.rotation_z()-lhs.rotation_z(), 2);
//    //    vectorrotation3d(f3, rhs.rotation_z()-lhs.rotation_z(), 2);

//    //    value_type phi;
//    //    value_type theta;

//    //    for (int i=1; i<= 360; i++){
//    //        for (int j=-90; j<= 90; j++) {

//    //            phi = i*2*M_PI/360;
//    //            theta = j*M_PI_4/90;

//    //            x = ((f0[0]+f1[0]*cos(theta)*cos(phi)+f2[0]*cos(theta)*sin(phi)+f3[0]*sin(theta))*(f0[0]+f1[0]*cos(theta)*cos(phi)+f2[0]*cos(theta)*sin(phi)+f3[0]*sin(theta)));
//    //            a = lhs.radius_a()*lhs.radius_a();
//    //            y = ((f0[1]+f1[1]*cos(theta)*cos(phi)+f2[1]*cos(theta)*sin(phi)+f3[1]*sin(theta))*(f0[1]+f1[1]*cos(theta)*cos(phi)+f2[1]*cos(theta)*sin(phi)+f3[1]*sin(theta)));
//    //            b = lhs.radius_b()*lhs.radius_b();
//    //            z = ((f0[2]+f1[2]*cos(theta)*cos(phi)+f2[2]*cos(theta)*sin(phi)+f3[2]*sin(theta))*(f0[2]+f1[2]*cos(theta)*cos(phi)+f2[2]*cos(theta)*sin(phi)+f3[2]*sin(theta)));
//    //            c = lhs.radius_c()*lhs.radius_c();

//    //            if (((x/a)+(y/b)+(z/c)) <= 1){
//    //                std::cout<<"Kollision"<<std::endl;
//    //                return true;
//    //            };

//    //        }
//    //    }
//    //    return false;
//}
}
#endif // CHECK_DISTANCE_H
