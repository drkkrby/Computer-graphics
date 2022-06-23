#include "ray_tracing.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <iostream>
#include <limits>

bool pointInTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n, const glm::vec3& p)
{
    glm::vec3 v0v1 = glm::cross(n, v1 - v0);
    v0v1 = v0v1 / (glm::dot(v2 - v0, v0v1));

    glm::vec3 v0v2 = glm::cross(n, v0 - v2);
    v0v2 = v0v2 / (glm::dot(v1 - v0, v0v2));

    float gamma = glm::dot(p - v2, v0v2);
    float beta = glm::dot(p - v1, v0v1);
    float alpha = 1 - (gamma + beta);

    if (alpha < 0 || alpha >1)
        return false;
    if (beta < 0 || beta >1)
        return false;
    if (gamma < 0 || gamma >1)
        return false;

    return true;
}

bool intersectRayWithPlane(const Plane& plane, Ray& ray)
{
    float t_n = plane.D - glm::dot(ray.origin, plane.normal);
    float t_d = glm::dot(ray.direction, plane.normal);
    float t = t_n / t_d;

    if (t <= 0)
        return false;
    else
    {
        if (t < ray.t)
        {
            ray.t = t;
            return true;
        }
        else
            return false;
    }
}

Plane trianglePlane(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2)
{
    Plane plane;

    plane.normal = glm::normalize(glm::cross(v0 - v2, v1 - v2));
    plane.D = glm::dot(plane.normal, v0);

    return plane;
}

/// Input: the three vertices of the triangle
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithTriangle(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, Ray& ray, HitInfo& hitInfo)
{

    bool final = false;

    Plane planes = trianglePlane(v0, v1, v2);
    float o_t = ray.t;
    bool intersect = intersectRayWithPlane(planes, ray);
    glm::vec3 p = ray.origin + ray.direction * ray.t;
    if (intersect)
    {
        final = pointInTriangle(v0, v1, v2, planes.normal, p);
    }
    if (final == false && intersect == true)
    {
        ray.t = o_t;
    }

    return final;

}

/// Input: a sphere with the following attributes: sphere.radius, sphere.center
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const Sphere& sphere, Ray& ray, HitInfo& hitInfo)
{
    glm::vec3 oc = ray.origin - sphere.center;
    glm::vec3 dir = glm::normalize(ray.direction);

    float a = glm::dot(dir, dir);
    float b = 2.0f * glm::dot(oc, dir);
    float c = glm::dot(oc, oc) - sphere.radius * sphere.radius;

    float discriminant = b * b - 4.0f * a * c;
    if (discriminant < 0)
    {
        return false;
    }
    else
    {
        float t0 = (-b - sqrt(discriminant)) / (2.0f * a);
        float t1 = (-b + sqrt(discriminant)) / (2.0f * a);

        float t = 0.0f;

        if (t0 > t1)
            std::swap(t0, t1);

        if (t0 < 0)
        {
            t0 = t1;
            if (t0 < 0)
                return false;
        }

        t = t0;

        if (t <= 0)
            return false;
        else
        {
            if (t < ray.t)
            {
                ray.t = t;
                glm::vec3 normal = glm::normalize((ray.origin + ray.direction * ray.t) - sphere.center);
                hitInfo.normal = normal;
                return true;
            }
            else
                return false;
        }

    }
}

/// Input: an axis-aligned bounding box with the following parameters: minimum coordinates box.lower and maximum coordinates box.upper
/// Output: if intersects then modify the hit parameter ray.t and return true, otherwise return false
bool intersectRayWithShape(const AxisAlignedBox& box, Ray& ray)
{
    float tx_min = (box.lower.x - ray.origin.x) / ray.direction.x;
    float tx_max = (box.upper.x - ray.origin.x) / ray.direction.x;

    float ty_min = (box.lower.y - ray.origin.y) / ray.direction.y;
    float ty_max = (box.upper.y - ray.origin.y) / ray.direction.y;

    float tz_min = (box.lower.z - ray.origin.z) / ray.direction.z;
    float tz_max = (box.upper.z - ray.origin.z) / ray.direction.z;

    float tx_in = glm::min(tx_min, tx_max);
    float tx_out = glm::max(tx_min, tx_max);

    float ty_in = glm::min(ty_min, ty_max);
    float ty_out = glm::max(ty_min, ty_max);

    float tz_in = glm::min(tz_min, tz_max);
    float tz_out = glm::max(tz_min, tz_max);

    float t_in1 = glm::max(tx_in, ty_in);
    float t_in = glm::max(t_in1, tz_in);

    float t_out1 = glm::min(tx_out, ty_out);
    float t_out = glm::min(t_out1, tz_out);

    if (t_in > t_out || t_out < 0)
        return false;
    else
    {
        if (t_in < ray.t)
        {
            if (t_in < 0)
            {
                ray.t = t_out;
                return true;
            }
            ray.t = t_in;
            return true;
        }
        else
            return false;

    }
}
