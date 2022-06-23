#include "bounding_volume_hierarchy.h"
#include "draw.h"

#include <glm/geometric.hpp>
#include <glm/gtx/component_wise.hpp>
#include <glm/vector_relational.hpp>
#include <iostream>
#ifdef __APPLE__
#include <OpenGL/GLU.h>
#else
#ifdef WIN32
// Windows.h includes a ton of stuff we don't need, this macro tells it to include less junk.
#define WIN32_LEAN_AND_MEAN
// Disable legacy macro of min/max which breaks completely valid C++ code (std::min/std::max won't work).
#define NOMINMAX
// GLU requires Windows.h on Windows :-(.
#include <Windows.h>
#endif
#include <GL/glu.h>
#endif



BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    // 1) get all the vertices 
    Node all_nodes;

    for (int meshIndex = 0; meshIndex < pScene->meshes.size(); meshIndex++) {
        Mesh m = pScene->meshes[meshIndex];
        for (int tIndex = 0; tIndex < m.triangles.size(); tIndex++) {
            glm::vec3 t = m.triangles[tIndex];
            TriangleMeshIndex tMeshIndex;
            tMeshIndex.meshIndex = meshIndex;
            tMeshIndex.triangleIndex = tIndex;
            allTrianglesIndex.push_back(tMeshIndex);
            Triangle tri;
            tri.triangles = {
                    m.vertices[t[0]].position,
                    m.vertices[t[1]].position,
                    m.vertices[t[2]].position
            };
            tri.triangleMeshPointer = allTrianglesIndex.size() - 1;
            all_nodes.triangles.push_back(tri);

           
           
        }
    }
    this->nodes.push_back(all_nodes);
    hierarchy(all_nodes, 0, 0, 0);
    // 3) select the greatest and lowest vertex in any direction for box
    // 4) set up box
    // 5) select half of vertices for both sides. 
    // 6) recursive call

}

void BoundingVolumeHierarchy::hierarchy(Node node, int level, int direction, int position) {


    if(level == this->numLevels() || node.triangles.size() <= 8) {
        node.triangles = mergeSort(node.triangles, direction);
        node.isLeaf = true;
        if (maxLevel < level) {
            maxLevel = level;
        }
        node.depth = maxLevel;
        for (Triangle t : node.triangles)
        {
            node.childrenTriangles.push_back(t.triangleMeshPointer);
        }
        node.box = getBox(node);
        this->nodes[position] = node;
        return;
    }

    node.triangles = mergeSort(node.triangles, direction);
    node.box = getBox(node);
    node.depth = level;
    
    Node left;
    Node right;
            
    for (int i = 0; i < node.triangles.size(); i++) {
        if (i < floor(node.triangles.size() / 2)) {
            left.triangles.push_back(node.triangles[i]);
        }
        else {
            right.triangles.push_back(node.triangles[i]);
        }
    }
    direction = (direction + 1) % 3;

    // 6) recursive call
    //std::cout << "Start recursion for layer " << level + 1 << "(left)" << std::endl;
    node.index_right_child = this->nodes.size();
    this->nodes.push_back(right);
    node.index_left_child = this->nodes.size();
    this->nodes.push_back(left);
    if (level == 0) {
        rootNode = position;
    }
    this->nodes[position] = node;
    
    hierarchy(left, level + 1, direction, node.index_left_child);
    hierarchy(right, level + 1, direction, node.index_right_child);  
    return;
}
// Creates the AABB for the node.
AxisAlignedBox BoundingVolumeHierarchy::getBox(BoundingVolumeHierarchy::Node node)
{
    float min_x = node.triangles[0].triangles[0].x;
    float min_y = node.triangles[0].triangles[0].y;
    float min_z = node.triangles[0].triangles[0].z;
    float max_x = node.triangles[0].triangles[0].x;
    float max_y = node.triangles[0].triangles[0].y;
    float max_z = node.triangles[0].triangles[0].z;
    //std::cout << "stil going 1 " << level << std::endl;
    for (Triangle t : node.triangles) {
        for (glm::vec3 v : t.triangles) {
            if (min_x > v.x) {
                min_x = v.x;
            }
            if (min_y > v.y) {
                min_y = v.y;
            }
            if (min_z > v.z) {
                min_z = v.z;
            }
            if (max_x < v.x) {
                max_x = v.x;
            }
            if (max_y < v.y) {
                max_y = v.y;
            }
            if (max_z < v.z) {
                max_z = v.z;
            }
        }
    }

    // 4) set up box and save in tree;

    return { { min_x, min_y, min_z },{ max_x, max_y, max_z } };
}

std::vector<Triangle> BoundingVolumeHierarchy::mergeSort(std::vector<Triangle> triangles, int direction){
    
    if (triangles.size() <= 1) {
        return triangles;
    }
    else {
        std::vector<Triangle> left = std::vector<Triangle>();
        std::vector<Triangle> right = std::vector<Triangle>();
        
        for (int i = 0; i < triangles.size(); i++) {
            if (i < floor(triangles.size() / 2)) {
                left.push_back(triangles[i]);
            }
            else {
                right.push_back(triangles[i]);
            }
        }
        left = mergeSort(left, direction);
        right = mergeSort(right, direction);


        return merge(left, right, direction);

    }
}
// merge method for mergeSort
std::vector<Triangle> BoundingVolumeHierarchy::merge(std::vector<Triangle> left, std::vector<Triangle> right, int direction) {
    int l_index = 0;
    int r_index = 0;

    std::vector<Triangle> sorted_trianges = std::vector<Triangle>();

    if (direction == 2) {
        while (l_index < left.size() && r_index < right.size()) {
            if ((left[l_index].triangles[0].x + left[l_index].triangles[1].x + left[l_index].triangles[2].x) / 3 < (right[r_index].triangles[0].x + right[r_index].triangles[1].x + right[r_index].triangles[2].x) / 3) {

                sorted_trianges.push_back(left[l_index]);
                l_index++;
            }
            else {

                sorted_trianges.push_back(right[r_index]);
                r_index++;
            }
        }
        while (l_index < left.size()) {
            sorted_trianges.push_back(left[l_index]);
            l_index++;
        }
        while (r_index < right.size()) {
            sorted_trianges.push_back(right[r_index]);
            r_index++;
        }
    }
    if (direction == 1) {
        while (l_index < left.size() && r_index < right.size()) {
            if ((left[l_index].triangles[0].y + left[l_index].triangles[1].y + left[l_index].triangles[2].y) / 3 < (right[r_index].triangles[0].y + right[r_index].triangles[1].y + right[r_index].triangles[2].y) / 3) {

                sorted_trianges.push_back(left[l_index]);
                l_index++;
            }
            else {

                sorted_trianges.push_back(right[r_index]);
                r_index++;
            }
        }
        while (l_index < left.size()) {
            sorted_trianges.push_back(left[l_index]);
            l_index++;
        }
        while (r_index < right.size()) {
            sorted_trianges.push_back(right[r_index]);
            r_index++;
        }
    }
    if (direction == 0) {
        while (l_index < left.size() && r_index < right.size()) {
            if ((left[l_index].triangles[0].z + left[l_index].triangles[1].z + left[l_index].triangles[2].z) / 3 < (right[r_index].triangles[0].z + right[r_index].triangles[1].z + right[r_index].triangles[2].z) / 3) {

                sorted_trianges.push_back(left[l_index]);
                l_index++;
            }
            else {

                sorted_trianges.push_back(right[r_index]);
                r_index++;
            }
        }
        while (l_index < left.size()) {
            sorted_trianges.push_back(left[l_index]);
            l_index++;
        }
        while (r_index < right.size()) {
            sorted_trianges.push_back(right[r_index]);
            r_index++;
        }
    }
    return sorted_trianges;
}
// Return the depth of the tree that is constructed. This is used to tell the
// slider in the UI how many steps it should display.
int BoundingVolumeHierarchy::numLevels() const
{
    return 15;
}

// Visualizes BV hierarchy. Creates different color values for each iteration in the loop that makes 
// all the triangles inside one leaf node in the BVH the same color, while other leaf nodes are different color.
void BoundingVolumeHierarchy::debugDraw(int level)
{
    // Draw the AABB as a (white) wireframe box.
    float R = 1.0f;
    float G = 0.0f;
    float B = 0.8f;
    for (Node n : this->nodes) 
    {
        if (n.depth == level)
        {
            if (n.isLeaf)
            {
                R = fmod(R + 0.1f, 1.0f);
                G = fmod(G + 0.2f, 1.0f);
                B = fmod(B + 0.3f, 1.0f);

                for (Triangle t : n.triangles)
                {
                    //Gets the triangles from the mesh and then gets the vertices of the triangle
                    TriangleMeshIndex tIndex = allTrianglesIndex[t.triangleMeshPointer];
                    glm::vec3 tr = m_pScene->meshes[tIndex.meshIndex].triangles[tIndex.triangleIndex];
                    Vertex v0 = m_pScene->meshes[tIndex.meshIndex].vertices[tr.x];
                    Vertex v1 = m_pScene->meshes[tIndex.meshIndex].vertices[tr.y];
                    Vertex v2 = m_pScene->meshes[tIndex.meshIndex].vertices[tr.z];
                    //Draws a triangle
                    glColor3f(R,G,B);
                    glBegin(GL_TRIANGLES);
                    glVertex3f(v0.position.x, v0.position.y, v0.position.z);
                    glVertex3f(v1.position.x, v1.position.y, v1.position.z);
                    glVertex3f(v2.position.x, v2.position.y, v2.position.z);
                    glEnd();
                   
                    drawAABB(n.box, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.2f);

                }
            }
            else {
                drawAABB(n.box, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 0.4f);
            }
        }
    }
    
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h
bool BoundingVolumeHierarchy::intersectHelper(Ray& ray, HitInfo& hitInfo, Node node) const
{
    Ray temp = ray;
    float closestBoxDistance = std::numeric_limits<float>::max();
    if (intersectRayWithShape(node.box, temp))
    {
        glm::vec3 boxIntersectionPoint = temp.origin + temp.direction * temp.t;
        float distance_to_ray = glm::distance(temp.origin, boxIntersectionPoint);
        if (distance_to_ray > closestBoxDistance)
            return true;
        closestBoxDistance = distance_to_ray;
        //Leaf node, enter and search
        if (node.isLeaf)
        {
            bool hit = false;
            Ray min_ray; min_ray.t = 0.0f;
            int min_t_Index;
            for (int tIndex : node.childrenTriangles)
            {
                TriangleMeshIndex tMesh = allTrianglesIndex[tIndex];
                glm::vec3 t = m_pScene->meshes[tMesh.meshIndex].triangles[tMesh.triangleIndex];

                const auto v0 = m_pScene->meshes[tMesh.meshIndex].vertices[t[0]];
                const auto v1 = m_pScene->meshes[tMesh.meshIndex].vertices[t[1]];
                const auto v2 = m_pScene->meshes[tMesh.meshIndex].vertices[t[2]];
                //Ray intersects with triangle, so update values.
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    if (min_ray.t < ray.t)
                    {
                        min_ray = ray;
                        min_t_Index = tIndex;
                        hit = true;
                    }
                }
            }
            if (hit)
            {
                TriangleMeshIndex tMesh = allTrianglesIndex[min_t_Index];
                glm::vec3 t = m_pScene->meshes[tMesh.meshIndex].triangles[tMesh.triangleIndex];
                //Vertices of triangle
                const auto v0 = m_pScene->meshes[tMesh.meshIndex].vertices[t[0]];
                const auto v1 = m_pScene->meshes[tMesh.meshIndex].vertices[t[1]];
                const auto v2 = m_pScene->meshes[tMesh.meshIndex].vertices[t[2]];
                //The shortest ray
                glm::vec3 p = min_ray.origin + min_ray.direction * min_ray.t;
                //Vectors of the plane of the triangle
                glm::vec3 a = v1.position - v0.position;
                glm::vec3 b = v2.position - v0.position;
                //Vector for the point and vertex
                glm::vec3 c = p - v0.position;
                //Calculates the u, v, w values for barycentric coordinates
                float d00 = glm::dot(a, a);
                float d01 = glm::dot(a, b);
                float d11 = glm::dot(b, b);
                float d20 = glm::dot(c, a);
                float d21 = glm::dot(c, b);
                float denom = 1 / (d00 * d11 - d01 * d01);
                float v = (d11 * d20 - d01 * d21) * denom;
                float w = (d00 * d21 - d01 * d20) * denom;
                float u = 1.0f - v - w;
                //Texture coordinate calculation and interpolation
                glm::vec2 v0T = v0.texCoord;
                glm::vec2 v1T = v1.texCoord;
                glm::vec2 v2T = v2.texCoord;
                hitInfo.texel = u * v0T + v * v1T + w * v2T;
                //Normal interpolation
                glm::vec3 v0N = v0.normal;
                glm::vec3 v1N = v1.normal;
                glm::vec3 v2N = v2.normal;
                hitInfo.normal = u * v0N + v * v1N + w * v2N;
                //Draws the vertex normals and the interpolated normal
                Ray v0R;
                Ray v1R;
                Ray v2R;
                Ray vIR;
                v0R.origin = v0.position;
                v1R.origin = v1.position;
                v2R.origin = v2.position;
                vIR.origin = p;
                v0R.direction = v0N;
                v1R.direction = v1N;
                v2R.direction = v2N;
                vIR.direction = hitInfo.normal;
                drawRay(v0R, glm::vec3(1.0f, 0.0f, 1.0f));
                drawRay(v1R, glm::vec3(1.0f, 0.0f, 1.0f));
                drawRay(v2R, glm::vec3(1.0f, 0.0f, 1.0f));
                drawRay(vIR, glm::vec3(0.0f, 1.0f, 0.0f));
                //Update the material
                hitInfo.material = m_pScene->meshes[tMesh.meshIndex].material;
                return true;
            }
        }
        else {
            //Interior node, continue
            intersectHelper(ray, hitInfo, this->nodes[node.index_left_child]);
            intersectHelper(ray, hitInfo, this->nodes[node.index_right_child]);
        }
    }
    else {
        return false;
    }
}
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo) const
{
    bool original_traversal = true;
    if (original_traversal) {
        return this->original_intersect(ray, hitInfo);
    }
    bool hit = false;
    // Intersect with all triangles of all meshes.
    Node node = this->nodes[this->rootNode];
    if (intersectHelper(ray, hitInfo, node))
        hit = true;
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres) {
        // hit |= intersectRayWithShape(sphere, ray, hitInfo);
        if (intersectRayWithShape(sphere, ray, hitInfo)) {
            hitInfo.material = sphere.material;
            hit = true;
        }
    }
    return hit;
}
bool BoundingVolumeHierarchy::original_intersect(Ray& ray, HitInfo& hitInfo) const
{
    bool hit = false;
    // Intersect with all triangles of all meshes.
    for (const auto& mesh : m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            const auto v0 = mesh.vertices[tri[0]];
            const auto v1 = mesh.vertices[tri[1]];
            const auto v2 = mesh.vertices[tri[2]];
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                glm::vec3 p = ray.origin + ray.direction * ray.t;
                glm::vec3 a = v1.position - v0.position;
                glm::vec3 b = v2.position - v0.position;
                glm::vec3 c = p - v0.position;
                float d00 = glm::dot(a, a);
                float d01 = glm::dot(a, b);
                float d11 = glm::dot(b, b);
                float d20 = glm::dot(c, a);
                float d21 = glm::dot(c, b);
                float denom = 1 / (d00 * d11 - d01 * d01);
                float v = (d11 * d20 - d01 * d21) * denom;
                float w = (d00 * d21 - d01 * d20) * denom;
                float u = 1.0f - v - w;

                glm::vec2 v0T = v0.texCoord;
                glm::vec2 v1T = v1.texCoord;
                glm::vec2 v2T = v2.texCoord;
                hitInfo.texel = u * v0T + v * v1T + w * v2T;

                glm::vec3 v0N = v0.normal;
                glm::vec3 v1N = v1.normal;
                glm::vec3 v2N = v2.normal;
                hitInfo.normal = u * v0N + v * v1N + w * v2N;
                Ray v0R;
                Ray v1R;
                Ray v2R;
                Ray vIR;
                v0R.origin = v0.position;
                v1R.origin = v1.position;
                v2R.origin = v2.position;
                vIR.origin = p;
                v0R.direction = v0N;
                v1R.direction = v1N;
                v2R.direction = v2N;
                vIR.direction = hitInfo.normal;
                drawRay(v0R, glm::vec3(1.0f, 0.0f, 1.0f));
                drawRay(v1R, glm::vec3(1.0f, 0.0f, 1.0f));
                drawRay(v2R, glm::vec3(1.0f, 0.0f, 1.0f));
                drawRay(vIR, glm::vec3(0.0f, 1.0f, 0.0f));
                hitInfo.material = mesh.material;
                hit = true;
            }
        }
    }
    // Intersect with spheres.
    for (const auto& sphere : m_pScene->spheres) {
        // hit |= intersectRayWithShape(sphere, ray, hitInfo);
        if (intersectRayWithShape(sphere, ray, hitInfo)) {
            hitInfo.material = sphere.material;
            hit = true;
        }
    }
    return hit;
}