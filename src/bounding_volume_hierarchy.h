#pragma once
#include "ray_tracing.h"
#include "bounding_volume_hierarchy.h"
#include "scene.h"
#include <array>
#include <span>

struct Triangle {
    std::vector<glm::vec3> triangles;
    int triangleMeshPointer;
};

class BoundingVolumeHierarchy {    

    struct TriangleMeshIndex {
        int meshIndex;
        int triangleIndex;
    };

    struct Node {
        bool isLeaf = false;
        std::vector<Triangle> triangles;
        int depth;
        AxisAlignedBox box;
        int index_left_child = -1;
        int index_right_child = -1;
        std::vector<int> childrenTriangles;
    };
    
public:
    BoundingVolumeHierarchy(Scene* pScene);
    //The constructor for the BVH hierarchy
    void hierarchy(Node all_triangles, int level, int direction, int position);
    //AABB generator for a node
    AxisAlignedBox getBox(BoundingVolumeHierarchy::Node node);
    //mergeSort for triangles
    std::vector<Triangle> mergeSort(std::vector<Triangle> triangles, int direction);
    std::vector<Triangle> merge(std::vector<Triangle> left, std::vector<Triangle> right, int direction);
    //Vector for keeping nodes
    std::vector<Node> nodes;
    //Vector for keeping the indices of the triangle in the mesh
    std::vector<TriangleMeshIndex> allTrianglesIndex;
    //Root Node index
    int rootNode;
    //Caps the max level for the tree depending on the object
    int maxLevel = 0;
    // The first function should return how many levels there are in the tree that you have constructed.
    // The second function should draw the bounding boxes of the nodes at the selected level.
    int numLevels() const;
    void debugDraw(int level);
    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool original_intersect(Ray& ray, HitInfo& hitInfo) const;
    bool intersect(Ray& ray, HitInfo& hitInfo) const;
    bool intersectNew(Ray& ray, HitInfo& hitInfo) const;
    bool intersectHelper(Ray& ray, HitInfo& hitInfo, Node node) const;
private:
    Scene* m_pScene;
};
