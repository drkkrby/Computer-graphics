#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "ray_tracing.h"
#include "screen.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <imgui.h>
#include <nfd.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
DISABLE_WARNINGS_POP()
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <framework/image.h>
#include <framework/imguizmo.h>
#include <framework/trackball.h>
#include <framework/variant_helper.h>
#include <framework/window.h>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <type_traits>
#include <variant>

// This is the main application. The code in here does not need to be modified.
constexpr glm::ivec2 windowResolution{ 800, 800 };
const std::filesystem::path dataPath{ DATA_DIR };

// parameters that we use for debugging purposses. 
bool recursive = true;
int GR_number_rays = 25;
float focal_length = 2;
float aperture = 0.2;
bool dof_off = true;

enum class ViewMode {
    Rasterization = 0,
    RayTracing = 1
};

glm::vec3 phongShading(const HitInfo& hitInfo, const glm::vec3& vertexPos,
    const glm::vec3& lightPos, const glm::vec3& cameraPos, const glm::vec3& color)
{

    glm::vec3 lightDir = lightPos - vertexPos;  // Direction of light

    glm::vec3 normally = glm::normalize(hitInfo.normal);
    lightDir = glm::normalize(lightDir);

    float diff = glm::dot(normally, lightDir);
    diff = glm::max(diff, 0.0f);
    glm::vec3 diffuse;
    if (hitInfo.material.kdTexture) {
        diffuse = diff * hitInfo.material.kdTexture->getTexel(hitInfo.texel);
    }
    else {
        diffuse = diff * hitInfo.material.kd; // Diffuse component of Phong
    }

    glm::vec3 reflectDir = glm::normalize(glm::reflect((-lightDir), normally));
    glm::vec3 posToViewDir = glm::normalize(cameraPos - vertexPos);

    float specularConstant = glm::pow(glm::max(glm::dot(posToViewDir, reflectDir), 0.0f), hitInfo.material.shininess);
    glm::vec3 specular = hitInfo.material.ks * specularConstant; // Specular component of phong

    glm::vec3 phong = color * (diffuse + specular); // Final Phong model

    return phong;

}

bool hard_shadows(const BoundingVolumeHierarchy& bvh, const Ray ray,
    const glm::vec3& lightPos, const glm::vec3& cameraPos, HitInfo hit)
{
    glm::vec3 intersection = ray.t * ray.direction + ray.origin;
    Ray shadow_ray;
    HitInfo hitInfo;
    float epsilon = 0.0001f;
    glm::vec3 final_colour{ 0.0f };

    if (glm::dot(hit.normal, ray.direction) > 0)
    {
        hit.normal = -hit.normal;
    }

    shadow_ray.direction = glm::normalize(lightPos - intersection);
    shadow_ray.origin = intersection + epsilon * shadow_ray.direction;
    bool result = bvh.intersect(shadow_ray, hitInfo);

    float distIL = glm::distance(intersection, lightPos);
    float distIS = glm::distance(intersection, shadow_ray.t * shadow_ray.direction + shadow_ray.origin);
    float dott = glm::dot(glm::normalize(shadow_ray.direction), glm::normalize(hit.normal));

    bool shade;


    if ((result == true && distIS <= distIL) || (dott < 0)) // Check for occulsion before the lights
    {
        //drawRay(shadow_ray, glm::vec3(0.0f, 0.0f, 1.0f));  
        shade = true;
        return shade;

    }
    else
    {
        shadow_ray.t = distIL;
        //drawRay(shadow_ray, glm::vec3(0.0f, 1.0f, 1.0f)); // FOr debugging
        shade = false;
        return shade;

    }
}


bool soft_shadows(const BoundingVolumeHierarchy& bvh, const Ray ray,
    const glm::vec3& lightPos, const glm::vec3& cameraPos, HitInfo hit, const glm::vec3& lightPosColor)
{
    glm::vec3 intersection = ray.t * ray.direction + ray.origin;
    Ray shadow_ray;
    HitInfo hitInfo;
    float epsilon = 0.0001f;
    glm::vec3 final_colour{ 0.0f };

    if (glm::dot(hit.normal, ray.direction) > 0)
    {
        hit.normal = -hit.normal;
    }

    shadow_ray.direction = glm::normalize(lightPos - intersection);
    shadow_ray.origin = intersection + epsilon * shadow_ray.direction;
    bool result = bvh.intersect(shadow_ray, hitInfo);

    float distIL = glm::distance(intersection, lightPos);
    float distIS = glm::distance(intersection, shadow_ray.t * shadow_ray.direction + shadow_ray.origin);
    float dott = glm::dot(glm::normalize(shadow_ray.direction), glm::normalize(hit.normal));

    bool shade;


    if ((result == true && distIS <= distIL) || (dott < 0))
    {
        //drawRay(shadow_ray, glm::vec3{ 1.0f, 1.0f, 1.0f });
        shade = true;
        return shade;

    }
    else
    {
        shadow_ray.t = distIL;
        //drawRay(shadow_ray, lightPosColor); // For debugging
        shade = false;
        return shade;

    }
}

static glm::vec3 getFinalColor(const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray, int depth)
{

    HitInfo hitInfo;
    glm::vec3 final_output{ 0.0f };

    if (bvh.intersect(ray, hitInfo))
    {
        for (const auto& light : scene.lights) // Iterate over all the lights in the scene
        {
            if (std::holds_alternative<PointLight>(light))
            {
                // Calculating color in case light is a point.
                const PointLight pointLight = std::get<PointLight>(light);

                glm::vec3 lightPos = pointLight.position;
                glm::vec3 camPos = ray.origin;
                glm::vec3 vertexPos = ray.origin + ray.t * ray.direction;

                if (hard_shadows(bvh, ray, lightPos, camPos, hitInfo) == false) // Check if region is under shadow
                {
                    if (glm::dot(hitInfo.normal, ray.direction) > 0)
                    {
                        hitInfo.normal = -hitInfo.normal;
                    }
                    final_output += phongShading(hitInfo, vertexPos,
                        lightPos, camPos, pointLight.color); // Shade where there are no shadows
                }
            }
            else if (std::holds_alternative<SegmentLight>(light))
            {
                // Calculating color in case light is a segment.

                const SegmentLight segmentLight = std::get<SegmentLight>(light);

                glm::vec3 endpoint0 = segmentLight.endpoint0;
                glm::vec3 endpoint1 = segmentLight.endpoint1;
                glm::vec3 color0 = segmentLight.color0;
                glm::vec3 color1 = segmentLight.color1;

                glm::vec3 final_color = { 0, 0, 0 };
                glm::vec3 camPos = ray.origin;
                glm::vec3 vertexPos = ray.origin + ray.t * ray.direction;
                const float intervals = 100;
                glm::vec3 color{ 0.0f };
                glm::vec3 seg_color{ 0.0f };

                for (float i = 0.0; i <= intervals; i++) 
                {

                    glm::vec3 point = (endpoint0 + (i / (intervals)) *
                        (endpoint1 - endpoint0)); // Each point acts as a individual point source

                    color = color0 * (1 - (i / intervals)) + color1 * (i / intervals); // Linear Interpolation

                    if (soft_shadows(bvh, ray, point, camPos, hitInfo, color) == false)
                    {

                        if (glm::dot(hitInfo.normal, ray.direction) > 0)
                        {
                            hitInfo.normal = -hitInfo.normal;
                        }
                        seg_color = phongShading(hitInfo, vertexPos, point, camPos, color);
                    }
                    glm::vec3 c_pointColor = seg_color / glm::vec3((intervals + 1));
                    final_output += c_pointColor;
                }
            }
            else if (std::holds_alternative<ParallelogramLight>(light))
            {
                // Calculating color in case light is a parallelogram.

                const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);

                glm::vec3 c = glm::vec3(0, 0, 0);


                glm::vec3 camPos = ray.origin;
                glm::vec3 vertexPos = ray.origin + ray.t * ray.direction;
                float intervals = 10.0f;
                // Construct 100 rays in the parallelogram uniformly distributed.
                for (float i = 0.0f; i <= 10; i++)
                {
                    for (float j = 0.0f; j <= 10; j++)
                    {
                        glm::vec3 light_pos =
                            parallelogramLight.v0 +
                            (parallelogramLight.edge01) * glm::vec3(i / intervals) +
                            (parallelogramLight.edge02) * glm::vec3(j / intervals); // Calculating points on the paralellogram


                        // calculate the color based on the distances between the origin and the corners. 
                        // The closer to the corner the more that color is present in that ray's color.

                        glm::vec3 light_pos_color = parallelogramLight.color0 * (1 - (i / intervals)) * (1 - (j / intervals))
                            + parallelogramLight.color1 * (i / intervals) * (1 - (j / intervals))
                            + parallelogramLight.color2 * (1 - (i / intervals)) * (j / intervals)
                            + parallelogramLight.color3 * (i / intervals) * (j / intervals); // Bilinear Interpolation

                        if (soft_shadows(bvh, ray, light_pos, camPos, hitInfo, light_pos_color) == false)
                        {

                            if (glm::dot(hitInfo.normal, ray.direction) > 0)
                            {
                                hitInfo.normal = -hitInfo.normal;
                            }
                            c += phongShading(hitInfo, vertexPos, light_pos, camPos, light_pos_color);
                        }
                    }

                }
                glm::vec3 c_pointColor = c / glm::vec3(glm::pow(intervals + 1, 2));
                final_output += c_pointColor;

            }
        }
        // reflection 
        // takes place if specular reflection is at least higher than 0
        if (glm::length(hitInfo.material.ks) > 0 && depth > 0)
        {


            // recursive ray 
            Ray reflectionRay;

            reflectionRay.direction = -glm::normalize(glm::reflect((-glm::normalize(ray.direction)), hitInfo.normal));
            reflectionRay.origin = ray.t * ray.direction + ray.origin + reflectionRay.direction * glm::vec3(0.0001);


            if (recursive) { // To switch between Reflections and glossy reflections
                // recursive call
                final_output += hitInfo.material.ks * getFinalColor(scene, bvh, reflectionRay, depth - 1);
            }
            else {
                // Glossy reflection
                glm::vec3 u = glm::normalize(glm::vec3{ 1, 0, 0 } + 
                    reflectionRay.direction * -glm::dot({ 1, 0, 0 }, reflectionRay.direction));
                glm::vec3 v = glm::normalize(glm::vec3{ 0, 1, 0 } + 
                    reflectionRay.direction * -glm::dot({ 0, 1, 0 }, reflectionRay.direction));

                float blur_degree = hitInfo.material.shininess;

                Ray glossyReflectionRay;
                glm::vec3 gloss_color = glm::vec3(0);

                // take average of multiple rays to create the blur.
                for (int i = 0; i < GR_number_rays; i++) {
                    glossyReflectionRay.direction = reflectionRay.direction + ((-blur_degree / 2) + 
                        ((float)rand() / RAND_MAX) * blur_degree) * u + ((-blur_degree / 2) + 
                            ((float)rand() / RAND_MAX) * blur_degree) * v;
                    glossyReflectionRay.origin = ray.t * ray.direction + ray.origin + 
                        glossyReflectionRay.direction * glm::vec3(0.0001);
                    gloss_color += getFinalColor(scene, bvh, glossyReflectionRay, depth - 1);
                }
                final_output += hitInfo.material.ks * gloss_color / glm::vec3(GR_number_rays);
            }

           
            //std::cout << std::endl;
            //std::cout << glm::dot(u, reflectionRay.direction) << std::endl;
            //std::cout << glm::dot(v, reflectionRay.direction) << std::endl;
            //std::cout << glm::dot(v, u) << std::endl;
            
            // glossy reflection
            
        }

        // Transparency 
        if (hitInfo.material.transparency < 1 && depth > 0) {
            Ray transparency_ray;
            transparency_ray.origin = ray.direction * glm::vec3(ray.t + 0.0001) + ray.origin;
            transparency_ray.direction = ray.direction;

            // Get color at next intersection point
            glm::vec3 transparency_color = getFinalColor(scene, bvh, transparency_ray, depth);

            final_output = final_output * hitInfo.material.transparency +
                transparency_color * (1 - hitInfo.material.transparency);
            //drawRay(transparency_ray, glm::vec3(0.0f, 1.0f, 0.0f)); // Disable the other VDs to see this
        }

        // depth of field
        //drawRay(ray, final_output);
        //return final_output;
        if ((ray.t <= focal_length + aperture && ray.t >= focal_length - aperture) || depth <= 0 || dof_off) {
            // case 1: inside focal range -> standard shading
            return final_output;
        }
        else {
            // case 2: outside focal range -> blurred shading

            glm::vec3 dof_color = { 0.0f, 0.0f, 0.0f };
            int num_rays = 5;
            // take average of multiple rays to get blur.
            for (int i = 0; i < num_rays; i++) {
                // randomize origin
                glm::vec3 shift = { (((double)rand() / (RAND_MAX)) - 0.5), 
                    (((double)rand() / (RAND_MAX)) - 0.5), 
                    (((double)rand() / (RAND_MAX)) - 0.5) };

                Ray convergence_ray;

                convergence_ray.origin = ray.origin + (shift * aperture);
                convergence_ray.direction = glm::normalize((ray.direction * focal_length) + ray.origin - convergence_ray.origin);

                //drawRay(convergence_ray, glm::vec3{ 1, 1, 1 });

                dof_color += getFinalColor(scene, bvh, convergence_ray, depth - 1);

            }
            return dof_color / glm::vec3(num_rays);

        }

        // Draw a white debug ray if the ray hits.
        drawRay(ray, final_output);
        // Set the color of the pixel to white if the ray hits.
        //return (final_output + dof_color);
    }
    else
    {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}
static glm::vec3 getFinalColor(const Scene& scene, const BoundingVolumeHierarchy& bvh, Ray ray) {
    // introduce depth to the getFinalColor.
    return getFinalColor(scene, bvh, ray, 5);
}

static void setOpenGLMatrices(const Trackball& camera);
static void drawLightsOpenGL(const Scene& scene, const Trackball& camera, int selectedLight);
static void drawSceneOpenGL(const Scene& scene);

// This is the main rendering function. You are free to change this function in any way (including the function signature).
static void renderRayTracing(const Scene& scene, const Trackball& camera, const BoundingVolumeHierarchy& bvh, Screen& screen)
{
#ifndef NDEBUG
    // Single threaded in debug mode
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos{
                float(x) / windowResolution.x * 2.0f - 1.0f,
                float(y) / windowResolution.y * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay));
        }
    }
#else
    // Multi-threaded in release mode
    const tbb::blocked_range2d<int, int> windowRange{ 0, windowResolution.y, 0, windowResolution.x };
    tbb::parallel_for(windowRange, [&](tbb::blocked_range2d<int, int> localRange) {
        for (int y = std::begin(localRange.rows()); y != std::end(localRange.rows()); y++) {
            for (int x = std::begin(localRange.cols()); x != std::end(localRange.cols()); x++) {
                // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
                const glm::vec2 normalizedPixelPos{
                    float(x) / windowResolution.x * 2.0f - 1.0f,
                    float(y) / windowResolution.y * 2.0f - 1.0f
                };
                const Ray cameraRay = camera.generateRay(normalizedPixelPos);
                screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay));
            }
        }
        });
#endif
}

int main(int argc, char** argv)
{
    Trackball::printHelp();
    std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
        << std::endl;

    Window window{ "Final Project", windowResolution, OpenGLVersion::GL2 };
    Screen screen{ windowResolution };
    Trackball camera{ &window, glm::radians(50.0f), 3.0f };
    camera.setCamera(glm::vec3(0.0f, 0.0f, 0.0f), glm::radians(glm::vec3(20.0f, 20.0f, 0.0f)), 3.0f);

    SceneType sceneType{ SceneType::SingleTriangle };
    std::optional<Ray> optDebugRay;
    Scene scene = loadScene(sceneType, dataPath);
    BoundingVolumeHierarchy bvh{ &scene };

    int screen_threshold = 68;
    int screen_scale = 19;
    int bvhDebugLevel = 0;
    bool debugBVH{ false };
    ViewMode viewMode{ ViewMode::Rasterization };

    window.registerKeyCallback([&](int key, int /* scancode */, int action, int /* mods */) {
        if (action == GLFW_PRESS) {
            switch (key) {
            case GLFW_KEY_R: {
                // Shoot a ray. Produce a ray from camera to the far plane.
                const auto tmp = window.getNormalizedCursorPos();
                optDebugRay = camera.generateRay(tmp * 2.0f - 1.0f);
            } break;
            case GLFW_KEY_ESCAPE: {
                window.close();
            } break;
            };
        }
        });

    int selectedLightIdx = scene.lights.empty() ? -1 : 0;
    int index = 0;

    while (!window.shouldClose()) {
        window.updateInput();

        // === Setup the UI ===
        ImGui::Begin("Final Project");
        {
            constexpr std::array items{ "SingleTriangle", "Cube (segment light)", "Cornell Box (with mirror)", "Cornell Box (parallelogram light and mirror)", "Monkey", "Teapot", "Dragon", /* "AABBs",*/ "Spheres", /*"Mixed",*/ "Custom", "Texture Quad" };
            if (ImGui::Combo("Scenes", reinterpret_cast<int*>(&sceneType), items.data(), int(items.size()))) {
                optDebugRay.reset();
                scene = loadScene(sceneType, dataPath);
                selectedLightIdx = scene.lights.empty() ? -1 : 0;
                bvh = BoundingVolumeHierarchy(&scene);
                if (optDebugRay) {
                    HitInfo dummy{};
                    bvh.intersect(*optDebugRay, dummy);
                }
            }
        }
        {
            constexpr std::array items{ "Rasterization", "Ray Traced" };
            ImGui::Combo("View mode", reinterpret_cast<int*>(&viewMode), items.data(), int(items.size()));
        }
        if (ImGui::Button("Render to file")) {
            // Show a file picker.
            nfdchar_t* pOutPath = nullptr;
            const nfdresult_t result = NFD_SaveDialog("bmp", nullptr, &pOutPath);
            if (result == NFD_OKAY) {
                std::filesystem::path outPath{ pOutPath };
                free(pOutPath); // NFD is a C API so we have to manually free the memory it allocated.
                outPath.replace_extension("bmp"); // Make sure that the file extension is *.bmp

                // Perform a new render and measure the time it took to generate the image.
                using clock = std::chrono::high_resolution_clock;
                const auto start = clock::now();
                renderRayTracing(scene, camera, bvh, screen);
                screen.bloom_filter();
                const auto end = clock::now();
                std::cout << "Time to render image: " << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds" << std::endl;

                // Store the new image.
                screen.writeBitmapToFile(outPath);
            }
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Debugging");
        if (viewMode == ViewMode::Rasterization) {
            ImGui::Checkbox("Draw BVH", &debugBVH);
            if (debugBVH)
                ImGui::SliderInt("BVH Level", &bvhDebugLevel, 0, bvh.maxLevel);
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("reflection");
        if (ImGui::Button("Switch")) {          
            recursive = !recursive;
        }
        if (recursive) {
            ImGui::Text("Recursive reflection");

        }
        else {
            ImGui::Text("Glossy reflection");
            ImGui::SliderInt("Number of rays", &GR_number_rays, 5, 50);
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("depth of field");
        if (ImGui::Button("Switch depth of field")) {
            
            dof_off = !dof_off;
        }
        if (dof_off) {
            ImGui::Text("off");

        }
        else {
            ImGui::Text("on");
        }
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Bloom filter");
        ImGui::SliderInt("Threshold (* 0.01)", &screen_threshold, 60, 100);
        ImGui::SliderInt("Scale (* 0.1)", &screen_scale, 12, 21);

        screen.set_brightness_threshold(((float)screen_threshold) / 100.0f);
        screen.set_scale(((float)screen_scale) / 10.0f);

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Text("Lights");
        {
            std::vector<std::string> options;
            options.push_back("None");
            for (size_t i = 0; i < scene.lights.size(); i++) {
                options.push_back("Light " + std::to_string(i));
            }
            std::vector<const char*> optionsPointers;
            std::transform(std::begin(options), std::end(options), std::back_inserter(optionsPointers), [](const auto& str) { return str.c_str(); });

            // Offset such that selectedLightIdx=-1 becomes item 0 (None).
            ++selectedLightIdx;
            ImGui::Combo("Selected light", &selectedLightIdx, optionsPointers.data(), static_cast<int>(optionsPointers.size()));
            --selectedLightIdx;

            if (selectedLightIdx >= 0) {
                setOpenGLMatrices(camera);
                std::visit(
                    make_visitor(
                        [&](PointLight& light) {
                            showImGuizmoTranslation(window, camera, light.position); // 3D controls to translate light source.
                            ImGui::DragFloat3("Light position", glm::value_ptr(light.position), 0.01f, -3.0f, 3.0f);
                            ImGui::ColorEdit3("Light color", glm::value_ptr(light.color));
                        },
                        [&](SegmentLight& light) {
                            static int selectedEndpoint = 0;
                            // 3D controls to translate light source.
                            if (selectedEndpoint == 0)
                                showImGuizmoTranslation(window, camera, light.endpoint0);
                            else
                                showImGuizmoTranslation(window, camera, light.endpoint1);

                            const std::array<const char*, 2> endpointOptions{ "Endpoint 0", "Endpoint 1" };
                            ImGui::Combo("Selected endpoint", &selectedEndpoint, endpointOptions.data(), (int)endpointOptions.size());
                            ImGui::DragFloat3("Endpoint 0", glm::value_ptr(light.endpoint0), 0.01f, -3.0f, 3.0f);
                            ImGui::DragFloat3("Endpoint 1", glm::value_ptr(light.endpoint1), 0.01f, -3.0f, 3.0f);
                            ImGui::ColorEdit3("Color 0", glm::value_ptr(light.color0));
                            ImGui::ColorEdit3("Color 1", glm::value_ptr(light.color1));
                        },
                            [&](ParallelogramLight& light) {
                            glm::vec3 vertex1 = light.v0 + light.edge01;
                            glm::vec3 vertex2 = light.v0 + light.edge02;

                            static int selectedVertex = 0;
                            // 3D controls to translate light source.
                            if (selectedVertex == 0)
                                showImGuizmoTranslation(window, camera, light.v0);
                            else if (selectedVertex == 1)
                                showImGuizmoTranslation(window, camera, vertex1);
                            else
                                showImGuizmoTranslation(window, camera, vertex2);

                            const std::array<const char*, 3> vertexOptions{ "Vertex 0", "Vertex 1", "Vertex 2" };
                            ImGui::Combo("Selected vertex", &selectedVertex, vertexOptions.data(), (int)vertexOptions.size());
                            ImGui::DragFloat3("Vertex 0", glm::value_ptr(light.v0), 0.01f, -3.0f, 3.0f);
                            ImGui::DragFloat3("Vertex 1", glm::value_ptr(vertex1), 0.01f, -3.0f, 3.0f);
                            light.edge01 = vertex1 - light.v0;
                            ImGui::DragFloat3("Vertex 2", glm::value_ptr(vertex2), 0.01f, -3.0f, 3.0f);
                            light.edge02 = vertex2 - light.v0;

                            ImGui::ColorEdit3("Color 0", glm::value_ptr(light.color0));
                            ImGui::ColorEdit3("Color 1", glm::value_ptr(light.color1));
                            ImGui::ColorEdit3("Color 2", glm::value_ptr(light.color2));
                            ImGui::ColorEdit3("Color 3", glm::value_ptr(light.color3));
                        },
                            [](auto) { /* any other type of light */ }),
                    scene.lights[selectedLightIdx]);
            }
        }

        if (ImGui::Button("Add point light")) {
            selectedLightIdx = int(scene.lights.size());
            scene.lights.push_back(PointLight{ .position = glm::vec3(0.0f), .color = glm::vec3(1.0f) });
        }
        if (ImGui::Button("Add segment light")) {
            selectedLightIdx = int(scene.lights.size());
            scene.lights.push_back(SegmentLight{ .endpoint0 = glm::vec3(0.0f), .endpoint1 = glm::vec3(1.0f), .color0 = glm::vec3(1, 0, 0), .color1 = glm::vec3(0, 0, 1) });
        }
        if (ImGui::Button("Add parallelogram light")) {
            selectedLightIdx = int(scene.lights.size());
            scene.lights.push_back(ParallelogramLight{
                .v0 = glm::vec3(0.0f),
                .edge01 = glm::vec3(1, 0, 0),
                .edge02 = glm::vec3(0, 1, 0),
                .color0 = glm::vec3(1, 0, 0), // red
                .color1 = glm::vec3(0, 1, 0), // green
                .color2 = glm::vec3(0, 0, 1), // blue
                .color3 = glm::vec3(1, 1, 1) // white
                });
        }
        if (selectedLightIdx >= 0 && ImGui::Button("Remove selected light")) {
            scene.lights.erase(std::begin(scene.lights) + selectedLightIdx);
            selectedLightIdx = -1;
        }

        // Clear screen.
        glViewport(0, 0, window.getFrameBufferSize().x, window.getFrameBufferSize().y);
        glClearDepth(1.0f);
        glClearColor(0.0, 0.0, 0.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        setOpenGLMatrices(camera);

        // Draw either using OpenGL (rasterization) or the ray tracing function.
        switch (viewMode) {
        case ViewMode::Rasterization: {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            drawSceneOpenGL(scene);
           
            if (optDebugRay) {
                // Call getFinalColor for the debug ray. Ignore the result but tell the function that it should
                // draw the rays instead.
                enableDrawRay = true;
                glDepthFunc(GL_LEQUAL);
                (void)getFinalColor(scene, bvh, *optDebugRay);
                enableDrawRay = false;
            }
            glPopAttrib();
        } break;
        case ViewMode::RayTracing: {
            screen.clear(glm::vec3(0.0f));
            renderRayTracing(scene, camera, bvh, screen);
            screen.setPixel(0, 0, glm::vec3(1.0f));
            screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
        } break;
        default:
            break;
        };

        drawLightsOpenGL(scene, camera, selectedLightIdx);

        if (debugBVH) {
            glPushAttrib(GL_ALL_ATTRIB_BITS);
            setOpenGLMatrices(camera);
            glDisable(GL_LIGHTING);
            glEnable(GL_DEPTH_TEST);
            glDepthFunc(GL_LEQUAL);

            // Enable alpha blending. More info at:
            // https://learnopengl.com/Advanced-OpenGL/Blending
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            enableDrawRay = true;
            bvh.debugDraw(bvhDebugLevel);
            enableDrawRay = false;
            glPopAttrib();
        }

        ImGui::End();
        window.swapBuffers();
    }

    return 0;
}

static void setOpenGLMatrices(const Trackball& camera)
{
    // Load view matrix.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const glm::mat4 viewMatrix = camera.viewMatrix();
    glMultMatrixf(glm::value_ptr(viewMatrix));

    // Load projection matrix.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const glm::mat4 projectionMatrix = camera.projectionMatrix();
    glMultMatrixf(glm::value_ptr(projectionMatrix));
}

static void drawLightsOpenGL(const Scene& scene, const Trackball& camera, int selectedLight)
{
    // Normals will be normalized in the graphics pipeline.
    glEnable(GL_NORMALIZE);
    // Activate rendering modes.
    glEnable(GL_DEPTH_TEST);
    // Draw front and back facing triangles filled.
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    // Interpolate vertex colors over the triangles.
    glShadeModel(GL_SMOOTH);

    glDisable(GL_LIGHTING);
    // Draw all non-selected lights.
    for (size_t i = 0; i < scene.lights.size(); i++) {
        std::visit(
            make_visitor(
                [](const PointLight& light) { drawSphere(light.position, 0.01f, light.color); },
                [](const SegmentLight& light) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glBegin(GL_LINES);
                    glColor3fv(glm::value_ptr(light.color0));
                    glVertex3fv(glm::value_ptr(light.endpoint0));
                    glColor3fv(glm::value_ptr(light.color1));
                    glVertex3fv(glm::value_ptr(light.endpoint1));
                    glEnd();
                    glPopAttrib();
                    drawSphere(light.endpoint0, 0.01f, light.color0);
                    drawSphere(light.endpoint1, 0.01f, light.color1);
                },
                [](const ParallelogramLight& light) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glBegin(GL_QUADS);
                    glColor3fv(glm::value_ptr(light.color0));
                    glVertex3fv(glm::value_ptr(light.v0));
                    glColor3fv(glm::value_ptr(light.color1));
                    glVertex3fv(glm::value_ptr(light.v0 + light.edge01));
                    glColor3fv(glm::value_ptr(light.color3));
                    glVertex3fv(glm::value_ptr(light.v0 + light.edge01 + light.edge02));
                    glColor3fv(glm::value_ptr(light.color2));
                    glVertex3fv(glm::value_ptr(light.v0 + light.edge02));
                    glEnd();
                    glPopAttrib();
                },
                    [](auto) { /* any other type of light */ }),
            scene.lights[i]);
    }

    // Draw a colored sphere at the location at which the trackball is looking/rotating around.
    glDisable(GL_LIGHTING);
    drawSphere(camera.lookAt(), 0.01f, glm::vec3(0.2f, 0.2f, 1.0f));
}

void drawSceneOpenGL(const Scene& scene)
{
    // Activate the light in the legacy OpenGL mode.
    glEnable(GL_LIGHTING);

    // Tell OpenGL where the lights are (so it nows how to shade surfaces in the scene).
    // This is only used in the rasterization view. OpenGL only supports point lights so
    // we replace segment/parallelogram lights by point lights.
    int i = 0;
    const auto enableLight = [&](const glm::vec3& position, const glm::vec3 color) {
        glEnable(GL_LIGHT0 + i);
        const glm::vec4 position4{ position, 1 };
        glLightfv(GL_LIGHT0 + i, GL_POSITION, glm::value_ptr(position4));
        const glm::vec4 color4{ glm::clamp(color, 0.0f, 1.0f), 1.0f };
        const glm::vec4 zero4{ 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0 + i, GL_AMBIENT, glm::value_ptr(zero4));
        glLightfv(GL_LIGHT0 + i, GL_DIFFUSE, glm::value_ptr(color4));
        glLightfv(GL_LIGHT0 + i, GL_SPECULAR, glm::value_ptr(zero4));
        // NOTE: quadratic attenuation doesn't work like you think it would in legacy OpenGL.
        // The distance is not in world space but in NDC space!
        glLightf(GL_LIGHT0 + i, GL_CONSTANT_ATTENUATION, 1.0f);
        glLightf(GL_LIGHT0 + i, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(GL_LIGHT0 + i, GL_QUADRATIC_ATTENUATION, 0.0f);
        i++;
    };
    for (const auto& light : scene.lights) {
        std::visit(
            make_visitor(
                [&](const PointLight& light) {
                    enableLight(light.position, light.color);
                },
                [&](const SegmentLight& light) {
                    // Approximate with two point lights: one at each endpoint.
                    enableLight(light.endpoint0, 0.5f * light.color0);
                    enableLight(light.endpoint1, 0.5f * light.color1);
                },
                    [&](const ParallelogramLight& light) {
                    enableLight(light.v0, 0.25f * light.color0);
                    enableLight(light.v0 + light.edge01, 0.25f * light.color1);
                    enableLight(light.v0 + light.edge02, 0.25f * light.color2);
                    enableLight(light.v0 + light.edge01 + light.edge02, 0.25f * light.color3);
                },
                    [](auto) { /* any other type of light */ }),
            light);
    }

    // Draw the scene and the ray (if any).
    drawScene(scene);
}
