#pragma once
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <vector>

class Screen {
public:
    Screen(const glm::ivec2& resolution);

    void clear(const glm::vec3& color);
    void setPixel(int x, int y, const glm::vec3& color);
    glm::vec3 getPixel(int x, int y);
    void bloom_filter();
    glm::vec3 filter(int i, int j);
    void writeBitmapToFile(const std::filesystem::path& filePath);
    void draw();
    float get_brightness_threshold();
    void set_brightness_threshold(float threshold);
    
    float get_scale();
    void set_scale(float scale);


private:
    glm::ivec2 m_resolution;
    std::vector<glm::vec3> m_textureData;

    float brightness_threshold = 0.9;
    float scale = 1.5;
    uint32_t m_texture;
};
