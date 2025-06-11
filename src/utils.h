#ifndef UTILS_H
#define UTILS_H

#include "GL/glew.h"
#include <glm/glm.hpp>
#include <vector>

struct CPU_Geometry
{
    std::vector<glm::vec3> position;
    std::vector<glm::vec3> color;
};

class GPU_Geometry
{
private:
    GLuint vao_;
    GLuint vbo_position_;
    GLuint vbo_color_;

public:
    GPU_Geometry();
    ~GPU_Geometry();

    void initOpenGLResources(); // Initialize OpenGL resources

    void uploadGeometry(const CPU_Geometry& cpu_geom);  // Upload CPU geometry data to GPU
    void bindVAO() const;   // Bind the VAO for rendering
    void unbindVAO() const; // Unbind the VAO
};

#endif