#include "utils.h"

#include <iostream>

GPU_Geometry::GPU_Geometry() : vao_(0), vbo_position_(0), vbo_color_(0)
{}

GPU_Geometry::~GPU_Geometry()
{
    // Delete VAO and VBOs
    glDeleteBuffers(1, &vbo_position_);
    glDeleteBuffers(1, &vbo_color_);
    glDeleteVertexArrays(1, &vao_);
}

void GPU_Geometry::initOpenGLResources()
{
    // Initialize OpenGL buffers
    glGenVertexArrays(1, &vao_);
    bindVAO();

    // Create VBOs for positions and colors
    glGenBuffers(1, &vbo_position_);
    glGenBuffers(1, &vbo_color_);
}

void GPU_Geometry::uploadGeometry(const CPU_Geometry& cpu_geom)
{
    // Upload position data to GPU
    glBindBuffer(GL_ARRAY_BUFFER, vbo_position_);
    glBufferData(
        GL_ARRAY_BUFFER, 
        cpu_geom.position.size() * sizeof(glm::vec3), 
        cpu_geom.position.data(), 
        GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), nullptr);
    glEnableVertexAttribArray(0);

    // Upload color data to GPU
    glBindBuffer(GL_ARRAY_BUFFER, vbo_color_);
    glBufferData(
        GL_ARRAY_BUFFER, 
        cpu_geom.color.size() * sizeof(glm::vec3), 
        cpu_geom.color.data(), 
        GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), nullptr);
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    unbindVAO();
}

void GPU_Geometry::bindVAO() const
{
    glBindVertexArray(vao_);
}

void GPU_Geometry::unbindVAO() const
{
    glBindVertexArray(0);
}