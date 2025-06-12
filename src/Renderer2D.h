#ifndef RENDERER_2D_H
#define RENDERER_2D_H

#include "MACGrid2D.h"
#include "Solver.h"
#include "utils.h"
#include "GLFW/glfw3.h"
#include <string>

class Renderer2D
{
private:
    GLFWwindow *window_; // Pointer to the GLFW window

    int window_width_   = 800;  // Width of the rendering window
    int window_height_  = 600;  // Height of the rendering window

    MACGrid2D *grid_;           // Pointer to the MACGrid2D instance
    Solver  *solver_;           // Pointer to the Solver instance
    
    CPU_Geometry grid_cpu_geom_;    // Grid CPU geometry data
    GPU_Geometry grid_gpu_geom_;    // Grid GPU geometry data
    GLuint shader_program_;         // OpenGL shader program ID

    std::string readShaderFile(const std::string &file_path);

public:
    Renderer2D(int width, int height, MACGrid2D *grid, Solver *glNamedRenderbufferStorageMultisampleCoverageEXT);
    ~Renderer2D();

    void initRenderer();
    void initShader();
    void render();
    void cleanup();

    // Getter for the GLFW window
    GLFWwindow* getWindow();
};

#endif