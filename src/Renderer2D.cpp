#include "Renderer2D.h"
#include <iostream>

Renderer2D::Renderer2D(int width, int height, MACGrid2D* grid) : 
    window_width_(width), window_height_(height), grid_(grid), window_(nullptr) 
{}

Renderer2D::~Renderer2D()
{
    cleanup();
}

void Renderer2D::initialize()
{
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        exit(-1);
    }

    window_ = glfwCreateWindow(window_width_, window_height_, "Stable Fluids 2D", nullptr, nullptr);
    if (!window_)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        exit(-1);
    }

    glfwMakeContextCurrent(window_);

    if (glewInit() != GLEW_OK)
    {
        std::cerr << "Failed to initialize GLEW" << std::endl;
        exit(-1);
    }

    // Initialize OpenGL resources (VAO, VBOs, etc.)
    grid_gpu_geom_.initOpenGLResources();

    // Check whether the grid pointer is valid
    if (!grid_)
    {
        std::cerr << "Error: MACGrid2D pointer is null in Renderer2D::initialize()" << std::endl;
        glfwDestroyWindow(window_);
        glfwTerminate();
        exit(-1);
    }

    // Prepare CPU geometry data from the grid
    const std::vector<glm::vec2> &cell_coords = grid_->getCellCoord();
    grid_cpu_geom_.position.clear();
    grid_cpu_geom_.color.clear();
    for (const glm::vec2 &coord : cell_coords)
    {
        grid_cpu_geom_.position.push_back(glm::vec3(coord, 0.0f)); // z = 0 for 2D
        grid_cpu_geom_.color.push_back(glm::vec3(0.0f, 0.5f, 1.0f));
    }

    // Upload CPU geometry data to GPU
    grid_gpu_geom_.uploadGeometry(grid_cpu_geom_);
}

void Renderer2D::render()
{
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT);

    // Bind the VAO and draw the grid, and then unbind it
    grid_gpu_geom_.bindVAO();
    glDrawArrays(GL_POINTS, 0, grid_cpu_geom_.position.size());
    grid_gpu_geom_.unbindVAO();

    glfwSwapBuffers(window_);
    glfwPollEvents();
}

void Renderer2D::cleanup()
{
    if (window_)
    {
        glfwDestroyWindow(window_);
        window_ = nullptr;
    
        glfwTerminate();
    }
}

GLFWwindow* Renderer2D::getWindow()
{
    return window_;
}