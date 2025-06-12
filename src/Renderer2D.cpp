#include "Renderer2D.h"
#include <iostream>
#include <fstream>
#include <sstream>

Renderer2D::Renderer2D(int width, int height, MACGrid2D* grid, Solver *solver) : 
    window_width_(width), window_height_(height), 
    grid_(grid), solver_(solver), window_(nullptr) 
{}

Renderer2D::~Renderer2D()
{
    cleanup();
}

void Renderer2D::initRenderer()
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

    // Set the background color
    glClearColor(0.2f, 0.2f, 0.2f, 1.0f);

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
    grid_cpu_geom_.position.clear();
    grid_cpu_geom_.color.clear();
    const std::vector<glm::vec2>& cell_coords = grid_->getCellCoord();
    for (const glm::vec2 &coord : cell_coords)
    {
        grid_cpu_geom_.position.push_back(glm::vec3(coord, 0.0f)); // z = 0 for 2D
        grid_cpu_geom_.color.push_back(glm::vec3(0.0f, 0.5f, 1.0f));
    }

    // Upload CPU geometry data to GPU
    grid_gpu_geom_.uploadGeometry(grid_cpu_geom_);
}

void Renderer2D::initShader()
{
    // Read shader files
    std::string vertex_shader_src = readShaderFile("../shader/vertex_shader.vert");
    std::string fragment_shader_src = readShaderFile("../shader/fragment_shader.frag");

    // Compile Vertex Shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    const char* vertex_src = vertex_shader_src.c_str();
    glShaderSource(vertex_shader, 1, &vertex_src, nullptr);
    glCompileShader(vertex_shader);

    GLint success;

    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char info_log[512];
        glGetShaderInfoLog(vertex_shader, 512, nullptr, info_log);
        std::cerr << "Vertex Shader Compilation Failed:\n" << info_log << std::endl;
        exit(-1);
    }

    // Compile Fragment Shader
    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    const char* fragment_src = fragment_shader_src.c_str();
    glShaderSource(fragment_shader, 1, &fragment_src, nullptr);
    glCompileShader(fragment_shader);

    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char info_log[512];
        glGetShaderInfoLog(fragment_shader, 512, nullptr, info_log);
        std::cerr << "Fragment Shader Compilation Failed:\n" << info_log << std::endl;
        exit(-1);
    }

    // Link Shaders into a Program
    shader_program_ = glCreateProgram();
    glAttachShader(shader_program_, vertex_shader);
    glAttachShader(shader_program_, fragment_shader);
    glLinkProgram(shader_program_);
    
    glGetProgramiv(shader_program_, GL_LINK_STATUS, &success);
    if (!success)
    {
        char info_log[512];
        glGetProgramInfoLog(shader_program_, 512, nullptr, info_log);
        std::cerr << "Shader Program Linking Failed:\n" << info_log << std::endl;
        exit(-1);
    }

    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);
}

void Renderer2D::render()
{
    // Clear the screen
    glClear(GL_COLOR_BUFFER_BIT);

    const std::vector<float> &smoke_density     = solver_->getSmokeDensityVector();
    const std::vector<glm::vec2> &cell_coords   = grid_->getCellCoord();
    float cell_size = 0.01;

    // for (size_t i = 0; i < smoke_density.size(); ++i)
    // {
    //     if (smoke_density[i] != 0.0f)
    //     {
    //         std::cout << "Smoke Density[" << i << "]: " << smoke_density[i] << std::endl;
    //     }
    // }

    grid_cpu_geom_.position.clear();
    grid_cpu_geom_.color.clear();

    for (int i = 0; i < cell_coords.size(); i++)
    {
        glm::vec2 coord = cell_coords[i];
        float density   = smoke_density[i];
        glm::vec3 color = glm::vec3(density, density, density); // Grayscale based on density

        grid_cpu_geom_.position.push_back(glm::vec3(coord.x - cell_size / 2, coord.y - cell_size / 2, 0.0f));
        grid_cpu_geom_.position.push_back(glm::vec3(coord.x + cell_size / 2, coord.y - cell_size / 2, 0.0f));
        grid_cpu_geom_.position.push_back(glm::vec3(coord.x + cell_size / 2, coord.y + cell_size / 2, 0.0f));
        grid_cpu_geom_.position.push_back(glm::vec3(coord.x - cell_size / 2, coord.y + cell_size / 2, 0.0f));

        for (int j = 0; j < 4; j++)
        {
            grid_cpu_geom_.color.push_back(color);
        }
    }

    grid_gpu_geom_.uploadGeometry(grid_cpu_geom_);

    glUseProgram(shader_program_);

    // Bind the VAO and draw the grid, and then unbind it
    grid_gpu_geom_.bindVAO();
    glDrawArrays(GL_QUADS, 0, grid_cpu_geom_.position.size());
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

std::string Renderer2D::readShaderFile(const std::string &file_path)
{
    std::ifstream shader_file(file_path);
    if (!shader_file.is_open())
    {
        std::cerr << "Failed to open shader file: " << file_path << std::endl;
        exit(-1);
    }

    std::stringstream shader_stream;
    shader_stream << shader_file.rdbuf();
    shader_file.close();

    return shader_stream.str();
}