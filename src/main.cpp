#include <GLFW/glfw3.h>
#include <iostream>

// Test
int main()
{
    // 1. Initialize GLFW
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    // 2. Create a windowed mode window and its OpenGL context
    GLFWwindow* window = glfwCreateWindow(800, 600, "Stable Fluids 2D", nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    // 3. Make the window's context current
    glfwMakeContextCurrent(window);

    // 4. Main loop
    while (!glfwWindowShouldClose(window))
    {
        // Render here

        // Swap buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    // 5. Clean up and exit
    glfwTerminate();
    return 0;
}