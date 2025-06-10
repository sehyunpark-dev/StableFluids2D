#include "MACGrid2D.h"
#include "Renderer2D.h"

// Test
int main()
{
    MACGrid2D *grid     = new MACGrid2D(128, 0.01);
    Renderer2D renderer = Renderer2D(800, 600, grid);
    
    renderer.initialize();
    GLFWwindow* window = renderer.getWindow();
    
    while (!glfwWindowShouldClose(window))
    {
        renderer.render();
    }

    renderer.cleanup();
    delete grid;
    return 0;
}