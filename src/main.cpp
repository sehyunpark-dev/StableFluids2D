#include "MACGrid2D.h"
#include "Renderer2D.h"

// Test
int main()
{
    MACGrid2D *grid         = new MACGrid2D(128, 0.01);
    Renderer2D *renderer    = new Renderer2D(1024, 768, grid);
    
    renderer->initialize();
    GLFWwindow* window = renderer->getWindow();
    
    while (!glfwWindowShouldClose(window))
    {
        renderer->render();
    }

    delete renderer;
    delete grid;
    return 0;
}