#include "MACGrid2D.h"
#include "Renderer2D.h"

// Test
int main()
{
    MACGrid2D *grid         = new MACGrid2D(128, 0.01f);
    Solver *solver          = new Solver(grid, grid->getRes(), grid->getCellSize(), 0.03f);
    Renderer2D *renderer    = new Renderer2D(800, 800, grid, solver);
    
    renderer->initRenderer();
    renderer->initShader();
    GLFWwindow* window = renderer->getWindow();
    
    while (!glfwWindowShouldClose(window))
    {
        renderer->render();
    }

    delete renderer;
    delete grid;
    return 0;
}