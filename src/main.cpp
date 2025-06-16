#include "MACGrid2D.h"
#include "Renderer2D.h"

// Test
int main()
{
    MACGrid2D *grid         = new MACGrid2D(256, 0.005f);
    Solver *solver          = new Solver(grid, grid->getRes(), grid->getCellSize(), 0.03f);
    Renderer2D *renderer    = new Renderer2D(800, 800, grid, solver);
    
    renderer->initRenderer();
    renderer->initShader();
    GLFWwindow* window = renderer->getWindow();

    bool isSpacePressedLastFrame = false;
    bool isRPressedLastFrame = false;
    
    while (!glfwWindowShouldClose(window))
    {
        bool isSpacePressedNow = (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS);
        bool isRPressedNow = (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS);

        if (isSpacePressedNow && !isSpacePressedLastFrame)
        {
            if (solver->getIsSimulating()) 
            {
                solver->deactiveSimulation();
            }
            else
            {
                solver->activeSimulation();
            }
        }
        if (isRPressedNow && !isRPressedLastFrame)
        {
            solver->reset();
        }

        isSpacePressedLastFrame = isSpacePressedNow;
        isRPressedLastFrame = isRPressedNow;

        if (solver->getIsSimulating())
        {
            solver->step();
        }

        renderer->renderSmoke();
    }

    delete renderer;
    delete grid;
    return 0;
}