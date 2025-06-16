#include "MACGrid2D.h"
#include "Renderer2D.h"
#include <iostream>
#include <string>
#include <sstream>

// Test
int main()
{
    MACGrid2D *grid         = new MACGrid2D(128, 0.01f);
    Solver *solver          = new Solver(grid, grid->getRes(), grid->getCellSize(), 0.03f);
    Renderer2D *renderer    = new Renderer2D(800, 800, grid, solver);
    
    renderer->initRenderer();
    renderer->initShader();
    GLFWwindow* window = renderer->getWindow();

    bool isSpacePressedLastFrame = false;
    bool isRPressedLastFrame = false;

    double lastTime = glfwGetTime();
    int nbFrames = 0;
    
    while (!glfwWindowShouldClose(window))
    {
        double currentTime = glfwGetTime();
        nbFrames++;

        if (currentTime - lastTime >= 1.0)
        {
            double msPerFrame = 1000.0 / double(nbFrames);
            double fps = double(nbFrames) / (currentTime - lastTime);

            std::stringstream ss;
            ss << "2D Fluid Solver | " << msPerFrame << " ms/frame (" << fps << " FPS)";
            
            glfwSetWindowTitle(window, ss.str().c_str());

            nbFrames = 0;
            lastTime = currentTime;
        }

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