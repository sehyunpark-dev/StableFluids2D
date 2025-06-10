#ifndef MAC_GRID_2D_H
#define MAC_GRID_2D_H

#include "glm/glm.hpp"

#include <iostream>
#include <vector>

class MACGrid2D
{
private:
    int res_            = 128; // Resolution of the grid in x and y directions (e.g. 128x128)
    float cell_size_    = 0.01; // Size of each cell in the grid (e.g. 0.01))

    glm::vec2 grid_center_ = {0.0f, 0.0f};

    std::vector<glm::vec2> cell_coord_; // Coordinates of each cell center (size : height * width)
    std::vector<double> u_;             // x-component of velocity  (size : height * (width + 1)))
    std::vector<double> v_;             // y-component of velocity  (size : (height + 1) * width))
    std::vector<double> p_;             // pressure of each cell    (size : height * width)

public:
    MACGrid2D(int res, float cell_size);

    int getRes() const;
    float getCellSize() const;
    glm::vec2 getGridCenter() const;
    std::vector<glm::vec2> getCellCoord() const;
    glm::vec2 getCellCoord(int i, int j) const;

    double getU(int i, int j) const;
    double getV(int i, int j) const;
    double getP(int i, int j) const;

    void setU(int i, int j, double val);
    void setV(int i, int j, double val);
    void setP(int i, int j, double val);
};

#endif