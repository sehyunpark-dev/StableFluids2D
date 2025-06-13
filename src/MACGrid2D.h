#ifndef MAC_GRID_2D_H
#define MAC_GRID_2D_H

#include "glm/glm.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>

class MACGrid2D
{
private:
    // The grid is considered as a square grid, which means height = width.
    int res_            = 128;  // Resolution of the grid in x and y directions (e.g. 128x128)
    float cell_size_    = 0.01f; // Size of each cell in the grid (e.g. 0.01))

    glm::vec2 grid_center_ = {0.0f, 0.0f};
    glm::vec2 grid_min_ = {0.0f, 0.0f};
    glm::vec2 grid_max_ = {0.0f, 0.0f};

    std::vector<glm::vec2> cell_coord_; // Coordinates of each cell center (size : height * width)
    std::vector<glm::vec2> u_coord_;    // Coordinates of each u face center (size : height * (width + 1))
    std::vector<glm::vec2> v_coord_;    // Coordinates of each u face center (size : (height + 1) * width)

public:
    MACGrid2D(int res, float cell_size);
    ~MACGrid2D() = default;

    // Getters for grid properties

    int getRes() const;
    float getCellSize() const;
    glm::vec2 getGridCenter() const;
    
    // Get the coordinate of all cell centers
    std::vector<glm::vec2> getCellCoord() const;

    // Get the coordinate of a specific cell center at given indices (x, y)
    glm::vec2 getCellCoord(int x, int y) const;

    // Get the coordinate of all u grid face centers
    std::vector<glm::vec2> getUCoord() const;

    // Get the coordinate of a specific u grid face center at given indices (x, y)
    glm::vec2 getUCoord(int x, int y) const;

    // Get the coordinate of all v grid face centers
    std::vector<glm::vec2> getVCoord() const;

    // Get the coordinate of a specific v grid face center at given indices (x, y)
    glm::vec2 getVCoord(int x, int y) const;
    
    // Get the index of the cell that contains the given coordinate
    int getCellIndex(const glm::vec2 &cell_coord) const;

    // Get the index of the u grid that contains the given coordinate
    int getUIndex(const glm::vec2 &u_coord) const;

    // Get the index of the v grid that contains the given coordinate
    int getVIndex(const glm::vec2 &v_coord) const;
};

#endif