#include "MACGrid2D.h"

MACGrid2D::MACGrid2D(int res, float cell_size) : res_(res), cell_size_(cell_size)
{
    // Initialize the grid parameters
    cell_coord_.resize(res_ * res_, glm::vec2(0.0f, 0.0f));
    
    float center_idx = (res_ - 1) / 2.0f;
    for (int i = 0; i < res_; i++)
    {
        for (int j = 0; j < res_; j++)
        {
            float y = (center_idx - i) * cell_size;
            float x = (j - center_idx) * cell_size;
            cell_coord_[i * res_ + j] = glm::vec2(x, y);
        }
    }
}

int MACGrid2D::getRes() const
{
    return res_;
}

float MACGrid2D::getCellSize() const
{
    return cell_size_;
}

glm::vec2 MACGrid2D::getGridCenter() const
{
    return grid_center_;
}

std::vector<glm::vec2> MACGrid2D::getCellCoord() const
{
    return cell_coord_;
}

glm::vec2 MACGrid2D::getCellCoord(int x, int y) const
{
    if (0 <= y && y < res_ && 0 <= x && x < res_)
    {
        return cell_coord_[y * res_ + x];
    }
    else
    {
        std::cerr << "[Error] Invalid cell coordinate index" << std::endl;
        return glm::vec2(-1.0f, -1.0f);
    }
}

int MACGrid2D::getCellIndex(const glm::vec2 &coord) const
{
    glm::vec2 temp_coord = coord + glm::vec2(res_ / 2.0f * cell_size_, res_ / 2.0f * cell_size_);
    
    int x = static_cast<int>(temp_coord.x / cell_size_);
    int y = static_cast<int>(temp_coord.y / cell_size_);

    if (0 <= y && y < res_ && 0 <= x && x < res_)
    {
        return y * res_ + x;
    }
    else
    {
        std::cerr << "[Error] Invalid cell coordinate" << std::endl;
        return -1; // Invalid index
    }
}