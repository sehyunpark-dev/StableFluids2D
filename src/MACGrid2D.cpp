#include "MACGrid2D.h"
#include <cmath>

MACGrid2D::MACGrid2D(int res, float cell_size) : res_(res), cell_size_(cell_size)
{
    // Initialize the grid parameters
    cell_coord_.resize(res_ * res_, glm::vec2(0.0f, 0.0f));
    u_coord_.resize(res_ * (res_ + 1), glm::vec2(0.0f, 0.0f));
    v_coord_.resize((res_ + 1) * res_, glm::vec2(0.0f, 0.0f));
    
    // Set the coordinate of each cell center
    float center_idx = (res_ - 1) / 2.0f;
    for (int i = 0; i < res_; i++)
    {
        for (int j = 0; j < res_; j++)
        {
            float x = (j - center_idx) * cell_size;
            float y = (center_idx - i) * cell_size;
            cell_coord_[i * res_ + j] = glm::vec2(x, y);
        }
    }

    // Set the coordinate of u face center
    float center_idx_x = res_ / 2.0f;
    float center_idx_y = (res_ - 1) / 2.0f;
    for (int i = 0; i < res_; i++)
    {
        for (int j = 0; j < res_ + 1; j++)
        {
            float x = (j - center_idx_x) * cell_size;
            float y = (center_idx_y - i) * cell_size;
            u_coord_[i * (res_ + 1) + j] = glm::vec2(x, y);
        }
    }

    // Set the coordinate of v face center
    center_idx_x = (res_ - 1) / 2.0f;
    center_idx_y = res_ / 2.0f;
    for (int i = 0; i < res_ + 1; i++)
    {
        for (int j = 0; j < res_; j++)
        {
            float x = (j - center_idx_x) * cell_size;
            float y = (center_idx_y - i) * cell_size;
            v_coord_[i * res_ + j] = glm::vec2(x, y);
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
        std::cerr << "[Error] Invalid cell coordinate index : input idx (" << x << ", " << y << ")" << std::endl;
        return glm::vec2(-999.0f, -999.0f);
    }
}

std::vector<glm::vec2> MACGrid2D::getUCoord() const
{
    return u_coord_;
}

glm::vec2 MACGrid2D::getUCoord(int x, int y) const
{
    if (0 <= y && y < res_ && 0 <= x && x < (res_ + 1))
    {
        return u_coord_[y * (res_ + 1) + x];
    }
    else
    {
        std::cerr << "[Error] Invalid u coordinate index : input idx (" << x << ", " << y << ")" << std::endl;
        return glm::vec2(-999.0f, -999.0f);
    }
}

std::vector<glm::vec2> MACGrid2D::getVCoord() const
{
    return v_coord_;
}

glm::vec2 MACGrid2D::getVCoord(int x, int y) const
{
    if (0 <= y && y < (res_ + 1) && 0 <= x && x < res_)
    {
        return v_coord_[y * res_ + x];
    }
    else
    {
        std::cerr << "[Error] Invalid v coordinate index : input idx (" << x << ", " << y << ")" << std::endl;
        return glm::vec2(-999.0f, -999.0f);
    }
}

int MACGrid2D::getCellIndex(const glm::vec2 &cell_coord) const
{
    glm::vec2 temp_coord = cell_coord + glm::vec2(res_ / 2.0f * cell_size_, -(res_ / 2.0f * cell_size_));
    
    int x = static_cast<int>(temp_coord.x / cell_size_);
    int y = static_cast<int>(-temp_coord.y / cell_size_);

    x = glm::clamp(x, 0, res_ - 1);
    y = glm::clamp(y, 0, res_ - 1);

    if (0 <= y && y < res_ && 0 <= x && x < res_)
    {
        return y * res_ + x;
    }
    else
    {
        std::cerr << "[Error] Invalid cell coordinate : input coord (" << cell_coord.x << ", " << cell_coord.y << ")" << std::endl;
        std::cerr << "Calculated cell idx : (" << x << ", " << y << ")" << std::endl;
        return -1; // Invalid index
    }
}

int MACGrid2D::getUIndex(const glm::vec2 &u_coord) const
{
    glm::vec2 temp_coord = u_coord + glm::vec2(res_ / 2.0f * cell_size_, -(res_ / 2.0f * cell_size_));

    int x = static_cast<int>(std::round(temp_coord.x / cell_size_));
    int y = static_cast<int>(-temp_coord.y / cell_size_);

    x = glm::clamp(x, 0, res_);
    y = glm::clamp(y, 0, res_ - 1);

    if (0 <= y && y < res_ && 0 <= x && x < (res_ + 1))
    {
        return y * (res_ + 1) + x;
    }
    else
    {
        std::cerr << "[Error] Invalid U coordinate : input coord (" << u_coord.x << ", " << u_coord.y << ")" << std::endl;
        std::cerr << "Calculated U idx : (" << x << ", " << y << ")" << std::endl;
        return -1; // Invalid index
    }
}

int MACGrid2D::getVIndex(const glm::vec2 &v_coord) const
{
    glm::vec2 temp_coord = v_coord + glm::vec2(res_ / 2.0f * cell_size_, -(res_ / 2.0f * cell_size_));

    int x = static_cast<int>(temp_coord.x / cell_size_);
    int y = static_cast<int>(std::round(-temp_coord.y / cell_size_));

    x = glm::clamp(x, 0, res_ - 1);
    y = glm::clamp(y, 0, res_);

    if (0 <= y && y < (res_ + 1) && 0 <= x && x < res_)
    {
        return y * res_ + x;
    }
    else
    {
        std::cerr << "[Error] Invalid V coordinate : input coord (" << v_coord.x << ", " << v_coord.y << ")" << std::endl;
        std::cerr << "Calculated V idx : (" << x << ", " << y << ")" << std::endl;
        return -1; // Invalid index
    }
}