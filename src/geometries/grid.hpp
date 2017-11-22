#pragma once

#include "geometry.hpp"



struct GridGeometry
    :public Geometry
{    
    GridGeometry();
    GridGeometry(std::size_t, std::size_t, std::size_t);

    virtual void generate() override;
    virtual void scale(const cartesian&) override;
    virtual void shift(const cartesian&) override;
    cartesian& matrixView(const std::size_t&, const std::size_t&, const std::size_t&);
    const cartesian& matrixView(const std::size_t&, const std::size_t&, const std::size_t&) const;

    std::size_t x {0};
    std::size_t y {0};
    std::size_t z {0};
};