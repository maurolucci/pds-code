//
// Created by max on 09.08.22.
//

#ifndef PDS_DRAW_GRID_HPP
#define PDS_DRAW_GRID_HPP

#include "pds.hpp"

namespace pds {
namespace style {
enum class Shape {
    Rectangle,
    Circle,
    Triangle,
    Hexagon,
    Octagon
};

struct Color {
    uint8_t r, g, b, a;
    static Color fromRGB(uint8_t r, uint8_t g, uint8_t b, uint8_t a=255);
    static Color fromHSL(float h, float s, float l, float a=1.0);
};

struct NodeStyle {
    Shape shape;
    float size;
    Color fillColor;
    Color drawColor;
    float thickness;
    std::string label;
};

struct EdgeStyle {
    Color color;
    float thickness;
};

struct NodeState {
    long id;
    PmuState pmu;
    bool zeroInjection;
    bool observed;
};

class DrawingOptions {
public:
    virtual NodeStyle nodeStyle(NodeState node) const = 0;
    virtual EdgeStyle edgeStyle(NodeState source, NodeState target) const = 0;
};

class DefaultDrawingOptions : public virtual DrawingOptions {
public:
    Color unobservedColor;
    Color observedColor;
    Color activeColor;
    float nodeSize;
    float nodeThickness;
    float lineThickness;
    Shape nonZeroInjectionShape;
    Shape zeroInjectionShape;

    DefaultDrawingOptions();
    DefaultDrawingOptions(const DefaultDrawingOptions&) = default;
    DefaultDrawingOptions(DefaultDrawingOptions&&) = default;

    DefaultDrawingOptions& operator=(const DefaultDrawingOptions&) = default;
    DefaultDrawingOptions& operator=(DefaultDrawingOptions&&) = default;

    virtual NodeStyle nodeStyle(pds::style::NodeState node) const override;
    virtual EdgeStyle edgeStyle(pds::style::NodeState source, pds::style::NodeState target) const override;
};
}

struct Coordinate {
    double x, y;
};
using Layout = pds::map<long, pds::Coordinate>;

Layout layoutGraph(const PowerGrid& graph);

void drawGrid(const PdsState &state,
              const std::string &filename,
              const Layout& layout,
              const style::DrawingOptions &style);

inline void drawGrid(const PdsState &state,
                     const std::string &filename,
                     const Layout& layout
) {
    return drawGrid(state, filename, layout, style::DefaultDrawingOptions{});
}

} // namespace pds

#endif //PDS_DRAW_GRID_HPP
