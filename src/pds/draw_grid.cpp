//
// Created by max on 09.08.22.
//
#include "draw_grid.hpp"

#include <ogdf/basic/Graph.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/fmmm/FMMMOptions.h>
#include <ogdf/fileformats/GraphIO.h>

namespace pds {
namespace style {
Color Color::fromRGB(uint8_t r, uint8_t g, uint8_t b, uint8_t a) { return {r, g, b, a}; }
Color Color::fromHSL(float h, float s, float l, float a) {
    float r, g, b;

    if(s == 0){
        r = g = b = l * 255; // achromatic
    }else{
        auto hue2rgb = [] (float p, float q, float t) -> float {
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        };

        float q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        float p = 2 * l - q;
        r = hue2rgb(p, q, h + 1/3);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1/3);
    }

    return {static_cast<uint8_t>(r * 255), static_cast<uint8_t>(g * 255), static_cast<uint8_t>(b * 255), static_cast<uint8_t>(a * 255)};
}

DefaultDrawingOptions::DefaultDrawingOptions() :
        unobservedColor(Color::fromRGB(0,0,0)),
        observedColor(Color::fromRGB(53, 221, 95)),
        activeColor(Color::fromRGB(221, 53, 95)),
        nodeSize(8),
        nodeThickness(2),
        lineThickness(2),
        nonZeroInjectionShape(Shape::Triangle),
        zeroInjectionShape(Shape::Circle)
{ }

NodeStyle DefaultDrawingOptions::nodeStyle(pds::style::NodeState node) const {
    Color nodeColor = node.observed ? observedColor : unobservedColor;
    Color fillColor = nodeColor;
    switch (node.pmu) {
        case PmuState::Inactive:
            fillColor = Color::fromRGB(255, 255, 255, 0);
            break;
        case PmuState::Active:
            nodeColor = activeColor;
            break;
        case PmuState::Blank:
            break;
    }
    Shape shape = node.zeroInjection ? zeroInjectionShape : nonZeroInjectionShape;
    return NodeStyle {
            .shape=shape,
            .size=nodeSize,
            .fillColor=fillColor,
            .drawColor=nodeColor,
            .thickness=nodeThickness,
            .label={}
    };
}

EdgeStyle DefaultDrawingOptions::edgeStyle(pds::style::NodeState, pds::style::NodeState) const {
    return {.color=Color::fromRGB(0, 0, 0, 255), .thickness=lineThickness};
}

ogdf::Shape toOgdf(Shape shape) {
    switch (shape) {
        case Shape::Rectangle:
            return ogdf::Shape::Rect;
        case Shape::Circle:
            return ogdf::Shape::Ellipse;
        case Shape::Hexagon:
            return ogdf::Shape::Hexagon;
        case Shape::Octagon:
            return ogdf::Shape::Octagon;
        case Shape::Triangle:
            return ogdf::Shape::Triangle;
    }
}

ogdf::Color toOgdf(Color color) {
    return ogdf::Color(color.r, color.g, color.b, color.a);
}

} // namespace style

struct FixedBBGraphAttributes : public virtual ogdf::GraphAttributes {
    ogdf::DRect bounds;
    FixedBBGraphAttributes(const ogdf::GraphAttributes& attr, ogdf::DRect bounds) : ogdf::GraphAttributes(attr), bounds(bounds) { }

    inline virtual ogdf::DRect boundingBox() const override {
        return bounds;
    }
};

Layout layoutGraph(const PowerGrid &graph) {
    ogdf::Graph G;
    map<long, ogdf::node> id_to_node;
    for (auto v: graph.vertices()) {
        id_to_node.emplace(graph[v].id, G.newNode());
    }
    for (auto e: graph.edges()) {
        auto [s, t] = graph.endpoints(e);
        G.newEdge(id_to_node.at(graph[s].id), id_to_node.at(graph[t].id));
    }
    ogdf::GraphAttributes GA(G);
    ogdf::FMMMLayout layout;
    layout.useHighLevelOptions(true);
    layout.unitEdgeLength(30.0);
    layout.newInitialPlacement(true);
    layout.qualityVersusSpeed(ogdf::FMMMOptions::QualityVsSpeed::GorgeousAndEfficient);
    layout.randSeed(1234);
    layout.call(GA);

    Layout coordinates;
    for (auto v: graph.vertices()) {
        auto node = id_to_node.at(graph[v].id);
        coordinates.emplace(v, Coordinate{GA.x(node), GA.y(node)});
    }
    return coordinates;
}

void drawGrid(const PdsState& state,
              const std::string &filename,
              const Layout& layout,
              const style::DrawingOptions &style) {
    ogdf::Graph G;
    ogdf::GraphAttributes GA(
            G,
            ogdf::GraphAttributes::nodeGraphics
            | ogdf::GraphAttributes::edgeGraphics
            | ogdf::GraphAttributes::nodeStyle
            | ogdf::GraphAttributes::edgeStyle
            | ogdf::GraphAttributes::edgeArrow
            | ogdf::GraphAttributes::nodeLabel
    );
    map<PowerGrid::VertexDescriptor, ogdf::node> id_to_node;
    const auto& graph = state.graph();
    auto nodeState = [&state, &graph](const PowerGrid::VertexDescriptor& v) -> style::NodeState {
        return {
                .id=graph[v].id,
                .pmu=state.activeState(v),
                .zeroInjection=graph[v].zero_injection,
                .observed=bool(state.isObserved(v))
        };
    };
    double minx = std::numeric_limits<double>::max(), maxx = 0, miny = std::numeric_limits<double>::max(), maxy = 0;
    double centerx = 0.0, centery = 0.0;
    for (auto [v, coord]: layout) {
        centerx += coord.x;
        centery += coord.y;
    }
    centerx /= layout.size();
    centery /= layout.size();
    double maxSize = 0;
    for (auto v: graph.vertices()) {
        auto node = G.newNode();
        id_to_node.emplace(v, node);
        if (layout.contains(graph[v].id)) {
            GA.x(node) = layout.at(graph[v].id).x;
            GA.y(node) = layout.at(graph[v].id).y;
        } else {
            GA.x(node) = centerx;
            GA.y(node) = centery;
        }
        auto nodeStyle = style.nodeStyle(nodeState(v));
        GA.fillColor(node) = style::toOgdf(nodeStyle.fillColor);
        GA.strokeColor(node) = style::toOgdf(nodeStyle.drawColor);
        GA.strokeWidth(node) = nodeStyle.thickness;
        GA.shape(node) = style::toOgdf(nodeStyle.shape);
        GA.label(node) = nodeStyle.label;
        GA.width(node) = nodeStyle.size;
        GA.height(node) = nodeStyle.size;
        maxSize = std::max(maxSize, double(nodeStyle.size + nodeStyle.thickness));
    }
    for (auto [v, coord]: layout) {
        minx = std::min(minx, coord.x - maxSize);
        maxx = std::max(maxx, coord.x + maxSize);
        miny = std::min(miny, coord.y - maxSize);
        maxy = std::max(maxy, coord.y + maxSize);
    }

    for (auto e: graph.edges()) {
        auto [s, t] = graph.endpoints(e);
        auto edge = G.newEdge(id_to_node.at(s), id_to_node.at(t));
        auto edgeStyle = style.edgeStyle(nodeState(s), nodeState(t));
        GA.strokeColor(edge) = style::toOgdf(edgeStyle.color);
        GA.strokeWidth(edge) = edgeStyle.thickness;
        if (state.isObservedEdge(s, t)) {
            if (state.isObservedEdge(t, s)) {
                GA.arrowType(edge) = ogdf::EdgeArrow::Both;
            } else {
                GA.arrowType(edge) = ogdf::EdgeArrow::Last;
            }
        } else if (state.isObservedEdge(t, s)) {
            GA.arrowType(edge) = ogdf::EdgeArrow::First;
        } else {
            GA.arrowType(edge) = ogdf::EdgeArrow::None;
        }
    }
    minx = std::min(minx, maxx);
    miny = std::min(miny, maxy);
    ogdf::DRect bounds{minx, miny, maxx, maxy};
    ogdf::GraphIO::drawSVG(FixedBBGraphAttributes(GA, bounds), filename);
}
} // namespace pds
