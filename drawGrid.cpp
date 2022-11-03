#include "pds.hpp"
#include "graphio.hpp"
#include "draw_grid.hpp"

#include <boost/program_options.hpp>
#include <fmt/format.h>

#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/fmmm/FMMMOptions.h>

using Layout = pds::map<pds::PowerGrid::VertexDescriptor, pds::Coordinate>;

void writeLayout(const Layout& layout, const std::string& filename) {
    FILE* outfile = fopen(filename.c_str(), "wb");
    if (outfile == nullptr) throw std::system_error(std::error_code(errno, std::system_category()));
    fmt::print(outfile, "Layout V1.0\n");
    for (auto [k, v]: layout) {
        fmt::print(outfile, "{} {} {}\n", k, v.x, v.y);
    }
    fclose(outfile);
}

struct ParseError : std::exception {
    std::string reason;
    ParseError(const std::string& reason, size_t line) : reason(fmt::format("{}: {}", line, reason)) { }
    const char * what() const noexcept override { return reason.c_str(); }
};

std::vector<std::string_view> split(const std::string_view& in, char delim) {
    std::vector<std::string_view> pieces;
    size_t start = 0;
    size_t end;
    while (start != std::string::npos) {
        end = in.find(delim, start);
        pieces.emplace_back(in.substr(start, end - start));
        if (end == std::string::npos) {
            start = std::string::npos;
        } else {
            start = end + 1;
        }

    }
    return pieces;
}

Layout readLayout(const std::string& filename) {
    std::ifstream infile(filename);
    std::string line;
    size_t lineno = 0;
    Layout layout;
    std::getline(infile, line);
    ++lineno;
    if (line != "Layout V1.0") throw ParseError { fmt::format("unknown format: {}", line), lineno };
    while (infile) {
        std::getline(infile, line);
        ++lineno;
        auto pieces = split(line, ' ');
        if (pieces.size() == 1 && pieces[0] == "") continue; // empty line
        if (pieces.size() != 3) throw ParseError { "not a triple", lineno };
        try {
            unsigned long v = std::stoi(std::string{pieces[0]});
            double x = std::stoi(std::string{pieces[1]});
            double y = std::stoi(std::string{pieces[2]});
            layout[v] = {x, y};
        } catch(const std::exception& ex) {
            throw ParseError { fmt::format("could not read number: {}", ex.what()), lineno };
        }
    }
    return layout;
}

void restrictLayout(Layout& layout, const pds::PowerGrid& graph) {
    std::vector<unsigned long> keys = layout | ranges::views::transform([](const auto& kv) { return kv.first;}) | ranges::to<std::vector>;
    for (auto k: keys) {
        erase_if(layout,[&graph](const auto& kv) { return !graph.hasVertex(kv.first);});
    }
}

int main(int argc, const char** argv) {
    namespace po = boost::program_options;
    using std::string;
    using namespace pds;

    po::options_description desc;
    desc.add_options()
            ("layout,l", po::value<string>(), "load layout from file")
            ("save-layout,s", po::value<string>(), "save layout to file")
            ("graph,g", po::value<string>()->required(), "graph file")
            ("draw,d", po::value<string>(), "create graph svg file")
            ("layout-bounds,b", "use bounds from layout for drawing")
            ("help,h", "show this help")
    ;
    po::positional_options_description pos;
    pos.add("graph", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);
    if (vm.count("help")) {
        desc.print(std::cout);
        return EXIT_FAILURE;
    }

    auto graph = readAutoGraph(vm["graph"].as<string>());

    Layout layout;
    if (vm.count("layout")) {
        layout = readLayout(vm["layout"].as<string>());
        if (!vm.count("layout-bounds")) {
            restrictLayout(layout, graph);
        }
    } else {
        layout = layoutGraph(graph);
    }
    if (vm.count("draw")) {
        PdsState state(graph);
        drawGrid(state, vm["draw"].as<string>(), layout);
    }
    if (vm.count("save-layout")) {
        writeLayout(layout, vm["save-layout"].as<string>());
    }
}