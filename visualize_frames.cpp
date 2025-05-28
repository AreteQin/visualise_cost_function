#include <matplot/matplot.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <array>

using namespace matplot;

// Helper: read "x,y,z" CSV into three vectors
static std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
read_csv(const std::string &filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open " + filename);
    std::string line;
    std::getline(in, line);  // skip header

    std::vector<double> xs, ys, zs;
    while (std::getline(in, line)) {
        std::stringstream ss(line);
        double x, y, z;
        char comma;
        ss >> x >> comma >> y >> comma >> z;
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
    }
    return {xs, ys, zs};
}

int main() {
    // 1) Load Frame A & B points (as generated previously)
    auto [xA, yA, zA] = read_csv("points_frame_a.csv");
    auto [xB, yB, zB] = read_csv("points_frame_b.csv");

    // 2) Transformation parameters (must match how points were generated)
    double scale = 0.5;
    double angle_rad = 45.0 * M_PI / 180.0;
    double cosA = std::cos(angle_rad);
    double sinA = std::sin(angle_rad);
    double tx = 1.0, ty = 2.0, tz = 1.0;        // Frame B origin in Frame A coords

    // 3) Axis lengths for display
    double L_A = 5.0;                          // Frame A axes length
    double L_B = L_A * scale;                  // Frame B axes length

    // 4) Compute Frame B axes endpoints (in Frame A / plot coords)
    //    Rotate each B‐axis unit‐vector by R_x(45°), then translate by (tx,ty,tz)
    std::array<double,3> Bx_end = {
        tx +  L_B,
        ty +   0.0,
        tz +   0.0
    };
    std::array<double,3> By_end = {
        tx +   0.0,
        ty + cosA * L_B,
        tz + sinA * L_B
    };
    std::array<double,3> Bz_end = {
        tx +   0.0,
        ty - sinA * L_B,
        tz + cosA * L_B
    };

    // 5) Plot everything
    auto fig = figure(true);
    hold(on);

    // Point clouds
    plot3(xA, yA, zA, ".b")->display_name("Frame A pts");
    plot3(xB, yB, zB, ".r")->display_name("Frame B pts");

    // --- before plotting axes, build vectors explicitly ---
    std::vector<double> Ax{0.0, L_A},  Ay_x{0.0, 0.0},  Az_x{0.0, 0.0};
    std::vector<double> Ax_b{tx, Bx_end[0]}, Ay_b{ty, Bx_end[1]}, Az_b{tz, Bx_end[2]};

    std::vector<double> Ay_a{0.0, 0.0},     By_a{0.0, L_A},  Zy_a{0.0, 0.0};
    std::vector<double> Ay_b2{tx, By_end[0]}, By_b{ty, By_end[1]}, Zy_b{tz, By_end[2]};

    std::vector<double> Az_a{0.0, 0.0},     Zy2_a{0.0, 0.0}, Za_a{0.0, L_A};
    std::vector<double> Az_b2{tx, Bz_end[0]}, Zy_b2{ty, Bz_end[1]}, Za_b{tz, Bz_end[2]};

    // Frame A axes (length L_A)
    plot3(Ax,   Ay_x,   Az_x,   "r-")->display_name("A X-axis");
    plot3(Ay_a, By_a,   Zy_a,   "g-")->display_name("A Y-axis");
    plot3(Az_a, Zy2_a,  Za_a,   "b-")->display_name("A Z-axis");

    // Frame B axes (length L_B), rotated+translated
    plot3(Ax_b, Ay_b,   Az_b,   "r-")->display_name("B X-axis");
    plot3(Ay_b2,By_b,   Zy_b,   "g-")->display_name("B Y-axis");
    plot3(Az_b2,Zy_b2,  Za_b,   "b-")->display_name("B Z-axis");

    // Labels & legend
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    title("Frame A vs Frame B (points + axes)");
    legend();

    // Set axis limits: adjust the min and max values as needed.
    xlim({-10, 10});
    ylim({-10, 10});
    zlim({-10, 10});
    show();
    return 0;
}
