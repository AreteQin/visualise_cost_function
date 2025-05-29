// combined.cpp
// – Generates N noisy points along a line, writes points_frame_a.csv and points_frame_b.csv,
//   then visualizes both clouds and their coordinate axes in 3D.
// Requires matplot-cpp: https://github.com/alandefreitas/matplotplusplus

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>
#include <vector>
#include <array>
#include <matplot/matplot.h>

using namespace std;
using namespace matplot;

int main()
{
    // --- 1) Generation parameters ---
    const int N = 100;
    double scale = 0.5;
    double angle_deg = 45.0; // MUST match generation & visualization
    double angle_rad = angle_deg * M_PI / 180.0;
    double cosA = cos(angle_rad);
    double sinA = sin(angle_rad);
    double tx = 1.0, ty = 2.0, tz = 1.0; // translation for Frame B

    // Line direction (unit)
    double dx = 1.0 / sqrt(3.0),
           dy = 1.0 / sqrt(3.0),
           dz = 1.0 / sqrt(3.0);

    // RNGs
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> t_dist(-10.0, 10.0);
    normal_distribution<double> noise(0.0, 0.2);

    // Prepare output files
    ofstream fileA("points_frame_a.csv"), fileB("points_frame_b.csv");
    fileA << "x,y,z\n";
    fileB << "x,y,z\n";
    fileA << fixed << setprecision(4);
    fileB << fixed << setprecision(4);

    // Containers for plotting
    vector<double> xA, yA, zA, xB, yB, zB;
    xA.reserve(N);
    yA.reserve(N);
    zA.reserve(N);
    xB.reserve(N);
    yB.reserve(N);
    zB.reserve(N);

    // Rotation matrix about X axis by angle_deg
    array<array<double, 3>, 3> R = {
        {
            {1, 0, 0},
            {0, cosA, -sinA},
            {0, sinA, cosA}
        }
    };
    // print rotation matrix
    cout << "Rotation matrix R:\n";
    for (const auto& row : R)
    {
        for (double val : row)
            cout << setw(8) << val << " ";
        cout << "\n";
    }

    // Generate
    for (int i = 0; i < N; ++i)
    {
        double t = t_dist(gen);
        // point on the ideal line
        double x_line = dx * t;
        double y_line = dy * t;
        double z_line = dz * t;
        // add noise
        double x_noisy = x_line + noise(gen);
        double y_noisy = y_line + noise(gen);
        double z_noisy = z_line + noise(gen);
        // Frame A coords
        xA.push_back(x_noisy);
        yA.push_back(y_noisy);
        zA.push_back(z_noisy);

        // scale for Frame B
        double xs = x_noisy * scale;
        double ys = y_noisy * scale;
        double zs = z_noisy * scale;
        // rotate by R
        double xr = R[0][0] * xs + R[0][1] * ys + R[0][2] * zs;
        double yr = R[1][0] * xs + R[1][1] * ys + R[1][2] * zs;
        double zr = R[2][0] * xs + R[2][1] * ys + R[2][2] * zs;
        // translate
        double xb = xr + tx;
        double yb = yr + ty;
        double zb = zr + tz;
        xB.push_back(xb);
        yB.push_back(yb);
        zB.push_back(zb);

        // write CSV
        fileA << x_noisy << "," << y_noisy << "," << z_noisy << "\n";
        fileB << xb << "," << yb << "," << zb << "\n";
    }

    fileA.close();
    fileB.close();
    cout << "Wrote " << N << " points to points_frame_a.csv and points_frame_b.csv\n";

    // --- 2) Visualization ---
    double L_A = 5.0; // axis length for Frame A
    double L_B = L_A * scale; // axis length for Frame B

    // Compute Frame B axes endpoints in Frame A coords
    array<double, 3> Bx_end = {tx + L_B, ty, tz};
    array<double, 3> By_end = {tx, ty + cosA * L_B, tz + sinA * L_B};
    array<double, 3> Bz_end = {tx, ty - sinA * L_B, tz + cosA * L_B};

    auto fig = figure(true);
    hold(on);

    // scatter points
    plot3(xA, yA, zA, ".b")->display_name("Frame A pts");
    plot3(xB, yB, zB, ".r")->display_name("Frame B pts");

    // --- before plotting axes, build vectors explicitly ---
    std::vector<double> Ax{0.0, L_A}, Ay_x{0.0, 0.0}, Az_x{0.0, 0.0};
    std::vector<double> Ax_b{tx, Bx_end[0]}, Ay_b{ty, Bx_end[1]}, Az_b{tz, Bx_end[2]};

    std::vector<double> Ay_a{0.0, 0.0}, By_a{0.0, L_A}, Zy_a{0.0, 0.0};
    std::vector<double> Ay_b2{tx, By_end[0]}, By_b{ty, By_end[1]}, Zy_b{tz, By_end[2]};

    std::vector<double> Az_a{0.0, 0.0}, Zy2_a{0.0, 0.0}, Za_a{0.0, L_A};
    std::vector<double> Az_b2{tx, Bz_end[0]}, Zy_b2{ty, Bz_end[1]}, Za_b{tz, Bz_end[2]};

    // Frame A axes (length L_A)
    plot3(Ax, Ay_x, Az_x, "r-")->display_name("A X-axis");
    plot3(Ay_a, By_a, Zy_a, "g-")->display_name("A Y-axis");
    plot3(Az_a, Zy2_a, Za_a, "b-")->display_name("A Z-axis");

    // Frame B axes (length L_B), rotated+translated
    plot3(Ax_b, Ay_b, Az_b, "r-")->display_name("B X-axis");
    plot3(Ay_b2, By_b, Zy_b, "g-")->display_name("B Y-axis");
    plot3(Az_b2, Zy_b2, Za_b, "b-")->display_name("B Z-axis");

    // Labels and legend
    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    title("Frame A vs Frame B (points + axes)");
    legend();

    // Axes limits
    xlim({-10, 10});
    ylim({-10, 10});
    zlim({-10, 10});

    show();
    return 0;
}
