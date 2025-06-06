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
    double angle_deg = 45.0; // angle for both rotations
    double angle_rad = angle_deg * M_PI / 180.0;
    double cx = cos(angle_rad), sx = sin(angle_rad); // for X rotation
    double cy = cos(angle_rad), sy = sin(angle_rad); // for Y rotation
    double tx = 0, ty = 0, tz = 0; // translation for Frame B

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
    xA.reserve(N); yA.reserve(N); zA.reserve(N);
    xB.reserve(N); yB.reserve(N); zB.reserve(N);

    // --- Create the combined rotation matrix ---
    // Rotation about X axis:
    // R_x = [ 1,   0,    0 ]
    //       [ 0,  cx, -sx ]
    //       [ 0,  sx,  cx ]
    //
    // Rotation about Y axis:
    // R_y = [ cy, 0, sy ]
    //       [  0, 1,  0 ]
    //       [ -sy,0, cy ]
    //
    // Combined (R = R_y * R_x):
    // R[0][0] = cy;       R[0][1] = sy * sx;       R[0][2] = sy * cx;
    // R[1][0] = 0;        R[1][1] = cx;            R[1][2] = -sx;
    // R[2][0] = -sy;      R[2][1] = cy * sx;       R[2][2] = cy * cx;
    array<array<double, 3>, 3> R = {{
        { cy,       sy * sx,   sy * cx },
        { 0,        cx,        -sx },
        { -sy,      cy * sx,   cy * cx }
    }};

    cout << "Combined rotation matrix R:\n";
    for (const auto& row : R)
    {
        for (double val : row)
            cout << setw(8) << val << " ";
        cout << "\n";
    }

    // Generate points and transform for Frame B
    for (int i = 0; i < N; ++i)
    {
        double t = t_dist(gen);
        // ideal point on the line
        double x_line = dx * t;
        double y_line = dy * t;
        double z_line = dz * t;
        // add noise for Frame A
        double x_noisy = x_line + noise(gen);
        double y_noisy = y_line + noise(gen);
        double z_noisy = z_line + noise(gen);
        xA.push_back(x_noisy);
        yA.push_back(y_noisy);
        zA.push_back(z_noisy);

        // Scale noisy point for Frame B
        double xs = x_noisy * scale;
        double ys = y_noisy * scale;
        double zs = z_noisy * scale;
        // Apply combined rotation matrix (R = R_y * R_x)
        double xr = R[0][0] * xs + R[0][1] * ys + R[0][2] * zs;
        double yr = R[1][0] * xs + R[1][1] * ys + R[1][2] * zs;
        double zr = R[2][0] * xs + R[2][1] * ys + R[2][2] * zs;
        // Apply translation for Frame B
        double xb = xr + tx;
        double yb = yr + ty;
        double zb = zr + tz;
        xB.push_back(xb);
        yB.push_back(yb);
        zB.push_back(zb);

        fileA << x_noisy << "," << y_noisy << "," << z_noisy << "\n";
        fileB << xb << "," << yb << "," << zb << "\n";
    }

    fileA.close();
    fileB.close();
    cout << "Wrote " << N << " points to points_frame_a.csv and points_frame_b.csv\n";

    // --- 2) Visualization ---
    double L_A = 5.0;              // axis length for Frame A
    double L_B = L_A * scale;      // axis length for Frame B

    // Compute Frame B axes endpoints by rotating local axis vectors and adding translation.
    // Local x-axis: [L_B, 0, 0]
    // Local y-axis: [0, L_B, 0]
    // Local z-axis: [0, 0, L_B]
    array<double, 3> Bx_end = { tx + R[0][0] * L_B,
                                ty + R[1][0] * L_B,
                                tz + R[2][0] * L_B };
    array<double, 3> By_end = { tx + R[0][1] * L_B,
                                ty + R[1][1] * L_B,
                                tz + R[2][1] * L_B };
    array<double, 3> Bz_end = { tx + R[0][2] * L_B,
                                ty + R[1][2] * L_B,
                                tz + R[2][2] * L_B };

    auto fig = figure(true);
    hold(on);

    // scatter points from both frames
    plot3(xA, yA, zA, ".b")->display_name("Frame ECEF positions");
    plot3(xB, yB, zB, ".r")->display_name("Frame SLAM positions");

    // Frame A axes (length L_A)
    vector<double> Ax{0.0, L_A}, Ay{0.0, 0.0}, Az{0.0, 0.0};
    vector<double> Bx{0.0, 0.0}, By{0.0, L_A}, Bz{0.0, 0.0};
    vector<double> Cx{0.0, 0.0}, Cy{0.0, 0.0}, Cz{0.0, L_A};

    plot3(Ax, Ay, Az, "r-")->display_name("ECEF X-axis");
    plot3(Bx, By, Bz, "g-")->display_name("ECEF Y-axis");
    plot3(Cx, Cy, Cz, "b-")->display_name("ECEF Z-axis");

    // Frame B axes (rotated+translated)
    vector<double> Ax_b{tx, Bx_end[0]},
                   Ay_b{ty, Bx_end[1]},
                   Az_b{tz, Bx_end[2]};
    vector<double> Bx_b{tx, By_end[0]},
                   By_b{ty, By_end[1]},
                   Bz_b{tz, By_end[2]};
    vector<double> Cx_b{tx, Bz_end[0]},
                   Cy_b{ty, Bz_end[1]},
                   Cz_b{tz, Bz_end[2]};

    plot3(Ax_b, Ay_b, Az_b, "-")->color({0.5, 0.0, 0.5}).display_name("SLAM X-axis");      // purple
    plot3(Bx_b, By_b, Bz_b, "-")->color({0.0, 0.39, 0.0}).display_name("SLAM Y-axis");     // dark green
    plot3(Cx_b, Cy_b, Cz_b, "-")->color({0.0, 0.0, 0.55}).display_name("SLAM Z-axis");    // dark blue

    xlabel("X");
    ylabel("Y");
    zlabel("Z");
    title("Frame ECEF vs Frame SLAM");
    legend();

    // Axes limits
    xlim({-10, 10});
    ylim({-10, 10});
    zlim({-10, 10});

    show();
    return 0;
}