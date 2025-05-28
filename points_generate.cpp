#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>

int main() {
    // Transformation params (unchanged)
    double scale = 0.5;
    double angle_rad = 45.0 * M_PI / 180.0;
    double cosA = std::cos(angle_rad), sinA = std::sin(angle_rad);
    double tx=1.0, ty=2.0, tz=1.0;

    const int N = 100;

    // --- NEW: set up line + noise distributions ---
    // Direction vector (unit)
    double dx = 1.0/std::sqrt(3.0),
           dy = 1.0/std::sqrt(3.0),
           dz = 1.0/std::sqrt(3.0);
    // Parameter t in [-10,10]
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> t_dist(-10.0, 10.0);
    // Gaussian noise with sigma = 0.2
    std::normal_distribution<double> noise(0.0, 0.2);

    std::ofstream fileA("points_frame_a.csv"),
                  fileB("points_frame_b.csv");
    fileA << "x,y,z\n"; fileB << "x,y,z\n";
    fileA << std::fixed << std::setprecision(4);
    fileB << std::fixed << std::setprecision(4);

    for (int i = 0; i < N; ++i) {
        // 1) pick a random t along the line
        double t = t_dist(gen);

        // 2) base point on the line
        double x_line = dx * t;
        double y_line = dy * t;
        double z_line = dz * t;

        // 3) add small noise
        double xA = x_line + noise(gen);
        double yA = y_line + noise(gen);
        double zA = z_line + noise(gen);

        // 4) transform to Frame B as before
        //    scale
        double xs = xA * scale, ys = yA * scale, zs = zA * scale;
        //    rotate 45° about X
        double xr = xs;
        double yr = cosA*ys - sinA*zs;
        double zr = sinA*ys + cosA*zs;
        //    translate
        double xB = xr + tx, yB = yr + ty, zB = zr + tz;

        // 5) write to CSV
        fileA << xA << "," << yA << "," << zA << "\n";
        fileB << xB << "," << yB << "," << zB << "\n";
    }

    fileA.close();  fileB.close();
    std::cout << "Wrote " << N << " nearly-collinear points." << std::endl;
    return 0;
}
