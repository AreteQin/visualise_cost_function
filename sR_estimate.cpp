//
// Created by qin on 27/05/25.
//
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <iomanip>   // for std::setprecision

int main() {
    // Define transformation parameters
    double scale = 0.5;
    double angle_degrees = 45.0;
    double angle_rad = angle_degrees * M_PI / 180.0;  // Convert degrees to radians
    double cosA = std::cos(angle_rad);
    double sinA = std::sin(angle_rad);
    // Translation vector components
    double tx = 1.0, ty = 2.0, tz = 1.0;

    // Number of points
    const int N = 100;

    // Set up random number generation for Frame A points
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    // Open output CSV files
    std::ofstream fileA("points_frame_a.csv");
    std::ofstream fileB("points_frame_b.csv");
    if (!fileA.is_open() || !fileB.is_open()) {
        std::cerr << "Error: Unable to open output file(s)." << std::endl;
        return 1;
    }

    // Write CSV headers
    fileA << "x,y,z\n";
    fileB << "x,y,z\n";

    // Set numeric format for output (fixed with 4 decimal places)
    fileA << std::fixed << std::setprecision(4);
    fileB << std::fixed << std::setprecision(4);

    // Generate points in Frame A, apply transform to get Frame B points
    for (int i = 0; i < N; ++i) {
        // Random coordinates in Frame A
        double xA = dist(gen);
        double yA = dist(gen);
        double zA = dist(gen);

        // Apply transformation from Frame A to Frame B:

        // 1. Scale (Frame B has scale factor 0.5 relative to Frame A)
        double x_scaled = xA * scale;
        double y_scaled = yA * scale;
        double z_scaled = zA * scale;

        // 2. Rotate 45 degrees around X-axis (using rotation matrix for X-axis)
        // X coordinate remains the same under rotation about X-axis
        double x_rot = x_scaled;
        // Y and Z coordinates are rotated in the Y-Z plane
        double y_rot = cosA * y_scaled - sinA * z_scaled;
        double z_rot = sinA * y_scaled + cosA * z_scaled;

        // 3. Translate by (1, 2, 1)
        double xB = x_rot + tx;
        double yB = y_rot + ty;
        double zB = z_rot + tz;

        // Write the original and transformed point to the CSV files
        fileA << xA << "," << yA << "," << zA << "\n";
        fileB << xB << "," << yB << "," << zB << "\n";
    }

    // Close files
    fileA.close();
    fileB.close();

    std::cout << "Generated " << N << " points and saved to CSV files." << std::endl;
    return 0;
}
