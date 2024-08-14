#include <iostream>
#include <memory>
//#include <g2o/core/g2o_core_api.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
//#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/dense/linear_solver_dense.h>
#include <Eigen/Core>
#include <opencv2/core/core.hpp>
#include <cmath>
#include <chrono>
#include <matplot/matplot.h>
#include <glog/logging.h>
using namespace std;

// 曲线模型的顶点，模板参数：优化变量维度和数据类型
class CurveFittingVertex : public g2o::BaseVertex<2, Eigen::Vector2d>
{
public:
    //    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // 重置
    virtual void setToOriginImpl() override
    {
        _estimate << 0, 0;
    }

    //  update estimation
    virtual void oplusImpl(const double* update) override
    {
        _estimate += Eigen::Vector2d(update);
    }

    // 存盘和读盘：留空
    virtual bool read(istream& in)
    {
    }

    virtual bool write(ostream& out) const
    {
    }
};

// 误差模型 模板参数：观测值维度，类型，连接顶点类型
class CurveFittingEdge : public g2o::BaseUnaryEdge<1, double, CurveFittingVertex>
{
public:
    // constructor
    CurveFittingEdge(double x) : BaseUnaryEdge(), _x(x)
    {
    }

    // 计算曲线模型误差
    virtual void computeError() override
    {
        const CurveFittingVertex* v = static_cast<const CurveFittingVertex*>(_vertices[0]);
        const Eigen::Vector2d ab = v->estimate(); // get current estimated curve parameters.
        _error(0, 0) = _measurement - std::exp(ab(0) * _x * _x + ab(1) * _x);
    }

    // 计算雅可比矩阵
    virtual void linearizeOplus() override
    {
        const CurveFittingVertex* v = static_cast<const CurveFittingVertex*>(_vertices[0]);
        const Eigen::Vector2d ab = v->estimate();
        double y = exp(ab[0] * _x * _x + ab[1] * _x);
        _jacobianOplusXi[0] = -_x * _x * y; // d(error)/d(a)
        _jacobianOplusXi[1] = -_x * y; // d(error)/d(b)
    }

    virtual bool read(istream& in)
    {
    }

    virtual bool write(ostream& out) const
    {
    }

public:
    double _x; // x 值， y 值为 _measurement
};

int main(int argc, char** argv)
{
    double ar = 1.0, br = 2.0; // 真实参数值
    double ae = 3.0, be = 3.0; // 估计参数值
    int N = 100; // 数据点
    double w_sigma = 1.0; // 噪声Sigma值
    double inv_sigma = 1.0 / w_sigma;
    cv::RNG rng; // OpenCV随机数产生器

    std::vector<double> x_real, y_real, x_data, y_data;
    for (int i = 0; i < N; i++)
    {
        double x = i / 100.0;
        x_real.push_back(x);
        y_real.push_back(exp(ar * x * x + br * x));
        x_data.push_back(x);
        y_data.push_back(exp(ar * x * x + br * x) + rng.gaussian(w_sigma * w_sigma));
    }

    // Plot given data and real plot
    matplot::plot(x_data, y_data, "o");
    matplot::hold(matplot::on);
    matplot::plot(x_real, y_real, "-");
    matplot::show();
    matplot::hold(matplot::off);

    // 构建图优化，先设定g2o
    typedef g2o::BlockSolver<g2o::BlockSolverTraits<2, 1>> BlockSolverType; // 每个误差项优化变量维度为3，误差值维度为1
    typedef g2o::LinearSolverDense<BlockSolverType::PoseMatrixType> LinearSolverType; // 线性求解器类型

    // 梯度下降方法，可以从GN, LM, DogLeg 中选
    auto solver = new g2o::OptimizationAlgorithmLevenberg(
        std::make_unique<BlockSolverType>(std::make_unique<LinearSolverType>()));
    g2o::SparseOptimizer optimizer; // 图模型
    optimizer.setAlgorithm(solver); // 设置求解器
    optimizer.setVerbose(false); // 打开调试输出

    // 往图中增加顶点
    CurveFittingVertex* v = new CurveFittingVertex();
    v->setEstimate(Eigen::Vector2d(ae, be));
    v->setId(0);
    optimizer.addVertex(v);

    // 往图中增加边
    for (int i = 0; i < N; i++)
    {
        CurveFittingEdge* edge = new CurveFittingEdge(x_data[i]);
        edge->setId(i);
        edge->setVertex(0, v); // 设置连接的顶点
        edge->setMeasurement(y_data[i]); // 观测数值
        edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity() * 1 / (w_sigma * w_sigma)); // 信息矩阵：协方差矩阵之逆
        optimizer.addEdge(edge);
    }

    // 执行优化
    vector<double> a_estimated, b_estimated, error_estimated;
    a_estimated.push_back(ae);
    b_estimated.push_back(be);
    error_estimated.push_back(0);
    for (int i = 0; i < N; i++)
    {
        error_estimated[0] += 0.5*pow((y_data[i]-exp(ae * x_data[i] * x_data[i] + be * x_data[i])),2);
    }
    cout << "start optimization" << endl;
    chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
    optimizer.initializeOptimization();
    for (int i = 0; i < 10; i++)
    {
        optimizer.optimize(1);
        // 输出优化值
        Eigen::Vector2d ab_estimate = v->estimate();
        LOG(INFO) << "estimated model: " << ab_estimate.transpose();
        a_estimated.push_back(ab_estimate(0));
        b_estimated.push_back(ab_estimate(1));
        error_estimated.push_back(optimizer.activeChi2());
    }
    chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << "solve time cost = " << time_used.count() << " seconds. " << endl;

    // visualize cost function
    auto [A, B] = matplot::meshgrid(matplot::linspace(ar - 1, ar + 3, N),
                                    matplot::linspace(br - 2, br + 1, N));
    auto Error_sum = matplot::zeros(N, N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double error_sum = 0;
            for (int k = 0; k < N; k++)
            {
                double error = y_data[k] - exp(A[i][j] * x_data[k] * x_data[k] + B[i][j] * x_data[k]);
                error_sum += error * error * 0.5;
            }
            Error_sum[i][j] = error_sum;
        }
    }
    auto s = matplot::surf(A, B, Error_sum);
    // matplot::show();
    matplot::hold(matplot::on);
    matplot::plot3(a_estimated, b_estimated, error_estimated, "-or");

    // matplot::fplot3(
    matplot::hold(matplot::off);
    matplot::show();

    return 0;
}
