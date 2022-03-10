#include <cmath>
#include <math.h>
#include <memory>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/array.hpp>

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"

using namespace LeptonInjector;

std::vector<double> logspace(double x_min, double x_max, unsigned int n_edges) {
    double log_x_min = log10(x_min);
    double log_x_max = log10(x_max);

    double delta_log_x = (log_x_max - log_x_min) / double(int(n_edges) - 1);
    unsigned int half = (unsigned int)(double(n_edges) / 2.0);

    std::vector<double> edges;
    edges.reserve(n_edges);
    double edge;

    for(unsigned int i=0; i < half; ++i) {
        edge = log_x_min + i*delta_log_x;
        edge = std::pow(10.0, edge);
        edges.push_back(edge);
    }
    for(unsigned int i=half; i < n_edges; ++i) {
        edge = log_x_max + (i - int(n_edges) + 1)*delta_log_x;
        edge = std::pow(10.0, edge);
        edges.push_back(edge);
    }
    return edges;
}

std::vector<double> linspace(double x_min, double x_max, unsigned int n_edges) {

    double delta_x = (x_max - x_min) / double(int(n_edges) - 1);
    unsigned int half = (unsigned int)(double(n_edges) / 2.0);

    std::vector<double> edges;
    edges.reserve(n_edges);
    double edge;

    for(unsigned int i=0; i < half; ++i) {
        edge = x_min + i*delta_x;
        edges.push_back(edge);
    }
    for(unsigned int i=half; i < n_edges; ++i) {
        edge = x_max + (i - int(n_edges) + 1)*delta_x;
        edges.push_back(edge);
    }
    return edges;
}

TEST(Table1D, Constructor)
{
    double x_min = 1e-5;
    double x_max = 1.0;
    unsigned int n_divisions = 10;

    LeptonInjector::TableData1D<double> table_data;
    table_data.x = logspace(x_min, x_max, n_divisions+1);
    table_data.f = std::vector<double>(n_divisions+1, 1.0);
    std::cerr << table_data.f[0] << std::endl;
    for(unsigned int i=1; i<table_data.f.size(); ++i) {
        table_data.f[i] = table_data.f[i-1] + 2;
        std::cerr << table_data.f[i] << std::endl;
    }

    LeptonInjector::Interpolator1D<double> interp(table_data);

    std::cerr << "IsLog " << interp.IsLog() << std::endl;
    std::cerr << "MinX " << interp.MinX() << std::endl;
    std::cerr << "MaxX " << interp.MaxX() << std::endl;
    std::cerr << "RangeX " << interp.RangeX() << std::endl;

    std::vector<double> sample_points = logspace(1e-6, 2.0, 10+1);
    for(auto const & i : sample_points) {
        std::cerr << i << ": " << interp(i) << std::endl;
    }

}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

