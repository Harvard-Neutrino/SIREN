#ifndef LI_TEST_FakeEarthModel_H
#define LI_TEST_FakeEarthModel_H

#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include "FakeMaterialModel.h"

using namespace earthmodel;

template <typename T>
std::string str(std::vector<T> v) {
    std::stringstream ss;
    for(auto const & s : v) {
        ss << s << " ";
    }
    return ss.str();
}

class FakeLegacyEarthModelFile {
protected:
    bool file_exists;
    std::mt19937 rng_;
    std::string blank_chars;
    std::string comment_str;
    std::string line_break;
    std::string name_char_set;
    std::string cruft_char_set;
    std::uniform_real_distribution<double> uniform_distribution;
    std::string model_file;

    std::string file_contents;

    std::vector<std::string> layer_names;
    std::vector<double> layer_thicknesses;
    std::vector<double> layer_radii;
    std::vector<Polynom> layer_densities;
    std::vector<std::string> layer_materials;

    std::vector<std::string> material_names;

    void set_material_names(std::vector<std::string> const & names) {
        material_names = names;
    }

    void remove_file() {
        if(file_exists) {
            std::remove(model_file.c_str());
            clear();
        }
        file_exists = false;
    }

    void create_file(std::vector<std::string> mat_names) {
        unsigned int n_layers = RandomDouble()*23+1;
        unsigned int n_lines = RandomDouble()*40+1;
        int poly_max = 6;
        create_file(mat_names, n_layers, n_lines, poly_max);
    }

    void create_file(std::vector<std::string> mat_names, unsigned int n_layers, unsigned int n_lines, int poly_max) {
        if(file_exists) {
            remove_file();
        }
        set_material_names(mat_names);
        uniform_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
        blank_chars = " \t";
        comment_str = "#";
        line_break = "\n";
        name_char_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        cruft_char_set = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
        model_file = std::tmpnam(nullptr);

        std::ofstream out(model_file.c_str(), std::ofstream::out);
        std::stringstream ss;

        unsigned int lines = 0;
        unsigned int layers = 0;
        double max_radius = 0;
        while((lines < n_lines) or (layers < n_layers)) {
            unsigned int type = RandomDouble()*3;
            if(type == 0) {
                ss << blank_line();
            }
            else if(type == 1) {
                ss << comment_line();
            }
            else if(type == 2) {
                ss << add_layer(max_radius, poly_max);
                layers += 1;
            }
            lines += 1;
        }
        file_contents = ss.str();
        out << file_contents;
        file_exists = true;
    }

    std::string blank_line() {
        unsigned int n_chars = RandomDouble()*5;
        std::stringstream s;
        for(unsigned int i=0; i<n_chars; ++i) {
            s << blank_chars[(unsigned int)(RandomDouble()*blank_chars.size())];;
        }
        s << line_break;
        return s.str();
    }

    std::string comment_line() {
        unsigned int n_chars = RandomDouble()*1024;
        std::stringstream s;
        s << comment_str;
        for(unsigned int i=0; i<n_chars; ++i) {
            s << cruft_char_set[(unsigned int)(RandomDouble()*cruft_char_set.size())];;
        }
        s << line_break;
        return s.str();
    }

    std::string random_name() {
        unsigned int n_chars = RandomDouble()*46 + 1;
        std::stringstream s;
        for(unsigned int i=0; i<n_chars; ++i)
            s << name_char_set[(unsigned int)(RandomDouble()*name_char_set.size())];
        return s.str();
    }

    std::string random_cruft() {
        unsigned int n_chars = RandomDouble()*46;
        std::stringstream s;
        for(unsigned int i=0; i<n_chars; ++i)
            s << cruft_char_set[(unsigned int)(RandomDouble()*cruft_char_set.size())];
        return s.str();
    }

    std::string random_blank_cruft() {
        unsigned int n_chars = RandomDouble()*5 + 1;
        std::stringstream s;
        s << ' ';
        for(unsigned int i=0; i<n_chars; ++i) {
            s << blank_chars[(unsigned int)(RandomDouble()*blank_chars.size())];;
        }
        return s.str();
    }

    std::string add_layer(double & prev_radius, int poly_max) {
        std::stringstream out;

        double radius;
        std::string name;
        std::string material_name;
        unsigned int n_parameters;
        std::vector<double> parameters;

        name = random_name();
        layer_names.push_back(name);

        double layer_thickness = RandomDouble()*1e6;
        layer_thicknesses.push_back(layer_thickness);
        radius = prev_radius + layer_thickness;
        layer_radii.push_back(radius);

        assert(material_names.size() > 0);
        int material_id = RandomDouble()*material_names.size();
        material_name = material_names[material_id];
        layer_materials.push_back(material_name);

        int is_homogenous = RandomDouble()*2;
        if(is_homogenous or poly_max <= 1) {
            n_parameters = 1;
        } else {
            n_parameters = RandomDouble()*std::max(0, poly_max-1)+1;
        }

        double parameter_sum = 0.0;
        double starting_max = 15;
        double starting_min = 0.0;
        double upper_bound;
        double lower_bound;
        for(unsigned int i=0; i<n_parameters; ++i) {
            if(i>0) {
                upper_bound = std::min(std::abs(parameter_sum), std::pow(starting_max, 1.0/(10*i))/(std::pow(i, 10)));
                lower_bound = std::max(std::min(parameter_sum, starting_min), -upper_bound);
            } else {
                upper_bound = starting_max;
                lower_bound = starting_min;
            }
            double parameter_value = (RandomDouble()*(upper_bound - lower_bound) + lower_bound);
            parameters.push_back(parameter_value);
            parameter_sum += parameter_value;
            if(parameter_sum < 0)
                parameter_sum = 0;
        }
        Polynom poly(parameters);
        poly.scale(1.0/layer_thickness);
        poly.shift(-prev_radius);
        bool poly_positive;
        poly_positive = poly.evaluate(prev_radius) >= 0;
        assert(poly_positive);
        for(unsigned int i=0; i<1000; ++i) {
            poly_positive = poly.evaluate(RandomDouble()*layer_thickness + prev_radius) >= 0;
            assert(poly_positive);
        }
        poly_positive = poly.evaluate(radius) >= 0;
        assert(poly_positive);
        layer_densities.push_back(poly);

        parameters = poly.GetCoefficient();

        out << random_blank_cruft();
        out << std::setprecision(16) << radius;
        out << random_blank_cruft();
        out << name;
        out << random_blank_cruft();
        out << material_name;
        out << random_blank_cruft();
        out << n_parameters;
        for(unsigned int i=0; i<n_parameters; ++i) {
            out << random_blank_cruft();
            out << std::setprecision(16) << parameters[i];
        }

        out << line_break;

        prev_radius = radius;
        return out.str();
    }

    void clear() {
        file_contents.clear();

        layer_names.clear();
        layer_thicknesses.clear();
        layer_radii.clear();
        layer_densities.clear();
        layer_materials.clear();

        material_names.clear();
        model_file = "";
    }

    double RandomDouble() {
        return uniform_distribution(rng_);
    }

    Vector3D RandomVector(double max_radius, double min_radius=0.0) {
        double radius = RandomDouble()*(max_radius-min_radius) + min_radius;
        double theta = RandomDouble()*M_PI;
        double phi = RandomDouble()*2.0*M_PI;
        Vector3D result;
        result.SetSphericalCoordinates(radius, phi, theta);
        return result;
    }
};

class FakeLegacyEarthModelTest : public FakeLegacyEarthModelFile, public FakeMaterialModelFile, public ::testing::Test {
protected:
    void setup() {
        FakeMaterialModelFile::create_file(std::vector<std::string>{"air", "atmosphere", "ice"}, 6);
        FakeLegacyEarthModelFile::create_file(FakeMaterialModelFile::material_names_);
    }
    void setup(unsigned int n_layers, int poly_max) {
        FakeMaterialModelFile::create_file(std::vector<std::string>{"air", "atmosphere", "ice"}, 6);
        FakeLegacyEarthModelFile::create_file(FakeMaterialModelFile::material_names_, n_layers, 0, poly_max);
    }
    void reset() {
        setup();
    }
    void reset(unsigned int n_layers, int poly_max) {
        setup(n_layers, poly_max);
    }
    void SetUp() override {
        FakeMaterialModelFile::file_exists = false;
        FakeLegacyEarthModelFile::file_exists = false;
        setup();
    }
    void TearDown() override {
        //FakeMaterialModelFile::remove_file();
        //FakeLegacyEarthModelFile::remove_file();
        FakeMaterialModelFile::clear();
        FakeLegacyEarthModelFile::clear();
    }
};

#endif // LI_TEST_FakeEarthModel_H

