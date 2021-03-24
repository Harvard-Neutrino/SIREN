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

    std::vector<std::string> material_names;

    void set_material_names(std::vector<std::string> const & names) {
        material_names = names;
    }

    void remove_file() {
        if(file_exists) {
            std::remove(model_file.c_str());
        }
    }

    void create_file() {
        if(file_exists) {
            remove_file();
        }
        uniform_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
        blank_chars = " \t";
        comment_str = "#";
        line_break = "\n";
        name_char_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        cruft_char_set = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
        model_file = std::tmpnam(nullptr);

        unsigned int n_layers = RandomDouble()*23+1;
        unsigned int n_lines = RandomDouble()*40+1;
        std::ofstream out(model_file.c_str(), std::ofstream::out);
        std::stringstream ss;

        unsigned int lines = 0;
        unsigned int layers = 0;
        double max_radius = 0;
        while(lines < n_lines and layers < n_layers) {
            unsigned int type = RandomDouble()*3;
            if(type == 0) {
                ss << blank_line();
            }
            else if(type == 1) {
                ss << comment_line();
            }
            else if(type == 2) {
                ss << add_layer(max_radius);
                layers += 1;
            }
            lines += 1;
        }
        file_contents = ss.str();
        out << file_contents;
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

    std::string add_layer(double & prev_radius) {
        std::stringstream out;

        double radius;
        std::string name;
        std::string material_name;
        unsigned int n_parameters;
        std::vector<double> parameters;

        name = random_name();

        double layer_thickness = RandomDouble()*1e6;
        radius = prev_radius + layer_thickness;

        int material_id = RandomDouble()*material_names.size();
        material_name = material_names[material_id];

        bool is_homogenous = RandomDouble()*2;
        if(is_homogenous) {
            n_parameters = 1;
        } else {
            n_parameters = RandomDouble()*6+1;
        }

        double parameter_sum = 0.0;
        for(unsigned int i=0; i<n_parameters; ++i) {
            double parameter_value = (RandomDouble()*(parameter_sum + 15) - parameter_sum);
            parameters.push_back(parameter_value);
            parameter_sum += parameter_value;
            if(parameter_sum < 0)
                parameter_sum = 0;
        }
        Polynom poly(parameters);
        poly.scale(1.0/layer_thickness);
        poly.shift(-prev_radius/layer_thickness);
        EXPECT_TRUE(poly.evaluate(prev_radius) >= 0);
        EXPECT_TRUE(poly.evaluate(radius) >= 0);

        parameters = poly.GetCoefficient();

        out << random_blank_cruft();
        out << radius;
        out << random_blank_cruft();
        out << name;
        out << random_blank_cruft();
        out << material_name;
        out << random_blank_cruft();
        out << n_parameters;
        for(unsigned int i=0; i<n_parameters; ++i) {
            out << random_blank_cruft();
            out << parameters[i];
        }

        out << line_break;

        prev_radius = radius;
        return out.str();
    }

    double RandomDouble() {
        return uniform_distribution(rng_);
    }
};

class FakeLegacyEarthModelTest : public FakeLegacyEarthModelFile, public FakeMaterialModelFile, public ::testing::Test {
protected:
    void SetUp() override {
        FakeMaterialModelFile::file_exists = false;
        FakeMaterialModelFile::create_file();
        FakeLegacyEarthModelFile::set_material_names(FakeMaterialModelFile::material_names_);
        FakeLegacyEarthModelFile::file_exists = false;
        FakeLegacyEarthModelFile::create_file();
    }
    void TearDown() override {
        FakeMaterialModelFile::remove_file();
        FakeLegacyEarthModelFile::remove_file();
    }
};

#endif // LI_TEST_FakeEarthModel_H

