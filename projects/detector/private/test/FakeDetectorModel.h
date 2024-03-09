#ifndef LI_TEST_FakeDetectorModel_H
#define LI_TEST_FakeDetectorModel_H

#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include "FakeMaterialModel.h"

#include "SIREN/geometry/Geometry.h"
#include "SIREN/detector/DetectorModel.h"
#include "SIREN/math/Quaternion.h"
#include "SIREN/math/Polynomial.h"

using namespace siren::detector;
using namespace siren::math;

template <typename T>
std::string str(std::vector<T> v) {
    std::stringstream ss;
    for(auto const & s : v) {
        ss << s << " ";
    }
    return ss.str();
}

struct ModelToken {
    enum Token {none, line, newline, comment, entry, string_comment, object, detector, density_distribution, material_name, label, shape, sphere, box, shape_coords, shape_angles, constant, radial_polynomial, polynomial, axis_coords};
    Token tok = none;
    int state = 0;

    ModelToken(Token tok) : tok(tok), state(0) {
    }

    ModelToken & operator++() {
        state += 1;
        return *this;
    }
};


class FakeDetectorModelFile {
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
    std::vector<Polynom> layer_densities;
    std::vector<std::string> layer_materials;

    std::vector<std::string> material_names;

    std::vector<ModelToken> stack;

    double max_radius;

    inline
    std::string token_line(ModelToken & token) {
        stack.push_back(ModelToken(ModelToken::newline));
        stack.push_back(ModelToken(ModelToken::comment));
        stack.push_back(ModelToken(ModelToken::entry));
        return "";
    }

    inline
    std::string token_comment(ModelToken & token) {
        int choice = int(uniform_distribution(rng_) * 2);
        switch(choice) {
            case 0:
                break;
            case 1:
                stack.push_back(ModelToken(ModelToken::string_comment));
                break;
            default:
                break;
        }
        return "";
    }

    inline
    std::string token_string_comment(ModelToken & token) {
        unsigned int n_chars = RandomDouble()*1024;
        std::stringstream s;
        s << comment_str;
        for(unsigned int i=0; i<n_chars; ++i) {
            s << cruft_char_set[(unsigned int)(RandomDouble()*cruft_char_set.size())];;
        }
        return s.str();
    }

    inline
    std::string token_break(ModelToken & token) {
        return "\n";
    }

    inline
    std::string token_entry(ModelToken & token) {
        int choice = int(uniform_distribution(rng_) * 3);
        switch(choice) {
            case 0:
                break;
            case 1:
                stack.push_back(ModelToken(ModelToken::object));
                break;
            case 2:
                stack.push_back(ModelToken(ModelToken::detector));
                break;
            default:
                break;
        }
        return "";
    }

    inline
    std::string token_detector(ModelToken & token) {
        Vector3D v = RandomVector(max_radius);
        std::stringstream ss;
        ss << "detector " << v.GetX() << random_blank_cruft() << v.GetY() << random_blank_cruft() << v.GetZ();
        return ss.str();
    }

    inline
    std::string token_object(ModelToken & token) {
        stack.push_back(ModelToken(ModelToken::density_distribution));
        stack.push_back(ModelToken(ModelToken::material_name));
        stack.push_back(ModelToken(ModelToken::label));
        stack.push_back(ModelToken(ModelToken::shape));
        return "object" + random_blank_cruft();
    }

    inline
    std::string token_shape(ModelToken & token) {
        int choice = int(uniform_distribution(rng_) * 2);
        switch(choice) {
            case 0:
                stack.push_back(ModelToken(ModelToken::sphere));
                break;
            case 1:
                stack.push_back(ModelToken(ModelToken::box));
                break;
            default:
                break;
        }
        return "";
    }

    inline
    std::string token_sphere(ModelToken & token) {
        double r;
        std::stringstream ss;
        switch(token.state) {
            case 0:
                stack.push_back(++token);
                stack.push_back(ModelToken(ModelToken::shape_coords));
                stack.push_back(ModelToken(ModelToken::shape_angles));
                return "sphere" + random_blank_cruft();
            case 1:
                r = RandomDouble() * max_radius;
                ss << r << random_blank_cruft();
                return ss.str();
            default:
                return "";
        }
    }

    inline
    std::string token_box(ModelToken & token) {
        Vector3D v;
        std::stringstream ss;
        switch(token.state) {
            case 0:
                stack.push_back(++token);
                stack.push_back(ModelToken(ModelToken::shape_coords));
                stack.push_back(ModelToken(ModelToken::shape_angles));
                return "box" + random_blank_cruft();
            case 1:
                v = RandomVector(max_radius);
                ss << v.GetX() << random_blank_cruft() << v.GetY() << random_blank_cruft() << v.GetZ() << random_blank_cruft();
                return ss.str();
            default:
                return "";
        }
    }

    inline
    std::string token_shape_coords(ModelToken & token) {
        Vector3D v = RandomVector(max_radius);
        std::stringstream ss;
        ss << v.GetX() << random_blank_cruft() << v.GetY() << random_blank_cruft() << v.GetZ() << random_blank_cruft();
        return ss.str();
    }

    std::string token_shape_angles(ModelToken & token) {
        Vector3D v = RandomVector(1.0);
        double theta = RandomDouble() * 2.0 * M_PI;
        Quaternion q;
        q.SetAxisAngle(v, theta);
        double alpha, beta, gamma;
        q.GetEulerAnglesZXZr(alpha, beta, gamma);
        std::stringstream ss;
        ss << alpha << random_blank_cruft() << beta << random_blank_cruft() << gamma << random_blank_cruft();
        return ss.str();
    }

    inline
    std::string token_label(ModelToken & token) {
        std::string name = random_name();
        layer_names.push_back(name);
        return name + random_blank_cruft();
    }

    inline
    std::string token_material_name(ModelToken & token) {
        assert(material_names.size() > 0);
        int material_id = RandomDouble()*material_names.size();
        std::string material_name = material_names[material_id];
        layer_materials.push_back(material_name);
        return material_name + random_blank_cruft();
    }

    inline
    std::string token_density_distribution(ModelToken & token) {
        int choice = int(uniform_distribution(rng_) * 2);
        switch(choice) {
            case 0:
                stack.push_back(ModelToken(ModelToken::constant));
                break;
            case 1:
                stack.push_back(ModelToken(ModelToken::radial_polynomial));
                break;
            default:
                break;
        }
        return "";
    }

    inline
    std::string token_constant(ModelToken & token) {
        std::stringstream ss;
        ss << "constant" << random_blank_cruft() << RandomDouble() * 15 << " ";
        return ss.str();
    }

    inline
    std::string token_radial_polynomial(ModelToken & token) {
        stack.push_back(ModelToken(ModelToken::polynomial));
        stack.push_back(ModelToken(ModelToken::axis_coords));
        return "radial_polynomial" + random_blank_cruft();
    }

    inline
    std::string token_polynomial(ModelToken & token, unsigned int poly_max) {
        double n_parameters = RandomDouble()*std::max((unsigned int)0, poly_max-1)+1;

        double parameter_sum = 0.0;
        double starting_max = 15;
        double starting_min = 0.0;
        double upper_bound;
        double lower_bound;

        std::vector<double> parameters;

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
        poly.scale(1.0/(max_radius*2));
        // poly.shift(-prev_radius);
        bool poly_positive;
        poly_positive = poly.evaluate(0) >= 0;
        assert(poly_positive);
        for(unsigned int i=0; i<1000; ++i) {
            poly_positive = poly.evaluate(RandomDouble()*max_radius*2) >= 0;
            assert(poly_positive);
        }
        poly_positive = poly.evaluate(max_radius) >= 0;
        layer_densities.push_back(poly);

        parameters = poly.GetCoefficient();

        std::stringstream ss;
        ss << n_parameters;
        ss << random_blank_cruft();
        for(unsigned int i=0; i<n_parameters; ++i) {
            ss << random_blank_cruft();
            ss << std::setprecision(16) << parameters[i];
        }
        return ss.str();
    }


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
        unsigned int poly_max = 6;
        max_radius = 1e7;
        create_file(mat_names, n_layers, n_lines, poly_max);
    }

    void create_file(std::vector<std::string> mat_names, unsigned int n_layers, unsigned int n_lines, unsigned int poly_max) {
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

        stack.clear();

        unsigned int lines = 0;
        unsigned int layers = 0;
        while((lines < n_lines) or (layers < n_layers) or stack.size() > 0) {
            if(stack.size() == 0) {
                stack.push_back(ModelToken(ModelToken::line));
            } else {
                ModelToken token = stack.back();
                stack.pop_back();
                ModelToken::Token tok = token.tok;
                std::string next;
                switch(tok) {
                    case ModelToken::line:
                        next = token_line(token); ++lines; break;
                    case ModelToken::newline:
                        next = token_break(token); break;
                    case ModelToken::comment:
                        next = token_comment(token); break;
                    case ModelToken::string_comment:
                        next = token_string_comment(token); break;
                    case ModelToken::entry:
                        next = token_entry(token); break;
                    case ModelToken::object:
                        next = token_object(token); ++layers; break;
                    case ModelToken::detector:
                        next = token_detector(token); break;
                    case ModelToken::shape:
                        next = token_shape(token); break;
                    case ModelToken::sphere:
                        next = token_sphere(token); break;
                    case ModelToken::box:
                        next = token_box(token); break;
                    case ModelToken::shape_coords:
                        next = token_shape_coords(token); break;
                    case ModelToken::shape_angles:
                        next = token_shape_angles(token); break;
                    case ModelToken::axis_coords:
                        next = token_shape_coords(token); break;
                    case ModelToken::label:
                        next = token_label(token); break;
                    case ModelToken::material_name:
                        next = token_material_name(token); break;
                    case ModelToken::density_distribution:
                        next = token_density_distribution(token); break;
                    case ModelToken::constant:
                        next = token_constant(token); break;
                    case ModelToken::radial_polynomial:
                        next = token_radial_polynomial(token); break;
                    case ModelToken::polynomial:
                        next = token_polynomial(token, poly_max); break;
                    default:
                        throw; break;
                }
                ss << next;
            }
        }
        file_contents = ss.str();
        out << file_contents;
        file_exists = true;
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
        unsigned int n_chars = (int)(RandomDouble()*5 + 1);
        std::string s(n_chars, ' ');
        for(unsigned int i=0; i<n_chars; ++i) {
            s[i] = blank_chars[std::min((unsigned int)(RandomDouble()*blank_chars.size()), (unsigned int)(blank_chars.size()-1))];;
        }
        return s;
    }


    void clear() {
        file_contents.clear();

        layer_names.clear();
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
        result.CalculateCartesianFromSpherical();
        return result;
    }
};

class FakeDetectorModelTest : public FakeDetectorModelFile, public FakeMaterialModelFile, public ::testing::Test {
protected:
    void setup() {
        FakeMaterialModelFile::create_file(std::vector<std::string>{"air", "atmosphere", "ice"}, 6);
        FakeDetectorModelFile::create_file(FakeMaterialModelFile::material_names_);
    }
    void setup(unsigned int n_layers, int poly_max) {
        FakeMaterialModelFile::create_file(std::vector<std::string>{"air", "atmosphere", "ice"}, 6);
        FakeDetectorModelFile::create_file(FakeMaterialModelFile::material_names_, n_layers, 0, poly_max);
    }
    void reset() {
        setup();
    }
    void reset(unsigned int n_layers, int poly_max) {
        setup(n_layers, poly_max);
    }
    void SetUp() override {
        FakeMaterialModelFile::file_exists = false;
        FakeDetectorModelFile::file_exists = false;
        setup();
    }
    void TearDown() override {
        //FakeMaterialModelFile::remove_file();
        //FakeLegacyDetectorModelFile::remove_file();
        FakeMaterialModelFile::clear();
        FakeDetectorModelFile::clear();
    }
};

class FakeLegacyDetectorModelFile {
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
        result.CalculateCartesianFromSpherical();
        return result;
    }
};

class FakeLegacyDetectorModelTest : public FakeLegacyDetectorModelFile, public FakeMaterialModelFile, public ::testing::Test {
protected:
    void setup() {
        FakeMaterialModelFile::create_file(std::vector<std::string>{"air", "atmosphere", "ice"}, 6);
        FakeLegacyDetectorModelFile::create_file(FakeMaterialModelFile::material_names_);
    }
    void setup(unsigned int n_layers, int poly_max) {
        FakeMaterialModelFile::create_file(std::vector<std::string>{"air", "atmosphere", "ice"}, 6);
        FakeLegacyDetectorModelFile::create_file(FakeMaterialModelFile::material_names_, n_layers, 0, poly_max);
    }
    void reset() {
        setup();
    }
    void reset(unsigned int n_layers, int poly_max) {
        setup(n_layers, poly_max);
    }
    void SetUp() override {
        FakeMaterialModelFile::file_exists = false;
        FakeLegacyDetectorModelFile::file_exists = false;
        setup();
    }
    void TearDown() override {
        //FakeMaterialModelFile::remove_file();
        //FakeLegacyDetectorModelFile::remove_file();
        FakeMaterialModelFile::clear();
        FakeLegacyDetectorModelFile::clear();
    }
};

#endif // LI_TEST_FakeDetectorModel_H

