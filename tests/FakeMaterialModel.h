#ifndef LI_TEST_FakeMaterialModel_H
#define LI_TEST_FakeMaterialModel_H

#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

using namespace earthmodel;


class FakeMaterialModelFile {
public:
    bool file_exists;
    std::mt19937 rng_;
    std::string blank_chars;
    std::string comment_str;
    std::string line_break;
    std::string name_char_set;
    std::string cruft_char_set;
    std::uniform_real_distribution<double> uniform_distribution;
    std::string materials_file;
    std::vector<std::string> possible_material_names;

    std::string file_contents;

    std::vector<std::string> material_names_;
    std::vector<double> pne_ratios_;
    std::map<std::string, int> material_ids_;
    std::map<std::string, std::vector<std::pair<std::string, double>>> material_components_;
    unsigned int n_materials_added_ = 0;

    void set_possible_material_names(std::vector<std::string> names, unsigned int n_random_names) {
        possible_material_names = names;
        for(unsigned int i=0; i<n_random_names; ++i) {
            possible_material_names.push_back(random_name());
        }
    }

    void remove_file() {
        if(file_exists) {
            std::remove(materials_file.c_str());
            clear();
        }
        file_exists = false;
    }
    void create_file(std::vector<std::string> names, unsigned int n_random_names) {
        if(file_exists) {
            remove_file();
        }
        set_possible_material_names(names, n_random_names);
        assert(material_names_.size() == 0);
        uniform_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
        blank_chars = " \t";
        comment_str = "#";
        line_break = "\n";
        name_char_set = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
        cruft_char_set = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
        materials_file = std::tmpnam(nullptr);

        unsigned int n_materials = RandomDouble()*23+1;
        unsigned int n_lines = RandomDouble()*40+1;
        std::ofstream out(materials_file.c_str(), std::ofstream::out);
        std::stringstream ss;

        unsigned int lines = 0;
        unsigned int materials = 0;
        while((lines < n_lines) or (materials < n_materials)) {
            unsigned int type = RandomDouble()*3;
            if(type == 0) {
                ss << blank_line();
            }
            else if(type == 1) {
                ss << comment_line();
            }
            else if(type == 2) {
                ss << add_material();
                materials += 1;
            }
            lines += 1;
        }
        file_contents = ss.str();
        out << file_contents;
        file_exists = true;
        assert(n_materials_added_ > 0);
        assert(material_names_.size() > 0);
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

    std::string pdg_code(unsigned int L, unsigned int Z, unsigned int A, unsigned int I) {
        const unsigned int buf_size = 64;
        char buf[buf_size];
        std::snprintf(buf, buf_size, "10%01u%03u%03u%01u", L, Z, A, I);
        //return std::format("10{:01d}{:03d}{:03d}{:01d}", L, Z, A, I);
        return std::string(buf);
    }

    std::string random_pdg_code() {
        unsigned int L = 0;
        unsigned int A = RandomDouble()*998 + 2;
        unsigned int Z = RandomDouble()*A;
        unsigned int I = 0;
        return pdg_code(L, Z, A, I);
    }

    std::string random_name() {
        unsigned int n_chars = RandomDouble()*46 + 1;
        std::stringstream s;
        for(unsigned int i=0; i<n_chars; ++i)
            s << name_char_set[(int)(RandomDouble()*name_char_set.size())];
        return s.str();
    }

    std::string random_cruft() {
        unsigned int n_chars = RandomDouble()*46;
        std::stringstream s;
        for(unsigned int i=0; i<n_chars; ++i)
            s << cruft_char_set[(int)(RandomDouble()*cruft_char_set.size())];
        return s.str();
    }

    std::string random_blank_cruft() {
        unsigned int n_chars = RandomDouble()*5 + 1;
        std::stringstream s;
        s << ' ';
        for(unsigned int i=0; i<n_chars; ++i) {
            s << blank_chars[(int)(RandomDouble()*blank_chars.size())];;
        }
        return s.str();
    }

    std::string add_material() {
        n_materials_added_ += 1;
        std::stringstream out;
        std::string name;
        if(possible_material_names.size() > 0) {
            int material_name_index = RandomDouble()*possible_material_names.size();
            name = possible_material_names[material_name_index];
            possible_material_names.erase(possible_material_names.begin()+material_name_index);
        } else {
            name = random_name();
        }
        int id;
        if(material_ids_.find(name) == material_ids_.end()) {
            id = material_names_.size();
            material_ids_.insert({name, id});
            material_names_.push_back(name);
        } else {
            assert(material_names_.size() > 0);
            id = material_ids_[name];
        }
        unsigned int n_components = RandomDouble()*11;
        double component_total = 0.0;
        std::vector<std::pair<std::string, double>> components;
        double Z_tot = 0.0;
        double A_tot = 0.0;
        std::set<std::string> pdg_codes;
        for(unsigned int i=0; i<n_components; ++i) {
            double amount = RandomDouble() * ((unsigned int)(RandomDouble()*46));
            std::string pdg = random_pdg_code();
            while(pdg_codes.find(pdg) != pdg_codes.end()) {
                pdg = random_pdg_code();
            }
            pdg_codes.insert(pdg);
            std::stringstream ss;
            std::string s;
            ss << pdg[3] << pdg[4] << pdg[5];
            s = ss.str();
            ss = std::stringstream(s);
            int Z;
            ss >> Z;
            ss.clear();
            ss = std::stringstream();
            ss << pdg[6] << pdg[7] << pdg[8];
            s = ss.str();
            ss = std::stringstream(s);
            int A;
            ss >> A;
            ss.clear();
            Z_tot += Z*amount;
            A_tot += A*amount;
            component_total += amount;
            components.push_back({pdg, amount});
        }
        if(A_tot == 0) {
            A_tot = 1;
        }
        double pne_ratio = Z_tot / A_tot;
        pne_ratios_.push_back(pne_ratio);

        bool normalize = RandomDouble() > 0.5;
        bool add_cruft = RandomDouble() > 0.5;
        out << name << random_blank_cruft() << n_components;
        if(add_cruft)
            out << random_blank_cruft() << comment_str << random_blank_cruft() << random_cruft();
        out << line_break;

        for(unsigned int i=0; i<n_components; ++i) {
            add_cruft = RandomDouble() > 0.5;
            std::string pdg = components[i].first;
            EXPECT_EQ(pdg[0], '1');
            EXPECT_EQ(pdg[1], '0');
            if(normalize) {
                components[i].second /= component_total;
            }
            double amount = components[i].second;
            out << pdg << random_blank_cruft() << std::setprecision(16) << amount;
            if(add_cruft)
                out << random_blank_cruft() << comment_str << random_blank_cruft() << random_cruft();
            out << line_break;
        }

        assert(n_components == components.size());

        if(material_components_.find(name) == material_components_.end()) {
            material_components_.insert({name, components});
        } else {
            material_components_[name] = components;
        }
        assert(material_components_[name].size() == n_components);
        return out.str();
    }

    double RandomDouble() {
        return uniform_distribution(rng_);
    }

    void clear() {
        material_ids_.clear();
        material_names_.clear();
        material_components_.clear();
        pne_ratios_.clear();
        materials_file = "";
        n_materials_added_ = 0;
    }

};

class FakeMaterialModelTest : public FakeMaterialModelFile, public ::testing::Test {
protected:
    void SetUp() override {
        FakeMaterialModelFile::file_exists = false;
        FakeMaterialModelFile::create_file({}, 12);
    }
    void TearDown() override {
        FakeMaterialModelFile::remove_file();
    }
};

#endif // LI_TEST_FakeMaterialModel_H

