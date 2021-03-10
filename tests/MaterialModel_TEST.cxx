
#include <cmath>
#include <math.h>
#include <cstdio>
#include <string>
#include <random>
#include <fstream>
#include <iostream>

#include <gtest/gtest.h>

#include <earthmodel-service/MaterialModel.h>

using namespace earthmodel;

class MaterialTest : public ::testing::Test {
protected:
    void SetUp() override {
        file_exists = false;
        create_file();
    }
    void TearDown() override {
        remove_file();
    }
    void remove_file() {
        if(file_exists) {
            std::remove(materials_file.c_str());
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
        materials_file = std::tmpnam(nullptr);

        unsigned int n_materials = RandomDouble()*23+1;
        unsigned int n_lines = RandomDouble()*40+1;
        std::ofstream out(materials_file.c_str(), std::ofstream::out);
        std::stringstream ss;

        unsigned int lines = 0;
        unsigned int materials = 0;
        while(lines < n_lines and materials < n_materials) {
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
        unsigned int Z = RandomDouble()*1000;
        unsigned int A = RandomDouble()*1000;
        unsigned int I = 0;
        return pdg_code(L, Z, A, I);
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

    std::string add_material() {
        std::stringstream out;
        std::string name = random_name();
        int id;
        if(material_ids_.find(name) == material_ids_.end()) {
            id = material_names_.size();
            material_ids_.insert({name, id});
            material_names_.push_back(name);
        } else {
            id = material_ids_[name];
        }
        unsigned int n_components = RandomDouble()*11;
        double component_total = 0.0;
        std::vector<std::pair<std::string, double>> components;
        for(unsigned int i=0; i<n_components; ++i) {
            double amount = RandomDouble() * ((unsigned int)(RandomDouble()*46));
            std::string pdg = random_pdg_code();
            component_total += amount;
            components.push_back({pdg, amount});
        }

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
            out << pdg << random_blank_cruft() << amount;
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
    }

    bool file_exists;
    std::mt19937 rng_;
    std::string blank_chars;
    std::string comment_str;
    std::string line_break;
    std::string name_char_set;
    std::string cruft_char_set;
    std::uniform_real_distribution<double> uniform_distribution;
    std::string materials_file;

    std::string file_contents;

    std::vector<std::string> material_names_;
    std::map<std::string, int> material_ids_;
    std::map<std::string, std::vector<std::pair<std::string, double>>> material_components_;
};

TEST(Constructor, Default)
{
    MaterialModel A;
}

TEST_F(MaterialTest, FileConstructor)
{
    ASSERT_NO_THROW(MaterialModel("", materials_file));
    MaterialModel A("", materials_file);
}

TEST_F(MaterialTest, DuplicateFile)
{
    ASSERT_NO_THROW(MaterialModel("", materials_file));
    MaterialModel A("", materials_file);
    ASSERT_NO_THROW(A.AddModelFile(materials_file));
    ASSERT_NO_THROW(MaterialModel("", {materials_file, materials_file}));
    ASSERT_NO_THROW(A.HasMaterial(0));
}

std::string print_v(std::vector<std::string> v) {
    std::stringstream s;
    for(auto ss : v) {
        s << ss << " ";
    }
    return s.str();
}

TEST_F(MaterialTest, MaterialCount)
{
    ASSERT_NO_THROW(MaterialModel("", materials_file));

    unsigned int N_RAND = 1000;
    for(unsigned int i=0; i<N_RAND; ++i) {
        clear();
        create_file();
        MaterialModel A("", materials_file);

        int material_count = 0;
        std::vector<std::string> names;
        while(A.HasMaterial(material_count)) {
            names.push_back(A.GetMaterialName(material_count));
            EXPECT_EQ(std::string(names.back()), std::string(material_names_[material_count])) << file_contents;
            EXPECT_EQ(A.GetMaterialMap(material_count).size(), material_components_[material_names_[material_count]].size()) << file_contents;
            material_count += 1;
        }

        ASSERT_EQ(material_count, material_names_.size()) << file_contents;
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
