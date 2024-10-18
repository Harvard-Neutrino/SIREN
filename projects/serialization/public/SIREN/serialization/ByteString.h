#ifndef SIREN_serialization_ByteString_H
#define SIREN_serialization_ByteString_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <cereal/cereal.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

namespace siren {
namespace serialization {

template<typename T>
std::vector<char> to_byte_string(T const & obj) {
    std::ostringstream oss(std::ios::binary);
    {
        cereal::BinaryOutputArchive oarchive(oss);
        oarchive(obj);
    }
    size_t n_bytes = oss.tellp();
    std::string str = oss.str();
    if(not (str.size() == n_bytes))
        throw std::runtime_error("str.size() != n_bytes");
    std::vector<char> vec(str.begin(), str.end());
    if(not (vec.size() == n_bytes))
        throw std::runtime_error("vec.size() != n_bytes");
    return vec;
}

std::string bytes_to_hex_string(std::vector<char> const & bytes) {
    std::cout << "bytes_to_hex_string" << std::endl;
    std::cout << "bytes.size(): " << bytes.size() << std::endl;
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (char const & byte : bytes) {
        oss << std::setw(2) << (unsigned int)std::uint8_t(byte);
    }
    std::string s = oss.str();
    std::cout << "s.size(): " << s.size() << std::endl;
    if(not (s.size() % 2 == 0))
        throw std::runtime_error("s.size() % 2 != 0");
    if(not (s.size() == 2 * bytes.size()))
        throw std::runtime_error("s.size() != 2 * bytes.size()");
    return s;
}

template<typename T>
T from_byte_string(std::vector<char> const & byte_string) {
    std::cout << "from_byte_string" << std::endl;
    std::string str(byte_string.begin(), byte_string.end());
    if(not (str.size() == byte_string.size()))
        throw std::runtime_error("str.size() != byte_string.size()");
    std::cout << "str.size(): " << str.size() << std::endl;
    std::istringstream iss(str, std::ios::binary);
    std::shared_ptr<T> obj(cereal::access::construct<T>());
    {
        cereal::BinaryInputArchive iarchive(iss);
        iarchive(*obj);
    }
    size_t n_bytes = iss.tellg();
    if(not (n_bytes == str.size()))
        throw std::runtime_error("n_bytes != str.size()");
    std::cout << "obj: " << obj << std::endl;
    return *obj;
}

std::vector<char> hex_string_to_bytes(std::string const & hex_string) {
    std::cout << "hex_string_to_bytes" << std::endl;
    std::cout << "hex_string.size(): " << hex_string.size() << std::endl;
    std::vector<char> bytes;
    for (size_t i = 0; i < hex_string.size(); i += 2) {
        std::string byte_string = hex_string.substr(i, 2);
        char byte = static_cast<char>(std::stoi(byte_string, nullptr, 16));
        bytes.push_back(byte);
    }
    if(not (bytes.size() == hex_string.size() / 2))
        throw std::runtime_error("bytes.size() != hex_string.size() / 2");
    std::cout << "bytes.size(): " << bytes.size() << std::endl;
    return bytes;
}

template<typename T>
pybind11::tuple pickle_save(T const & cpp_obj) {
    std::vector<char> byte_string = to_byte_string(cpp_obj);
    std::string hex_string = bytes_to_hex_string(byte_string);
    std::cout << "Pickle save" << std::endl;
    std::cout << "hex_string.size(): " << hex_string.size() << std::endl;
    return pybind11::make_tuple(hex_string);
}

template<typename T>
T pickle_load(pybind11::tuple t) {
    std::cout << "pickle_load" << std::endl;
    std::cout << "t.size(): " << t.size() << std::endl;
    if (t.size() != 1) {
        throw std::runtime_error("Invalid state!");
    }
    std::string hex_string = t[0].cast<std::string>();
    std::cout << "hex_string.size(): " << hex_string.size() << std::endl;
    std::vector<char> byte_string = hex_string_to_bytes(hex_string);
    std::cout << "byte_string.size(): " << byte_string.size() << std::endl;
    T res = from_byte_string<T>(byte_string);
    std::cout << "res: " << res << std::endl;
    return res;
}

} // namespace serialization
} // namespace siren

#endif // SIREN_serialization_ByteString_H
