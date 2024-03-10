#ifndef SIREN_InjectionFailure_H
#define SIREN_InjectionFailure_H

namespace siren {
namespace utilities {

class InjectionFailure : public std::runtime_error {
public:
    InjectionFailure() : std::runtime_error("") {};
    InjectionFailure(const std::string& s) : std::runtime_error(s) {};
    InjectionFailure(const char * s) : std::runtime_error(s) {};
};

class AddProcessFailure : public std::runtime_error {
public:
    AddProcessFailure() : std::runtime_error("") {};
    AddProcessFailure(const std::string& s) : std::runtime_error(s) {};
    AddProcessFailure(const char * s) : std::runtime_error(s) {};
};

class SecondaryProcessFailure : public std::runtime_error {
public:
    SecondaryProcessFailure() : std::runtime_error("") {};
    SecondaryProcessFailure(const std::string& s) : std::runtime_error(s) {};
    SecondaryProcessFailure(const char * s) : std::runtime_error(s) {};
};

class PythonImplementationError : public std::runtime_error {
public:
    PythonImplementationError() : std::runtime_error("") {};
    PythonImplementationError(const std::string& s) : std::runtime_error(s) {};
    PythonImplementationError(const char * s) : std::runtime_error(s) {};
};

}
}

#endif
