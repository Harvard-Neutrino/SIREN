#ifndef SIREN_InjectionFailure_H
#define SIREN_InjectionFailure_H

#include <string>
#include <stdexcept>

// Exception classes are thrown from libSIREN and caught by exception
// translators registered in the pybind11 extension modules, which are
// compiled with hidden symbol visibility.  On Itanium-ABI platforms
// (macOS/Linux) RTTI equality is pointer-based, so default visibility on
// the class is required for the typeinfo to unify across shared-library
// boundaries; without it the typed catch fails and Python sees a plain
// RuntimeError.  MSVC compares RTTI by name, so no export attribute is
// needed on Windows.
#if defined(_WIN32)
#define SIREN_EXCEPTION_EXPORT
#else
#define SIREN_EXCEPTION_EXPORT __attribute__((visibility("default")))
#endif

namespace siren {
namespace utilities {

class SIREN_EXCEPTION_EXPORT InjectionFailure : public std::runtime_error {
public:
    InjectionFailure() : std::runtime_error("") {};
    InjectionFailure(const std::string& s) : std::runtime_error(s) {};
    InjectionFailure(const char * s) : std::runtime_error(s) {};
};

class SIREN_EXCEPTION_EXPORT AddProcessFailure : public std::runtime_error {
public:
    AddProcessFailure() : std::runtime_error("") {};
    AddProcessFailure(const std::string& s) : std::runtime_error(s) {};
    AddProcessFailure(const char * s) : std::runtime_error(s) {};
};

class SIREN_EXCEPTION_EXPORT SecondaryProcessFailure : public std::runtime_error {
public:
    SecondaryProcessFailure() : std::runtime_error("") {};
    SecondaryProcessFailure(const std::string& s) : std::runtime_error(s) {};
    SecondaryProcessFailure(const char * s) : std::runtime_error(s) {};
};

class SIREN_EXCEPTION_EXPORT PythonImplementationError : public std::runtime_error {
public:
    PythonImplementationError() : std::runtime_error("") {};
    PythonImplementationError(const std::string& s) : std::runtime_error(s) {};
    PythonImplementationError(const char * s) : std::runtime_error(s) {};
};

class SIREN_EXCEPTION_EXPORT ConfigurationError : public std::runtime_error {
public:
    ConfigurationError(const std::string& s) : std::runtime_error(s) {};
};

class SIREN_EXCEPTION_EXPORT MeasureCompatibilityError : public std::runtime_error {
public:
    MeasureCompatibilityError(const std::string& s) : std::runtime_error(s) {};
};

class SIREN_EXCEPTION_EXPORT WeightCalculationError : public std::runtime_error {
public:
    WeightCalculationError(const std::string& s) : std::runtime_error(s) {};
};

}
}

#endif
