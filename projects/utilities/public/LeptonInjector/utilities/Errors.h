#ifndef LI_InjectionFailure_H
#define LI_InjectionFailure_H

namespace LI {
namespace utilities {

class InjectionFailure : public std::runtime_error {
public:
    InjectionFailure() : std::runtime_error("") {};
    InjectionFailure(const std::string& s) : std::runtime_error(s) {};
    InjectionFailure(const char * s) : std::runtime_error(s) {};
};

}
}

#endif
