#ifndef SIREN_Pybind11Trampoline_H
#define SIREN_Pybind11Trampoline_H

#include <pybind11/pybind11.h>

#define SELF_OVERRIDE_PURE(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        const BaseType * ref; \
        if(selfname) { \
            ref = selfname.cast<BaseType *>(); \
        } else { \
            ref = this; \
        } \
        do { \
            do { \
                pybind11::gil_scoped_acquire gil; \
                pybind11::function override \
                    = pybind11::get_override(static_cast<const BaseType *>(ref), pyfuncname); \
                if (override) { \
                    auto o = override(__VA_ARGS__); \
                    if (pybind11::detail::cast_is_temporary_value_reference<returnType>::value) { \
                        static pybind11::detail::override_caster_t<returnType> caster; \
                        return pybind11::detail::cast_ref<returnType>(std::move(o), caster); \
                    } \
                    return pybind11::detail::cast_safe<returnType>(std::move(o)); \
                } \
            } while (false); \
            pybind11::pybind11_fail( \
                "Tried to call pure virtual function \"" PYBIND11_STRINGIFY(BaseType) "::" #cfuncname "\""); \
        } while (false);

#define SELF_OVERRIDE(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        const BaseType * ref; \
        if(selfname) { \
            ref = selfname.cast<BaseType *>(); \
        } else { \
            ref = this; \
        } \
        do { \
            do { \
                pybind11::gil_scoped_acquire gil; \
                pybind11::function override \
                    = pybind11::get_override(static_cast<const BaseType *>(ref), pyfuncname); \
                if (override) { \
                    auto o = override(__VA_ARGS__); \
                    if (pybind11::detail::cast_is_temporary_value_reference<returnType>::value) { \
                        static pybind11::detail::override_caster_t<returnType> caster; \
                        return pybind11::detail::cast_ref<returnType>(std::move(o), caster); \
                    } \
                    return pybind11::detail::cast_safe<returnType>(std::move(o)); \
                } \
            } while (false); \
            return BaseType::cfuncname(__VA_ARGS__); \
        } while (false);

#endif // SIREN_Pybind11Trampoline_H
