#pragma once
#ifndef SIREN_Pybind11Trampoline_H
#define SIREN_Pybind11Trampoline_H

#include <sstream>
#include <functional>
#include <cxxabi.h>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/utility.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>

#define SELF_OVERRIDE_PURE(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            if(!selfname) { \
                auto *_tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
                if(_tinfo) { \
                    pybind11::handle _h = pybind11::detail::get_object_handle(static_cast<const BaseType *>(this), _tinfo); \
                    if(_h) { \
                        selfname = pybind11::reinterpret_borrow<pybind11::object>(_h); \
                    } \
                } \
            } \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
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
            pybind11::pybind11_fail( \
                "Tried to call pure virtual function \"" PYBIND11_STRINGIFY(BaseType) "::" #cfuncname "\""); \
        } while (false);

#define SELF_OVERRIDE(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            if(!selfname) { \
                auto *_tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
                if(_tinfo) { \
                    pybind11::handle _h = pybind11::detail::get_object_handle(static_cast<const BaseType *>(this), _tinfo); \
                    if(_h) { \
                        selfname = pybind11::reinterpret_borrow<pybind11::object>(_h); \
                    } \
                } \
            } \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
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
        return BaseType::cfuncname(__VA_ARGS__);

#define SELF_OVERRIDE_PURE_REF(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            if(!selfname) { \
                auto *_tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
                if(_tinfo) { \
                    pybind11::handle _h = pybind11::detail::get_object_handle(static_cast<const BaseType *>(this), _tinfo); \
                    if(_h) { \
                        selfname = pybind11::reinterpret_borrow<pybind11::object>(_h); \
                    } \
                } \
            } \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
            pybind11::function override \
                = pybind11::get_override(static_cast<const BaseType *>(ref), pyfuncname); \
            if (override) { \
                auto o = override.operator()<pybind11::return_value_policy::reference>(__VA_ARGS__); \
                if (pybind11::detail::cast_is_temporary_value_reference<returnType>::value) { \
                    static pybind11::detail::override_caster_t<returnType> caster; \
                    return pybind11::detail::cast_ref<returnType>(std::move(o), caster); \
                } \
                return pybind11::detail::cast_safe<returnType>(std::move(o)); \
            } \
            pybind11::pybind11_fail( \
                "Tried to call pure virtual function \"" PYBIND11_STRINGIFY(BaseType) "::" #cfuncname "\""); \
        } while (false);

#define SELF_OVERRIDE_REF(selfname, BaseType, returnType, cfuncname, pyfuncname, ...) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            if(!selfname) { \
                auto *_tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
                if(_tinfo) { \
                    pybind11::handle _h = pybind11::detail::get_object_handle(static_cast<const BaseType *>(this), _tinfo); \
                    if(_h) { \
                        selfname = pybind11::reinterpret_borrow<pybind11::object>(_h); \
                    } \
                } \
            } \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
            pybind11::function override \
                = pybind11::get_override(static_cast<const BaseType *>(ref), pyfuncname); \
            if (override) { \
                auto o = override.operator()<pybind11::return_value_policy::reference>(__VA_ARGS__); \
                if (pybind11::detail::cast_is_temporary_value_reference<returnType>::value) { \
                    static pybind11::detail::override_caster_t<returnType> caster; \
                    return pybind11::detail::cast_ref<returnType>(std::move(o), caster); \
                } \
                return pybind11::detail::cast_safe<returnType>(std::move(o)); \
            } \
        } while (false); \
        return BaseType::cfuncname(__VA_ARGS__);

// Recover the python half of a trampoline instance into `selfname` when it
// has not been cached yet. Used by the SELF_OVERRIDE family and by trampoline
// methods that provide their own C++ defaults.
#define SELF_RECOVER(selfname, BaseType) \
        if(!selfname) { \
            auto *_tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
            if(_tinfo) { \
                pybind11::handle _h = pybind11::detail::get_object_handle(static_cast<const BaseType *>(this), _tinfo); \
                if(_h) { \
                    selfname = pybind11::reinterpret_borrow<pybind11::object>(_h); \
                } \
            } \
        }

// Dispatch Name() to a python override when one exists; otherwise report the
// python class name, falling back to the C++ base type name when the
// instance has no python half.
#define SELF_OVERRIDE_NAME_CLASSNAME_DEFAULT(selfname, BaseType) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            SELF_RECOVER(selfname, BaseType) \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
            pybind11::function override \
                = pybind11::get_override(static_cast<const BaseType *>(ref), "Name"); \
            if (override) { \
                return override().cast<std::string>(); \
            } \
            if(selfname) { \
                return selfname.attr("__class__").attr("__name__").cast<std::string>(); \
            } \
            return std::string(PYBIND11_STRINGIFY(BaseType)); \
        } while (false);

// Dispatch equal() to a python override when one exists; otherwise compare
// object identity. The framework deduplicates distributions by shared
// pointer identity, so identity is the safe default for python subclasses.
#define SELF_OVERRIDE_EQUAL_IDENTITY_DEFAULT(selfname, BaseType, CompareType, other) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            SELF_RECOVER(selfname, BaseType) \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
            pybind11::function override \
                = pybind11::get_override(static_cast<const BaseType *>(ref), "equal"); \
            if (override) { \
                return override(other).cast<bool>(); \
            } \
            return static_cast<CompareType const *>(this) == &other; \
        } while (false);

// Dispatch less() to a python override when one exists; otherwise order by
// object address so that python subclasses have a consistent total order
// within a process.
#define SELF_OVERRIDE_LESS_IDENTITY_DEFAULT(selfname, BaseType, CompareType, other) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            SELF_RECOVER(selfname, BaseType) \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
            pybind11::function override \
                = pybind11::get_override(static_cast<const BaseType *>(ref), "less"); \
            if (override) { \
                return override(other).cast<bool>(); \
            } \
            return std::less<CompareType const *>()(static_cast<CompareType const *>(this), &other); \
        } while (false);

// Dispatch clone() to a python override when one exists; otherwise fail with
// an instructive error. Nothing in the framework calls clone(), so python
// subclasses only need to define it when they use cloning themselves.
#define SELF_OVERRIDE_CLONE_REQUIRED(selfname, BaseType, ReturnElementType) \
        do { \
            pybind11::gil_scoped_acquire gil; \
            SELF_RECOVER(selfname, BaseType) \
            const BaseType * ref; \
            if(selfname) { \
                ref = selfname.cast<BaseType *>(); \
            } else { \
                ref = this; \
            } \
            pybind11::function override \
                = pybind11::get_override(static_cast<const BaseType *>(ref), "clone"); \
            if (override) { \
                return override().cast<std::shared_ptr<ReturnElementType>>(); \
            } \
            throw std::runtime_error("clone() is not implemented for this python-defined distribution; define clone(self) returning a fresh instance to enable cloning"); \
        } while (false);

#define Pybind11TrampolineCerealMethods(BaseType, TrampolineType) \
public: \
    mutable pybind11::object self; \
    pybind11::object get_representation() { \
        const BaseType * ref; \
        if(self) { \
            ref = self.cast<BaseType *>(); \
        } else { \
            ref = dynamic_cast<const BaseType *>(dynamic_cast<const TrampolineType *>(this)); \
            if (!ref) { \
                int status; \
                char * realname; \
                const std::type_info  &ti = typeid(BaseType); \
                realname = abi::__cxa_demangle(ti.name(), 0, 0, &status); \
                std::stringstream msg; \
                msg << "Cannot cast this to " << realname; \
                free(realname); \
                throw std::runtime_error(msg.str()); \
            } \
        } \
        auto *tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
        pybind11::function override_func = \
            tinfo ? pybind11::detail::get_type_override(static_cast<const BaseType *>(ref), tinfo, "get_representation") : pybind11::function(); \
        if (override_func) { \
            pybind11::object o = override_func(); \
            if(not pybind11::isinstance<pybind11::dict>(o)) { \
                throw std::runtime_error("get_representation must return a dict"); \
            } \
            return o; \
        } \
        pybind11::object _self; \
        if(this->self) { \
            _self = pybind11::reinterpret_borrow<pybind11::object>(this->self); \
        } else { \
            auto *tinfo = pybind11::detail::get_type_info(typeid(BaseType)); \
            pybind11::handle self_handle = get_object_handle(dynamic_cast<const BaseType *>(this), tinfo); \
            _self = pybind11::reinterpret_borrow<pybind11::object>(self_handle); \
        } \
        pybind11::dict d; \
        if (pybind11::hasattr(_self, "__dict__")) { \
            d = _self.attr("__dict__"); \
        } \
        return d; \
    } \
    static pybind11::tuple pickle_save(BaseType & cpp_obj) { \
        TrampolineType * trampoline = dynamic_cast<TrampolineType*>(&cpp_obj); \
        if(!trampoline) { \
            throw std::runtime_error("Cannot pickle a C++-defined instance through the python trampoline; pickling is only supported for python-defined subclasses"); \
        } \
        pybind11::object x = trampoline->get_representation(); \
        return pybind11::make_tuple(x); \
    } \
    static std::pair<std::unique_ptr<BaseType>, pybind11::dict> pickle_load(const pybind11::tuple &t) { \
        if (t.size() != 1) { \
            throw std::runtime_error("Invalid state!"); \
        } \
        auto cpp_state = std::unique_ptr<BaseType>(new TrampolineType); \
        auto py_state = t[0].cast<pybind11::dict>(); \
        return std::make_pair(std::move(cpp_state), py_state); \
    } \
 \
    template<typename Archive> \
    void save(Archive & archive, std::uint32_t const version) const { \
        if(version == 0) { \
 \
            pybind11::object obj; \
            if(this->self) { \
                obj = this->self; \
            } else { \
                auto *tinfo = pybind11::detail::get_type_info(typeid(TrampolineType)); \
                pybind11::handle self_handle = get_object_handle(dynamic_cast<const TrampolineType *>(this), tinfo); \
                obj = pybind11::reinterpret_borrow<pybind11::object>(self_handle); \
            } \
 \
			pybind11::module pkl = pybind11::module::import("pickle"); \
			pybind11::bytes bytes = pkl.attr("dumps")(obj); \
			std::string str_repr = (std::string)(bytes.attr("hex")().cast<std::string>()); \
 \
			archive(::cereal::make_nvp("PythonPickleBytesRepresentation", str_repr)); \
 \
            archive(cereal::virtual_base_class<BaseType>(dynamic_cast<const BaseType*>(this))); \
 \
        } else { \
            throw std::runtime_error("BaseType only supports version <= 0!"); \
        } \
    } \
 \
    template<typename Archive> \
    void load(Archive & archive, std::uint32_t version) { \
        if(version == 0) { \
            std::string str_repr; \
			archive(::cereal::make_nvp("PythonPickleBytesRepresentation", str_repr)); \
 \
            pybind11::module pkl = pybind11::module::import("pickle"); \
 \
            pybind11::object bytes_module = pybind11::module::import("builtins").attr("bytes"); \
            pybind11::object bytes = bytes_module.attr("fromhex")(str_repr); \
 \
            pkl.attr("loads")(bytes); \
            this->self = pkl.attr("loads")(bytes); \
 \
            archive(cereal::virtual_base_class<BaseType>(dynamic_cast<const BaseType*>(this))); \
 \
        } else { \
            throw std::runtime_error("BaseType only supports version <= 0!"); \
        } \
    } \

#define RegisterTrampolinePickleMethods(object, TrampolineType) object.def(pybind11::pickle(&TrampolineType::pickle_save, &TrampolineType::pickle_load));
#define TrampolinePickleMethods(TrampolineType) .def(pybind11::pickle(&TrampolineType::pickle_save, &TrampolineType::pickle_load))

#define RegisterTrampolineCerealMethods(BaseType, TrampolineType, Pybind11TrampolineType) \
    CEREAL_CLASS_VERSION(TrampolineType,0); \
    CEREAL_REGISTER_TYPE(TrampolineType); \
    CEREAL_REGISTER_POLYMORPHIC_RELATION(BaseType, TrampolineType); \
    CEREAL_CLASS_VERSION(Pybind11TrampolineType, 0); \
    CEREAL_REGISTER_TYPE(Pybind11TrampolineType); \
    CEREAL_REGISTER_POLYMORPHIC_RELATION(Pybind11TrampolineType, TrampolineType); \

#endif // SIREN_Pybind11Trampoline_H
