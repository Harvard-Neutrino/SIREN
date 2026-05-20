// Macros for generating the boilerplate operator methods (swap, operator=,
// equal, less) that every Geometry subclass implements identically.
//
// Usage in a .cxx file:
//
//   SIREN_GEOMETRY_SWAP(Box, x_, y_, z_)
//   SIREN_GEOMETRY_ASSIGN(Box)
//   SIREN_GEOMETRY_EQUAL(Box, x_, y_, z_)
//   SIREN_GEOMETRY_LESS(Box, x_, y_, z_)
//
// These expand to the standard dynamic_cast + per-field pattern shared by
// every Geometry subclass.

#pragma once
#ifndef SIREN_GeometryMacros_H
#define SIREN_GeometryMacros_H

#include <tuple>     // std::tie used in EQUAL/LESS expansions
#include <utility>   // std::swap used in SWAP expansion

// --- Internal helpers (variadic field expansion) ---

// Argument counter (up to 12 fields)
#define SIREN_ARG_N_(_1,_2,_3,_4,_5,_6,_7,_8,_9,_10,_11,_12,N,...) N
#define SIREN_NARGS_(...) SIREN_ARG_N_(__VA_ARGS__,12,11,10,9,8,7,6,5,4,3,2,1)
#define SIREN_CAT_(a, b) a##b
#define SIREN_CAT(a, b) SIREN_CAT_(a, b)

// std::swap(field, other->field) for each field
#define SIREN_SWAP1_(o, a) std::swap(a, o->a);
#define SIREN_SWAP_1(o, a) SIREN_SWAP1_(o,a)
#define SIREN_SWAP_2(o, a,b) SIREN_SWAP1_(o,a) SIREN_SWAP1_(o,b)
#define SIREN_SWAP_3(o, a,b,c) SIREN_SWAP_2(o,a,b) SIREN_SWAP1_(o,c)
#define SIREN_SWAP_4(o, a,b,c,d) SIREN_SWAP_3(o,a,b,c) SIREN_SWAP1_(o,d)
#define SIREN_SWAP_5(o, a,b,c,d,e) SIREN_SWAP_4(o,a,b,c,d) SIREN_SWAP1_(o,e)
#define SIREN_SWAP_6(o, a,b,c,d,e,f) SIREN_SWAP_5(o,a,b,c,d,e) SIREN_SWAP1_(o,f)
#define SIREN_SWAP_7(o, a,b,c,d,e,f,g) SIREN_SWAP_6(o,a,b,c,d,e,f) SIREN_SWAP1_(o,g)
#define SIREN_SWAP_8(o, a,b,c,d,e,f,g,h) SIREN_SWAP_7(o,a,b,c,d,e,f,g) SIREN_SWAP1_(o,h)
#define SIREN_SWAP_9(o, a,b,c,d,e,f,g,h,i) SIREN_SWAP_8(o,a,b,c,d,e,f,g,h) SIREN_SWAP1_(o,i)
#define SIREN_SWAP_10(o, a,b,c,d,e,f,g,h,i,j) SIREN_SWAP_9(o,a,b,c,d,e,f,g,h,i) SIREN_SWAP1_(o,j)
#define SIREN_SWAP_11(o, a,b,c,d,e,f,g,h,i,j,k) SIREN_SWAP_10(o,a,b,c,d,e,f,g,h,i,j) SIREN_SWAP1_(o,k)
#define SIREN_SWAP_12(o, a,b,c,d,e,f,g,h,i,j,k,l) SIREN_SWAP_11(o,a,b,c,d,e,f,g,h,i,j,k) SIREN_SWAP1_(o,l)
#define SIREN_SWAP_ALL(o, ...) SIREN_CAT(SIREN_SWAP_, SIREN_NARGS_(__VA_ARGS__))(o, __VA_ARGS__)

// Build std::tie(a, b, c) from field list
#define SIREN_TIE_1(a) std::tie(a)
#define SIREN_TIE_2(a,b) std::tie(a,b)
#define SIREN_TIE_3(a,b,c) std::tie(a,b,c)
#define SIREN_TIE_4(a,b,c,d) std::tie(a,b,c,d)
#define SIREN_TIE_5(a,b,c,d,e) std::tie(a,b,c,d,e)
#define SIREN_TIE_6(a,b,c,d,e,f) std::tie(a,b,c,d,e,f)
#define SIREN_TIE_7(a,b,c,d,e,f,g) std::tie(a,b,c,d,e,f,g)
#define SIREN_TIE_8(a,b,c,d,e,f,g,h) std::tie(a,b,c,d,e,f,g,h)
#define SIREN_TIE_9(a,b,c,d,e,f,g,h,i) std::tie(a,b,c,d,e,f,g,h,i)
#define SIREN_TIE_10(a,b,c,d,e,f,g,h,i,j) std::tie(a,b,c,d,e,f,g,h,i,j)
#define SIREN_TIE_11(a,b,c,d,e,f,g,h,i,j,k) std::tie(a,b,c,d,e,f,g,h,i,j,k)
#define SIREN_TIE_12(a,b,c,d,e,f,g,h,i,j,k,l) std::tie(a,b,c,d,e,f,g,h,i,j,k,l)
#define SIREN_TIE(...) SIREN_CAT(SIREN_TIE_, SIREN_NARGS_(__VA_ARGS__))(__VA_ARGS__)

// Build std::tie(o->a, o->b, o->c) from field list
#define SIREN_OTIE_1(o, a) std::tie(o->a)
#define SIREN_OTIE_2(o, a,b) std::tie(o->a,o->b)
#define SIREN_OTIE_3(o, a,b,c) std::tie(o->a,o->b,o->c)
#define SIREN_OTIE_4(o, a,b,c,d) std::tie(o->a,o->b,o->c,o->d)
#define SIREN_OTIE_5(o, a,b,c,d,e) std::tie(o->a,o->b,o->c,o->d,o->e)
#define SIREN_OTIE_6(o, a,b,c,d,e,f) std::tie(o->a,o->b,o->c,o->d,o->e,o->f)
#define SIREN_OTIE_7(o, a,b,c,d,e,f,g) std::tie(o->a,o->b,o->c,o->d,o->e,o->f,o->g)
#define SIREN_OTIE_8(o, a,b,c,d,e,f,g,h) std::tie(o->a,o->b,o->c,o->d,o->e,o->f,o->g,o->h)
#define SIREN_OTIE_9(o, a,b,c,d,e,f,g,h,i) std::tie(o->a,o->b,o->c,o->d,o->e,o->f,o->g,o->h,o->i)
#define SIREN_OTIE_10(o, a,b,c,d,e,f,g,h,i,j) std::tie(o->a,o->b,o->c,o->d,o->e,o->f,o->g,o->h,o->i,o->j)
#define SIREN_OTIE_11(o, a,b,c,d,e,f,g,h,i,j,k) std::tie(o->a,o->b,o->c,o->d,o->e,o->f,o->g,o->h,o->i,o->j,o->k)
#define SIREN_OTIE_12(o, a,b,c,d,e,f,g,h,i,j,k,l) std::tie(o->a,o->b,o->c,o->d,o->e,o->f,o->g,o->h,o->i,o->j,o->k,o->l)
#define SIREN_OTIE(o, ...) SIREN_CAT(SIREN_OTIE_, SIREN_NARGS_(__VA_ARGS__))(o, __VA_ARGS__)

// --- Public macros ---

#define SIREN_GEOMETRY_SWAP(ClassName, ...) \
void ClassName::swap(Geometry& geometry) { \
    ClassName* o_ = dynamic_cast<ClassName*>(&geometry); \
    if(!o_) return; \
    Geometry::swap(*o_); \
    SIREN_SWAP_ALL(o_, __VA_ARGS__) \
}

#define SIREN_GEOMETRY_ASSIGN(ClassName) \
ClassName& ClassName::operator=(const Geometry& geometry) { \
    if(this != &geometry) { \
        const ClassName* o_ = dynamic_cast<const ClassName*>(&geometry); \
        if(!o_) return *this; \
        ClassName tmp(*o_); \
        swap(tmp); \
    } \
    return *this; \
}

#define SIREN_GEOMETRY_EQUAL(ClassName, ...) \
bool ClassName::equal(const Geometry& geometry) const { \
    const ClassName* o_ = dynamic_cast<const ClassName*>(&geometry); \
    if(!o_) return false; \
    return SIREN_TIE(__VA_ARGS__) == SIREN_OTIE(o_, __VA_ARGS__); \
}

#define SIREN_GEOMETRY_LESS(ClassName, ...) \
bool ClassName::less(const Geometry& geometry) const { \
    const ClassName* o_ = dynamic_cast<const ClassName*>(&geometry); \
    if(!o_) return false; \
    return SIREN_TIE(__VA_ARGS__) < SIREN_OTIE(o_, __VA_ARGS__); \
}

#endif // SIREN_GeometryMacros_H
