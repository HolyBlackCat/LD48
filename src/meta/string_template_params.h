#pragma once

#include <algorithm>
#include <cstddef>

#include "meta/basic.h"
#include "program/errors.h"

// This file offers compile-time strings that can be used as template parameters.
// Example:
//     template <char ...Name> void foo(Meta::ConstStringParam<Name...>) {}
//     foo("123"_c);


namespace Meta
{
    template <char ...C>
    struct ConstStringParam {};

    #ifdef __clang__
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wgnu-string-literal-operator-template"
    #endif

    template <typename T, T ...C>
    [[nodiscard]] ConstStringParam<C...> operator""_c()
    {
        return {};
    }

    #ifdef __clang__
    #pragma clang diagnostic pop
    #endif


    // A better, conformant implementation. Requires Clang 12+.
    #if 0
    // This file offers compile-time strings that can be used as template parameters.
    // Example 1:
    //     template <Meta::ConstString Name> void foo() {std::cout << Name.str << '\n';}
    //     foo<"123">();
    // Example 2:
    //     template <Meta::ConstString Name> void foo(Meta::ConstStringParam<Name>) {std::cout << Name.str << '\n';}
    //     foo("123"_c);

    // A string that can be used as a template parameter.
    template <std::size_t N>
    struct ConstString
    {
        char str[N]{};

        static constexpr std::size_t size = N - 1;

        [[nodiscard]] std::string_view view() const
        {
            return {str, str + size};
        }

        constexpr ConstString() {}
        constexpr ConstString(const char (&new_str)[N])
        {
            ASSERT(new_str[N-1] == '\0');
            std::copy_n(new_str, N, str);
        }
    };

    // A tag structure returned by `operator""_c` below.
    template <Meta::ConstString S>
    struct ConstStringParam {};

    // Returns a string encoded into a template parameter of a tag structure `ConstStringParam`.
    template <Meta::ConstString S>
    [[nodiscard]] constexpr ConstStringParam<S> operator""_c()
    {
        return {};
    }
    #endif
}

using Meta::operator""_c;
