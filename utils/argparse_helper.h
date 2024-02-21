#pragma once

// HINT: It's recommended to include this "heavy" header
// (which involves huge amount of metaprogramming operations)
// into some specific source file instead of using it as header-only
// to reduce compilation overhead.

#include "rfl/named_tuple_t.hpp"
#include "rfl/to_view.hpp"
#include "utils/demangle.h"
#include "utils/reflect.h"
#include "utils/utils.h"
#include <argparse/argparse.hpp>
#include <concepts>
#include <fmt/format.h>
#include <magic_enum.hpp>
#include <nlohmann/json.hpp>
#include <nwgraph/util/demangle.hpp>
#include <rfl/Result.hpp>
#include <rfl/always_false.hpp>

using argparse::Argument;
using argparse::ArgumentParser;

namespace details {
template <class T>
auto add_argument(ArgumentParser& parser, std::string_view arg_name) -> Argument& {
  if constexpr (std::is_same_v<T, std::string>) {
    // String: std::string only
    return parser.add_argument(arg_name);
  } else if constexpr (std::is_same_v<T, bool>) {
    // Bool: As flag argument
    return parser.add_argument(arg_name).implicit_value(true);
  } else if constexpr (std::is_enum_v<T>) {
    // Enumerators: Flag enumerators are not supported
    return parser.add_argument(arg_name);
  } else if constexpr (std::signed_integral<T>) {
    // Signed integer
    return parser.add_argument(arg_name).template scan<'i', T>();
  } else if constexpr (std::unsigned_integral<T>) {
    // Unsigned integer
    return parser.add_argument(arg_name).template scan<'u', T>();
  } else if constexpr (std::is_floating_point_v<T>) {
    // Floating-point
    return parser.add_argument(arg_name).template scan<'g', T>();
  } else if constexpr (is_std_vector_v<T>) {
    // List: std::vector only.
    using ValueType = ranges::range_value_t<T>;
    static_assert(!ranges::range<ValueType>, "Nested range is not supported.");
    using ValueReflectType = ReflectionType<ValueType>;
    return add_argument<ValueReflectType>(parser, arg_name).nargs(argparse::nargs_pattern::any);
  } else {
    static_assert(rfl::always_false_v<T>, "Unsupported value type.");
  }
}

template <class Field, class TView>
auto init_argument_parser_visitor(ArgumentParser& parser, const TView& default_view) -> void {
  constexpr auto name_literal = Field::name_;

  auto arg_name = "--" + name_literal.str();
  ranges::replace(arg_name, '_', '-'); // Formats '--some_member' to '--some-member'

  using Value = FieldReflectionType<Field>;
  add_argument<Value>(parser, arg_name);
}

template <class Field, class TView>
void read_argument_parser_visitor(ArgumentParser& parser, const TView& dest_view) {
  constexpr auto name_literal = Field::name_;
  auto& dest = *dest_view.template get<name_literal>();

  auto arg_name = "--" + name_literal.str();
  ranges::replace(arg_name, '_', '-'); // Formats '--some_member' to '--some-member'

  using Value = FieldReflectionType<Field>;
  if constexpr (std::is_enum_v<Value>) {
    auto opt = parser.present(arg_name);
    if (opt) {
      auto enum_value = magic_enum::enum_cast<Value>(*opt);
      if (!enum_value) {
        constexpr auto msg_pattern = "Invalid enumerator string '{}' for type {}.";
        throw std::invalid_argument{fmt::format(msg_pattern, *opt, demangle_type_name<Value>())};
      }
      dest = *enum_value;
    } // Otherwise, the field is simply ignored
  } else {
    auto opt = parser.template present<Value>(arg_name);
    if (opt) {
      dest = *opt; // Performs check during operator= for rfl::Validator types
    }
  }
}

template <class FieldTuple, class TView, size_t... Indices>
auto init_argument_parser_impl(ArgumentParser& parser, const TView& default_view, std::index_sequence<Indices...>)
    -> void {
  (init_argument_parser_visitor<std::tuple_element_t<Indices, FieldTuple>>(parser, default_view), ...);
}

template <class FieldTuple, class TView, size_t... Indices>
auto read_argument_parser_impl(ArgumentParser& parser, const TView& dest_view, std::index_sequence<Indices...>)
    -> void {
  (read_argument_parser_visitor<std::tuple_element_t<Indices, FieldTuple>>(parser, dest_view), ...);
}

template <class FieldTuple, std::default_initializable T>
auto init_argument_parser(ArgumentParser& parser) -> void {
  auto default_value = T{};
  init_argument_parser_impl<FieldTuple>(parser, rfl::to_view(default_value),
                                        std::make_index_sequence<std::tuple_size_v<FieldTuple>>{});
}

template <class FieldTuple, std::default_initializable T>
auto read_argument_parser(ArgumentParser& parser, T& dest) -> void {
  auto dest_view = rfl::to_view(dest);
  read_argument_parser_impl<FieldTuple>(parser, dest_view, std::make_index_sequence<std::tuple_size_v<FieldTuple>>{});
}
} // namespace details

template <std::default_initializable T, invocable_or_nullptr<ArgumentParser&> AppendExtraArgsFn = std::nullptr_t>
auto init_argument_parser_generic(ArgumentParser& parser, AppendExtraArgsFn&& append_extra_args_fn = nullptr) noexcept
    -> void {
  using FieldTuple = typename rfl::named_tuple_t<T>::Fields;
  // Step 1: Adds all the members of type T to the parser
  details::init_argument_parser<FieldTuple, T>(parser);
  // Step 2: (Optional) Adds customized additional arguments
  if constexpr (!std::is_null_pointer_v<AppendExtraArgsFn>) {
    std::invoke(append_extra_args_fn, parser);
  }
}

template <std::default_initializable T, invocable_or_nullptr<ArgumentParser&, T&> GetExtraArgsFn = std::nullptr_t>
auto read_argument_parser_generic(ArgumentParser& parser, T& dest,
                                  GetExtraArgsFn&& get_extra_args_fn = nullptr) noexcept -> void {
  using FieldTuple = typename rfl::named_tuple_t<T>::Fields;
  // Step 1: (Optional) Handles customized arguments first
  if constexpr (!std::is_null_pointer_v<GetExtraArgsFn>) {
    std::invoke(get_extra_args_fn, parser, dest);
  }
  // Step 2: Gets all the members of T
  details::read_argument_parser<FieldTuple>(parser, dest);
}

template <std::default_initializable T, invocable_or_nullptr<ArgumentParser&> AppendExtraArgsFn = std::nullptr_t,
          invocable_or_nullptr<ArgumentParser&, T&> GetExtraArgsFn = std::nullptr_t>
auto parse_from_args_generic(int argc, char** argv, AppendExtraArgsFn&& append_extra_args_fn = nullptr,
                             GetExtraArgsFn&& get_extra_args_fn = nullptr) noexcept -> rfl::Result<T> try {
  using FieldTuple = typename rfl::named_tuple_t<T>::Fields;
  // Step 1: Prepares the argument parser (see above)
  auto parser = ArgumentParser();
  init_argument_parser_generic<T>(parser, std::forward<AppendExtraArgsFn>(append_extra_args_fn));
  // Step 2: Performs parsing
  parser.parse_args(argc, argv);
  // Step 3: Writes each member from the parsing result to the destination value (see above)
  auto res = T{};
  read_argument_parser_generic(parser, res, std::forward<GetExtraArgsFn>(get_extra_args_fn));
  return res;
} catch (std::exception& e) {
  constexpr auto msg_pattern = "Exception of type '{}' caught during argument parsing: `{}'";
  return rfl::Error{fmt::format(msg_pattern, demangle_type_name(e), e.what())};
} catch (...) {
  return rfl::Error{"Unknown exception caught during argument parsing."};
}
