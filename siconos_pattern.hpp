#ifndef SICONOS_PATTERN_HPP
#define SICONOS_PATTERN_HPP

#include "siconos_util.hpp"

#include <initializer_list>
namespace siconos
{
  namespace some
  {
    struct tag { using tag_t = void; };

    struct t {};
    struct scalar : t {};
    struct indice : t {};

    template<typename T>
    struct vdescriptor : t {};

    template<std::size_t N, std::size_t M>
    struct matrix : t {};

    template<std::size_t N>
    struct vector : t {};

    template<typename T>
    concept type = std::derived_from<T, some::t>;
  }

  template<string_literal Text>
  struct text
  {
    static constexpr auto data = Text;
  };

  template<string_literal Symbol>
  struct symbol : text<Symbol>{};;

  template<string_literal Descr>
  struct description : text<Descr>{};

  template<typename T>
  struct tag
  {
    using type = T;
  };

  // concept item
  template<typename T>
  struct use
  {
    using use_t = void;
    using type = T;
  };

  // concept some:...
  template<typename T>
  struct structure
  {
    using structure_t = void;
    using type = T;
  };

  template<
    typename Tag,
    typename MathSymbol,
    typename Description,
    typename Structure>
  struct attribute
  {
    using tag = typename Tag::type;
    using symbol = MathSymbol;
    using description = Description;
    using structure = Structure;

    static constexpr symbol symbol_v = symbol{};
  };

  template<typename T>
  struct attribute_p
  {
    static constexpr auto value = concepts::attribute<T>;
  };

  template<typename T>
  struct tag_p
  {
    static constexpr auto value = concepts::tag<T>;
  };

  template<typename T>
  struct use_p
  {
    static constexpr auto value = concepts::use<T>;
  };

  template<typename T>
  struct keep_p
  {
    static constexpr auto value = concepts::keep<T>;
  };

  template<typename Tag, std::size_t N>
  struct keep
  {
    using tag = Tag;
    static constexpr std::size_t size = N;
  };


  template<typename Descr, typename ...Args>
  struct item
  {
    static constexpr auto args = std::make_tuple(Args{}...);

    using description = Descr;
    using tags = decltype(filter<tag_p>(args));
    using attributes = decltype(filter<attribute_p>(args));
    using uses = decltype(filter<use_p>(args));
    using keeps = decltype(filter<keep_p>(args));

    template<string_literal symb>
    static constexpr auto get = proj(attributes{})(get_m<symb>);

  };

  template<typename Descr, typename ...Args>
  struct vertex_item : item<Descr, Args...>
  {
    using vertex_item_t = void;

    static constexpr auto args = std::make_tuple(Args{}...);

    using description = Descr;
    using tags = decltype(filter<tag_p>(args));
    using attributes = decltype(filter<attribute_p>(args));
    using uses = decltype(filter<use_p>(args));
    using keeps = decltype(filter<keep_p>(args));

    template<string_literal symb>
    static constexpr auto get = proj(attributes{})(get_m<symb>);
  };

  template<typename... Ts>
  struct handle
  {
    using type = std::tuple<Ts...>;
    type value = {};
    handle(Ts... ts) : value{std::forward<Ts>(ts)...} {};
  };

  namespace concepts
  {
    // T is a tag
    template<typename T, typename Data>
    concept vertex_item_t = requires (T t) { { static_cast<typename Data::vertex_items>(t) }; };
  }
}

#endif
