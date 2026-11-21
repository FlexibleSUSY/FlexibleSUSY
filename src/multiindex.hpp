// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef FS_MULTIINDEX_H
#define FS_MULTIINDEX_H

#include <array>
#include <cstddef>
#include <iterator>
#include <utility>

#include <boost/mpl/at.hpp>
#include <boost/mpl/size.hpp>

namespace flexiblesusy
{
namespace detail
{
/// converts a boost::mpl sequence of integral constants to a std::array
template <class Sequence, std::size_t... I>
constexpr std::array<int, sizeof...(I)>
mpl_to_array(std::index_sequence<I...> /* indices */)
{
   return {static_cast<int>(boost::mpl::at_c<Sequence, I>::type::value)...};
}

template <class Sequence>
constexpr auto mpl_to_array()
{
   return mpl_to_array<Sequence>(
      std::make_index_sequence<boost::mpl::size<Sequence>::value>{});
}
} // namespace detail

/**
 * @class multiindex
 * @brief forward iterator over all index combinations in [Begin, End)
 *
 * Begin and End are boost::mpl sequences of integral constants with
 * the lower (inclusive) and upper (exclusive) bound of each index.
 * The first index runs fastest.  Without indices there is exactly one
 * (empty) combination; if any index has an empty range, there is none.
 */
template <class Begin, class End>
class multiindex
{
public:
   static constexpr std::size_t size = boost::mpl::size<Begin>::value;
   using data_type = std::array<int, size>;

   using difference_type = std::ptrdiff_t;
   using value_type = data_type;
   using pointer = const value_type*;
   using reference = const value_type&;
   using iterator_category = std::forward_iterator_tag;

   static_assert(boost::mpl::size<End>::value == size,
                 "multiindex: Begin and End have different lengths");

   multiindex() = default;

   static multiindex begin() { return multiindex{lower, is_empty()}; }
   static multiindex end() { return multiindex{lower, true}; }

   multiindex& operator++()
   {
      at_end = !increment(std::make_index_sequence<size>{});
      return *this;
   }

   multiindex operator++(int)
   {
      multiindex copy(*this);
      operator++();
      return copy;
   }

   reference operator*() const { return data; }
   pointer operator->() const { return &data; }

   bool operator==(const multiindex& other) const
   {
      return at_end == other.at_end && (at_end || data == other.data);
   }

private:
   static constexpr data_type lower = detail::mpl_to_array<Begin>();
   static constexpr data_type upper = detail::mpl_to_array<End>();

   data_type data{};
   bool at_end{true};

   multiindex(const data_type& indices, bool is_end)
      : data(indices), at_end(is_end)
   {
   }

   /// increments the first index that is not at its upper bound and
   /// resets the ones before it; returns false if there is none
   template <std::size_t... I>
   bool increment(std::index_sequence<I...> /* positions */)
   {
      return ((++data[I] < upper[I] ? true : (data[I] = lower[I], false)) ||
              ...);
   }

   static constexpr bool is_empty()
   {
      for (std::size_t i = 0; i < size; ++i) {
         if (lower[i] >= upper[i]) {
            return true;
         }
      }
      return false;
   }
};

/// range of all index combinations in [Begin, End)
template <class Begin, class End>
struct multiindex_range {
   using iterator = multiindex<Begin, End>;

   iterator begin() const { return iterator::begin(); }
   iterator end() const { return iterator::end(); }
};

namespace meta
{
template <class ObjectWithIndexBounds>
struct index_bounds {
   using type = typename ObjectWithIndexBounds::index_bounds;
};
} // namespace meta

/**
 * range of all index combinations of an object with index bounds,
 * used as
 * @code{.cpp}
 * for (const auto& indices : index_range<Vertex>()) { ... }
 * @endcode
 *
 * Objects with the same bounds give the same type.
 */
template <class ObjectWithIndexBounds>
using index_range = multiindex_range<
   typename meta::index_bounds<ObjectWithIndexBounds>::type::first,
   typename meta::index_bounds<ObjectWithIndexBounds>::type::second>;

} // namespace flexiblesusy

#endif
