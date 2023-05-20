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

#ifndef FS_FIND_IF_H
#define FS_FIND_IF_H

#include <type_traits>
#include <utility>

#include <boost/mpl/begin_end.hpp>
#include <boost/mpl/deref.hpp>
#include <boost/mpl/next.hpp>

namespace flexiblesusy
{
namespace meta
{
namespace detail
{
template <class Iterator, class End, template <typename> class Predicate,
          template <typename> class F, class State>
bool find_if_impl(State&& state, std::true_type)
{
   return false;
}

template <class Iterator, class End, template <typename> class Predicate,
          template <typename> class F, class State>
bool find_if_impl(State&& state, std::false_type)
{
   if (Predicate<typename boost::mpl::deref<Iterator>::type>{}(state)) {
      F<typename boost::mpl::deref<Iterator>::type>{}(state);
      return true;
   }

   using NextIterator = typename boost::mpl::next<Iterator>::type;

   return find_if_impl<NextIterator, End, Predicate, F>(
      std::forward<State>(state),
      typename std::is_same<NextIterator, End>::type{});
}
} // namespace detail

template <class Sequence, template <typename> class Predicate,
          template <typename> class F, class State>
bool find_if(State&& state)
{
   return detail::find_if_impl<typename boost::mpl::begin<Sequence>::type,
                               typename boost::mpl::end<Sequence>::type,
                               Predicate, F>(
      std::forward<State>(state),
      typename std::is_same<typename boost::mpl::begin<Sequence>::type,
                            typename boost::mpl::end<Sequence>::type>::type{});
}
} // namespace meta
} // namespace flexiblesusy

#endif
