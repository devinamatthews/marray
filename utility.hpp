#ifndef _MARRAY_UTILITY_HPP_
#define _MARRAY_UTILITY_HPP_

#ifndef MARRAY_TEST
#define MARRAY_TEST(...)
#endif

#include <type_traits>
#include <array>
#include <vector>
#include <utility>
#include <iterator>
#include <cassert>
#include <algorithm>

namespace MArray
{
    /*
     * Create a vector from the specified elements, where the type of the vector
     * is taken from the first element.
     */
    template <typename T, typename... Args>
    std::vector<typename std::decay<T>::type>
    make_vector(T&& t, Args&&... args)
    {
        return std::vector<typename std::decay<T>::type>
            {{std::forward<T>(t), std::forward<Args>(args)...}};
    }

    MARRAY_TEST
    (
        std::vector<int> x1 = {1, 2, 3, 4};
        auto x2 = make_vector(1, 2, 3, 4);
        assert(x1 == x2);
    )

    /*
     * Create an array from the specified elements, where the type of the array
     * is taken from the first element.
     */
    template <typename T, typename... Args>
    std::array<typename std::decay<T>::type, sizeof...(Args)+1>
    make_array(T&& t, Args&&... args)
    {
        return std::array<typename std::decay<T>::type, sizeof...(Args)+1>
            {{std::forward<T>(t), std::forward<Args>(args)...}};
    }

    MARRAY_TEST
    (
        std::array<int, 4> x1 = {1, 2, 3, 4};
        auto x2 = make_array(1, 2, 3, 4);
        assert(x1 == x2);
    )

    template <typename T>
    class range_t
    {
        static_assert(std::is_integral<T>::value, "The type must be integral.");

        protected:
            T from;
            T to;

            typedef T value_type;
            typedef T size_type;

        public:
            class iterator : std::iterator<std::random_access_iterator_tag,T>
            {
                protected:
                    T val;

                public:
                    using typename std::iterator<std::random_access_iterator_tag,T>::iterator_category;
                    using typename std::iterator<std::random_access_iterator_tag,T>::value_type;
                    using typename std::iterator<std::random_access_iterator_tag,T>::difference_type;
                    using typename std::iterator<std::random_access_iterator_tag,T>::pointer;
                    using typename std::iterator<std::random_access_iterator_tag,T>::reference;

                    constexpr iterator() : val(0) {}

                    constexpr iterator(T val) : val(val) {}

                    bool operator==(const iterator& other)
                    {
                        return val == other.val;
                    }

                    bool operator!=(const iterator& other)
                    {
                        return val != other.val;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        assert(r.begin() == r.end()-4);
                        assert(r.begin() != r.end());
                    )

                    value_type operator*() const
                    {
                        return val;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        assert(*r.begin() == 1);
                    )

                    iterator& operator++()
                    {
                        ++val;
                        return *this;
                    }

                    iterator operator++(int x)
                    {
                        iterator old(*this);
                        ++val;
                        return old;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        auto i = r.begin();
                        ++i;
                        auto j = i++;
                        assert(*i == 3 && *j == 2);
                    )

                    iterator& operator--()
                    {
                        --val;
                        return *this;
                    }

                    iterator operator--(int x)
                    {
                        iterator old(*this);
                        --val;
                        return old;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        auto i = r.begin();
                        --i;
                        auto j = i--;
                        assert(*i == -1 && *j == 0);
                    )

                    iterator& operator+=(difference_type n)
                    {
                        val += n;
                        return *this;
                    }

                    iterator operator+(difference_type n)
                    {
                        return iterator(val+n);
                    }

                    friend iterator operator+(difference_type n, const iterator& i)
                    {
                        return iterator(i.val+n);
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        auto i = r.begin();
                        assert(*(i+2) == 3);
                        assert(*(5+i) == 6);
                        i += 1;
                        assert(*i == 2);
                    )

                    iterator& operator-=(difference_type n)
                    {
                        val -= n;
                        return *this;
                    }

                    iterator operator-(difference_type n)
                    {
                        return iterator(val-n);
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        auto i = r.begin();
                        assert(*(i-2) == -1);
                        i -= 1;
                        assert(*i == 0);
                    )

                    difference_type operator-(const iterator& other)
                    {
                        return val-other.val;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        assert(r.end()-r.begin() == 4);
                    )

                    bool operator<(const iterator& other)
                    {
                        return val < other.val;
                    }

                    bool operator<=(const iterator& other)
                    {
                        return val <= other.val;
                    }

                    bool operator>(const iterator& other)
                    {
                        return val > other.val;
                    }

                    bool operator>=(const iterator& other)
                    {
                        return val >= other.val;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        assert(r.begin() < r.end());
                        assert(r.end() > r.begin());
                        assert(r.begin() <= r.end());
                        assert(r.end() >= r.begin());
                        assert(r.begin() <= r.begin());
                        assert(r.begin() >= r.begin());
                    )

                    value_type operator[](difference_type n) const
                    {
                        return val+n;
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        assert(r.begin()[7] == 8);
                    )

                    friend void swap(iterator& a, iterator& b)
                    {
                        using std::swap;
                        swap(a.val, b.val);
                    }

                    MARRAY_TEST
                    (
                        range_t<int> r(1, 5);
                        auto b = r.begin();
                        auto e = r.end();
                        swap(b, e);
                        assert(b == r.end());
                        assert(e == r.begin());
                    )
            };

            constexpr range_t(T from, T to) : from(from), to(to) {}

            size_type size() const
            {
                return to-from;
            }

            MARRAY_TEST
            (
                range_t<int> r(1, 5);
                assert(r.size() == 4);
            )

            iterator begin() const
            {
                return iterator(from);
            }

            iterator end() const
            {
                return iterator(to);
            }

            MARRAY_TEST
            (
                range_t<int> r(1, 5);
                assert(*r.begin() == 1);
                assert(*r.end() == 5);
            )

            value_type front() const
            {
                return from;
            }

            value_type back() const
            {
                return to;
            }

            MARRAY_TEST
            (
                range_t<int> r(1, 5);
                assert(r.front() == 1);
                assert(r.back() == 5);
            )

            value_type operator[](size_type n) const
            {
                return from+n;
            }

            MARRAY_TEST
            (
                range_t<int> r(1, 5);
                assert(r[7] == 8);
            )

            operator std::vector<T>() const
            {
                return std::vector<T>(begin(), end());
            }

            MARRAY_TEST
            (
                range_t<int> r(1, 5);
                std::vector<int> x1 = {1, 2, 3, 4};
                std::vector<int> x2 = r;
                assert(x1 == x2);
            )
    };

    template <typename T>
    range_t<T> range(T to)
    {
        return range_t<T>(T(), to);
    }

    MARRAY_TEST
    (
        std::vector<int> x1 = {0, 1, 2, 3};
        std::vector<int> x2 = range(4);
        assert(x1 == x2);
    )

    template <typename T>
    range_t<T> range(T from, T to)
    {
        return range_t<T>(from, to);
    }

    MARRAY_TEST
    (
        std::vector<int> x1 = {1, 2, 3, 4};
        std::vector<int> x2 = range(1, 5);
        assert(x1 == x2);
    )
}

#endif
