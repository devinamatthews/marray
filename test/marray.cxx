#include "marray.hpp"

#include "gtest/gtest.h"

using namespace std;
using namespace MArray;
using namespace MArray::detail;
using namespace MArray::slice;
using namespace MArray::transpose;

TEST(marray, aligned_allocator)
{
    constexpr size_t alignment = 1ull<<21; // 2MB
    vector<double,aligned_allocator<double,alignment>> x(100);
    EXPECT_EQ(0, uintptr_t(x.data())&(alignment-1));
}

TEST(marray, marray_ref_assign_to_ref)
{
    marray<double,3> a(2, 3, 2);
    marray<double,4> b(4, 3, 2);
    marray<double,4> c(4, 4, 3, 2);

    for (int i = 0;i < 4*3*2;i++)
    {
        b.data()[i] = i;
    }

    for (int i = 0;i < 4*4*3*2;i++)
    {
        c.data()[i] = i;
    }

    a[0] = b[1][3];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+3*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+3*2]);
    }

    a[0] = c[1][1][3];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+3*2+4*3*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+3*2]);
    }
}

TEST(marray, marray_ref_assign_to_slice)
{
    marray<double,3> a(2, 3, 2);
    marray<double,3> b(4, 4, 2);
    marray<double,3> c(4, 4, 4, 2);

    for (int i = 0;i < 4*4*2;i++)
    {
        b.data()[i] = i;
    }

    for (int i = 0;i < 4*4*4*2;i++)
    {
        c.data()[i] = i;
    }

    a[0] = b[1][range(3)];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+4*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+3*2]);
    }

    a[0] = c[1][1][range(3)];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+4*2+4*4*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+3*2]);
    }
}

TEST(marray, const_marray_ref_data)
{
    const marray<double,3> a(2, 3, 2);
    EXPECT_EQ(a.data()+1*3*2+2*2, a[1][2].data());
    EXPECT_TRUE((std::is_same<const double*,decltype(a[1][2].data())>::value));
}

TEST(marray, marray_ref_data)
{
    marray<double,3> a(2, 3, 2);
    EXPECT_EQ(a.data()+1*3*2+2*2, a[1][2].data());
    EXPECT_TRUE((std::is_same<double*,decltype(a[1][2].data())>::value));
}

TEST(marray, const_marray_ref_view)
{
    marray<double,3> a_(2, 3, 2);
    const marray<double,3>& a = a_;

    for (int i = 0;i < 2*3*2;i++)
    {
        a.data()[i] = i;
    }

    const_marray_view<double,1> b = a[1][2];

    for (int i = 0;i < 2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, b[i]+1*3*2+2*2);
    }

    auto c = view(a[1][2]);

    for (int i = 0;i < 2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, c[i]+1*3*2+2*2);
    }

    EXPECT_TRUE((is_same<decltype(b),decltype(c)>::value));
}

TEST(marray, marray_ref_view)
{
    marray<double,3> a(2, 3, 2);

    for (int i = 0;i < 2*3*2;i++)
    {
        a.data()[i] = i;
    }

    marray_view<double,1> b = a[1][2];

    for (int i = 0;i < 2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, b[i]+1*3*2+2*2);
    }

    auto c = view(a[1][2]);

    for (int i = 0;i < 2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, c[i]+1*3*2+2*2);
    }

    EXPECT_TRUE((is_same<decltype(b),decltype(c)>::value));
}

TEST(marray, marray_slice_assign_to_ref)
{
    marray<double,3> a(2, 4, 2);
    marray<double,4> b(4, 3, 2);
    marray<double,4> c(4, 4, 3, 2);

    for (int i = 0;i < 4*3*2;i++)
    {
        b.data()[i] = i;
    }

    for (int i = 0;i < 4*4*3*2;i++)
    {
        c.data()[i] = i;
    }

    a[0][range(3)] = b[1];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+3*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+4*2]);
    }

    a[0][range(3)] = c[1][1];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+3*2+4*3*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+4*2]);
    }
}

TEST(marray, marray_slice_assign_to_slice)
{
    marray<double,3> a(2, 4, 2);
    marray<double,4> b(4, 4, 2);
    marray<double,4> c(4, 4, 4, 2);

    for (int i = 0;i < 4*4*2;i++)
    {
        b.data()[i] = i;
    }

    for (int i = 0;i < 4*4*4*2;i++)
    {
        c.data()[i] = i;
    }

    a[0][range(3)] = b[1][range(3)];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+4*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+4*2]);
    }

    a[0][range(3)] = c[1][1][range(3)];

    for (int i = 0;i < 3*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i+4*2+4*4*2, a.data()[i]);
        EXPECT_EQ(0, a.data()[i+4*2]);
    }
}

TEST(marray, const_marray_slice_data)
{
    const marray<double,3> a(2, 3, 2);
    EXPECT_EQ(a.data()+1*3*2+1*2, a[1][range(1,3)].data());
    EXPECT_TRUE((std::is_same<const double*,decltype(a[1][range(1,3)].data())>::value));
}

TEST(marray, marray_slice_data)
{
    marray<double,3> a(2, 3, 2);
    EXPECT_EQ(a.data()+1*3*2+1*2, a[1][range(1,3)].data());
    EXPECT_TRUE((std::is_same<double*,decltype(a[1][range(1,3)].data())>::value));
}

TEST(marray, const_marray_slice_view)
{
    marray<double,3> a_(2, 3, 2);
    const marray<double,3>& a = a_;

    for (int i = 0;i < 2*3*2;i++)
    {
        a.data()[i] = i;
    }

    const_marray_view<double,1> b = a[1][slice(1,3)];

    for (int i = 0;i < 2*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, b[i]+1*3*2+1*2);
    }

    auto c = view(a[1][slice(1,3)]);

    for (int i = 0;i < 2*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, c[i]+1*3*2+1*2);
    }

    EXPECT_TRUE((is_same<decltype(b),decltype(c)>::value));
}

TEST(marray, marray_slice_view)
{
    marray<double,3> a(2, 3, 2);

    for (int i = 0;i < 2*3*2;i++)
    {
        a.data()[i] = i;
    }

    const_marray_view<double,1> b = a[1][slice(1,3)];

    for (int i = 0;i < 2*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, b[i]+1*3*2+1*2);
    }

    auto c = view(a[1][slice(1,3)]);

    for (int i = 0;i < 2*2;i++)
    {
        SCOPED_TRACE(i);
        EXPECT_EQ(i, c[i]+1*3*2+1*2);
    }

    EXPECT_TRUE((is_same<decltype(b),decltype(c)>::value));
}

TEST(marray, align)
{
    EXPECT_EQ(0, align(0,8));
    EXPECT_EQ(8, align(4,8));
    EXPECT_EQ(8, align(1,8));
    EXPECT_EQ(8, align(7,8));
    EXPECT_EQ(8, align(8,8));
    EXPECT_EQ(280, align(277,8));
    EXPECT_EQ(160, align(160,8));
}

TEST(marray, marray_type)
{
    EXPECT_TRUE((is_same<double*,typename marray_type<false, double, 0>::type>::value));
    EXPECT_TRUE((is_same<marray_view<double,1>,typename marray_type<false, double, 1>::type>::value));
    EXPECT_TRUE((is_same<marray_view<double,2>,typename marray_type<false, double, 2>::type>::value));
    EXPECT_TRUE((is_same<marray_view<double,3>,typename marray_type<false, double, 3>::type>::value));
    EXPECT_TRUE((is_same<marray_view<double,4>,typename marray_type<false, double, 4>::type>::value));
    EXPECT_TRUE((is_same<const double*,typename marray_type<true, double, 0>::type>::value));
    EXPECT_TRUE((is_same<const_marray_view<double,1>,typename marray_type<true, double, 1>::type>::value));
    EXPECT_TRUE((is_same<const_marray_view<double,2>,typename marray_type<true, double, 2>::type>::value));
    EXPECT_TRUE((is_same<const_marray_view<double,3>,typename marray_type<true, double, 3>::type>::value));
    EXPECT_TRUE((is_same<const_marray_view<double,4>,typename marray_type<true, double, 4>::type>::value));
}

TEST(marray, is_slice)
{
    EXPECT_TRUE((is_slice<decltype(range(4))>::value));
    EXPECT_TRUE((is_slice<decltype(all)>::value));
    EXPECT_FALSE((is_slice<int>::value));
}

TEST(marray, num_slices)
{
    EXPECT_EQ(3, (num_slices<int,decltype(range(3)),double,decltype(all),decltype(range(1,5)),long>::N));
}

TEST(marray, return_type)
{
    EXPECT_TRUE((is_same<marray_view<double,3>,typename return_type<false,double,int,decltype(range(3)),double,decltype(all),decltype(range(1,5)),long>::type>::value));
    EXPECT_TRUE((is_same<const_marray_view<double,3>,typename return_type<true,double,int,decltype(range(3)),double,decltype(all),decltype(range(1,5)),long>::type>::value));
    EXPECT_TRUE((is_same<double*,typename return_type<false,double,int,long>::type>::value));
    EXPECT_TRUE((is_same<const double*,typename return_type<false,double,int,long>::type>::value));
}

#if 0

/*
 * These helper structs apply the given set of indices to an array
 * to generate either a reference to a specific element or a
 * lower-dimensional view of the array.
 */
template <typename RT, typename Index, typename... Indices>
struct get_slice
{
    template <typename ArrayOrSlice>
    RT operator()(ArrayOrSlice&& x, Index idx0, Indices... idx)
    {
        return get_slice<RT, Indices...>()(x[idx0], idx...);
    }
};

template <typename RT, typename Index>
struct get_slice<RT, Index>
{
    template <typename ArrayOrSlice>
    RT operator()(ArrayOrSlice&& x, Index idx0)
    {
        return x[idx0];
    }
};

#endif

TEST(marray, are_integral)
{
    EXPECT_TRUE((are_integral<0, are_empty>::value));
    EXPECT_TRUE((are_integral<1, are_empty, int>::value));
    EXPECT_TRUE((are_integral<2, are_empty, int, long>::value));
    EXPECT_TRUE((are_integral<3, are_empty, int, long, bool>::value));
    EXPECT_TRUE((are_integral<4, are_empty, int, long, bool, char>::value));

    EXPECT_FALSE((are_integral<1, are_empty>::value));
    EXPECT_FALSE((are_integral<2, are_empty>::value));
    EXPECT_FALSE((are_integral<3, are_empty>::value));
    EXPECT_FALSE((are_integral<4, are_empty>::value));
    EXPECT_FALSE((are_integral<2, are_empty, int>::value));
    EXPECT_FALSE((are_integral<3, are_empty, int>::value));
    EXPECT_FALSE((are_integral<4, are_empty, int>::value));
    EXPECT_FALSE((are_integral<3, are_empty, int, int>::value));
    EXPECT_FALSE((are_integral<4, are_empty, int, int>::value));
    EXPECT_FALSE((are_integral<4, are_empty, int, int, int>::value));

    EXPECT_FALSE((are_integral<0, are_empty, int>::value));
    EXPECT_FALSE((are_integral<0, are_empty, int, long>::value));
    EXPECT_FALSE((are_integral<1, are_empty, int, long>::value));
    EXPECT_FALSE((are_integral<0, are_empty, int, long, bool>::value));
    EXPECT_FALSE((are_integral<1, are_empty, int, long, bool>::value));
    EXPECT_FALSE((are_integral<2, are_empty, int, long, bool>::value));
    EXPECT_FALSE((are_integral<0, are_empty, int, long, bool, char>::value));
    EXPECT_FALSE((are_integral<1, are_empty, int, long, bool, char>::value));
    EXPECT_FALSE((are_integral<2, are_empty, int, long, bool, char>::value));
    EXPECT_FALSE((are_integral<3, are_empty, int, long, bool, char>::value));

    EXPECT_FALSE((are_integral<1, are_empty, double>::value));
    EXPECT_FALSE((are_integral<2, are_empty, double, int*>::value));
    EXPECT_FALSE((are_integral<3, are_empty, int, long, vector<char>>::value));
    EXPECT_FALSE((are_integral<4, are_empty, int, int*, bool, double>::value));
}

TEST(marray, are_empty)
{
    EXPECT_TRUE((are_empty<>::value));
    EXPECT_FALSE((are_empty<int>::value));
    EXPECT_FALSE((are_empty<int,double>::value));
}

struct do_nothing
{
    template <typename... Args> do_nothing(Args&&...) {}
};

TEST(marray, set_len)
{
    array<size_t, 4> len;
    set_len<4, do_nothing, int, int, long, char>(nullptr, len, 1, 2, 3, 4);
    EXPECT_EQ(len[0], 1);
    EXPECT_EQ(len[1], 2);
    EXPECT_EQ(len[2], 3);
    EXPECT_EQ(len[3], 4);
}

TEST(marray, are_reset_args_direct)
{
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, double&&>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, int>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, uninitialized_t>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, double&&, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, int, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template direct<int, long, char, bool, uninitialized_t, Layout>::value));

    EXPECT_TRUE((are_reset_args<double,1>::template direct<int>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, double&&>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, int>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, uninitialized_t>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, double&&, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, int, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template direct<int, uninitialized_t, Layout>::value));

    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, double, bool>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, double, char, bool, double&&>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, double, char, bool, int>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, double, bool, uninitialized_t>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, double, bool, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, bool, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, double, bool, double&&, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, double, bool, int, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, double, uninitialized_t, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, bool, vector<double>, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, bool, double&&, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, bool, int, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, bool, uninitialized_t, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template direct<int, long, char, bool, vector<double>, vector<double>>::value));

    EXPECT_FALSE((are_reset_args<double,1>::template direct<double>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, double&&>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, int>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, uninitialized_t>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<char, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, double&&, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, int, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<double, uninitialized_t, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<char,vector<double>, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<char,double&&, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<char,int, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<char,uninitialized_t, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template direct<char,vector<double>, vector<double>>::value));
}

TEST(marray, are_reset_args_const_view)
{
    EXPECT_TRUE((are_reset_args<double,4>::template const_view<int, long, char, bool, const double*>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template const_view<int, long, char, bool, const double*, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template const_view<int, long, char, bool, const double*, int, long, char, bool>::value));

    EXPECT_TRUE((are_reset_args<double,1>::template const_view<int, const double*>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template const_view<int, const double*, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template const_view<int, const double*, int>::value));

    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, double, char, bool, const double*>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, double, char, bool, const double*, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, double, char, bool, const double*, int, long, char, bool>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, const double*, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>, int, long, char, bool>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, const double*, int, double, char, bool>::value));

    EXPECT_FALSE((are_reset_args<double,1>::template const_view<double, const double*>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<double, const double*, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<double, const double*, int>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, const double*, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>, int>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, const double*, double>::value));
}

TEST(marray, are_reset_args_view)
{
    EXPECT_TRUE((are_reset_args<double,4>::template const_view<int, long, char, bool, double*>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template const_view<int, long, char, bool, double*, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,4>::template const_view<int, long, char, bool, double*, int, long, char, bool>::value));

    EXPECT_TRUE((are_reset_args<double,1>::template const_view<int, double*>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template const_view<int, double*, Layout>::value));
    EXPECT_TRUE((are_reset_args<double,1>::template const_view<int, double*, int>::value));

    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, double, char, bool, double*>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, double, char, bool, double*, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, double, char, bool, double*, int, long, char, bool>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, double*, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, vector<double>, int, long, char, bool>::value));
    EXPECT_FALSE((are_reset_args<double,4>::template const_view<int, long, char, bool, double*, int, double, char, bool>::value));

    EXPECT_FALSE((are_reset_args<double,1>::template const_view<double, double*>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<double, double*, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<double, double*, int>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>, Layout>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, double*, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>, vector<double>>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, vector<double>, int>::value));
    EXPECT_FALSE((are_reset_args<double,1>::template const_view<long, double*, double>::value));
}

#if 0

struct reset_helper
{
    template <typename T, unsigned ndim>
    reset_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
    {
        array.reset(len, val, layout);
    }

    template <typename T, unsigned ndim>
    reset_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const uninitialized_t& u, Layout layout = DEFAULT)
    {
        array.reset(len, u, layout);
    }

    template <typename T, unsigned ndim>
    reset_helper(const_marray_view<T, ndim>& array, std::array<size_t, ndim>& len, const T* ptr, Layout layout = DEFAULT)
    {
        array.reset(len, const_cast<T*>(ptr), layout);
    }

    template <typename T, unsigned ndim, typename... Args>
    reset_helper(const_marray_view<T, ndim>& array, std::array<size_t, ndim>& len, const T* ptr, Args... args)
    {
        array.reset(len, const_cast<T*>(ptr), make_array(args...));
    }
};

template <typename T, typename ndim, typename... Args>
void reset(const_marray_view<T, ndim>& array, const Args&... args)
{
    std::array<size_t, ndim> len;
    set_len<T, ndim, reset_helper, Args...> tmp(array, len, args...);
}

template <typename T, typename ndim, typename... Args>
void reset(marray_view<T, ndim>& array, const Args&... args)
{
    std::array<size_t, ndim> len;
    set_len<T, ndim, reset_helper, Args...> tmp(array, len, args...);
}

template <typename T, typename ndim, typename... Args>
void reset(marray<T, ndim>& array, const Args&... args)
{
    std::array<size_t, ndim> len;
    set_len<T, ndim, reset_helper, Args...> tmp(array, len, args...);
}

/*
 * Helper class to determine if the arguments Args... are of one of the
 * forms which are acceptable for passing to marray<T, ndim>::resize.
 */
template <typename T, unsigned ndim>
struct are_resize_args
{
    template <typename... Args>
    struct direct_helper : std::false_type {};

    template <>
    struct direct_helper<> : std::true_type {};

    template <typename Arg>
    struct direct_helper<Arg>
    : std::integral_constant<bool, std::is_convertible<Arg,T>::value> {};

    template <typename... Args>
    using direct =  std::integral_constant<bool,
         are_integral<ndim, direct_helper,
                      typename std::decay<Args>::type...>::value>;
};

struct resize_helper
{
    template <typename T, unsigned ndim>
    resize_helper(marray<T, ndim>& array, std::array<size_t, ndim>& len, const T& val=T())
    {
        array.resize(len, val);
    }
};

template <typename T, typename ndim, typename... Args>
void resize(marray<T, ndim>& array, const Args&... args)
{
    std::array<size_t, ndim> len;
    set_len<T, ndim, resize_helper, Args...> tmp(array, len, args...);
}

template <typename T, unsigned ndim, typename Allocator=aligned_allocator<T,MARRAY_BASE_ALIGNMENT>>
class marray_base
{
    static_assert(ndim > 0, "0-dimensional marrays are not allowed.");

    template <typename T_, unsigned ndim_> friend class detail::marray_base;
    template <typename T_, unsigned ndim_> friend class const_marray_view;
    template <typename T_, unsigned ndim_> friend class marray_view;
    template <typename T_, unsigned ndim_> friend class marray;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

    public:
        typedef unsigned idx_type;
        typedef size_t size_type;
        typedef ptrdiff_t stride_type;
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T& reference;
        typedef const T& const_reference;

    protected:
        struct mem_t : Allocator
        {
            typedef std::allocator_traits<Allocator> alloc_traits;

            mem_t(Allocator alloc) : Allocator(std::move(alloc)) {}

            ~mem_t() { clear(); }

            void allocate(size_type size, unitialized_t u)
            {
                clear();
                size_ = size;
                is_alloced_ = true;
                data_ = alloc_traits::allocate(*this, size);
            }

            void allocate(size_type size, const T& value=T())
            {
                allocate(size, uninitialized);
                for (size_type i = 0;i < size;i++)
                {
                    alloc_traits::construct(*this, data_+i, value);
                }
            }

            void set(pointer data)
            {
                clear();
                data_ = data;
            }

            void set(mem_t&& other)
            {
                using std::swap;
                clear();
                swap(other.is_alloced_, is_alloced_);
                swap(other.data_, data_);
                swap(other.size_, size_);
            }

            void set(const mem_t& other)
            {
                clear();
                size_ = size;
                is_alloced_ = true;
                data_ = alloc.allocate(other.size_);
            }

            void clear()
            {
                if (is_alloced_ && data_)
                {
                    for (size_type i = 0;i < size;i++)
                    {
                        alloc_traits::destroy(*this, data_+i);
                    }
                    alloc_traits::deallocate(*this, data_, size_);
                }
                data_ = nullptr;
                is_alloced_ = false;
            }

            pointer data_ = nullptr;
            size_t size_ = 0;
            bool is_alloced_ = false;
        } mem_;

        std::array<idx_type,ndim> len_ = {};
        std::array<stride_type,ndim> stride_ = {};

        marray_base() {}

        marray_base(const marray_base& other) = delete;

        marray_base(marray_base&& other) = delete;

        marray_base& operator=(const marray_base& other) = delete;

        marray_base& operator=(marray_base&& other) = delete;

        marray_base& operator=(const marray_view<T,ndim>& other)
        {
            copy(other, static_cast<marray_view<T,ndim>&>(*this));
            return *this;
        }

        void view(const marray_base& other)
        {
            mem_.set(other.mem_.data_);
            len_ = other.len_;
            stride_ = other.stride_;
        }

        void copy(const marray_base& other, Layout layout=DEFAULT)
        {
            if (std::is_scalar<T>::value)
            {
                reset(other.len_, uninitialized, layout);
            }
            else
            {
                reset(other.len_, T(), layout);
            }

            *this = other;
        }

        void copy(marray_base&& other)
        {
            other.mem_.set(std::move(other.mem_));
            len_ = other.len_;
            stride_ = other.stride_;
        }

        template <typename U>
        void reset(const std::array<U, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
        {
            reset(len, uninitialized, layout);
            std::uninitialized_fill_n(mem_.data_, mem_.size_, val);
        }

        template <typename U>
        void reset(const std::array<U, ndim>& len, uninitialized_t u, Layout layout = DEFAULT)
        {
            std::copy_n(len.begin(), ndim, len_.begin());

            size_t size;
            if (layout == ROW_MAJOR)
            {
                stride_[ndim-1] = 1;
                for (unsigned i = ndim-1;i > 0;i--)
                {
                    stride_[i-1] = stride_[i]*len_[i];
                }
                size = stride_[0]*len_[0];
            }
            else
            {
                stride_[0] = 1;
                for (unsigned i = 1;i < ndim;i++)
                {
                    stride_[i] = stride_[i-1]*len_[i-1];
                }
                size = stride_[ndim-1]*len_[ndim-1];
            }

            mem_.allocate(size, u);
        }

        template <typename U>
        marray_view<T, ndim> permute(const std::array<U, ndim>& perm)
        {
            marray_view<T, ndim> view;

            view.mem_.set(mem_.data_);

            for (unsigned i = 0;i < ndim;i++)
            {
                assert(0 <= perm[i] && perm[i] < ndim);
                for (unsigned j = 0;j < i;j++) assert(perm[i] != perm[j]);
            }

            for (unsigned i = 0;i < ndim;i++)
            {
                view.len_[i] = len_[perm[i]];
                view.stride_[i] = stride_[perm[i]];
            }

            return view;
        }

        template <typename... Args>
        typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, marray_view<T, ndim>>::type
        permute(Args&&... args)
        {
            return permute(make_array(std::forward<Args>(args)...));
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        reference front()
        {
            assert(length() > 0);
            return mem_.data_[0];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        reference front(unsigned dim)
        {
            assert(dim == 0);
            return front();
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
        marray_view<T, ndim-1> front(unsigned dim)
        {
            assert(dim < ndim);
            assert(len_[dim] > 0);

            marray_view<T, ndim-1> view;
            view.mem_.set(mem_.data_);
            copy_n(stride_.begin(), dim, view.stride_.begin());
            copy_n(stride_.begin()+dim+1, ndim-dim-1, view.stride_.begin()+dim);
            copy_n(len_.begin(), dim, view.len_.begin());
            copy_n(len_.begin()+dim+1, ndim-dim-1, view.len_.begin()+dim);

            return view;
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        reference back()
        {
            assert(length() > 0);
            return mem_.data_[(length()-1)*stride()];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        reference back(unsigned dim)
        {
            assert(dim == 0);
            return back();
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
        marray_view<T, ndim-1> back(unsigned dim)
        {
            marray_view<T, ndim-1> view = front(dim);
            view.data_ += (len_[dim]-1)*stride_[dim];
            return view;
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        reference operator[](idx_type i)
        {
            assert(i < length());
            return mem_.data_[i*stride()];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
        marray_ref<T, ndim, 2> operator[](idx_type i)
        {
            assert(i < len_[0]);
            return const_marray_ref<T, ndim, 2>(*this, (stride_type)0, i);
        }

        template <typename I>
        marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x)
        {
            assert(x.front() >= 0 && x.back() <= len_[0]);
            return const_marray_slice<T, ndim, 2, 1>(*this, (stride_type)0, {}, {}, x);
        }

        marray_slice<T, ndim, 2, 1> operator[](slice::all_t x)
        {
            return const_marray_slice<T, ndim, 2, 1>(*this, (stride_type)0, {}, {}, range(idx_type(), len_[0]));
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        reference operator()(idx_type i)
        {
            return (*this)[i];
        }

        template <typename... Indices>
        typename std::enable_if<(sizeof...(Indices) == ndim) && (sizeof...(Indices) > 1),
                                typename detail::return_type<false, T, Indices...>::type>::type
        operator()(Indices... idx)
        {
            return detail::get_slice<typename detail::return_type<true, T, Indices...>::type, Indices...>()(*this, idx...);
        }

        pointer data()
        {
            return mem_.data_;
        }

        void rotate_dim(unsigned dim, idx_type shift)
        {
            assert(dim < ndim);

            idx_type n = len_[dim];
            stride_type s = stride_[dim];

            shift = shift%n;
            if (shift < 0) shift += n;

            if (shift == 0) return;

            std::array<idx_type, ndim-1> sublen;
            std::array<stride_type, ndim-1> substride;

            std::copy_n(len_.begin(), dim, sublen.begin());
            std::copy_n(len_.begin()+dim+1, ndim-dim-1, sublen.begin()+dim);

            std::copy_n(stride_.begin(), dim, substride.begin());
            std::copy_n(stride_.begin()+dim+1, ndim-dim-1, substride.begin()+dim);

            pointer p = mem_.data_;
            miterator<idx_type, stride_type, ndim> it(sublen, substride);
            while (it.next(p))
            {
                pointer a = p;
                pointer b = p+(shift-1)*s;
                for (idx_type i = 0;i < shift/2;i++)
                {
                    std::swap(*a, *b);
                    a += s;
                    b -= s;
                }

                a = p+shift*s;
                b = p+(n-1)*s;
                for (idx_type i = 0;i < (n-shift)/2;i++)
                {
                    std::swap(*a, *b);
                    a += s;
                    b -= s;
                }

                a = p;
                b = p+(n-1)*s;
                for (idx_type i = 0;i < n/2;i++)
                {
                    std::swap(*a, *b);
                    a += s;
                    b -= s;
                }
            }
        }

        template <typename U>
        void rotate(const std::array<U, ndim>& shift)
        {
            for (unsigned dim = 0;dim < ndim;dim++)
            {
                rotate_dim(dim, shift[dim]);
            }
        }

        template <typename... Args>
        typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value>::type
        rotate(Args&&... args)
        {
            rotate(make_array(std::forward<Args>(args)...));
        }

    public:
        void reset()
        {
            mem_.clear();
            len_.fill(0);
            stride_.fill(0);
        }

        template <typename U>
        void reset(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
        {
            if (layout == ROW_MAJOR)
            {
                stride_[ndim-1] = 1;
                for (unsigned i = ndim-1;i > 0;i--)
                {
                    stride_[i-1] = stride_[i]*len[i];
                }
            }
            else
            {
                stride_[0] = 1;
                for (unsigned i = 1;i < ndim;i++)
                {
                    stride_[i] = stride_[i-1]*len[i-1];
                }
            }

            reset(len, ptr, stride_);
        }

        template <typename U, typename V>
        void reset(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
        {
            mem_.set(ptr);
            std::copy_n(len.begin(), ndim, len_.begin());
            std::copy_n(stride.begin(), ndim, stride_.begin());
        }

        template <typename U>
        const_marray_view<T, ndim> permute(const std::array<U, ndim>& perm) const
        {
            return const_cast<marray_base&>(*this).permute(perm);
        }

        template <typename... Args>
        typename std::enable_if<detail::are_integral<ndim, detail::are_empty, Args...>::value, const_marray_view<T, ndim>>::type
        permute(Args&&... args) const
        {
            return permute(make_array(std::forward<Args>(args)...));
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        const_reference front() const
        {
            assert(length() > 0);
            return mem_.data_[0];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        const_reference front(unsigned dim) const
        {
            assert(dim == 0);
            return front();
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
        const_marray_view<T, ndim-1> front(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).front(dim);
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        const_reference back() const
        {
            assert(length() > 0);
            return mem_.data_[(length()-1)*stride()];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        const_reference back(unsigned dim) const
        {
            assert(dim == 0);
            return back();
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
        const_marray_view<T, ndim-1> back(unsigned dim) const
        {
            return const_cast<marray_base&>(*this).back(dim);
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        const_reference operator[](idx_type i) const
        {
            assert(i < length());
            return mem_.data_[i*stride()];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_!=1>::type>
        const_marray_ref<T, ndim, 2> operator[](idx_type i) const
        {
            return const_cast<marray_base&>(*this)[i];
        }

        template <typename I>
        const_marray_slice<T, ndim, 2, 1> operator[](const range_t<I>& x) const
        {
            return const_cast<marray_base&>(*this)[x];
        }

        const_marray_slice<T, ndim, 2, 1> operator[](slice::all_t x) const
        {
            return const_cast<marray_base&>(*this)[x];
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        const_reference operator()(idx_type i) const
        {
            return (*this)[i];
        }

        template <typename... Indices>
        typename std::enable_if<(sizeof...(Indices) == ndim) && (sizeof...(Indices) > 1),
                                typename detail::return_type<true, T, Indices...>::type>::type
        operator()(Indices... idx) const
        {
            return const_cast<marray_base&>(*this)(idx...);
        }

        const_pointer data() const
        {
            return mem_.data_;
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        idx_type length() const
        {
            return len_[0];
        }

        idx_type length(unsigned dim) const
        {
            return len_[dim];
        }

        const std::array<idx_type, ndim>& lengths() const
        {
            return len_;
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        stride_type stride() const
        {
            return stride_[0];
        }

        stride_type stride(unsigned dim) const
        {
            return stride_[dim];
        }

        const std::array<stride_type, ndim>& strides() const
        {
            return stride_;
        }
};
}

/*
 * This class represents a multidimensional array of elements of type T,
 * with ndim > 0 dimensions. An array may either own its own data (in
 * which case it is responsible for allocation and deallocation), or
 * represent a "view" of all or part of another array. It is the user's
 * responsibility to ensure that views are not referenced once the original
 * array has been destructed or reset with the resize() or reset() methods,
 * since no internal reference counting is performed.
 *
 * Arrays are constructed using one of the following forms:
 *
 * 1) const_marray<T, ndim>()
 *
 *  Constructs an array with all dimensions of zero length.
 *
 * 2) const_marray<T, ndim>(const const_marray<T, ndim>& other)
 *
 *  Creates a copy of the given array. If the array is a view, then so is
 *  the copy.
 *
 * 3) const_marray<T, ndim>(const_marray<T, ndim>&& other)
 *
 *  Move-constructs a copy of the given array. If the array is a view, then
 *  this behaves identically to the copy constructor. Otherwise, the given
 *  array is reinitialized to a state as if it were constructed using the
 *  empty constructor.
 *
 * 4) const_marray<T, ndim>(const const_marray<T, ndim>& other,
 *                          const construct_view_t& cv)
 *
 *  Constructs a view of the given array, even if it is not a view itself.
 *
 * 5) const_marray<T, ndim>(const const_marray<T, ndim>& other,
 *                          const construct_copy_t& cp)
 *
 *  Constructs a new (freshly allocated) copy of the given array, even if
 *  it is a view.
 *
 * 6.1) template <typename U>
 *      explicit const_marray<T, ndim>(const std::array<U, ndim>& len,
 *                                     Layout layout = DEFAULT)
 *
 * 6.2) template <typename U>
 *      const_marray<T, ndim>(const std::array<U, ndim>& len,
 *                            const T& val, Layout layout = DEFAULT)
 *
 * 6.3) template <typename U>
 *      const_marray<T, ndim>(const std::array<U, ndim>& len,
 *                            const uninitialized_t& u, Layout layout = DEFAULT)
 *
 *  Constructs an array with the given lengths and data layout. In the first
 *  form the elements are default-initialized. In the second form they are
 *  copy-initialized from val. In the third form they are not initialized.
 *
 * 7.4) template <typename U>
 *      explicit const_marray<T, ndim>(const std::array<U, ndim>& len,
 *                                     unsigned align_at, Layout layout = DEFAULT)
 *
 * 7.5) template <typename U>
 *      const_marray<T, ndim>(const std::array<U, ndim>& len,
 *                            const T& val, unsigned align_at, Layout layout = DEFAULT)
 *
 * 7.6) template <typename U>
 *      const_marray<T, ndim>(const std::array<U, ndim>& len, const uninitialized_t& u,
 *                            unsigned align_at, Layout layout = DEFAULT)
 *
 *  Constructs an array with the given lengths and data layout. In the first
 *  form the elements are default-initialized. In the second form they are
 *  copy-initialized from val. In the third form they are not initialized.
 *  The strides of the indices are adjusted such that stride (align_at-1)
 *  for row-major or stride (align_at) for column-major storage is a
 *  multiple of the processor's vector length. align_at = 0 specifies no
 *  additional alignment. In all cases (even in constructors without an
 *  align_at parameter), the base pointer of the array is aligned to a
 *  multiple of the cache line size.
 *
 * 8.1) template <typename I1, ..., typename IN>
 *      explicit const_marray<T, ndim>(I1 len1, ..., IN lenN, Layout layout = DEFAULT)
 *
 * 8.2) template <typename I1, ..., typename IN>
 *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const T& val, Layout layout = DEFAULT)
 *
 * 8.3) template <typename I1, ..., typename IN>
 *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const unsigned_t& u,
 *                            Layout layout = DEFAULT)
 *
 * 9.1) template <typename I1, ..., typename IN>
 *      const_marray<T, ndim>(I1 len1, ..., IN lenN, unsigned align_at,
 *                            Layout layout = DEFAULT)
 *
 * 9.1) template <typename I1, ..., typename IN>
 *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const T& val,
 *                            unsigned align_at, Layout layout = DEFAULT)
 *
 * 9.1) template <typename I1, ..., typename IN>
 *      const_marray<T, ndim>(I1 len1, ..., IN lenN, const unsigned_t& u,
 *                            unsigned align_at, Layout layout = DEFAULT)
 *
 *  The same as 6.1-6.3 and 7.1-7.3, except that the lengths are given as
 *  distinct parameters with possibly difference integral types. There must
 *  be exactly ndim length parameters.
 *
 * 10.1) template <typename U>
 *       const_marray<T, ndim>(const std::array<U, ndim>& len, const T* ptr,
 *                             Layout layout = DEFAULT)
 *
 * 10.2) template <typename U, typename V>
 *       const_marray<T, ndim>(const std::array<U, ndim>& len, const T* ptr,
 *                             const std::array<V, ndim>& stride)
 *
 *  Constructs an array with the given lengths and data layout, using ptr
 *  as the array data. The given pointer is referenced directly by the
 *  array, providing a "view" of an externally-allocated data segment.
 *  No internal allocation or copying is performed. In the second form the
 *  strides are specified directly, instead of inferred from the data
 *  layout. It is the user's reposnsibility to ensure that the information
 *  given results in legal accesses of ptr.
 *
 * 11.1) template <typename I1, ..., typename IN>
 *       const_marray<T, ndim>(I1 len1, ..., IN lenN, const T* ptr, Layout layout = DEFAULT)
 *
 * 11.2) template <typename I1, ..., typename IN, typename S1, ..., typename SN>
 *       const_marray<T, ndim>(I1 len1, ..., IN lenN, const T* ptr,
 *                             S1 stride1, ..., SN strideN)
 *
 *  The same as 10.1 and 10.2, except that the lengths and, in the second form
 *  the strides, are given as distinct parameters with possibly different
 *  integral types. There must be exactly ndim length and stride parameters.
 *
 * 12.1) template <...>
 *       const_marray<T, ndim>(array[][][]...)
 *
 * 12.2) template <...>
 *       const_marray<T, ndim>(array(...))
 *
 *  Creates a view from the given array reference as if it were repeatedly
 *  indexed with slice::all.
 *
 * Individual elements of the array may be accessed using one of two forms:
 *
 * 1) C-style access:
 *
 *  const_marray<double, 4> arr(4, 3, 9, 5);
 *  double el = arr[1][1][6][0];
 *
 * 2) FORTRAN-style access (using 0-based indices):
 *
 *  const_marray<double, 4> arr(4, 3, 9, 5);
 *  double el = arr(1, 1, 6, 0);
 *
 * These two methods are entirely equivalent. In addition, a (possibly
 * lower-dimensional) view of all or part of the array may be created
 * using this same notation. Instead of an integral index, a range_t
 * or slice::all_t argument denotes that the given range or all elements,
 * respectively, of a particular dimension are requested. For example:
 *
 *  using slice::all;
 *  const_marray<double, 4> arr(4, 3, 9, 5);
 *  const_marray<double, 2> sub_arr = arr[all][2][range(1,5)][1];
 *  // sub_arr has lengths (4, 4)
 *  double el = sub_arr[2][2];
 *
 * The last two constructors detailed above allow expressions of the type:
 *
 *  using slice::all;
 *  const_marray<double, 4> arr(4, 3, 9, 5);
 *  const_marray<double, 3> sub_arr = arr[all][2][range(1,5)];
 *  // this is identical to:
 *  // const_marray<double, 3> sub_arr = arr[all][2][range(1,5)][all];
 *  double el = sub_arr[2][2][0];
 *
 * Note that only the trailing dimensions may be implicitly sliced using
 * this notation.
 *
 * The type marray<T, ndim> provides non-const access to elements and for
 * the creation of non-const views.
 */

template <typename T, unsigned ndim>
class const_marray_view : detail::marray_base<T,ndim>
{
    static_assert(ndim > 0, "0-dimensional marrays are not allowed.");

    template <typename T_, unsigned ndim_> friend class detail::marray_base;
    template <typename T_, unsigned ndim_> friend class const_marray_view;
    template <typename T_, unsigned ndim_> friend class marray_view;
    template <typename T_, unsigned ndim_> friend class marray;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

    protected:
        typedef detail::marray_base<T,ndim> base;

    public:
        using typename base::idx_type;
        using typename base::stride_type;
        using typename base::value_type;
        using typename base::pointer;
        using typename base::const_pointer;
        using typename base::reference;
        using typename base::const_reference;

        const_marray_view() {}

        const_marray_view(const const_marray_view& other)
        {
            base::reset(other);
        }

        template <typename U>
        const_marray_view(const std::array<U, ndim>& len, const_pointer ptr, Layout layout = DEFAULT)
        {
            reset(len, ptr, layout);
        }

        template <typename U, typename V>
        const_marray_view(const std::array<U, ndim>& len, const_pointer ptr, const std::array<V, ndim>& stride)
        {
            reset(len, ptr, stride);
        }

        template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim, Args...>::value>::type>
        explicit const_marray_view(Args&&... args)
        {
            reset(std::forward<Args>(args)...);
        }

        template <typename U>
        void reset(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
        {
            base::reset(len, const_cast<pointer>(ptr), layout);
        }

        template <typename U, typename V>
        void reset(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
        {
            base::reset(len, const_cast<pointer>(ptr), stride);
        }

        template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim>::const_view<Args...>::value>::type>
        void reset(Args&&... args)
        {
            detail::reset(*this, std::forward<Args>(args)...);
        }

        operator const marray_view<T,ndim>&() const { return static_cast<const marray_view<T,ndim>&>(*this); }
};

template <typename T, unsigned ndim>
class marray_view : public const_marray_view<T, ndim>
{
    template <typename T_, unsigned ndim_> friend class detail::marray_base;
    template <typename T_, unsigned ndim_> friend class const_marray_view;
    template <typename T_, unsigned ndim_> friend class marray_view;
    template <typename T_, unsigned ndim_> friend class marray;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

    protected:
        typedef detail::marray_base<T,ndim> base;

    public:
        using typename base::idx_type;
        using typename base::stride_type;
        using typename base::value_type;
        using typename base::pointer;
        using typename base::const_pointer;
        using typename base::reference;
        using typename base::const_reference;

    public:
        marray_view() {}

        marray_view(marray_view& other)
        : const_marray_view<T, ndim>(other) {}

        marray_view(marray_view&& other)
        : const_marray_view<T, ndim>(other) {}

        template <typename U>
        marray_view(const std::array<U, ndim>& len, pointer ptr, Layout layout = DEFAULT)
        : const_marray_view<T, ndim>(len, ptr, layout) {}

        template <typename U, typename V>
        marray_view(const std::array<U, ndim>& len, pointer ptr, const std::array<V, ndim>& stride)
        : const_marray_view<T, ndim>(len, ptr, stride) {}

        template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim, Args...>::value>::type>
        marray_view(Args&&... args)
        : const_marray_view<T, ndim>(std::forward<Args>(args)...) {}

        template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim>::view<Args...>::value>::type>
        void reset(Args&&... args)
        {
            detail::reset(*this, std::forward<Args>(args)...);
        }

        using base::operator=;

        using base::permute;

        using base::front;

        using base::back;

        using base::operator[];

        using base::operator();

        using base::data;

        using base::rotate_dim;

        using base::rotate;
};

template <typename T, unsigned ndim>
class marray : public marray_view<T, ndim>
{
    template <typename T_, unsigned ndim_> friend class detail::marray_base;
    template <typename T_, unsigned ndim_> friend class const_marray_view;
    template <typename T_, unsigned ndim_> friend class marray_view;
    template <typename T_, unsigned ndim_> friend class marray;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_> friend class const_marray_ref;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class marray_slice;
    template <typename T_, unsigned ndim_, unsigned dim_, unsigned newdim_> friend class const_marray_slice;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_reset_args;
    template <typename T_, unsigned ndim_, typename... Args> friend class detail::are_resize_args;

    protected:
        typedef detail::marray_base<T,ndim> base;

        Layout layout_;

    public:
        using typename base::idx_type;
        using typename base::stride_type;
        using typename base::value_type;
        using typename base::pointer;
        using typename base::const_pointer;
        using typename base::reference;
        using typename base::const_reference;

    public:
        marray() : layout_(DEFAULT) {}

        marray(const marray_view<T, ndim>& other, Layout layout=DEFAULT)
        : layout_(layout)
        {
            base::copy(other, layout);
        }

        marray(const marray& other, Layout layout=DEFAULT)
        : layout_(layout)
        {
            base::copy(other, layout);
        }

        marray(marray&& other, Layout layout=DEFAULT)
        : layout_(layout)
        {
            base::copy(std::move(other), layout);
        }

        template <typename U>
        marray(const std::array<U, ndim>& len, const T& val=T(), Layout layout = DEFAULT)
        : layout_(layout)
        {
            base::reset(len, val, layout);
        }

        template <typename U>
        marray(const std::array<U, ndim>& len, uninitialized_t u, Layout layout = DEFAULT)
        : layout_(layout)
        {
            base::reset(len, u, layout);
        }

        template <typename... Args, typename std::enable_if<detail::are_reset_args<T, ndim, Args...>::value>::type>
        marray(Args&&... args)
        {
            reset(std::forward<Args>(args)...);
            layout_ = stride_[0] == 1 ? COLUMN_MAJOR : ROW_MAJOR;
        }

        using base::reset;

        void clear()
        {
            reset();
        }

        template <typename... Args, typename=typename std::enable_if<detail::are_reset_args<T, ndim>::direct<Args...>::value>::type>
        void reset(Args&&... args)
        {
            detail::reset(*this, std::forward<Args>(args)...);
        }

        template <typename U>
        void resize(const std::array<U, ndim>& len, const T& val=T())
        {
            marray a(std::move(*this));
            reset(len, val, layout_);
            marray_view<T,ndim> b(*this);

            /*
             * It is OK to change the geometry of 'a' even if it is not
             * a view since it is about to go out of scope.
             */
            for (unsigned i = 0;i < ndim;i++)
            {
                a.len_[i] = b.len_[i] = std::min(a.len_[i], b.len_[i]);
            }

            b = a;
        }

        template <typename... Args, typename=typename std::enable_if<detail::are_resize_args<T, ndim>::direct<Args...>::value>::type>
        void resize(Args&&... args)
        {
            detail::resize(*this, std::forward<Args>(args)...);
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        void push_back(const T& x)
        {
            resize(length()+1);
            back() = x;
        }

        template <unsigned ndim_=ndim, typename=typename std::enable_if<ndim_==1>::type>
        void pop_back()
        {
            resize(length()-1);
        }

        void push_back(unsigned dim, const marray_view<T, ndim-1>& x)
        {
            assert(dim < ndim);

            for (unsigned i = 0, j = 0;i < ndim;i++)
            {
                if (i != dim)
                {
                    assert(len_[i] == x.len_[j++]);
                }
            }

            std::array<idx_type, ndim> len = len_;
            len[dim]++;
            resize(len);
            back() = x;
        }

        void pop_back(unsigned dim)
        {
            assert(dim < ndim);
            assert(len_[dim] > 0);

            std::array<idx_type, ndim> len = len_;
            len[dim]--;
            resize(len);
        }

        void swap(marray&& other)
        {
            swap(other);
        }

        void swap(marray& other)
        {
            using std::swap;
            swap(this->mem_.data_, other.mem_.data_);
            swap(this->mem_.size_, other.mem_.size_);
            swap(this->len_,       other.len_);
            swap(this->stride_,    other.stride_);
            swap(this->layout_,    other.layout_);
        }

        friend void swap(marray& a, marray& b)
        {
            a.swap(b);
        }
};

template <typename T, int N>
void copy(const marray_view<T,N>& a, marray_view<T,N>& b)
{
    assert(a.lengths() == b.lengths());

    auto it = make_iterator(a.lengths(), a.strides(), b.strides());
    auto a_ = a.data();
    auto b_ = b.data();
    while (it.nextIteration(a_, b_)) *b_ = *a_;
}

template <typename T, int N>
void copy(const marray_view<T,N>& a, marray_view<T,N>&& b)
{
    copy(a, b);
}

template <typename T>
matrix_view<T> operator^(matrix_view<T>& m, transpose_t t)
{
    return matrix_view<T>(m.length(1), m.length(0), m.data(),
                          m.stride(1), m.stride(0));
}

template <typename T>
matrix_view<T> operator^(matrix_view<T>&& m, transpose_t t)
{
    return m^transpose::T;
}

template <typename T>
const_matrix_view<T> operator^(const matrix_view<T>& m, transpose_t t)
{
    return const_cast<matrix_view<T>&>(m)^transpose::T;
}

template <typename T>
void gemm(const T& alpha, const matrix_view<T>& a,
                          const matrix_view<T>& b,
          const T&  beta,       matrix_view<T>& c)
{
    //TODO
}

template <typename T>
void gemm(const T& alpha, const matrix_view<T>&  a,
                          const matrix_view<T>&  b,
          const T&  beta,       matrix_view<T>&& c)
{
    gemm(alpha, a, b, beta, c);
}

template <typename U>
void gemm(char transa, char transb,
          const U& alpha, const matrix_view<U>& a,
                          const matrix_view<U>& b,
          const U&  beta,       matrix_view<U>& c)
{
    using transpose::T;

    if (toupper(transa) == 'T')
    {
        if (toupper(transb) == 'T') gemm(alpha, a^T, b^T, beta, c);
        else                        gemm(alpha, a^T, b  , beta, c);
    }
    else
    {
        if (toupper(transb) == 'T') gemm(alpha, a  , b^T, beta, c);
        else                        gemm(alpha, a  , b  , beta, c);
    }
}

template <typename T>
void gemm(char transa, char transb,
          const T& alpha, const matrix_view<T>&  a,
                          const matrix_view<T>&  b,
          const T&  beta,       matrix_view<T>&& c)
{
    gemm(transa, transb, alpha, a, b, beta, c);
}
}

#endif
