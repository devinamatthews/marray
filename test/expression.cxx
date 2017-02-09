#include "gtest/gtest.h"
#include "expression.hpp"

using namespace std;
using namespace MArray;

TEST(expression, assign)
{
    using namespace slice;

    double data1[12] = { 1, 2, 3, 4, 5, 6,
                         7, 8, 9,10,11,12};
    double data2[12] = {12,11,10, 9, 8, 7,
                         6, 5, 4, 3, 2, 1};
    double data3[6] = {-1,-1,-1,-1,-1,-1};

    marray<double,3> v1({2, 3, 2});
    marray_view<double,3> v2({2, 3, 2}, data1);
    marray<double,4> v3({4, 3, 2, 2}, 1.0);
    marray<double,3> v4({2, 3, 2}, 4.0);
    marray_view<double,3> v5({2, 3, 2}, data2);
    marray_view<double,2> v6({3, 2}, data3);
    marray<double,2> v7({3, 2}, 5.0);

    v1 = 2.0;
    EXPECT_EQ((array<double,12>{2, 2, 2, 2, 2, 2,
                                2, 2, 2, 2, 2, 2}), *(array<double,12>*)v1.data());

    v1 = 3;
    EXPECT_EQ((array<double,12>{3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3}), *(array<double,12>*)v1.data());

    v1 = v2;
    EXPECT_EQ((array<double,12>{ 1, 2, 3, 4, 5, 6,
                                 7, 8, 9,10,11,12}), *(array<double,12>*)v1.data());

    v1 = v3[range(2)][all][1][all];
    EXPECT_EQ((array<double,12>{ 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1}), *(array<double,12>*)v1.data());

    v1 = v4;
    EXPECT_EQ((array<double,12>{ 4, 4, 4, 4, 4, 4,
                                 4, 4, 4, 4, 4, 4}), *(array<double,12>*)v1.data());

    v2 = 2.0;
    EXPECT_EQ((array<double,12>{2, 2, 2, 2, 2, 2,
                                2, 2, 2, 2, 2, 2}), *(array<double,12>*)v2.data());

    v2 = 3;
    EXPECT_EQ((array<double,12>{3, 3, 3, 3, 3, 3,
                                3, 3, 3, 3, 3, 3}), *(array<double,12>*)v2.data());

    v2 = v5;
    EXPECT_EQ((array<double,12>{12,11,10, 9, 8, 7,
                                 6, 5, 4, 3, 2, 1}), *(array<double,12>*)v2.data());

    v2 = v3[range(2)][all][1][all];
    EXPECT_EQ((array<double,12>{ 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1}), *(array<double,12>*)v2.data());

    v2 = v4;
    EXPECT_EQ((array<double,12>{ 4, 4, 4, 4, 4, 4,
                                 4, 4, 4, 4, 4, 4}), *(array<double,12>*)v2.data());

    v1[1][all][all] = 2.0;
    EXPECT_EQ((array<double,12>{4, 4, 4, 4, 4, 4,
                                2, 2, 2, 2, 2, 2}), *(array<double,12>*)v1.data());

    v1[1][all][all] = 3;
    EXPECT_EQ((array<double,12>{4, 4, 4, 4, 4, 4,
                                3, 3, 3, 3, 3, 3}), *(array<double,12>*)v1.data());

    v1[1][all][all] = v6;
    EXPECT_EQ((array<double,12>{ 4, 4, 4, 4, 4, 4,
                                -1,-1,-1,-1,-1,-1}), *(array<double,12>*)v1.data());

    v1[1][all][all] = v3[1][all][1][all];
    EXPECT_EQ((array<double,12>{ 4, 4, 4, 4, 4, 4,
                                 1, 1, 1, 1, 1, 1}), *(array<double,12>*)v1.data());

    v1[1][all][all] = v7;
    EXPECT_EQ((array<double,12>{ 4, 4, 4, 4, 4, 4,
                                 5, 5, 5, 5, 5, 5}), *(array<double,12>*)v1.data());
}

TEST(expression, bcast)
{
    using namespace slice;

    double data[3] = {1, 2, 3};

    marray<double,3> v1({3, 2, 3});
    marray_view<double,1> v2({3}, data);

    v1 = v2;
    EXPECT_EQ((array<double,18>{1, 2, 3, 1, 2, 3,
                                1, 2, 3, 1, 2, 3,
                                1, 2, 3, 1, 2, 3}), *(array<double,18>*)v1.data());

    v1 = 0;
    EXPECT_EQ((array<double,18>{0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0}), *(array<double,18>*)v1.data());

    v1 = v2[bcast][bcast];
    EXPECT_EQ((array<double,18>{1, 2, 3, 1, 2, 3,
                                1, 2, 3, 1, 2, 3,
                                1, 2, 3, 1, 2, 3}), *(array<double,18>*)v1.data());

    v1 = v2[all][bcast][bcast];
    EXPECT_EQ((array<double,18>{1, 1, 1, 1, 1, 1,
                                2, 2, 2, 2, 2, 2,
                                3, 3, 3, 3, 3, 3}), *(array<double,18>*)v1.data());
}

TEST(expression, add)
{
    double data1[3] = {1, 2, 3};
    double data2[3] = {3, 2, 1};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);
    marray_view<double,1> v3({3}, data2);

    v1 = v2 + v3;
    EXPECT_EQ((array<double,3>{4, 4, 4}), *(array<double,3>*)v1.data());

    v1 = v2 + 1;
    EXPECT_EQ((array<double,3>{2, 3, 4}), *(array<double,3>*)v1.data());

    v1 = 2.0 + v3;
    EXPECT_EQ((array<double,3>{5, 4, 3}), *(array<double,3>*)v1.data());

    v1 += v2;
    EXPECT_EQ((array<double,3>{6, 6, 6}), *(array<double,3>*)v1.data());

    v1 += 1;
    EXPECT_EQ((array<double,3>{7, 7, 7}), *(array<double,3>*)v1.data());
}

TEST(expression, sub)
{
    double data1[3] = {1, 2, 3};
    double data2[3] = {3, 2, 1};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);
    marray_view<double,1> v3({3}, data2);

    v1 = v2 - v3;
    EXPECT_EQ((array<double,3>{-2, 0, 2}), *(array<double,3>*)v1.data());

    v1 = v2 - 1;
    EXPECT_EQ((array<double,3>{0, 1, 2}), *(array<double,3>*)v1.data());

    v1 = 2.0 - v3;
    EXPECT_EQ((array<double,3>{-1, 0, 1}), *(array<double,3>*)v1.data());

    v1 -= v2;
    EXPECT_EQ((array<double,3>{-2, -2, -2}), *(array<double,3>*)v1.data());

    v1 -= 1;
    EXPECT_EQ((array<double,3>{-3, -3, -3}), *(array<double,3>*)v1.data());
}

TEST(expression, mul)
{
    double data1[3] = {1, 2, 3};
    double data2[3] = {3, 2, 1};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);
    marray_view<double,1> v3({3}, data2);

    v1 = v2 * v3;
    EXPECT_EQ((array<double,3>{3, 4, 3}), *(array<double,3>*)v1.data());

    v1 = v2 * 1;
    EXPECT_EQ((array<double,3>{1, 2, 3}), *(array<double,3>*)v1.data());

    v1 = 2.0 * v3;
    EXPECT_EQ((array<double,3>{6, 4, 2}), *(array<double,3>*)v1.data());

    v1 *= v2;
    EXPECT_EQ((array<double,3>{6, 8, 6}), *(array<double,3>*)v1.data());

    v1 *= 2;
    EXPECT_EQ((array<double,3>{12, 16, 12}), *(array<double,3>*)v1.data());
}

TEST(expression, div)
{
    double data1[3] = {1, 2, 3};
    double data2[3] = {3, 2, 1};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);
    marray_view<double,1> v3({3}, data2);

    v1 = v2 / v3;
    EXPECT_EQ((array<double,3>{1.0/3, 1, 3}), *(array<double,3>*)v1.data());

    v1 = v2 / 1;
    EXPECT_EQ((array<double,3>{1, 2, 3}), *(array<double,3>*)v1.data());

    v1 = 2.0 / v3;
    EXPECT_EQ((array<double,3>{2.0/3, 1, 2}), *(array<double,3>*)v1.data());

    v1 /= v2;
    EXPECT_EQ((array<double,3>{2.0/3, 0.5, 2.0/3}), *(array<double,3>*)v1.data());

    v1 /= 2;
    EXPECT_EQ((array<double,3>{1.0/3, 0.25, 1.0/3}), *(array<double,3>*)v1.data());
}

TEST(expression, pow)
{
    double data1[3] = {1, 2, 3};
    double data2[3] = {3, 2, 1};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);
    marray_view<double,1> v3({3}, data2);

    v1 = pow(v2, v3);
    EXPECT_EQ((array<double,3>{1, 4, 3}), *(array<double,3>*)v1.data());

    v1 = pow(v2, 2);
    EXPECT_EQ((array<double,3>{1, 4, 9}), *(array<double,3>*)v1.data());

    v1 = pow(2.0, v3);
    EXPECT_EQ((array<double,3>{8, 4, 2}), *(array<double,3>*)v1.data());
}

TEST(expression, negate)
{
    double data1[3] = {1, 2, 3};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);

    v1 = -v2;
    EXPECT_EQ((array<double,3>{-1, -2, -3}), *(array<double,3>*)v1.data());
}

TEST(expression, exp)
{
    double data1[3] = {1, 2, 3};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);

    v1 = exp(v2);
    EXPECT_EQ((array<double,3>{exp(1), exp(2), exp(3)}), *(array<double,3>*)v1.data());
}

TEST(expression, sqrt)
{
    double data1[3] = {4, 9, 16};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);

    v1 = sqrt(v2);
    EXPECT_EQ((array<double,3>{2, 3, 4}), *(array<double,3>*)v1.data());
}

TEST(expression, compound)
{
    double data1[3] = {1, 2, 3};
    double data2[3] = {3, 2, 1};
    double data3[3] = {4, 7, 2};

    marray<double,1> v1({3});
    marray_view<double,1> v2({3}, data1);
    marray_view<double,1> v3({3}, data2);
    marray_view<double,1> v4({3}, data3);

    v1 = (pow(v2, 2) * v3 + 1)/4 + sqrt(v4);
    EXPECT_DOUBLE_EQ(3, v1[0]);
    EXPECT_DOUBLE_EQ(9.0/4 + sqrt(7), v1[1]);
    EXPECT_DOUBLE_EQ(5.0/2 + sqrt(2), v1[2]);
}

TEST(expression, mixed_rank)
{
    double data1[12] = {1, 2, 3, 4, 5, 6,
                        7, 8, 9,10,11,12};
    double data2[6] = {1, 2, 3, 4, 5, 6};
    double data3[3] = {3, 2, 1};

    marray<double,3> v1({2, 2, 3});
    marray_view<double,3> v2({2, 2, 3}, data1);
    marray_view<double,2> v3({2, 3}, data2);
    marray_view<double,1> v4({3}, data3);

    v1 = v2 * 2 + v3 / v4;
    EXPECT_DOUBLE_EQ( 2 + 1.0/3, v1.data()[ 0]);
    EXPECT_DOUBLE_EQ( 4 + 2.0/2, v1.data()[ 1]);
    EXPECT_DOUBLE_EQ( 6 + 3.0/1, v1.data()[ 2]);
    EXPECT_DOUBLE_EQ( 8 + 4.0/3, v1.data()[ 3]);
    EXPECT_DOUBLE_EQ(10 + 5.0/2, v1.data()[ 4]);
    EXPECT_DOUBLE_EQ(12 + 6.0/1, v1.data()[ 5]);
    EXPECT_DOUBLE_EQ(14 + 1.0/3, v1.data()[ 6]);
    EXPECT_DOUBLE_EQ(16 + 2.0/2, v1.data()[ 7]);
    EXPECT_DOUBLE_EQ(18 + 3.0/1, v1.data()[ 8]);
    EXPECT_DOUBLE_EQ(20 + 4.0/3, v1.data()[ 9]);
    EXPECT_DOUBLE_EQ(22 + 5.0/2, v1.data()[10]);
    EXPECT_DOUBLE_EQ(24 + 6.0/1, v1.data()[11]);
}
