#ifndef _MARRAY_ITERATOR_HPP_
#define _MARRAY_ITERATOR_HPP_

#ifndef MARRAY_TEST
#define MARRAY_TEST(...)
#endif

#include <vector>
#include <array>
#include <cstddef>

namespace MArray
{

    template <typename idx_type = int, typename size_type = size_t>
    class Iterator
    {
        public:
            template <size_t ndim>
            Iterator(const std::array<idx_type, ndim>& len,
                     const std::array<size_type, ndim>& stride_0)
            : first(true), pos(ndim), len(len.begin(), len.end()), stride(1)
            {
                stride[0].assign(stride_0.begin(), stride_0.end());
                check();
            }

            template <size_t ndim>
            Iterator(const std::array<idx_type, ndim>& len,
                     const std::array<size_type, ndim>& stride_0,
                     const std::array<size_type, ndim>& stride_1)
            : first(true), pos(ndim), len(len.begin(), len.end()), stride(2)
            {
                stride[0].assign(stride_0.begin(), stride_0.end());
                stride[1].assign(stride_1.begin(), stride_1.end());
                check();
            }

            template <size_t ndim>
            Iterator(const std::array<idx_type, ndim>& len,
                     const std::array<size_type, ndim>& stride_0,
                     const std::array<size_type, ndim>& stride_1,
                     const std::array<size_type, ndim>& stride_2)
            : first(true), pos(ndim), len(len.begin(), len.end()), stride(3)
            {
                stride[0].assign(stride_0.begin(), stride_0.end());
                stride[1].assign(stride_1.begin(), stride_1.end());
                stride[2].assign(stride_2.begin(), stride_2.end());
                check();
            }

            template <typename cv_ptr_0>
            bool nextIteration(cv_ptr_0& ptr_0)
            {
                if (stride.size() != 1) abort();

                if (first)
                {
                    first = false;
                }
                else
                {
                    if (len.size() == 0)
                    {
                        first = true;
                        return false;
                    }

                    for (idx_type i = 0;i < len.size();i++)
                    {
                        if (pos[i] == len[i]-1)
                        {
                            ptr_0 -= pos[i]*stride[0][i];
                            pos[i] = 0;

                            if (i == len.size()-1)
                            {
                                first = true;
                                return false;
                            }
                        }
                        else
                        {
                            ptr_0 += stride[0][i];
                            pos[i]++;

                            return true;
                        }
                    }
                }

                return true;
            }

            template <typename cv_ptr_0, typename cv_ptr_1>
            bool nextIteration(cv_ptr_0& ptr_0, cv_ptr_1& ptr_1)
            {
                if (stride.size() != 2) abort();

                if (first)
                {
                    first = false;
                }
                else
                {
                    if (len.size() == 0)
                    {
                        first = true;
                        return false;
                    }

                    for (idx_type i = 0;i < len.size();i++)
                    {
                        if (pos[i] == len[i]-1)
                        {
                            ptr_0 -= pos[i]*stride[0][i];
                            ptr_1 -= pos[i]*stride[1][i];
                            pos[i] = 0;

                            if (i == len.size()-1)
                            {
                                first = true;
                                return false;
                            }
                        }
                        else
                        {
                            ptr_0 += stride[0][i];
                            ptr_1 += stride[1][i];
                            pos[i]++;

                            return true;
                        }
                    }
                }

                return true;
            }

            template <typename cv_ptr_0, typename cv_ptr_1, typename cv_ptr_2>
            bool nextIteration(cv_ptr_0& ptr_0, cv_ptr_1& ptr_1, cv_ptr_2& ptr_2)
            {
                if (stride.size() != 3) abort();

                if (first)
                {
                    first = false;
                }
                else
                {
                    if (len.size() == 0)
                    {
                        first = true;
                        return false;
                    }

                    for (idx_type i = 0;i < len.size();i++)
                    {
                        if (pos[i] == len[i]-1)
                        {
                            ptr_0 -= pos[i]*stride[0][i];
                            ptr_1 -= pos[i]*stride[1][i];
                            ptr_2 -= pos[i]*stride[2][i];
                            pos[i] = 0;

                            if (i == len.size()-1)
                            {
                                first = true;
                                return false;
                            }
                        }
                        else
                        {
                            ptr_0 += stride[0][i];
                            ptr_1 += stride[1][i];
                            ptr_2 += stride[2][i];
                            pos[i]++;

                            return true;
                        }
                    }
                }

                return true;
            }

        private:
            void check()
            {
                for (idx_type j = 0;j < stride.size();j++)
                {
                    if (stride[j].size() != len.size()) abort();
                }

                for (idx_type i = 0;i < len.size();i++)
                {
                    if (len[i] <= 0) abort();
                }
            }

            bool first;
            std::vector<idx_type> pos;
            std::vector<idx_type> len;
            std::vector<std::vector<size_type>> stride;
    };

}

#endif
