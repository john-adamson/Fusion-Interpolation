// MIT License
//
// Copyright (c) 2016 John Adamson, www.oxforddynamics.com
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
 
#ifndef FUSION_INTERPOLATION_0980078A_D871_4284_8FC6_15D2E4A59BCC
#define FUSION_INTERPOLATION_0980078A_D871_4284_8FC6_15D2E4A59BCC
 
//----------------------------------------------------------------------------
//
// Value Order Convention
// ----------------------
//
// 1D
// ~~
//
//   +----------+
//   x0         x1
//
//   values = [[x0],
//             [x1]
//            ]
//
//
// 2D
// ~~
//
//   x0,y1      x1,y1
//     +----------+
//     |          |
//     |          |
//     |          |
//     |          |
//     +----------+
//   x0,y0      x1,y0
//
//   values = [[x0, y0],
//             [x0, y1],
//             [x1, y0],
//             [x1, y1]
//            ]
//
//
// 3D
// ~~
//
//                x0,y1,z1            x1,y1,z1
//                   +-------------------+
//                  /|                  /|
//                 / |                 / |
//                /  |                /  |
//               /   |               /   |
//              /    |              /    |
//    x0,y0,z1 +-------------------+ x1,y0,z1
//             |     |             |     |
//             |     +-------------|-----+
//             |    / x0,y1,z0     |    / x1,y1,z0
//             |   /               |   /
//             |  /                |  /
//             | /                 | /
// z  y        |/                  |/
// | /         +-------------------+
// |/       x0,y0,z0            x1,y0,z0
// +--x
//
//   values = [[x0, y0, z0],
//             [x0, y0, z1],
//             [x0, y1, z0],
//             [x0, y1, z1],
//             [x1, y0, z0],
//             [x1, y0, z1],
//             [x1, y1, z0],
//             [x1, y1, z1]
//            ]
//
// ----
//
// And so on for higher dimensions.
//
//----------------------------------------------------------------------------
 
#include <array>
 
#include <armadillo> // v6.500+
 
namespace fusion_interpolation
{
 
//----------------------------------------------------------------------------
constexpr uint64_t power_of_two(const uint64_t exponent)
{
    return static_cast<uint64_t>(1) << exponent;
}
 
//----------------------------------------------------------------------------
template<uint64_t N>
using point_type = std::array<double, N>;
 
template<uint64_t N>
using limits_type = std::array<std::array<double, 2>, N>;
 
template<uint64_t N>
using point_array_type = std::array<point_type<N>, power_of_two(N)>;
 
template<uint64_t N>
using value_array_type = std::array<double, power_of_two(N)>;
 
template<uint64_t N>
using row_type = std::array<bool, N>;
 
template<uint64_t N>
using row_list_type = std::array<row_type<N>, power_of_two(N)>;
 
//----------------------------------------------------------------------------
template<uint64_t N>
row_list_type<N> multipliers()
{
    row_list_type<N> rows = {};
    for (uint64_t i(0); i < rows.size(); ++i)
    {
        for (uint64_t j(0); j < rows[i].size(); ++j)
        {
            rows[i][j] = i & power_of_two(j);
        }
    }
    return rows;
}
 
//----------------------------------------------------------------------------
template<uint64_t N>
point_array_type<N> points_from_limits(const limits_type<N> & limits)
{
    point_type<N> p0 = {};
    point_type<N> p1 = {};
    for (uint64_t i(0); i < N; ++i)
    {
        p0[i] = limits[i][0];
        p1[i] = limits[i][1];
    }
 
    point_array_type<N> pa;
    const auto two_to_N = power_of_two(N);
    for (uint64_t i(0); i < two_to_N; ++i)
    {
        for (uint64_t j(0); j < N; ++j)
        {
            pa[i][j] = (i / (two_to_N / power_of_two(j + 1)) % 2 == 0) ? p0[j] : p1[j];
        }
    }
    return pa;
}
//----------------------------------------------------------------------------
template<uint64_t N>
class regular_cell
{
public:
    regular_cell
    (
        const limits_type<N>      & limits
    ,   const value_array_type<N> & values
    )
    {
        const point_array_type<N> pa = points_from_limits(limits);
 
        arma::mat A(two_to_N_, two_to_N_);
        for (uint64_t i(0); i < two_to_N_; ++i)
        {
            for (uint64_t j(0); j < two_to_N_; ++j)
            {
                double v = 1.0;
                for (uint64_t k(0); k < N; ++k)
                {
                    v *= filter_[j][k] ? pa[i][k] : 1.0;
                }
                A(i, j) = v;
            }
        }
#ifdef _DEBUG
        A.print("A: ");
#endif
 
        arma::vec b(two_to_N_);
        for (uint64_t i(0); i < two_to_N_; ++i)
        {
            b(i) = values[i];
        }
#ifdef _DEBUG
        b.print("b: ");
#endif
 
        a_ = arma::solve(A, b, arma::solve_opts::equilibrate + arma::solve_opts::no_approx);
#ifdef _DEBUG
        a_.print("a: ");
#endif
    }
 
    double interpolate(const point_type<N> & p)
    {
        double result = 0.0;
        for (uint64_t i(0); i < two_to_N_; ++i)
        {
            double v = a_[i];
            for (uint64_t j(0); j < N; ++j)
            {
                v *= filter_[i][j] ? p[j] : 1.0;
            }
            result += v;
        }
 
        return result;
    }
 
private:
 
    arma::vec a_;
 
private:
 
    // These variables apply to an entire mesh so don't need to be stored on a
    // per cell basis.  If a mesh class were available they could be moved to
    // that.
    static const uint64_t two_to_N_;
    static const row_list_type<N> filter_;
};
 
//----------------------------------------------------------------------------
template<uint64_t N>
const uint64_t regular_cell<N>::two_to_N_ = power_of_two(N);
 
template<uint64_t N>
const row_list_type<N> regular_cell<N>::filter_ = multipliers<N>();
 
} // namespace fusion_interpolation
 
#endif // FUSION_INTERPOLATION_0980078A_D871_4284_8FC6_15D2E4A59BCC