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
 
#include "fusion_interpolation.hpp"
 
#ifdef USE_GOOGLE_BENCHMARK
#include "benchmark/benchmark.h"
 
#include <random>
 
using std::random_device;
using std::mt19937;
using std::uniform_real_distribution;
 
static void BM_RegularCellInterpolation(benchmark::State& state)
#else
#include <iostream>
 
using std::cout;
using std::endl;
 
int main()
#endif
{
    using namespace fusion_interpolation;
   
    const uint64_t N = 3;
   
    limits_type<N> limits = { { { {0.2, 1.0} } // x [min, max]
                              , { {0.4, 2.0} } // y [min, max]
                              , { {0.8, 3.0} } // z [min, max]
                              } };
   
    const value_array_type<N> values = { { 0 // [x0, y0, z0]
                                         , 1 // [x0, y0, z1]
                                         , 2 // [x0, y1, z0]
                                         , 3 // [x0, y1, z1]
                                         , 4 // [x1, y0, z0]
                                         , 5 // [x1, y0, z1]
                                         , 6 // [x1, y1, z0]
                                         , 7 // [x1, y1, z1]
                                         } };
       
    regular_cell<N> cell(limits, values);
   
#ifdef USE_GOOGLE_BENCHMARK
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 1);
    const point_type<N> point = { { dis(gen)
                                  , dis(gen)
                                  , dis(gen)
                                  } };
    while (state.KeepRunning())
    {
        cell.interpolate(point);
    }
#else
        const point_type<N> point = { { 0.5
                                      , 0.5
                                      , 0.5
                                      } };
        cout << cell.interpolate(point) << endl;
#endif
}
 
#ifdef USE_GOOGLE_BENCHMARK
// Register the function as a benchmark
BENCHMARK(BM_RegularCellInterpolation);
 
BENCHMARK_MAIN();
#endif