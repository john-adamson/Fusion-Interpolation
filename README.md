Fusion Interpolation
====================
 
Background
----------
 
Fusion Interpolation is a small single header containing a class for
linear interpolation of values within a regular cell aligned with the reference
axes.  The class can be of any dimension (well, any useful dimension anyway!).
 
The code is for a single cell only and not a mesh (but could easily be adapted
to be part of a mesh class).
 
An example is provided for a 3D cell.
 
 
Prerequisites
-------------
 
Fusion Interpolation relies upon the
[Armadillo C++ Linear Algebra Library](http://arma.sourceforge.net)
 
The example uses [Google Benchmark](https://github.com/google/benchmark) if
the user also defines `USE_GOOGLE_BENCHMARK`.
 
 
How it Works
------------
 
The interpolated values for the following dimensions can be calculated as:
 
* 1D: `f(x)     = a0 + a1.x`
* 2D: `f(x,y)   = a0 + a1.x + a2.y + a3.x.y`
* 3D: `f(x,y,z) = a0 + a1.x + a2.y + a3.z + a4.x.y + a5.x.z + a6.y.z + a7.x.y.z`
 
and so on, where `x`, `y` and `z` are the coordinates of the point
requiring an interpolated value.
 
By constructing a matrix of the equations given above for each cell corner the
`a` coefficients can be calculated.
 
On construction, the class `regular_cell` calculates the `a` coefficients
and stores them for future use.
 
When the interpolation function is called it applies the calculation given
above to the supplied point to determine the interpolated result.
 
 
Value Order Convention
----------------------
 
The order of the values must conform to the following pattern.  Dimensions 1 to
3 are shown below.  Higher dimensions follow the same pattern.
 
### 1D
 
    +----------+
    x0         x1
    
    values = [[x0],
              [x1]
             ]
 
### 2D
 
    x0,y1      x1,y1
      +----------+
      |          |
      |          |
      |          |
      |          |
      +----------+
    x0,y0      x1,y0
    
    values = [[x0, y0],
              [x0, y1],
              [x1, y0],
              [x1, y1]
             ]
 
### 3D
 
                    x0,y1,z1            x1,y1,z1
                       +-------------------+
                      /|                  /|
                     / |                 / |
                    /  |                /  |
                   /   |               /   |
                  /    |              /    |
        x0,y0,z1 +-------------------+ x1,y0,z1
                |     |             |     |
                |     +-------------|-----+
                |    / x0,y1,z0     |    / x1,y1,z0
                |   /               |   /
                |  /                |  /
                | /                 | /
    z  y        |/                  |/
    | /         +-------------------+
    |/       x0,y0,z0            x1,y0,z0
    +--x
    
      values = [[x0, y0, z0],
                [x0, y0, z1],
                [x0, y1, z0],
                [x0, y1, z1],
                [x1, y0, z0],
                [x1, y0, z1],
                [x1, y1, z0],
                [x1, y1, z1]
               ]
 
 
### 2D Example
 
    f(0,1)=3   f(1,1)=7
       +----------+
       |          |
       |          |
       |          |
       |          |
       +----------+
    f(0,0)=0   f(1,0)=5
    
    values = [0, 3, 5, 7]
 
 
Usage
-----
 
See `example_3d.cpp` for usage.
 
 
License
-------
 
MIT license.  If this code is used in any commercial and/or academic work please
notify the author.
