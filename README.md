# The Fourier Scan Matching (FSM) implementation in C++

FSM's code is located in a single header file for ease of use (`include/fsm.h`).
This repository provides testing code for FSM and PLICP+GPM (`src/sm{_node}.cpp`.

## Dependencies
`CGAL 4.7`
`FFTW3`
'boost/random'

## Building

As always
```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Run tests located in the `dataset` directory with


```
./sm_node A B C D E F G H I J K L M N O P Q R S
```

e.g.

```
./sm_node 2 1 0 778 0.2 0.786 0.1 0.0 0.0 0.0 360 360 FSM 200 0 3 10 0 0
```

where

- `A`: The number of iterations for the translational component (the larger the
       location displacement between scans the higher this value needs to be)
- `B`: How many times to iterate over all instances of `dataset`
- `C`: The start sample id; typically 0
- `D`: The end sample id; here `|dataset| = 778`
- `E`: The maximum location displacement. Sample locations are generated uniformly in [-E,+E]
- `F`: The maximum orientation displacement. Sample orientations are generated uniformly in [-F,+F]
- `G`: The standard deviation of measurement noise corrupting (real) scans
- `H`: The standard deviation of measurement noise corrupting virtual scans. Suitable in the scan--to--map-scan matching context; irrelevant in scan-matching
- `I`: Proportion of randomly invalidated rays; in [0,1]
- `J`: Proportion of sequentially invalidated rays; in [0,1]
- `K`: The size of scans
- `L`: The size of the map. Relevant for scan--to--map-scan-matching tests; irrelevant in scan-matching
- `M`: Identifier of the method used. Legitimate values are FSM, DBH, KU, CSM (PLICP + GPM)
- `N`: The maximum number of iterations over one sampling degree
- `O`: The minimum sampling degree
- `P`: The maximum sampling degree. Larger values result in extra accuracy and extra execution time
- `Q`: The maximum number of recoveries. Larger values result in extra accuracy and extra execution time
- `R`: Enforcement of a terminal constraint; attempts a recovery if not fulfilled, if true. Set to false for scan-matching
- `S`: Enforcement of an early gear-up feature, see code for details. Set to false for scan-matching; true for scan--to--map-scan matching
