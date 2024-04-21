<h1 align='center' style="text-align:center; font-weight:bold; font-size:2.0em;letter-spacing:2.0px;"> FSM: Correspondenceless scan-matching of panoramic 2D range scans </h1>

[![ieeexplore.ieee.org](https://img.shields.io/badge/IEEE/RSJ_IROS_2022_paper-00629B)](https://ieeexplore.ieee.org/abstract/document/9981228)
[![youtube.com](https://img.shields.io/badge/1'_presentation-YouTube-FF0000)](https://www.youtube.com/watch?v=hB4qsHCEXGI)
[![github.com](https://img.shields.io/badge/pdf_presentation-333333)](https://github.com/phd-li9i/fsm_presentation_iros22/blob/master/main.pdf)

This repository provides the source code of FSM and testing code.
FSM's code is located in a single header file for ease of use (`include/fsm.h`).

For the ROS 1 package of FSM visit the `fsm_lidom_ros` package [github page](https://github.com/li9i/fsm_lidom_ros).

## Why use FSM for lidar odometry from 2d panoramic scans

![Experimental results at a glance](https://i.imgur.com/GvFlHgF.png)


## Dependencies
`CGAL 4.7`
`FFTW3`
`boost/random`

## Building

As always
```console
mkdir build
cd build
cmake ..
make
```

Run tests located in the `dataset` directory with


```console
./sm_node A B C D E F G H 0 0 K L M N O P Q R S
```

e.g.

```console
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
- `K`: The size of scans
- `L`: The size of the map. Relevant for scan--to--map-scan-matching tests; irrelevant in scan-matching
- `M`: Identifier of the method used. Only FSM supported in this repository
- `N`: The maximum number of iterations over one sampling degree
- `O`: The minimum sampling degree
- `P`: The maximum sampling degree. Larger values result in extra accuracy and extra execution time
- `Q`: The maximum number of recoveries. Larger values result in extra accuracy and extra execution time
- `R`: Enforcement of a terminal constraint; attempts a recovery if not fulfilled, if true. Set to false for scan-matching
- `S`: Enforcement of an early gear-up feature, see code for details. Set to false for scan-matching; true for scan--to--map-scan matching

## System at a glance

![System at a glance](https://i.imgur.com/PHIThGJ.png)
