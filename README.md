# Agricultural Disease in Continuous Time

This is a simulation of the spread of disease among
farms, done with a stochastic, continuous-time framework.

## Installation

Prerequisites

- [Semimarkov Library](https://github.com/afidd/Semi-Markov) The library is header-only. Note the installation location of the headers for this package.

- Gnu Scientific Library (GSL)

- [HDF5](http://www.hdfgroup.org/HDF5/)

- [Boost Library](http://www.boost.org/) version 1.55 or later.

- A C++ compiler that uses -std=c++11, so a recent g++ or clang++.


This package uses GNU Autotools for installation, which
means the two steps are:

```bash
./configure
make
```

The typical options for `configure` apply and are described in the
INSTALL file, but two options will be most important:

```bash
./configure --with-boost=/home/ajd27/Documents/boost_1_57_0 \
            --with-semimarkov=/usr/local/include/semimarkov-0.1
```

These tell configure where to find boost libraries and
the Semi-Markov library.
