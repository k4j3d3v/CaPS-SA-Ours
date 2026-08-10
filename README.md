# CaPS-SA

CaPS-SA is a simple, parallel, and cache-friendly Suffix Array and LCP array construction algorithm.

## Installation

From source:

```bash
git clone https://github.com/k4j3d3v/CaPS-SA-Ours
cd CaPS-SA-mod/
mkdir build && cd build/
cmake -DCMAKE_INSTALL_PREFIX=../ ..
make install VERBOSE=1
cd ..
```

This installs `caps_sa` in a sub-directory named `bin`, inside the project root directory.
To test:

```bash
bin/caps_sa --help
```
should produce 
```
CaPS-SA driver 


bin/caps_sa [OPTIONS] input output


POSITIONALS:
  input TEXT REQUIRED         input path 
  output TEXT REQUIRED        output path 

OPTIONS:
  -h,     --help              Print this help message and exit 
          --data-type TEXT    type of input data [text: "t", genomic: "g", or integer: "i"] 
          --symbol-width TEXT:{32,64} 
                              Symbol width for integer inputs (32 or 64) 
          --ext-mem           pass this flag to use external memor construction 
          --output-lcp        pass this flag to output the LCP array along with the SA 
          --collate-extmem-result Needs: --ext-mem 
                              collate the external memory buckets into a single file 
          --subproblem-count UINT 
                              subproblem count to use 
          --bounded-context UINT 
                              bounded context to use (default: unlimited) 
```

## Usage

Set the number of threads to be used:
```bash
export PARLAY_NUM_THREADS=<thread-count>
```

To compute the SA and LCP array of a uint8 sequence:
```bash
bin/caps_sa  README.md README.md.sa --output-lcp 
```
this produces the file `README.md.sa` and `README.md.sa.lcp`

For uint8 input files smaller than $2**32$ the SA and LCP array are 4x the size of the input, for files larger than $2**32$  SA and LCP array are 8x the size of the input


### From the old usage instruction


Note that by default the subproblem count is set to 8000. If `caps_sa` is run on small datasets it may produce a segmentation fault if a given subproblem is of size 0. Future releases will dynamically set subproblem count, but as of the version 1 release, please use a small subproblem count for datasets significantly smaller than the human genome.
