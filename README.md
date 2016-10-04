# FNTRU Encryption Library
The library provides a homomorphic encryption scheme by Doröz and Sunar. Our scheme adopts the flattening technique proposed in GSW to the Stehlé and Steinfeld NTRU variant. Similar to GSW, the scheme does not require evaluation keys or key switching. Also, our scheme uses wide key distributions, and hence is immune to the Subfield Lattice Attack. In case of timings, we are compatitive compared to the existing schemes by achieving 5.8 msec for homomorphic multiplications (5 levels).    

The details of the scheme can be seen in : https://eprint.iacr.org/2016/315.pdf.

# Thread
User can set the number of thread that they want to use. The parallelism is achieved by using NTL thread settings. In order to set the number of threads, you can use the following function:

SetNumThreads(thread);

# Installation
This library uses NTL (v9.9.1 or later) with GMP support (with C++11 on). Tested on Linux environment. You can use the makefile to compile the project with he following command:

make
