/*
MIT License

Copyright (c) 2016 Vernam Group

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


//////////////////////////////////////////////////////////////////////////////
//	This is a strip down version of the earlier LTV code (DHS-LTV). The		//
// 	original code has class name ltv and source file ltv.cpp. Since we need	//
//	the basics for Stehle and Steinfelds NTRU, it was easier to strip down	//
//	it from the previously constructed LTV. Although, we get rid of			//
//	the unnecessary functions, the code might include some left overs.		//
//////////////////////////////////////////////////////////////////////////////

#ifndef SSNTRU_H_
#define SSNTRU_H_

#include <ctime>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>

#include "def.h"
#include "general.h"

using namespace std;
NTL_CLIENT

class ssntru{
public:

///////////////////////////////////////////////////////////////////////////////
	ssntru(int dm, int num, int max_bit, int tableS, ZZ qq,  int wordSize);
////// Primitive Functins
		ZZX 	Prim_Encrypt(ZZX m, int index);									//Encrypt
		ZZX		Prim_Decrypt(ZZX c, int index);									//Decrypt
////// Some Useful Functions
		ZZX 	Func_Sample();													// Create Sample poly with B-bound
		ZZX 	Func_Sample(long logq);											// Create Sample poly with given log(q)-bound
		void 	Func_SetModulus();												// Set Modulus
		//TODO: Change GSW
		void    Func_ModulusFind(int num, int max_bit);					// Create modulus q
		void 	Func_ComputeKeys();											// Compute Secret/Public Keys
//////	Arithmetic Operations
		void	Arith_PolyReduce(ZZX &out, ZZX &in);							// Polynomial Reduction
		void	Arith_CoeffReduce(ZZX &out, ZZX &in, int index);				// Coeff Modular Reduction
		ZZX 	Arith_PolyCoeffReduce(ZZX &in, int index);						// Poly and Coeff Reduction
		ZZX		Arith_MulModZZX(ZZX &a, ZZX &b, int index);						// Polynomial Modular Mult.
		ZZX		Arith_AddModZZX(ZZX &a, ZZX &b, int index);						// Polynomial Modular Add
///////	Pointers
		ZZ		*Pointer_Q(int i);												// Pointer to modulus q
		ZZX		*Pointer_PolyModulus();											// Pointer to polynomial modulus
		ZZX		*Pointer_PublicKey();											// Pointer to Public Key
		ZZX		*Pointer_SecretKey(int i);										// Pointer to Secret Key
		int 	*Pointer_N();													// Pointer to polynomial degree
		GlobalParam *Pointer_GP();
		void	SetSigma(long sig);												// Set Sigma for noise
private:
	int word;
	int N, degree_m;  		// poly degree, euler phi_m
	ZZ 	p, q, B;			// message modulus, encryption modulus, noise bound
	ZZX modulus;			// polynomial modulus
	long sigma;

	ZZ  	*q_list;		// list of q modulus for each level
	ZZX 	*pk, *sk;		// list of public and secret keys for each level

	GlobalParam gp;			// Global param struct to store some parameters
	myReduction myr;		// class to complete fast polynomial reduction
	bool reducsetflag;		//
	int rand_seed, message_rand, max_bitsize, init_bitsize, word_blocksize; //

/////// Private Functions
	ZZX 	ComputeFastCycModulus(int n);	//Compute Cyclotomic Polynomial
	void 	SetFastPolyModulus();
	int 	NumofWords(int bs);
};
///////////////////////////////////////////////////////////////////////////////
////// Functions used by ssntru member functions
	void 	clear(ZZX &x, int size);
	void 	find_inverse(ZZX &f_inv, ZZX &f, ZZ &q, int &degree, bool &isfound, ZZX modulus);
	void 	coeff_reduction(ZZX &out, ZZX &in, ZZ &q, int &degree);
	void 	coeff_reduction_q_2(ZZX &out, ZZX &in, ZZ &q, int N);
	void 	getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2);
	int 	MobuisFunction(int n);
	RR		randn (long mu, long sigmaLog, long qLog);
	double 	randn (double mu, double sigma);
///////////////////////////////////////////////////////////////////////////////

#endif /* SSNTRU_H_ */



