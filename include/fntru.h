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

///////////////////////////////////////////////////////////////////
// 	FNTRU class. Contains functions for basic operations:
// 	Encrypt, Decrypt, Flatten
// 	It also contains logic functions:
//	XOR, OR, AND
///////////////////////////////////////////////////////////////////

#ifndef FNTRU_H_
#define FNTRU_H_

#include "ssntru.h"
#include <NTL/BasicThreadPool.h>

struct mat_ZZX{
	vec_ZZX *row;
	int len;
};

struct fntru_cipher{
	vec_ZZX row;
	int len;
};

class fntru{
public:
	fntru(int cyclotM_, int bitsizeMod_, int radixSize_,  ZZ messagePrime_);
	void Flatten(fntru_cipher &c);			//Flattening function
	void Encrypt(fntru_cipher &c, ZZX m);
	void Decrypt(ZZX &m, fntru_cipher &c);
	void InitializeCipher(fntru_cipher &c);	// Initialize the ciphertext before any usage
	void Mult(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel);

	//Logic Functions
	void AND(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel);
	void XOR(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel);
	void OR(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel);
	void NOT(fntru_cipher &c, fntru_cipher &a, int sel);
	void NAND(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel);

	//TODO: Destructor
	//~fntru();
	ssntru 	*n;
private:
	void Flatten(mat_ZZX &c, int &degree, int &l);
	void Flatten(vec_ZZX &BD, int &degree, int &l);
	void Flatten(fntru_cipher &c, int &degree, int &l);
	void BitDecomp(vec_ZZX &BD, ZZX &c, int N, int bitlen);
	void BitDecompInv(ZZX &out, vec_ZZX &BD, int bitlen);

	int 	radixSize;	// radix-bitsize
	int		matrixDim; // Dimension of the matrix

	ZZ		*PowTwo;
	mat_ZZX mulMat;
};

void PolyReduce(fntru_cipher &c, int &vectorLength, ssntru &n);
void CoeffReduce(fntru_cipher &c, int &vectorLength, ssntru &n);
void MatrixVectorMul(vec_ZZX &c, mat_ZZX &a, vec_ZZX &b, int &vectorLength, ssntru &n);
void MatrixMul(mat_ZZX &c, mat_ZZX &a, mat_ZZX &b, int &vectorLength, ssntru &n);
void VectorDotProd(ZZX &c, vec_ZZX &a, vec_ZZX &b, int &vectorLength, ssntru &n);
#endif /* FNTRU_H_ */
