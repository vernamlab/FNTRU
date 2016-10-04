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

#include "../include/fntru.h"
#include <math.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// Compute sigma it is larger than sqrt(q)
long sigma(long n, long logq, long e){
	long s;

	RR n_ 		= to_RR(n);
	RR q_		= to_RR("2");	q_ = power(q_, logq);
	RR e_ 		= power(to_RR("2"), -e);//to_RR(e);
	RR s_;

	s_ = 2*n_*sqrt(log(8*n_*q_)/log(2.0));
	s_ = log(s_)/log(2.0);
	s_ = s_ + to_RR((0.5-2/e_)*logq);
	return s = to_long(s_);
}
// Create fntru
fntru::fntru(int cyclotM_, int bitsizeMod_, int radixSize_,  ZZ messagePrime_){

	n =  new ssntru(cyclotM_, 1, bitsizeMod_, 1, messagePrime_, 1);
	int sig = sigma(*n->Pointer_N(), bitsizeMod_, 80);	//Compute sigma and set

	n->SetSigma(sig);	// Set Sigma
	n->Func_SetModulus();
	n->Func_ModulusFind(1, bitsizeMod_);
	n->Func_ComputeKeys();

	radixSize	= radixSize_;	// Set radix size
	matrixDim 	= (bitsizeMod_/radixSize) + ((bitsizeMod_%radixSize == 0)? 0 : 1); // Compute matrix size = total bitsize / wordsize
/////////////////////////////////////////////////
	//Compute powers of radix
	PowTwo 		= new ZZ[matrixDim];
	ZZ c = to_ZZ("1");
	c = c << radixSize;
	PowTwo[0] = 1;
	for(int i=1; i<matrixDim; i++)
		PowTwo[i] = PowTwo[i-1]*c;
/////////////////////////////////////////////////
	//Set arrays in the matrix
	mulMat.row = new vec_ZZX[matrixDim];
	for(int i=0; i<matrixDim; i++)
		mulMat.row[i].SetLength(matrixDim);
}

//Initialize matrix
void fntru::InitializeCipher(fntru_cipher &m){
	m.row.SetLength(matrixDim);
	m.len = matrixDim;
}
//Compute the bit decomposition oof the given polynomial
void fntru::BitDecomp(vec_ZZX &BD, ZZX &c, int N, int bitlen){
	ZZX r;
	ZZ t;
	int b[radixSize], t_;
	for(int i=0; i<bitlen; i++){
		for(int j=0; j<N; j++){
			t_ = 0;
			for(int k=0; k<radixSize; k++)
				b[k] = bit(coeff(c, j), radixSize*i+k);

			for(int k=0; k<radixSize; k++)
				t_  = t_ + (b[k] << k);

			t = to_ZZ(t_);
			SetCoeff(r, j, t);
		}
		BD[i] = r;
	}
}
///////////////////////////////////////////////
// Create the polynomial by using the bit decomposition vector
void fntru::BitDecompInv(ZZX &out, vec_ZZX &BD, int bitlen){
	ZZX t;
	t = 0;
	for(int i=0; i<bitlen; i++){
		t = t + BD[i]*PowTwo[i];
	}
	out = t;
}

//////////////// FLATTEN ////////////////////////
//Take the ciphertext and apply flattening
void fntru::Flatten(fntru_cipher &c){
	Flatten(c, *n->Pointer_N(), matrixDim);
}

void fntru::Flatten(fntru_cipher &c, int &degree, int &l){
		Flatten(c.row, degree, l);
}
// Apply flattening for the vector
// First apply inverse bit decomposition and later apply bit decomposition
void fntru::Flatten(vec_ZZX &BD, int &degree, int &l){
	ZZX t;

	BitDecompInv(t, BD, l);
	for(int i=0; i<degree; i++)
		SetCoeff(t, i, coeff(t, i) %(*n->Pointer_Q(0)));

	BitDecomp(BD, t, degree, l);
}
// Apply flattening to the matrix (row by row)
void fntru::Flatten(mat_ZZX &c, int &degree, int &l){
	for(int i=0; i<l; i++){
		Flatten(c.row[i], degree, l);
	}
}
//////////////// FLATTEN ////////////////////////


void fntru::Encrypt(fntru_cipher &c, ZZX m){
	ZZX zero;
	ZZX E[matrixDim];
	mat_ZZX tempCipher;

	// Set Temp matrix
	tempCipher.row = new vec_ZZX[matrixDim];
	for(int i=0; i<matrixDim; i++)
		tempCipher.row[i].SetLength(matrixDim);

	// Encryption vector
	zero = 0;
	for(int i=0; i<matrixDim; i++)
		E[i] = n->Prim_Encrypt(zero, 0);

	// Binary Decomp of encryption vector into the temp matrix
	for(int i=0;i<matrixDim; i++)
		BitDecomp(tempCipher.row[i], E[i], *n->Pointer_N(), matrixDim);

	// add message
	for(int i=0; i<matrixDim; i++)
		tempCipher.row[i][i] = tempCipher.row[i][i]+m;

	// flatten
	Flatten(tempCipher, *n->Pointer_N(), matrixDim);

	// convert matrix into the encryption vector
	for(int i=0;i<matrixDim; i++)
		BitDecompInv(c.row[i], tempCipher.row[i], matrixDim);

	// free memory
	for(int i=0; i<matrixDim; i++)
		tempCipher.row[i].kill();
}

void fntru::Decrypt(ZZX &m, fntru_cipher &c){
	m = n->Prim_Decrypt(c.row[0], 0);
}

// if sel == 1: do l*l matrix multiplication by matrix*vector
// if sel == 0: do 1*l matrix multiplication by vector*vector
// When we are computing multiplication, we optimize it by performing
// vector*matrix multiplication. This way we do not need the flattening, since
// we do not do matrix*matrix multiplication.
//
// The function has 2 options. When we do vector=matrix*vector, it results in
// a ciphertext vector. All the elements (ntru ciphers) inside fntru_cipher is computed.
// In case of vector=vector*vector, we only compute the first row (first ntru cipher).
// This way we reduce the computation by a factor of matrixDim (number of ciphertexts inside fntru_cipher)
// If the second option is used it destroys the structure of fntru_cipher. Only first element
// protects its ciphertext form. Therefore, the user should use that as the output while evaluating all the circuit.
void fntru::Mult(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel){
	if(sel == 1){
		//Bit decompose all the rows of the operand to turn it into a matrix
		for(int i=0; i<matrixDim; i++)
			BitDecomp( mulMat.row[i], a.row[i], *n->Pointer_N(), matrixDim);
		// compute matrix multiplication
		MatrixVectorMul(c.row, mulMat, b.row, matrixDim, *n);

		PolyReduce(c, matrixDim, *n);
		CoeffReduce(c, matrixDim, *n);
		// clear temp matrix
		for(int i=0; i<matrixDim; i++)
			for(int j=0; j<matrixDim; j++)
				mulMat.row[i][j]=0;
	}
	else {
		// Bit decompose only first ciphertext
		BitDecomp( mulMat.row[0], a.row[0], *n->Pointer_N(), matrixDim);
		VectorDotProd(c.row[0], mulMat.row[0], b.row, matrixDim, *n);

		// !!!TODO: Add a function for row
		n->Arith_PolyReduce(c.row[0], c.row[0]);
		for(int j=0; j< (*n->Pointer_N()); j++)
			SetCoeff(c.row[0], j, coeff(c.row[0], j) % (*n->Pointer_Q(0)));
		// Clear Matrix
		for(int j=0; j<matrixDim; j++)
			mulMat.row[0][j]=0;
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//Compute matrix multiplication
void MatrixMul(mat_ZZX &c, mat_ZZX &a, mat_ZZX &b, int &vectorLength, ssntru &n){
	ZZX t;
	for(int i=0; i<vectorLength; i++){
		for(int j=0; j<vectorLength; j++){
			t = 0;
			for(int k=0; k<vectorLength; k++)
				t = t + a.row[i][k]*b.row[k][j];

			n.Arith_PolyReduce(t,t);
			c.row[i][j] = t;
		}
	}
}

// Compute the dot product of two vectors. In other words, we compute multiplication
// between the first row of the matrix and the input ciphertext vector.
// If thread is defined it uses NTL thread routine to compute the multiplications.
void VectorDotProd(ZZX &c, vec_ZZX &a, vec_ZZX &b, int &vectorLength, ssntru &n){
#ifdef Thread
	ZZX tt[vectorLength];
	PartitionInfo pinfo(vectorLength);
	long cnt = pinfo.NumIntervals();
	NTL_EXEC_INDEX(cnt, index)
			long first, last;
			pinfo.interval(first, last, index);
			for (int k = first; k<last; k++) {
				mul(tt[k], a[k], b[k]);
			}
	NTL_EXEC_INDEX_END

	c = 0;
	for (int k = 0; k<vectorLength; k++)
		c = c + tt[k];
#endif

#ifndef Thread
	c = 0;
	ZZX tt;
	for (int k = 0; k<vectorLength; k++){
		mul(tt, a[k], b[k]);
//		mul(tt, a[k], b[k]);
		c=c+tt;
	}
#endif
}


//Matrix vector multiplication. TODO: Check 3 options for which one is more optimized for speed and area.
void MatrixVectorMul(vec_ZZX &c, mat_ZZX &a, vec_ZZX &b, int &vectorLength, ssntru &n){

#ifdef MULT1
	ZZX acc, tmp;

	ZZX t, t2;
	ZZX tt[vectorLength];

	for (int i = 0; i<vectorLength; i++) {
		//clear(acc);
		acc = 0;

		PartitionInfo pinfo(vectorLength);
		long cnt = pinfo.NumIntervals();

		NTL_EXEC_INDEX(cnt, index)
			long first, last;
			pinfo.interval(first, last, index);
			for (int k = first; k<last; k++) {
				mul(tt[k], a.row[i][k], b[k]);
			}
		NTL_EXEC_INDEX_END

		acc = 0;
		for(int j=0; j<vectorLength; j++)
			acc = acc + tt[j];
		c[i] = acc;
	}
#endif
#ifdef MULT2
	ZZX acc, tmp;

	ZZX t, t2;
	ZZX tt[vectorLength];

	myTimerReal tr;
	double multtime = 0, addtime = 0;
	for (int i = 0; i<vectorLength; i++) {
		//clear(acc);
		acc = 0;
		start(tr);
		PartitionInfo pinfo(vectorLength);
		long cnt = pinfo.NumIntervals();

		NTL_EXEC_INDEX(cnt, index)
			long first, last;
			pinfo.interval(first, last, index);
			for (int k = first; k<last; k++) {
				mul(tt[k], a.row[i][k], b[k]);
			}
		NTL_EXEC_INDEX_END
		stop(tr);
		multtime += getseconds(tr);
		acc = 0;
		start(tr);
		for(int j=0; j<vectorLength; j++)
			acc = acc + tt[j];
		c[i] = acc;
		stop(tr);
		addtime += getseconds(tr);
	}

#endif
#ifdef MULT3
	ZZX tmp;

	ZZX t, t2;
	ZZX acc[vectorLength], tt[vectorLength];

	for(int j=0; j<vectorLength; j++){
		acc[j] = 0;
		c[j] = 0;
	}

	for (int i = 0; i<vectorLength; i++) {

		PartitionInfo pinfo(vectorLength);
		long cnt = pinfo.NumIntervals();

		NTL_EXEC_INDEX(cnt, index)
			long first, last;
			pinfo.interval(first, last, index);
			for (int k = first; k<last; k++) {
				mul(tt[k], a.row[k][i], b[i]);
			}
		NTL_EXEC_INDEX_END


		for(int j=0; j<vectorLength; j++)
			c[j] = c[j] + tt[j];

//		for(int j=0; j<vectorLength; j++)
//			c[j] = c[j] + acc[j];
	}
#endif
}
// Coefficient modulus reduction for all ciphertexts in fntru ciphertext
void CoeffReduce(fntru_cipher &c, int &vectorLength, ssntru &n){
	for(int i=0; i<vectorLength; i++)
		for(int j=0; j< (*n.Pointer_N()); j++)
			SetCoeff(c.row[i], j, coeff(c.row[i], j) % (*n.Pointer_Q(0)));
}
// Polynomial modulus reduction for all ciphertexts in fntru ciphertext
void PolyReduce(fntru_cipher &c, int &vectorLength, ssntru &n){
	for(int i=0; i<vectorLength; i++)
		n.Arith_PolyReduce(c.row[i], c.row[i]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

void fntru::AND(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel){
	Mult(c, a, b, sel);
}

void fntru::XOR(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel){

	if(sel == 1){
		ZZX t[a.len];
		for(int i=0; i<a.len; i++)
			t[i] = a.row[i] + b.row[i];

		Mult(c, a, b, 1);
		for(int i=0; i<a.len; i++)
			c.row[i] = t[i] - 2*c.row[i];

		PolyReduce(c, c.len, *n);
		CoeffReduce(c, c.len, *n);
	}
	else if(sel == 0){
		ZZX t;
		int vecSize = 1;
		t = a.row[0] + b.row[0];
		Mult(c, a, b, 0);
		c.row[0] = t - 2*c.row[0];
		PolyReduce(c, vecSize, *n);
		CoeffReduce(c, vecSize, *n);
	}
}

void fntru::OR(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel){

	if(sel == 1){
		ZZX t[a.len];
		for(int i=0; i<a.len; i++)
			t[i] = a.row[i] + b.row[i];

		Mult(c, a, b, 1);
		for(int i=0; i<a.len; i++)
			c.row[i] = t[i] - c.row[i];

		PolyReduce(c, c.len, *n);
		CoeffReduce(c, c.len, *n);
	}
	else if(sel == 0){
		ZZX t;
		int vecSize = 1;
		t = a.row[0] + b.row[0];
		Mult(c, a, b, 0);
		c.row[0] = t - c.row[0];
		PolyReduce(c, vecSize, *n);
		CoeffReduce(c, vecSize, *n);
	}
}

void fntru::NOT(fntru_cipher &c, fntru_cipher &a, int sel){

	if(sel == 1){
		ZZX two;
		two = to_ZZ("1");
		for(int i=0; i<a.len; i++){
			c.row[i] = two - a.row[i];
			two = two << 1;
		}
		CoeffReduce(c, c.len, *n);
	}
	else if(sel == 0){
		int vecSize = 1;
		c.row[0] = to_ZZ("1") - a.row[0];
		CoeffReduce(c, vecSize, *n);
	}
}


void fntru::NAND(fntru_cipher &c, fntru_cipher &a, fntru_cipher &b, int sel){

	Mult(c, a, b, sel);
	if(sel == 1){
		ZZX two;
		two = to_ZZ("1");
		for(int i=0; i<a.len; i++){
			c.row[i] = two - c.row[i];
			two = two << 1;
		}
		CoeffReduce(c, c.len, *n);
	}
	else if(sel == 0){
		int vecSize = 1;
		c.row[0] = to_ZZ("1") - c.row[0];
		CoeffReduce(c, vecSize, *n);
	}
}

































