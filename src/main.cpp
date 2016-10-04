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
#include <sys/time.h>
#include <NTL/GF2X.h>
#include <NTL/mat_ZZ.h>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//And and XOR operations for test preparation
int AND(int b0, int b1){
	return(b0*b1)%2;
}

int XOR(int b0, int b1){
	return(b0+b1)%2;
}

//Convert message in ZZX to int. If not bit print -1
int ConvMessPoly(ZZX in){
	ZZX zero, one;
	one = 1;
	zero = 0;

	if(in == 0)
		return 0;
	else if(in == 1)
		return 1;
	else
		return -1;
}

//Check if all the values of boolean circuit is equal
bool CheckBooleanCircuitResult(int midResult[], ZZX messmidResult[], int opSize){
	bool isEqual = true;

	for(int i=0; i<opSize+1; i++)
		if(midResult[i] != ConvMessPoly(messmidResult[i])){
			isEqual = false;
			break;
		}

	return isEqual;
}

//Print all the values of the Boolean Circuit (mid results)
void PrintBooleanCircuitValues(int midResult[], ZZX messmidResult[], int opSize){
	for(int i=0; i<opSize+1; i++){
		cout << "Op\t" << i << "\t" << midResult[i] << "\t" << ConvMessPoly(messmidResult[i]) << "\t" << endl;
	}
}

void RandomBooleanTest(){

	int opSize = 16; 	// Set boolean circuit operation size
	int b[opSize+1]; 	// Number of bits
	int op[opSize];		// Operations. 0 is AND, 1 is XOR

	//Create Random Bits and operations
	srand(time(NULL));
	for(int i=0; i<opSize; i++)
			op[i] = rand()%2;
	for(int i=0; i<opSize+1; i++)
			b[i] = rand()%2;

	//Compute the random boolean function and store the mid-results
	int opResult, midResult[opSize+1];
	opResult = b[0];
	for(int i=0; i<opSize; i++){
		midResult[i] = opResult;
		if(op[i] == 0)
			opResult = AND(opResult, b[i+1]);
		else
			opResult = XOR(opResult, b[i+1]);
	}
	//Store the last result
	midResult[opSize] = opResult;


	// Set fntru parameters and create ciphertexts to encrypt the bits of boolean function
	int M 			= M_;
	int radix		= 16;
	int bitSize 	= Dif_Prime;
	ZZ p 			= to_ZZ("2");

	fntru g(M, bitSize, radix,  p);
	ZZX message, messaget, m, r0, r1;
	fntru_cipher gResult, gb[opSize+1];

	// Initialize the ciphertexts
	g.InitializeCipher(gResult);
	for(int i=0; i<opSize+1; i++)
		g.InitializeCipher(gb[i]);

	// Encrypt the boolean bits
	for(int i=0; i<opSize+1; i++){
		message = b[i];
		g.Encrypt(gb[i], message);
	}

	// Compute the boolean function and store the mid-results by decrypting the result
	ZZX messResult, messb[opSize+1], messmidResult[opSize+1];
	gResult = gb[0];
	for(int i=0; i<opSize; i++){

		g.Decrypt(messb[i], gb[i]);
		g.Decrypt(messResult, gResult);

		messmidResult[i] = messResult;
		if(op[i] == 0)
			g.AND(gResult, gResult, gb[i+1], 0);
		else
			g.XOR(gResult, gResult, gb[i+1], 0);
	}
	// Decrypt the last result
	g.Decrypt(messResult, gResult);
	messmidResult[opSize] = messResult;

	cout << "Boolean Equal:\t" << CheckBooleanCircuitResult(midResult, messmidResult, opSize) << endl;
	PrintBooleanCircuitValues(midResult, messmidResult, opSize);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void TestAND(fntru &f){

	fntru_cipher a, b, c;

	f.InitializeCipher(a);
	f.InitializeCipher(b);
	f.InitializeCipher(c);

	ZZX one, zero, result;
	one 	= 1;
	zero 	= 0;

	cout << "\\\\\\\\\\\\ AND \\\\\\\\\\\\" << endl;
	cout << "a\t" << "b\t" << "c\t" << endl << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, zero);
	f.AND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, one);
	f.AND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "1\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, zero);
	f.AND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, one);
	f.AND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "1\t" << ConvMessPoly(result) << endl;

}

void TestNAND(fntru &f){

	fntru_cipher a, b, c;

	f.InitializeCipher(a);
	f.InitializeCipher(b);
	f.InitializeCipher(c);

	ZZX one, zero, result;
	one 	= 1;
	zero 	= 0;

	cout << "\\\\\\\\\\\\ NAND \\\\\\\\\\\\" << endl;
	cout << "a\t" << "b\t" << "c\t" << endl << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, zero);
	f.NAND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, one);
	f.NAND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "1\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, zero);
	f.NAND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, one);
	f.NAND(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "1\t" << ConvMessPoly(result) << endl;

}

void TestOR(fntru &f){

	fntru_cipher a, b, c;

	f.InitializeCipher(a);
	f.InitializeCipher(b);
	f.InitializeCipher(c);

	ZZX one, zero, result;
	one 	= 1;
	zero 	= 0;

	cout << "\\\\\\\\\\\\ OR \\\\\\\\\\\\" << endl;
	cout << "a\t" << "b\t" << "c\t" << endl << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, zero);
	f.OR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, one);
	f.OR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "1\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, zero);
	f.OR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, one);
	f.OR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "1\t" << ConvMessPoly(result) << endl;

}

void TestXOR(fntru &f){

	fntru_cipher a, b, c;

	f.InitializeCipher(a);
	f.InitializeCipher(b);
	f.InitializeCipher(c);

	ZZX one, zero, result;
	one 	= 1;
	zero 	= 0;

	cout << "\\\\\\\\\\\\ XOR \\\\\\\\\\\\" << endl;
	cout << "a\t" << "b\t" << "c\t" << endl << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, zero);
	f.XOR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, zero);
	f.Encrypt(b, one);
	f.XOR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "0\t" << "1\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, zero);
	f.XOR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.Encrypt(b, one);
	f.XOR(c, a, b, 1);
	f.Decrypt(result, c);
	cout << "1\t" << "1\t" << ConvMessPoly(result) << endl;
}

void TestNOT(fntru &f){

	fntru_cipher a, c;

	f.InitializeCipher(a);
	f.InitializeCipher(c);

	ZZX one, zero, result;
	one 	= 1;
	zero 	= 0;

	cout << "\\\\\\\\\\\\ NOT \\\\\\\\\\\\" << endl;
	cout << "a\t" << "c\t" << endl << endl;

	f.Encrypt(a, zero);
	f.NOT(c, a, 1);
	f.Decrypt(result, c);
	cout << "0\t" << ConvMessPoly(result) << endl;

	f.Encrypt(a, one);
	f.NOT(c, a, 1);
	f.Decrypt(result, c);
	cout << "1\t" << ConvMessPoly(result) << endl;
}

void GateTest(){
	int M 			= M_;
	int radix		= 16;
	int bitSize 	= Dif_Prime;
	ZZ p 			= to_ZZ("2");

	fntru f(M, bitSize, radix,  p);

	TestNOT(f);
	TestXOR(f);
	TestOR(f);
	TestAND(f);
	TestNAND(f);
}


int main(){

	RandomBooleanTest();
	GateTest();

	return 0;
}



















