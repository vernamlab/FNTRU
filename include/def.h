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

#ifndef DEF_H_
#define DEF_H_

#define Dif_Prime   (200)
#define M_			1285

// Define thread if you want parallelism
#define Thread
// Matrix vector multiplication. TODO: Check 3 options for which one is more optimized for speed and area.
#define MULT2

//0 is for random Poly modulus.
//1 is x^n-1 poly modulus.
//2 is x^n+1 poly modulus.
#define PolyModXN_1		2

//Set Randomness
#define Random 1
//Randomness
#if Random == 1
	#define TrueRandom
#else
	#define PseudoRandom
#endif


#endif /* DEF_H_ */






