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


#ifndef GENERAL_H_
#define GENERAL_H_

#include <time.h>
#include <iostream>
#include <sys/time.h>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

using namespace std;
NTL_CLIENT

//////////////////////////////////////////////////////////////////////////////////

typedef struct{
	ZZ q;
	ZZ B;
	int M_CycDegree;
	int N_PolyDegree;
	int D_FactDegree;
	int FactSize;
}GlobalParam;

void 	Set(GlobalParam &gp, ZZ q, ZZ B, int M);
ZZ 		euler_toient(ZZ x);
int		ComputeFactorDegree(int m);
void 	RandomPolyGen(ZZX &x, int N, int bs, ZZ q);

//////////////////////////////////////////////////////////////////////////////////

class myTimer{
public:
	void Start();
	void Stop();
	double GetTime();
	void ShowTime(string s);
private:
	clock_t startTime;
	clock_t stopTime;
};

//////////////////////////////////////////////////////////////////////////////////

class myNewTimer{
public:
	void Start();
	void Stop();
	double GetTime();
	void ShowTime(string s);
private:
	struct timespec startTime, stopTime;
};

//////////////////////////////////////////////////////////////////////////////////

class myReduction{
public:
	void Set_degree(int n);
	void SetModulus(ZZX mod);
	void ComputeTable();
	void Reduce(ZZX &out, ZZX &in);
	void Div_high_low(ZZX &high, ZZX &low, ZZX &t, int low_deg, int high_deg);
private:
	ZZX modulus;
	int degree_n;
	ZZX *x_n;
	int loop;
};

//////////////////////////////////////////////////////////////////////////////////


struct myTimerReal{
	struct timeval begin, end;
};

void start(myTimerReal &t);
void stop(myTimerReal &t);
double getseconds(myTimerReal &t);

#endif

















