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

#include "../include/general.h"

void Set(GlobalParam &gp, ZZ q, ZZ B, int M){
	gp.q = q;
	gp.B = B;
	gp.M_CycDegree = M;

	gp.N_PolyDegree = to_long(euler_toient(to_ZZ(gp.M_CycDegree)));
	gp.D_FactDegree = ComputeFactorDegree(gp.M_CycDegree);
	gp.FactSize = gp.N_PolyDegree/gp.D_FactDegree;
}

ZZ euler_toient(ZZ x){
	if(ProbPrime(x) == 1)
		return x-1;

	if(x < to_ZZ("3"))
		return x;

	ZZ res = x;
	ZZ t = to_ZZ("2");
	bool is = false;

	while(x != to_ZZ("1")){
		while(GCD(x,t) == t){
			x=x/t;
			is = true;
		}
		if(is)
			res = res*(t-1)/t;
		is = false;
		t = NextPrime(t+1);
	}
	return res;
}

int ComputeFactorDegree(int m){
	int p=1;
	bool loop = true;

	ZZ t;
	while(loop){
		t = power(to_ZZ("2"), p)-1;
		if(t%m == 0)
			loop = false;
		else
			p++;
	}
	return p;
}

void RandomPolyGen(ZZX &x, int N, int bs, ZZ q){
	clear(x);
	for(int i=0; i<N; i++)
		SetCoeff(x, i, RandomBits_ZZ(bs)%q);
}


//////////////////////////////////////////////////////////////////////////////////


void myTimer::Start(){
	startTime = clock();
}

void myTimer::Stop(){
	stopTime = clock();
}

double myTimer::GetTime(){
	return (double)(stopTime-startTime)/CLOCKS_PER_SEC;
}

void myTimer::ShowTime(string s){
	cout << s << ":\t" << GetTime() << endl;
}

///////////////////////////////////////////////////////////////////////////////////

void myNewTimer::Start(){
	clock_gettime(CLOCK_MONOTONIC, &startTime);
}

void myNewTimer::Stop(){
	clock_gettime(CLOCK_MONOTONIC, &stopTime);
}

double myNewTimer::GetTime(){
	double secs = (stopTime.tv_sec - startTime.tv_sec);
	secs += (stopTime.tv_nsec - startTime.tv_nsec) / 1000000000.0;
	return secs;
}

void myNewTimer::ShowTime(string s){
	cout << s << ":\t" << GetTime() << endl;
}

///////////////////////////////////////////////////////////////////////////////////

void myReduction::SetModulus(ZZX mod){
	modulus = mod;
}

void myReduction::ComputeTable(){

	ZZX temp;
	for(int i=0; i<loop; i++){
		clear(temp);
		SetCoeff(temp, degree_n+(degree_n>>(i+1)), 1);
		x_n[i] = temp % modulus;
	}
}

void myReduction::Set_degree(int n){
	degree_n = n;
	loop = 0;
	while(n>128){
		loop++;
		n=n/2;
	}
	x_n = new ZZX[loop];
}

void myReduction::Div_high_low(ZZX &high, ZZX &low, ZZX &t, int low_deg, int high_deg){
	for(int i=0; i<low_deg; i++)
		SetCoeff(low, i, coeff(t, i));

	for(int i=low_deg; i<high_deg+1; i++)
		SetCoeff(high, i-low_deg, coeff(t, i));
}

void myReduction::Reduce(ZZX &out, ZZX &in){

	int low_deg, high_deg;
	ZZX high, low, t;
	t=in;
	for(int i=0; i<loop; i++){
		high_deg = degree_n+(degree_n>>(i));
		low_deg = degree_n+(degree_n>>(i+1));

		Div_high_low(high, low, t, low_deg, high_deg);
		t = low + high*x_n[i];
		clear(low);
		clear(high);
	}
	out = t%modulus;
}

///////////////////////////////////////////////////////////////////////////////////


void start(myTimerReal &t){
	gettimeofday(&t.begin, NULL);
}

void stop(myTimerReal &t){
	gettimeofday(&t.end, NULL);
}

double getseconds(myTimerReal &t){
	double begin_time = (double)(t.begin.tv_sec) * 1000000 + (double)(t.begin.tv_usec);
	double end_time   = (double)(t.end.tv_sec) * 1000000 + (double)(t.end.tv_usec);

	double r = (end_time - begin_time)/1000000.0;
	return r;
}












