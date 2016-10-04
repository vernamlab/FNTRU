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

#include "../include/ssntru.h"

//////////////////////////////////////////////////////////////////////////////////////////////
// Create SSNTRU
ssntru::ssntru(int dm, int num, int max_bit, int tableS, ZZ qq, int wordSize){

	word 	= wordSize;
	ZZ BB 	= to_ZZ("1");
	Set(gp, qq, BB, dm);

	p 			= gp.q;
	B 			= gp.B;
	N 			= gp.N_PolyDegree;
	degree_m 	= gp.M_CycDegree;
	rand_seed 	 = 0;
	message_rand = 0;
	reducsetflag = false;

	clear(modulus, N+1);
	ZZ_p::init(p);

	word_blocksize = NumofWords(max_bit);
	max_bitsize = max_bit;

	pk = new ZZX[num];
	sk  = new ZZX[num];
	q_list = new ZZ[num];

	//Set Randomness
#ifdef TrueRandom
	srand(time(NULL));
	SetSeed(to_ZZ(time(NULL)));
#endif

#ifdef PseudoRandom
	srand(5);
	SetSeed(to_ZZ("5"));
#endif
}

int ssntru::NumofWords(int bs){
	int size;

	size = bs/word;
	if (bs%word != 0)
		size = size+1;

	return size;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set sigma
void ssntru::SetSigma(long sig){
	sigma = sig;
}
// Encryption function
ZZX ssntru::Prim_Encrypt(ZZX m, int i){
	//log(N*log(N)/log(2.0))/log(2.0)
	ZZX s, e, result;
	s = Func_Sample(log(sqrt(N*log(N)/log(2.0)))/log(2.0)+3);
	e = Func_Sample(log(sqrt(N*log(N)/log(2.0)))/log(2.0)+3);
	coeff_reduction(s, s, q_list[i], N);
	coeff_reduction(e, e, q_list[i], N);
	result = pk[i]*s + e*p + m;
	result = Arith_PolyCoeffReduce(result, i);

	return result;
}
//Compute secret and public keys
void ssntru::Func_ComputeKeys(){

	if(reducsetflag == false)
		SetFastPolyModulus();

////////////////////////First Key//////////////////////////////////////
	ZZX f, g, ft, f_inv;
	bool isfound = false;
//log(N*log(N)/log(2.0))/log(2.0)
//sigma
	ft  = Func_Sample(sigma+3);
	ZZX ft_ = ft;

	int i=0;
	while(!isfound){
		isfound = true;
		coeff_reduction(ft, ft_, q_list[i], N);
		f = ft*p + 1;
		coeff_reduction(f, f, q_list[i], N);
		find_inverse(f_inv, f, q_list[i], N, isfound, modulus);
		coeff_reduction(f_inv, f_inv, q_list[i], N);
	}

	isfound = false;
	g 	  = Func_Sample(sigma+3);

	coeff_reduction(g, g, q_list[i], N);
	pk[i] = g*f_inv;
	Arith_PolyReduce(pk[i], pk[i]);

	pk[i] = pk[i]*p;
	sk[i] = f;

	coeff_reduction(pk[i], pk[i], q_list[i], N);
	coeff_reduction(sk[i], sk[i], q_list[i], N);
}

void ssntru::Func_ModulusFind(int num, int max_bit){

	ZZ 	t, co, k;
	max_bitsize = max_bit;
	init_bitsize = max_bitsize;
	GenPrime(t, max_bitsize);

	k = t/p;
	t = k*p+1;
	while(ProbPrime(t, 80) != 1){
		k++;
		t = k*p+1;
	}

	q_list[0] = t;

}

ZZX ssntru::Prim_Decrypt(ZZX c, int index){

	ZZ x;
	ZZX result, t;

	t = sk[index]*c;
	t = Arith_PolyCoeffReduce(t, index);

	for(int i=0; i<N; i++){
		x = coeff(t, i);
		if(x>((q_list[index]-1)/2))
			x = x-q_list[index];
		SetCoeff(result, i, (x%p));
	}
	return result;
}

ZZ Max2(ZZX in, int deg){
	ZZ max;
	max = in[0];
	for(int i=1; i<deg; i++)
		if(max<in[i])
			max = in[i];

	return max;
}

ZZ Min2(ZZX in, int deg){
	ZZ min;
	min = in[0];
	for(int i=1; i<deg; i++)
		if(min>in[i])
			min = in[i];

	return min;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

ZZX ssntru::Func_Sample(){
	ZZX result;
	for(int i = 0; i < N;i ++)
	{
		ZZ tp;
		RandomBits(tp, to_long(B));
		if(rand()%2)
			SetCoeff(result,i, tp);
		else
			SetCoeff(result,i, 0-tp);
	}
	return result;
}

ZZX ssntru::Func_Sample(long logq){

	ZZX result;
	for(int i = 0; i < N;i ++){
		ZZ tp;
		RandomBits(tp, logq);
		if(rand()%2)SetCoeff(result,i, tp);
		else SetCoeff(result,i, 0-tp);
	}
	return result;
}

void ssntru::Func_SetModulus(){
#if PolyModXN_1 == 1
	SetCoeff(modulus, N, 1);
	SetCoeff(modulus, 0, -1);
#elif PolyModXN_1 == 2
	SetCoeff(modulus, N, 1);
	SetCoeff(modulus, 0, 1);
#else
	modulus = ComputeFastCycModulus(degree_m);
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void PolyRedXN_1(ZZX &out, ZZX &in, int N){
	ZZX res;
	for(int i=0; i<N; i++){
		SetCoeff(res, i, coeff(in,i)+coeff(in,i+N));
	}
	out = res;
}

void PolyRedXNP1(ZZX &out, ZZX &in, int N){
	ZZX res;
	for(int i=0; i<N; i++){
		SetCoeff(res, i, coeff(in,i)-coeff(in,i+N));
	}
	out = res;
}

void ssntru::Arith_PolyReduce(ZZX &out, ZZX &in){
#if PolyModXN_1 == 0
	myr.Reduce(out, in);
#elif PolyModXN_1 == 1
	PolyRedXN_1(out, in, N);
#elif PolyModXN_1 == 2
	PolyRedXNP1(out, in, N);
#endif
}

void ssntru::Arith_CoeffReduce(ZZX &out, ZZX &in, int index){
	int l = in.rep.length();
	coeff_reduction(out, in, q_list[index], l);
}

ZZX ssntru::Arith_PolyCoeffReduce(ZZX &in, int index){
	ZZX r;
	r = in;
	Arith_CoeffReduce(r, r, index);
	Arith_PolyReduce(r, r);
	Arith_CoeffReduce(r, r, index);
	return r;
}

ZZX ssntru::Arith_MulModZZX(ZZX &a, ZZX &b, int index){
	ZZX t;
	t = a*b;
	Arith_PolyReduce(t, t);
	coeff_reduction(t, t, q_list[index], N);
	return t;
}

ZZX ssntru::Arith_AddModZZX(ZZX &a, ZZX &b, int index){
	ZZX t = (a+b);
	Arith_PolyReduce(t, t);
	coeff_reduction(t, t, q_list[index], N);
	return t;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

ZZ *ssntru::Pointer_Q(int i){
	return &(q_list[i]);
}

ZZX	*ssntru::Pointer_PolyModulus(){
	return (&modulus);
}

GlobalParam *ssntru::Pointer_GP(){
	return (&gp);
}

ZZX	*ssntru::Pointer_PublicKey(){
	return &(pk[0]);
}

ZZX	*ssntru::Pointer_SecretKey(int i){
	return &(sk[i]);
}

int *ssntru::Pointer_N(){
	return &(N);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


ZZX ssntru::ComputeFastCycModulus(int n){
	ZZX modulus;
	int s;
	modulus = 1;

	ZZX t_vec[n];
	int s_vec[n];

	for(int i=0; i<n; i++)
		s_vec[i] = 0;

	for(int d=1; d<=n; d++){
			if(GCD(d, n) == d){

				ZZX t;
				SetCoeff(t, 0 , -1);
				SetCoeff(t, n/d, 1);

				s = MobuisFunction(d);

				t_vec[d-1] = t;
				s_vec[d-1] = s;
			}
	}

	for(int i=0; i<n; i++)
		if(s_vec[i] == 1)
			modulus = modulus * t_vec[i];

	for(int i=0; i<n; i++)
		if(s_vec[i] == -1)
			modulus = modulus /  t_vec[i];

	return modulus;
}

void ssntru::SetFastPolyModulus(){
	myr.SetModulus(modulus);
	myr.Set_degree(N);
	myr.ComputeTable();
	reducsetflag = true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

void clear(ZZX &x, int size){
	for(int i=0; i<size; i++)
		SetCoeff(x, i, to_ZZ("0"));
}

void coeff_reduction(ZZX &out, ZZX &in, ZZ &q, int &degree){
	for(int i=0; i<degree; i++)
		SetCoeff(out, i, (coeff(in,i)%q));
}

void find_inverse(ZZX &f_inv, ZZX &f, ZZ &q, int &degree, bool &isfound, ZZX modulus){

	ZZ_p::init(q);
	ZZ_pX phi;
	phi = to_ZZ_pX(modulus);
	ZZ_pE::init(phi);

	ZZ_pE f_, f_inv_;
	f_ = to_ZZ_pE(to_ZZ_pX(f));
	try{ f_inv_ = inv(f_); }
	catch(int e){ isfound = false; }
	ZZ_pX tp = rep(f_inv_);

	for(int i=0; i<degree; i++)
		SetCoeff(f_inv, i, rep(coeff(tp, i)));
}

void coeff_reduction_q_2(ZZX &out, ZZX &in, ZZ &q, int N){
	ZZ t;
	ZZ s = (q-1)/2;
	for(int i=0; i<N; i++){
		t = coeff(in, i);
		if(t>s)
			t = t-q;
		SetCoeff(out, i, t);
	}
}

void getClosestMod2(ZZ& out, ZZ& in, ZZ& p1, ZZ& p2)
{
	ZZ result;

	result = in  * p2 / p1;
	if(bit(result,0)!=bit(in,0)){

		if(result > 0)result ++;
		else result --;
	}
	out = result;
}

int MobuisFunction(int n){
	int t, primes;
	primes = 0;

	if(n == 1)
		return 1;
	else{
		for(int i=2; i<=n; i++){
			if(ProbPrime(i)){
				if(GCD(i,n) == i){
					t=n/i;
					primes++;
					if(GCD(i, t) == i)
						return 0;
				}
			}
		}
		if(primes%2 == 0)
			return 1;
		else
			return -1;
	}
}


double randn (double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;

  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }

  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (double) X1);
}


ZZ uniformRand(long qLog){

	ZZ q = power(to_ZZ("2"), qLog);
	int l = qLog/32+2;

	ZZ t;
	t = 0;
	for(int i=0; i<l; i++)
		t = (t<<32) + to_ZZ(rand());
	t = t%q;

	return t;
}

RR	randn (long mu, long sigmaLog, long qLog){

	RR sigma = pow(to_RR("2"), to_RR(sigmaLog)) / pow(to_RR("2"), to_RR(qLog));
	RR MAX = pow(to_RR("2"), to_RR(qLog));
	RR U1, U2, W, mult;
	static RR X1, X2;
	static int call = 0;

	if (call == 1){
		call = !call;
		return (mu + sigma * X2);
	}

	RR u1, u2;
	do{
//		U1 = -1 + (to_RR(RandomBits_ZZ(qLog)) / MAX) * 2;
//		U2 = -1 + (to_RR(RandomBits_ZZ(qLog)) / MAX) * 2;

//		U1 = -1 + (to_RR(RandomBnd(power(to_ZZ("2"), qLog))) / MAX) * 2;
//		U2 = -1 + (to_RR(RandomBnd(power(to_ZZ("2"), qLog))) / MAX) * 2;

		U1 = -1 + to_RR(uniformRand(qLog)) / MAX * 2;
		U2 = -1 + to_RR(uniformRand(qLog)) / MAX * 2;


		W = U1*U1+U2*U2;

	}while (W >= 1 || W == 0);

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return (mu + sigma * X1);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

