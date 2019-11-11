
// MFCApplication1Dlg.cpp : 实现文件
//

#include "stdafx.h"
#include "MFCApplication1.h"
#include "MFCApplication1Dlg.h"
#include "afxdialogex.h"
#include <pthread.h>
using namespace std;
uniform_int_distribution<unsigned int> distribution(0, MAXN - 1);// 左闭右闭 
default_random_engine rng(static_cast<unsigned int>(time(nullptr)));
LARGE_INTEGER nFreq, t1, t2;
double tt = 0;

//static vector<unsigned int>littlePrime;
int RSA::fac[8];
bool RSA::threadok[8];
BI RSA::threadPrimes[8];
bool RSA::stop1 = 0;
bool RSA::stop2 = 0;;
vector<unsigned int>RSA::littlePrime;
double RSA::keygenTime = 0;
int BI::fac = 0;
double BI::pi = acos(-1.0);
int BI::G = 3;
int BI::NTTP = (119 << 23) + 1;

RSA::RSA() {
	// init little primes with linear
	unsigned int n = MAXN; 
	//	unsigned int n=1000;
	littlePrime.clear();
	for (unsigned int i = 2; i <= n; i++) {
		if (!ok[i])littlePrime.push_back(i);
		for (unsigned int j = 0; j<littlePrime.size() && i*littlePrime[j] <= n; j++) {
			ok[i*littlePrime[j]] = 1;
			if (i%littlePrime[j] == 0)break;
		}
	}
	// optimize: littlePrime combine used for gcd
	priority_queue<unsigned int, vector<unsigned int>, greater<unsigned int> >q;
	for (int i = 0; i<littlePrime.size(); i++)q.push(littlePrime[i]);
	unsigned int a, b;
	while (1) {
		if (!q.empty()) { a = q.top(); q.pop(); }
		else break;
		if (!q.empty()) { b = q.top(); q.pop(); }
		else { q.push(a); break; }
		if (a*b<65536)q.push(a*b);
		else { q.push(a); q.push(b); break; }
	}
	littlePrime.clear();
	while (!q.empty()) { littlePrime.push_back(q.top()); q.pop(); }

	//	cout<<littlePrime.size()<<endl; // 65536内6530个。。有点多，用1000内的吧168个
}
string RSA::publish() {
	string s;
	s += "e:\r\n";
	s += e.tostring();
	s += "(65537)\r\n";
	s += "d:\r\n";
	s += d.tostring();
	s += "\r\n";
	s += "n:\r\n";
	s += n.tostring();
	s += "\r\n";
	return s;
}
string BI::tostring() {
	string s="";
	char c[4];
	for (int i = len - 1; i >= 0; i--) {
		unsigned int b = a[i];
		for (int j = 0; j < 4; j++) {
			c[j] = (b & 15) > 9 ? ((b & 15) - 10 + 'a') : ((b & 15) + '0');
			b /=16;
		}
		for(int j=3;j>=0;j--)s += c[j];
	}
	return s;
}
void RSA::keygen(int bit) {
	QueryPerformanceCounter(&t1);
	int threads = 8;
	pthread_t tids[8];
	for (int i = 0; i < 8; i++)threadok[i] = 0;
	stop1 = stop2 = 0;
	for (int i = 0; i < threads / 2; i++) {
		threadmessage* msg=(threadmessage*)malloc(sizeof(threadmessage));
		msg->bit = bit / 2;
		msg->tid= i;
		int ret = pthread_create(&tids[i], NULL, primeGenThread1, (void*)msg);
		if (ret != 0) cout << "pthread_create error: error_code=" << ret << endl;
	}
	for (int i = 4; i < threads; i++) {
		threadmessage* msg = (threadmessage*)malloc(sizeof(threadmessage));;
		msg->bit = bit / 2;
		msg->tid = i;
		int ret = pthread_create(&tids[i], NULL, primeGenThread2, (void*)msg);
		if (ret != 0) cout << "pthread_create error: error_code=" << ret << endl;
	}
	void * status;
	for(int i=0;i<threads;i++)pthread_join(tids[i], &status);
	for (int i = 0; i < 4; i++)if (threadok[i])p = threadPrimes[i];
	for (int i = 4; i < 8; i++)if (threadok[i])q = threadPrimes[i];
	//p = primeGen(bit / 2);
	//q = primeGen(bit / 2);
	n = p*q;
	phi = (p - 1)*(q - 1);
	e = BI(65537);
	BI temp;
	int f1, f2;
	//	QueryPerformanceCounter(&t1);
	BI G = BI::exgcd(e, d, phi, temp, f1, f2);
	//	QueryPerformanceCounter(&t2);
	//	cout<<(t2.QuadPart - t1.QuadPart) / (double)nFreq.QuadPart<<endl;
	G.print();// 1 is right
	if (f1 == -1) {
		f1 = 1;
		d = phi - d%phi;
	}
	//	e.print();
	//	d.print();
	(e*d%phi).print();// 1 is right
					  //	phi.print();

	QueryPerformanceCounter(&t2);
	double dt = (t2.QuadPart - t1.QuadPart) / (double)nFreq.QuadPart;
	keygenTime = dt;
}
string RSA::encode(string s) {
	// s^e
	BI M;
	M.len = 0;
	for (int i = s.length() - 1; i >= 0; i -= 2) {
		if (i == 0) {
			unsigned int a = s[i];
			M.a[M.len] = a;
			M.len++;
			break;
		}
		unsigned int a = s[i - 1];
		if (i != 0) {
			a <<= 8;
			a += s[i];
		}
		M.a[M.len] = a;
		M.len++;
	}
	//BI ni = BI::newton(n, n.len+1);
	//BI E = BI::powNewton(M, e, n, ni);
		BI E=BI::pow(M,e,n);
	//	E.print();
	//	E=hex2b(b2hex(E));
	//	E.print();
	return b2hex(E);
}
string RSA::decode(string s) {
	BI E = hex2b(s);
	BI M;
	{// CTR + newton
		BI ni = BI::newton(p, p.len + 1);
		BI a1 = BI::powNewton(E, d, p, ni);
		ni = BI::newton(q, q.len + 1);
		BI a2 = BI::powNewton(E, d, q, ni);
		//		M=CRT(a1,a2,p,q);
		BI qni, pni;
		int f1, f2;
		BI::exgcd(q, qni, p, pni, f1, f2);
		if (f1 == -1)qni = p - qni%p;
		if (f2 == -1)pni = q - pni%q;
		M = (a1*q*qni + a2*p*pni) % n;
	}
		//{// simple newton
		//	BI ni=BI::newton(n,n.len+1);
		//	M=BI::powNewton(E,d,n,ni);
		//}
	s = "";
	for (int i = 0; i<M.len; i++) {
		char c = M.a[i] & 255;
		s += c;
		if (i != M.len - 1) {
			c = (M.a[i] >> 8) & 255; s += c;
		}
		else {
			c = (M.a[i] >> 8) & 255;
			if (c != 0)s += c;
		}
	}
	reverse(s.begin(), s.end());
	return s;
}
string RSA::b2hex(BI A) {
	string s = "";
	for (int i = 0; i<A.len; i++) {
		char c = '0' + (A.a[i] & 15); s += (c>'9' ? (c - ':') + 'a' : c);
		c = '0' + ((A.a[i] >> 4) & 15); s += (c>'9' ? (c - ':') + 'a' : c);
		c = '0' + ((A.a[i] >> 8) & 15); s += (c>'9' ? (c - ':') + 'a' : c);
		c = '0' + ((A.a[i] >> 12) & 15); s += (c>'9' ? (c - ':') + 'a' : c);
	}
	return s;
}
BI RSA::hex2b(string s) {
	while (s.length() % 4 != 0)s = "0" + s;
	BI E; E.len = 0;
	unsigned a, b;
	for (int i = 0; i<s.length(); i += 4) {
		a = (s[i]>'9' ? (s[i] - 'a') + 10 : s[i] - '0');
		b = (s[i + 1]>'9' ? (s[i + 1] - 'a') + 10 : s[i + 1] - '0'); b <<= 4; a += b;
		b = (s[i + 2]>'9' ? (s[i + 2] - 'a') + 10 : s[i + 2] - '0'); b <<= 8; a += b;
		b = (s[i + 3]>'9' ? (s[i + 3] - 'a') + 10 : s[i + 3] - '0'); b <<= 12; a += b;
		E.a[E.len] = a;
		E.len++;
	}
	return E;
}
void* RSA::primeGenThread1(void* mmsg) {
	threadmessage* msg = (threadmessage*)mmsg;
	BI P;
	int bit = msg->bit;
	int tid = msg->tid;
	int k = bit / 16;
	P.len = k;

	int againTime = 0;
	P.a[k - 1] = 32768 + distribution(rng) % 32768;
	for (int i = 0; i<k - 1; i++)P.a[i] = distribution(rng);
again:
	if (stop1)return 0;
	//		P=P+littlePrime[littlePrime.size()-1];
	P = P + 19260817;// 玄学加成 
	P.len = k; P.a[k - 1] |= 32768;
	while (P.len>k) {
		P.a[P.len - 1] = 0;
		P.len--;
	}P.a[P.len - 1] |= MAXN / 2;
	for (int i = 0; i<littlePrime.size(); i++) {
		// 大数gcd小数快速检验 -- 优化版 
		unsigned int r = 0;
		for (int j = P.len - 1; j >= 0; j--)r = (r*MAXN + P.a[j]) % littlePrime[i];
		if (r == 0 || BI::__gcd(r, littlePrime[i]) != 1)goto again;
	}
	// (P-1)'s suffix 0s
	BI Q = P - BI(1);
	int cnt0 = 0;
	for (int i = 0; i<k; i++)
		for (unsigned int j = 1; j<MAXN; j <<= 1)
			if (Q.a[i] & j)goto out;
			else cnt0++;
	out:
		unsigned int shift = cnt0 % 16;// 16 is log_2(MAXN)
		unsigned int offset = cnt0 / 16;
		Q.len = k - offset;
		for (int i = 0; i<Q.len; i++)Q.a[i] = (Q.a[i + offset] >> shift) + (Q.a[i + offset + 1] << (16 - shift)) % MAXN;
	// miller rabin
	int times = 2;
	if (bit >= 1536 / 2)times = 1;
	againTime++;
	while (times--) {
		if (stop1)return 0;
		BI R; R.len = k;
		for (int i = 0; i<k; i++)R.a[i] = distribution(rng);
		if (R>P)R = R - P;

		{// newton method
//			BI ni = BI::newton(P, P.len + 1);
			BI ni = BI::newtonThread(P, P.len + 1,tid);
			R = BI::powNewtonThread(R, Q, P, ni,tid);
			if (R.len != 1 && R.a[0] != 1)goto again;
			while (cnt0-- && !(R.len == 1 && R.a[0] == 1)) {
				if (R.len != 1 || R.a[0] != MAXN - 1)goto again;
				R = R*R;
				if (R>P)R = BI::modNewtonThread(R, P, ni, tid);
			}
		}
	}
	threadok[tid] = 1;
	stop1 = 1;
	threadPrimes[tid] = P;
}
void* RSA::primeGenThread2(void* mmsg) {
	threadmessage* msg = (threadmessage*)mmsg;
	BI P;
	int bit = msg->bit;
	int tid = msg->tid;
	int k = bit / 16;
	P.len = k;

	int againTime = 0;
	P.a[k - 1] = 32768 + distribution(rng) % 32768;
	for (int i = 0; i<k - 1; i++)P.a[i] = distribution(rng);
again:
	if (stop2)return 0;
	//		P=P+littlePrime[littlePrime.size()-1];
	P = P + 19260817;// 玄学加成 
	P.len = k; P.a[k - 1] |= 32768;
	while (P.len>k) {
		P.a[P.len - 1] = 0;
		P.len--;
	}P.a[P.len - 1] |= MAXN / 2;
	for (int i = 0; i<littlePrime.size(); i++) {
		// 大数gcd小数快速检验 -- 优化版 
		unsigned int r = 0;
		for (int j = P.len - 1; j >= 0; j--)r = (r*MAXN + P.a[j]) % littlePrime[i];
		if (r == 0 || BI::__gcd(r, littlePrime[i]) != 1)goto again;
	}
	// (P-1)'s suffix 0s
	BI Q = P - BI(1);
	int cnt0 = 0;
	for (int i = 0; i<k; i++)
		for (unsigned int j = 1; j<MAXN; j <<= 1)
			if (Q.a[i] & j)goto out;
			else cnt0++;
	out:
		unsigned int shift = cnt0 % 16;// 16 is log_2(MAXN)
		unsigned int offset = cnt0 / 16;
		Q.len = k - offset;
		for (int i = 0; i<Q.len; i++)Q.a[i] = (Q.a[i + offset] >> shift) + (Q.a[i + offset + 1] << (16 - shift)) % MAXN;
	// miller rabin
	int times = 2;
	if (bit >= 1536 / 2)times = 1;
	againTime++;
	while (times--) {
		if (stop2)return 0;
		BI R; R.len = k;
		for (int i = 0; i<k; i++)R.a[i] = distribution(rng);
		if (R>P)R = R - P;

		{// newton method
			//			BI ni = BI::newton(P, P.len + 1);
			BI ni = BI::newtonThread(P, P.len + 1, tid);
			R = BI::powNewtonThread(R, Q, P, ni, tid);
			if (R.len != 1 && R.a[0] != 1)goto again;
			while (cnt0-- && !(R.len == 1 && R.a[0] == 1)) {
				if (R.len != 1 || R.a[0] != MAXN - 1)goto again;
				R = R*R;
				if (R>P)R = BI::modNewtonThread(R, P, ni, tid);
			}
		}
	}
	threadok[tid] = 1;
	stop2 = 1;
	threadPrimes[tid] = P;
}
BI RSA::primeGen(int bit) {
	BI P;
	int k = bit / 16;
	P.len = k;

	int againTime = 0;
	P.a[k - 1] = 32768 + distribution(rng) % 32768;
	for (int i = 0; i<k - 1; i++)P.a[i] = distribution(rng);
again:
	//		P=P+littlePrime[littlePrime.size()-1];
	P = P + 19260817;// 玄学加成 
	P.len = k; P.a[k - 1] |= 32768;
	while (P.len>k) {
		P.a[P.len - 1] = 0;
		P.len--;
	}P.a[P.len - 1] |= MAXN / 2;
	//		for(int i=0;i<littlePrime.size();i++){
	//			// 大数模小数快速检验
	//			unsigned int r=0;
	//			for(int j=P.len-1;j>=0;j--)r=(r*MAXN+P.a[j])%littlePrime[i];
	//			if(r==0)goto again;
	//		}
	for (int i = 0; i<littlePrime.size(); i++) {
		// 大数gcd小数快速检验 -- 优化版 
		unsigned int r = 0;
		for (int j = P.len - 1; j >= 0; j--)r = (r*MAXN + P.a[j]) % littlePrime[i];
		if (r == 0 || BI::__gcd(r, littlePrime[i]) != 1)goto again;
	}
	// (P-1)'s suffix 0s
	BI Q = P - BI(1);
	int cnt0 = 0;
	for (int i = 0; i<k; i++)
		for (unsigned int j = 1; j<MAXN; j <<= 1)
			if (Q.a[i] & j)goto out;
			else cnt0++;
		out:
			unsigned int shift = cnt0 % 16;// 16 is log_2(MAXN)
			unsigned int offset = cnt0 / 16;
			Q.len = k - offset;
			for (int i = 0; i<Q.len; i++)Q.a[i] = (Q.a[i + offset] >> shift) + (Q.a[i + offset + 1] << (16 - shift)) % MAXN;
			// miller rabin
			int times = 1;
			againTime++;
			while (times--) {
				BI R; R.len = k;
				for (int i = 0; i<k; i++)R.a[i] = distribution(rng);
				if (R>P)R = R - P;

				{// newton method
					BI ni = BI::newton(P, P.len + 1);
					R = BI::powNewton(R, Q, P, ni);
					if (R.len != 1 && R.a[0] != 1)goto again;
					while (cnt0-- && !(R.len == 1 && R.a[0] == 1)) {
						if (R.len != 1 || R.a[0] != MAXN - 1)goto again;
						R = R*R;
						if (R>P)R = BI::modNewton(R, P, ni);
					}
				}

				//			{// simple method
				//				R=BI::pow(R,Q,P);// R^Q mod P
				//	//			R=BI::pow(R,P-BI(1),P);
				//	//			if(R.len!=1||R.a[0]!=1)goto again;
				//	//			while(cnt0--)R=R*R%P;
				//				if(R.len!=1&&R.a[0]!=1)goto again;
				//				while(cnt0--&&!(R.len==1&&R.a[0]==1)){
				//					if(R.len!=1||R.a[0]!=MAXN-1)goto again;
				//					R=R*R;
				//					if(R>P)R=R%P;
				//				}
				//			}
			}
			cout << againTime << endl;
			cout << tt << endl;
			tt = 0;
			return P;
}
//BI::BI(const char*s){
//    int t,k,index,l,i;
//    const int DLEN=4;
//    memset(a,0,sizeof(a));
//    l=strlen(s);
//    len=l/DLEN;
//    if(l%DLEN)len++;
//    index=0;
//    for(i=l-1;i>=0;i-=DLEN){
//        t=0;
//        k=i-DLEN+1;
//        if(k<0)
//            k=0;
//        for(int j=k;j<=i;j++)
//            t=t*10+s[j]-'0';
//        a[index++]=t;
//    }
//}
BI BI::pow(BI A, const BI &B, const BI &Q) {
	BI res(1);
	int b = 0, k = B.len * 16;
	for (int i = B.len - 1; i >= 0; i--)
		for (int shift = 15; shift >= 0; shift--)
			if ((B.a[i] >> shift) & 1)goto over;
			else k--;
		over:
			;
			while (b <= k) {
				if ((B.a[b / 16] >> (b % 16)) & 1)res = res*A%Q;
				b++;
				A = A*A%Q;
			}
			return res;
}
BI BI::divNewton(const BI &A, const BI &ni) {
	BI R = A*ni;
	//R.fix999(R.len - BI::fac);
	R.fix(R.len - BI::fac);
	return R;
}
BI BI::modNewton(const BI &A, const BI &B, const BI &ni) {
	BI R = A*ni;
	//R.fix999(R.len - BI::fac);
	R.fix(R.len - BI::fac);
	R = A - R*B;
	return R;
}
BI BI::powNewton(BI A, const BI &B, const BI &Q, const BI& ni) {
	BI res(1);
	int b = 0, k = B.len * 16;
	for (int i = B.len - 1; i >= 0; i--) {
		for (int shift = 15; shift >= 0; shift--) {
			if ((B.a[i] >> shift) & 1)goto over;
			else k--;
		}
	}
	over:
		;
	while (b <= k) {
		if ((B.a[b / 16] >> (b % 16)) & 1) {
			res = res*A;
			if(res>Q)
				res = modNewton(res, Q, ni);
		}
		b++;
		A = A*A;
		A = modNewton(A, Q, ni);
	}
	return res;
}
BI::BI(const unsigned int b) {
	unsigned int d = b; len = 0; memset(a, 0, sizeof(a));
	while (d >= MAXN) { a[len++] = d%MAXN; d = d / MAXN; }a[len++] = d;
}
BI::BI(const BI &B) :len(B.len) {
	memset(a, 0, sizeof(a));
	for (int i = 0; i<len; i++)a[i] = B.a[i];
}
BI &BI::operator=(const BI &B) {
	memset(a, 0, sizeof(a)); len = B.len;
	for (int i = 0; i<len; i++)a[i] = B.a[i];
	return *this;
}
BI BI::operator+(const BI &B) const {
	BI A(*this);
	int llen = max(len, B.len);
	for (int i = 0; i<llen; i++) {
		A.a[i] += B.a[i];
		if (A.a[i] >= MAXN) { A.a[i + 1]++; A.a[i] -= MAXN; }
	}
	if (A.a[llen] != 0)A.len = llen + 1;
	else A.len = llen;
	return A;
}
// 一定保证大减小 
BI BI::operator-(const BI &B) const {
	BI A(*this);
	int llen = A.len;
	for (int i = 0; i<llen; i++) {
		if (A.a[i]<B.a[i]) {
			int j = i + 1;
			while (A.a[j] == 0)j++;
			A.a[j--]--;
			while (j>i)A.a[j--] += MAXN - 1;
			A.a[i] += MAXN - B.a[i];
		}
		else A.a[i] -= B.a[i];
	}
	A.len = llen;
	while (A.a[A.len - 1] == 0 && A.len>1)A.len--;
	return A;
}
int BI::rev(int id, int len) {
	int ret = 0;
	for (int i = 0; (1 << i)<len; i++) {
		ret <<= 1;
		if (id&(1 << i))ret |= 1;
	}
	return ret;
}
//void BI::FFT(complex<double> *x, int n, int type, int *R) {
//	for (int i = 0; i < n; i++) if (i < R[i]) swap(x[i], x[R[i]]);
//	for (int i = 1; i < n; i <<= 1) {
//		complex<double> wn(cos(pi / i), type*sin(pi / i));
//		for (int j = 0; j < n; j += i << 1) {
//			complex<double> w(1, 0);
//			for (int k = 0; k < i; k++, w *= wn) {
//				complex<double> X = x[j + k], Y = w*x[j + k + i];
//				x[j + k] = X + Y;
//				x[j + k + i] = X - Y;
//			}
//		}
//	}
//}
BI BI::operator*(const BI &B) const {
	BI A;
	BI A2;
	if (len >= 200 && B.len >= 200) {// NTT
		int size = 1;
		while (size <= len + B.len)size *= 2;
		size *= 2;
		//int fa[size];
		//int fb[size];
		int *fa = (int*)malloc(sizeof(int)*size);
		int *fb = (int*)malloc(sizeof(int)*size);
		memset(fa, 0, sizeof(fa));
		memset(fb, 0, sizeof(fb));
		for (int i = 0; i<len; i++)fa[i * 2] = a[i] & 255, fa[i * 2 + 1] = (a[i] >> 8) & 255;
		for (int i = 0; i<B.len; i++)fb[i * 2] = B.a[i] & 255, fb[i * 2 + 1] = (B.a[i] >> 8) & 255;
		NTT(fa, size, 1, NTTP);
		NTT(fb, size, 1, NTTP);
		for (int i = 0; i<size; i++)fa[i] = 1ll * fa[i] * fb[i] % NTTP;
		NTT(fa, size, -1, NTTP);
		A.len = size / 2;
		for (int i = 0; i<size / 2; i++) {
			A.a[i] += fa[i * 2] + fa[i * 2 + 1] * 256;
			A.a[i + 1] += A.a[i] / MAXN;
			A.a[i] %= MAXN;
		}
		while (A.a[A.len - 1] == 0 && A.len>1)A.len--;
	}
	else {// simple
		A.len = len + B.len;
		//	memset(A.a,0,sizeof(A.a));
		for (int i = 0; i<len; i++) {
			for (int j = 0; j<B.len; j++) {
				A.a[i + j] += a[i] * B.a[j];
				if (A.a[i + j] >= MAXN) {
					A.a[i + j + 1] += A.a[i + j] / MAXN;
					A.a[i + j] %= MAXN;
				}
			}
		}
		while (A.a[A.len - 1] == 0 && A.len>1)A.len--;
	}
	return A;
}
//BI BI::operator*(const BI &B) const{
//	BI A;
//	if(len>=10&&B.len>=10){// FFT
//		int size=1;
//		while(size<=len+B.len)size*=2;
//		cp fa[size];
//		cp fb[size];
//		for(int i=0;i<len;i++)fa[i].r=a[i];
//		for(int i=len;i<size;i++)fa[i].r=0;
//		for(int i=0;i<B.len;i++)fb[i].r=B.a[i];
//		for(int i=B.len;i<size;i++)fb[i].r=0;
//		for(int i=0;i<size;i++)fa[i].i=fb[i].i=0;
//		FFT(fa,size,1);
//		FFT(fb,size,1);
//		for(int i=0;i<size;i++)fa[i]=fa[i]*fb[i];
//		FFT(fa,size,-1);
//		A.len=size;
//		for(int i=0;i<size;i++){
//			A.a[i]=fa[i].r;
//			A.a[i+1]+=A.a[i]/MAXN;
//			A.a[i]%=MAXN;
//		}
//	}else{// simple
//	if(len>=10&&B.len>=10){// FFT
//		int size=1;
//		int L=0;
//		while(size<=len+B.len)size*=2,L++;
//		complex<double> fa[size];
//		complex<double> fb[size];
//		memset(fa,0,sizeof(fa));
//		memset(fb,0,sizeof(fb));
//		for(int i=0;i<len;i++)fa[i]=a[i];
//		for(int i=0;i<B.len;i++)fb[i]=B.a[i];
//		int R[size];
//        for (int i=0;i<size;i++)R[i]=(R[i>>1]>>1)|((i&1)<<L-1);
//		FFT(fa,size,1,R);
//		FFT(fb,size,1,R);
//		for(int i=0;i<size;i++)fa[i]=fa[i]*fb[i];
//		FFT(fa,size,-1,R);
//		A.len=size;
//		for(int i=0;i<size;i++){
//			A.a[i]=fa[i].real();
//			A.a[i+1]+=A.a[i]/MAXN;
//			A.a[i]%=MAXN;
//		}
//	A.len=len+B.len;
////	memset(A.a,0,sizeof(A.a));
//	for(int i=0;i<len;i++){
//		for(int j=0;j<B.len;j++){
//			A.a[i+j]+=a[i]*B.a[j];
//			if(A.a[i+j]>=MAXN){
//				A.a[i+j+1]+=A.a[i+j]/MAXN;
//				A.a[i+j]%=MAXN;
//			}
//		}
//	}
//	while(A.a[A.len-1]==0&&A.len>1)A.len--;
//	return A;
//}
// 这里只管除以2的 
int BI::power(int x, int n, int P) {
	int res = 1;
	for (; n; n >>= 1) {
		if (n & 1)res = 1ll * res*x%P;
		x = 1ll * x*x%P;
	}return res;
}
void BI::rader(int* y, int len) {
	for (int i = 1, j = len / 2; i < len - 1; i++) {
		if (i < j) swap(y[i], y[j]);
		int k = len / 2;
		while (j >= k) { j -= k; k /= 2; }
		if (j < k) j += k;
	}
}
void BI::NTT(int* y, int len, int op, int P) {
	rader(y, len);
	for (int h = 2; h <= len; h <<= 1) {
		int wn = power(G, (P - 1) / h, P);
		if (op == -1)wn = power(wn, P - 2, P);
		for (int j = 0; j<len; j += h) {
			int w = 1;
			for (int k = j; k<j + h / 2; k++) {
				int u = y[k];
				int t = 1ll * w*y[k + h / 2] % P;
				y[k] = (u + t) % P;
				y[k + h / 2] = (u - t + NTTP) % P;
				w = 1ll * w*wn%P;
			}
		}
	}
	if (op == -1) {
		int inv = power(len, P - 2, P);
		for (int i = 0; i<len; i++)y[i] = 1ll * y[i] * inv%P;
	}
}
BI BI::operator/(const unsigned int &) const {
	BI A(*this);
	for (int i = A.len - 1; i>0; i--) { if (A.a[i] & 1)A.a[i - 1] += MAXN; A.a[i] /= 2; }A.a[0] /= 2;
	while (A.a[A.len - 1] == 0 && A.len>1)A.len--;
	return A;
}
BI BI::operator/(const BI &B) const {
	//	QueryPerformanceCounter(&t1);

	BI A(*this);
	if (A<B)return BI(0);
	// 二分试商法
	BI L, R(A), one(1), two(2);
	BI M;
	int ttimes = 0;
	while (L<R) {
		ttimes++;
		M = (L + R + one) / 2;
		if (A<M*B)R = M - one;
		else L = M;
	}
	while (L.a[L.len - 1] == 0 && L.len>1)L.len--;

	//	QueryPerformanceCounter(&t2);
	//	tt+=(t2.QuadPart - t1.QuadPart) / (double)nFreq.QuadPart;
	return L;
}
BI BI::operator%(const BI &B) const {

	BI A(*this);
	A = A - (A / B)*B;
	while (A.a[A.len - 1] == 0 && A.len>1)A.len--;

	return A;
	//	return A-(A/B)*B;
}
BI BI::sub(const int l, const int r) {
	BI B; B.len = r - l + 1;
	for (int i = 0; i<B.len; i++)B.a[i] = a[l + i];
	return B;
}
BI BI::head(const int k) { return sub(len - k, len - 1); }
BI BI::tail(const int k) { return sub(0, k - 1); }
bool BI::operator>(const BI &B) const {
	if (len>B.len)return true;
	else if (len == B.len) {
		int i = len - 1;
		while (a[i] == B.a[i] && i >= 0)i--;
		if (i >= 0 && a[i]>B.a[i])return true;
		else return false;
	}
	else return false;
}
bool BI::operator<(const BI &B) const {
	if (len<B.len)return true;
	else if (len == B.len) {
		int i = len - 1;
		while (a[i] == B.a[i] && i >= 0)i--;
		if (i >= 0 && a[i]<B.a[i])return true;
		else return false;
	}
	else return false;
}
bool BI::operator>(const unsigned int &b) const { BI B(b); return *this>B; }
bool BI::operator<(const unsigned int &b) const { BI B(b); return *this<B; }
bool BI::operator==(const BI &B) const {
	if (len != B.len)return false;
	for (int i = 0; i<len; i++)if (a[i] != B.a[i])return false;
	return true;
}
void BI::print()const {
	BI A(*this); char c[17]; c[16] = '\0';
	for (int i = len - 1; i >= 0; i--) { for (int j = 0; j<16; j++) { c[16 - j - 1] = ((A.a[i] & 1) ? '1' : '0'); A.a[i] /= 2; }cout << c; }
	cout << endl;
}
void BI::printN()const {
	for (int i = len - 1; i >= 0; i--) { printf("%04d", a[i]); }
	cout << endl;
}
void BI::print10()const {
	for (int i = len - 1; i >= 0; i--) { printf("%01d", a[i]); }
	cout << endl;
}
// 万进制的
//ostream& operator<<(ostream& out, BI &b) {
//	int i;
//	out << b.a[b.len - 1];
//	for (i = b.len - 2; i >= 0; --i) {
//		out.width(4);
//		out.fill('0');
//		out << b.a[i];
//	}
//	return out;
//}
void BI::fix(int lim) {// 保留lim的int数组位数
	if (len<lim)return;
	int m = len - lim;
	m = min(m, len);
	for (int i = 0; i<lim; i++)a[i] = a[i + m];
	for (int i = lim; i<len; i++)a[i] = 0;
	len = lim;
}
void BI::fix999(int lim) {
	if (len<lim)return;
	int m = len - lim;
	m = min(m, len);
	int k9 = 0;
	for (int i = m - 1; i >= 0; i--)
		if (a[i] == MAXN - 1)k9++;
		else break;
		for (int i = 0; i<lim; i++)a[i] = a[i + m];
		for (int i = lim; i<len; i++)a[i] = 0;
		len = lim;
		if (k9 >= m / 2) {
			a[0] += 1;
			for (int i = 0; i<len; i++)if (a[i] >= MAXN) { a[i + 1]++; a[i] -= MAXN; }
			if (a[len] != 0)len = len + 1;
		}
}
BI BI::exgcd(const BI A, BI&x, const BI B, BI&y, int &f1, int &f2) {
	if (B == 0) { x = 1; y = 0; f1 = f2 = 1; return A; }
	BI D = exgcd(B, x, A%B, y, f1, f2);
	BI tem = x;
	int temf = f1;
	x = y; f1 = f2;
	BI R = A / B*y;
	if (temf == 1) {
		if (f2 == 1) {
			if (tem<R) { y = R - tem; f2 = -1; }
			else { y = tem - R; f2 = 1; }
		}
		else {
			y = tem + R; f2 = 1;
		}
	}
	else {
		if (f2 == 1) { y = tem + R; f2 = -1; }
		else {
			if (tem<R) { y = R - tem; f2 = 1; }
			else { y = tem - R; f2 = -1; }
		}
	}
	return D;
}
BI BI::newton(const BI &B, int bit) {
	// calc 1/B by x=x(2-Bx)
	//	QueryPerformanceCounter(&t1);

	BI ni(1);
	//int times=50;
	int lim = bit * 2;
	fac = B.len;
	BI temp = ni;

	//    while(times--){
	while (1) {
		BI two;
		two.len = fac + 1;
		two.a[fac] = 2;
		ni = ni*(two - B*ni);
		fac *= 2;
		if (ni.len>lim) {
			fac -= ni.len - lim;
			ni.fix(lim);
		}
		if (temp == ni) {
			//			cout<<times<<endl;
			break;
		}
		temp = ni;
	}
	//	QueryPerformanceCounter(&t2);
	//	tt+=(t2.QuadPart - t1.QuadPart) / (double)nFreq.QuadPart;
	return ni;
}
BI BI::divNewtonThread(const BI &A, const BI &ni,const int& tid) {
	BI R = A*ni;
	//R.fix999(R.len - BI::fac);
	R.fix(R.len - RSA::fac[tid]);
	return R;
}
BI BI::modNewtonThread(const BI &A, const BI &B, const BI &ni, const int& tid) {
	BI R = A*ni;
	//R.fix999(R.len - BI::fac);
	R.fix(R.len - RSA::fac[tid]);
	R = A - R*B;
	return R;
}
BI BI::powNewtonThread(BI A, const BI &B, const BI &Q, const BI& ni, const int& tid) {
	BI res(1);
	int b = 0, k = B.len * 16;
	for (int i = B.len - 1; i >= 0; i--) {
		for (int shift = 15; shift >= 0; shift--) {
			if ((B.a[i] >> shift) & 1)goto over;
			else k--;
		}
	}
over:
	;
	while (b <= k) {
		if ((B.a[b / 16] >> (b % 16)) & 1) {
			res = res*A;
			res = modNewtonThread(res, Q, ni, tid);
		}
		b++;
		A = A*A;
		A = modNewtonThread(A, Q, ni,tid);
	}
	return res;
}
BI BI::newtonThread(const BI &B, int bit, const int& tid) {
	// calc 1/B by x=x(2-Bx)
	//	QueryPerformanceCounter(&t1);

	BI ni(1);
	//    int times=30;
	int lim = bit * 2;
	RSA::fac[tid] = B.len;
	BI temp = ni;

	//    while(times--){
	while (1) {
		BI two;
		two.len = RSA::fac[tid] + 1;
		two.a[RSA::fac[tid]] = 2;
		ni = ni*(two - B*ni);
		RSA::fac[tid] *= 2;
		if (ni.len>lim) {
			RSA::fac[tid] -= ni.len - lim;
			ni.fix(lim);
		}
		if (temp == ni) {
			//			cout<<times<<endl;
			break;
		}
		temp = ni;
	}
	//	QueryPerformanceCounter(&t2);
	//	tt+=(t2.QuadPart - t1.QuadPart) / (double)nFreq.QuadPart;
	return ni;
}
#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CMFCApplication1Dlg 对话框



CMFCApplication1Dlg::CMFCApplication1Dlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_MFCAPPLICATION1_DIALOG, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMFCApplication1Dlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CMFCApplication1Dlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_BUTTON2, &CMFCApplication1Dlg::OnBnClickedButton2)
	ON_BN_CLICKED(IDC_BUTTON3, &CMFCApplication1Dlg::OnBnClickedButton3)
	ON_BN_CLICKED(IDC_BUTTON1, &CMFCApplication1Dlg::OnBnClickedButton1)
END_MESSAGE_MAP()


// CMFCApplication1Dlg 消息处理程序

BOOL CMFCApplication1Dlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: 在此添加额外的初始化代码
	CEdit* pBoxOne;
	pBoxOne = (CEdit*)GetDlgItem(EDIT1);
	pBoxOne->SetWindowText("input message");
	QueryPerformanceFrequency(&nFreq);
	CButton* radio;
	radio = (CButton*)GetDlgItem(radio768);
	radio->SetCheck(1);
	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CMFCApplication1Dlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CMFCApplication1Dlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CMFCApplication1Dlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CMFCApplication1Dlg::OnBnClickedButton2()
{
	// 加密
	CEdit* pBoxOne;
	pBoxOne = (CEdit*)GetDlgItem(EDIT1);
	CString cs;
	pBoxOne->GetWindowTextA(cs);
	pBoxOne->SetWindowText("");
	string s = cs;
	s=rsa.encode(s);
	pBoxOne = (CEdit*)GetDlgItem(EDIT2);
	pBoxOne->SetWindowText(s.data());
}


void CMFCApplication1Dlg::OnBnClickedButton3()
{
	// 解密
	CEdit* pBoxOne;
	pBoxOne = (CEdit*)GetDlgItem(EDIT2);
	CString cs;
	pBoxOne->GetWindowTextA(cs);
	pBoxOne->SetWindowText("");
	string s=cs;
	s=rsa.decode(s);
	pBoxOne = (CEdit*)GetDlgItem(EDIT1);
	pBoxOne->SetWindowText(s.data());
}
//int bb;
//void* add(void* a){
//	bb++;
//	return 0;
//}
void CMFCApplication1Dlg::OnBnClickedButton1()
{
	// 生成密钥

	
	int bit = 768;
	CButton* radio;
	radio = (CButton*)GetDlgItem(radio256);
	if (radio->GetState())bit = 256;
	radio = (CButton*)GetDlgItem(radio512);
	if (radio->GetState())bit = 512;
	radio = (CButton*)GetDlgItem(radio768);
	if (radio->GetState())bit = 768;
	radio = (CButton*)GetDlgItem(radio1024);
	if (radio->GetState())bit = 1024;
	radio = (CButton*)GetDlgItem(radio1536);
	if (radio->GetState())bit = 1536;
	radio = (CButton*)GetDlgItem(radio2048);
	if (radio->GetState())bit = 2048;
	rsa.keygen(bit);
	CEdit* pBoxOne;
	pBoxOne = (CEdit*)GetDlgItem(EDIT0);	
	pBoxOne->SetWindowText(rsa.publish().data());
	string s = "生成时间：" + to_string(rsa.keygenTime) + "秒";
	CStatic *cs = (CStatic*)GetDlgItem(static_time);
	cs->SetWindowText(s.data());
	
	//pthread_t tids[10];
	//int t;
	//for (int i = 0; i < 10; ++i)
	//{
	//	//参数依次是：创建的线程id，线程参数，调用的函数，传入的函数参数
	//	int ret = pthread_create(&tids[i], NULL, add, (void*)t);
	//	if (ret != 0)
	//	{
	//		cout << "pthread_create error: error_code=" << ret << endl;
	//	}
	//}
	////等各个线程退出后，进程才结束，否则进程强制结束了，线程可能还没反应过来；
	//pthread_exit(NULL);
}
