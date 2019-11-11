
// MFCApplication1Dlg.h : 头文件
//

#pragma once
#include<stdio.h>
#include<cmath>
#include<string.h>
#include<algorithm>
#include<vector>
#include<iostream>
#include<queue>
#include<windows.h>
#include<functional>
#include<complex>
#include<ctime>
#include<random>
using namespace std;
const unsigned int MAXN = 65536;
//const unsigned int MAXN=10000;
//const unsigned int MAXN=10;
const unsigned int MAXSIZE = 1025;// 1024:192 

class BI {
public:
	struct cp {
		double r, i;
		cp(double _r = 0.0, double _i = 0.0) { r = _r; i = _i; }
		cp operator+(const cp &b) { return cp(r + b.r, i + b.i); }
		cp operator-(const cp &b) { return cp(r - b.r, i - b.i); }
		cp operator*(const cp &b) { return cp(r*b.r - i*b.i, r*b.i + i*b.r); }
	};
	unsigned int a[MAXSIZE];
	int len;
	static int fac;
	static double pi;
	static int G, NTTP;

	BI() { len = 1, memset(a, 0, sizeof(a)); }
	BI(const unsigned int);
	BI(const BI &);
	//    BI(const char*);

	BI &operator =(const BI &);
	BI operator +(const BI &) const;
	BI operator -(const BI &) const;
	BI operator *(const BI &) const;
	BI operator /(const BI &) const;
	BI operator /(const unsigned int &) const;
	BI operator %(const BI &) const;
	BI sub(const int l, const int r);
	BI head(const int k);
	BI tail(const int k);
	bool operator <(const BI &) const;
	bool operator <(const unsigned int &b) const;
	bool operator >(const BI &) const;
	bool operator >(const unsigned int &b) const;
	bool operator ==(const BI &) const;
	void print()const;
	void printN()const;
	void print10()const;
	void fix(int lim);
	void fix999(int lim);// 999..进位成1 
	string tostring();
	static BI pow(BI, const BI &, const BI &);
	static BI powNewton(BI, const BI &, const BI &, const BI &);
	static BI newton(const BI &, int bit);
	static BI divNewton(const BI &, const BI &ni);
	static BI modNewton(const BI &, const BI &, const BI &ni);
	static BI powNewtonThread(BI, const BI &, const BI &, const BI &,const int&);
	static BI newtonThread(const BI &, int bit, const int&);
	static BI divNewtonThread(const BI &, const BI &ni, const int&);
	static BI modNewtonThread(const BI &, const BI &, const BI &ni, const int&);
	static BI exgcd(const BI, BI&, const BI, BI&, int &, int &);
	static int power(int x, int n, int P);
	static int rev(int id, int len);
	//    static void FFT(cp *s,int n,int inv);
	static void rader(int* y, int len);
	static void NTT(int* y, int len, int op, int P);
	//static void FFT(complex<double> *, int, int, int *);
	static int __gcd(int a, int b) { return (b) ? __gcd(b, a%b) : a; }
	//friend ostream& operator <<(ostream&, BI&);
};
class RSA {
public:
	struct threadmessage {
		int bit;
		int tid;
		threadmessage() { ; }
		threadmessage(int _bit,int _tid):bit(_bit),tid(_tid) { ; }
	};
	BI p, q;
	BI d, e;
	BI n, phi;
	//vector<unsigned int>littlePrime;
	static vector<unsigned int>littlePrime;
	bool ok[MAXN + 5];
	static int fac[8];
	static bool threadok[8];
	static BI threadPrimes[8];
	static bool stop1,stop2;
	static double keygenTime;
	RSA();
	void keygen(int bit);
	string encode(string s);
	string decode(string s);
	string b2hex(BI s);
	string publish();
	BI hex2b(string s);
	BI primeGen(int bit);// prime div & miller rabin
	static void* primeGenThread1(void* mmsg);
	static void* primeGenThread2(void* mmsg);
};
// CMFCApplication1Dlg 对话框
class CMFCApplication1Dlg : public CDialogEx
{
// 构造
public:
	RSA rsa;
	CMFCApplication1Dlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MFCAPPLICATION1_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
	BOOL CMFCApplication1Dlg::PreTranslateMessage(MSG* pMsg)
	{
		if (pMsg->message == WM_KEYDOWN && pMsg->wParam == VK_ESCAPE) return TRUE;
		if (pMsg->message == WM_KEYDOWN && pMsg->wParam == VK_RETURN) return TRUE;
		else
			return CDialog::PreTranslateMessage(pMsg);
	}
public:
	afx_msg void OnBnClickedButton2();
	afx_msg void OnBnClickedButton3();
	afx_msg void OnBnClickedButton1();
};
