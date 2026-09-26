////////////////////////////////////////////////////////////////
//
// template for ModInt (Mod isn't assign)
// version 1.0 (Last Update Jun 16th, 2026)
//
// usage:
//   mint::setM(mod); init_fact(n); mint x; cin >> x;
//
////////////////////////////////////////////////////////////////
struct mint {
	static unsigned M;
	static unsigned long long NEG_INV_M;
	static void setM(unsigned m) {
		assert(1<m&&m<=INT_MAX);M = m;
		NEG_INV_M = -1ULL / M;
	}
	unsigned x;
	mint() : x(0U) {}
	mint(unsigned x_) : x(x_ % M) {}
	mint(unsigned long long x_) : x(x_ % M) {}
	mint(int x_) : x(((x_ %= static_cast<int>(M)) < 0) ? (x_ + static_cast<int>(M)) : x_) {}
	mint(long long x_) : x(((x_ %= static_cast<long long>(M)) < 0) ? (x_ + static_cast<long long>(M)) : x_) {}
	mint& operator+=(const mint &a) {
		x = ((x += a.x) >= M) ? (x - M) : x;
		return *this;
	}
	mint& operator-=(const mint &a) {
		x = ((x -= a.x) >= M) ? (x + M) : x;
		return *this;
	}
	mint& operator*=(const mint &a) {
		const unsigned long long y = static_cast<unsigned long long>(x) * a.x;
		const unsigned long long q = static_cast<unsigned long long>(
			(static_cast<unsigned __int128>(NEG_INV_M) * y) >> 64
		);
		const unsigned long long r = y - M * q;
		x = r - M * (r >= M);
		return *this;
	}
	mint& operator/=(const mint &a) { return (*this *= a.inv()); }
	mint pow(long long e) const {
		unsigned long long k=e;
		mint a=*this,b=1;
		if(e<0){a=a.inv();k=0-k;}
		for(;k;k>>=1){if(k&1)b*=a;a*=a;}
		return b;
	}
	mint inv() const {
		unsigned a = M, b = x;
		ll y = 0, z = 1;
		for (; b; ) {
			const unsigned q = a / b;
			const unsigned c = a - q * b;
			a = b; b = c;
			const ll w = y - (ll)q * z;
			y = z; z = w;
		}
		assert(a == 1U);
		return mint(y);
	}
	mint operator+() const { return *this; }
	mint operator-() const { mint a; a.x = x ? (M - x) : 0U; return a;}
	mint operator+(const mint &a) const { return (mint(*this) += a); }
	mint operator-(const mint &a) const { return (mint(*this) -= a); }
	mint operator*(const mint &a) const { return (mint(*this) *= a); }
	mint operator/(const mint &a) const { return (mint(*this) /= a); }
	mint& operator++() { *this += 1; return *this; }
	mint& operator--() { *this -= 1; return *this; }
	bool operator!() const { return x == 0; }
	explicit operator bool() const { return x != 0; }
	bool operator==(const mint &a) const { return (x == a.x); }
	bool operator!=(const mint &a) const { return (x != a.x); }
	template<class T> friend mint operator+(T a, const mint &b) { return (mint(a) += b); }
	template<class T> friend mint operator-(T a, const mint &b) { return (mint(a) -= b); }
	template<class T> friend mint operator*(T a, const mint &b) { return (mint(a) *= b); }
	template<class T> friend mint operator/(T a, const mint &b) { return (mint(a) /= b); }
	friend istream& operator>>(istream &i, mint &a) { ll x;if(i>>x)a=mint(x);return i; }
	friend ostream& operator<<(ostream &o, const mint &a) { return o << a.x; }
};
unsigned mint::M = 998244353;
unsigned long long mint::NEG_INV_M = -1ULL / mint::M;


vector<mint> fac,ifac,inv;
unsigned fact_mod=0;
void init_fact(int n=0){
	assert(0<=n&&(unsigned)n<mint::M);
	fac.assign(n+1,1);ifac.resize(n+1);inv.assign(n+1,0);
	for(int i=1;i<=n;i++)fac[i]=fac[i-1]*i;
	ifac[n]=fac[n].inv();
	for(int i=n;i>=1;i--)ifac[i-1]=ifac[i]*i;
	for(int i=1;i<=n;i++)inv[i]=ifac[i]*fac[i-1];
	fact_mod=mint::M;
}
inline mint C(int x,int y){
	if(x<0||y<x)return 0;
	assert(fact_mod==mint::M&&y<(int)fac.size());
	return fac[y]*ifac[x]*ifac[y-x];
}
inline mint binom(int y,int x){return C(x,y);}
// end for math/dynamic-mod-int.cpp
/////////////////////////
// !!!!! Call mint::setM first; changing the modulus invalidates old values and
// factorial tables; choose only one mint definition. !!!!
