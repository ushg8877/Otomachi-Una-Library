////////////////////////////////////////////////////////////////
//
// template for Hash, prepared by Otomachi_Una
// version 1.2 (Last Update Jun 13th, 2026)
//
// usage:
//   init_hash(); auto H = Hash_of(" abc"); Hash sub = H[r]-H[l-1];
//   sub.v() for compare; random bases each run
//
////////////////////////////////////////////////////////////////
mt19937 rnd(time(0));
const int MOD1 = 998244353;
const int MOD2 = 1e9 + 7;
const int B1 = 256 + rnd() % (MOD1-257);
const int B2 = 256 + rnd() % (MOD2-257);

vector<ll> pw1{1}, pw2{1};
void init_hash(int n=0);

struct Hash {
	int l=0;
	ll h1=0,h2=0;
	ll v() const { return h1 << 30 | h2; }
};

inline Hash operator+(const Hash &x, const Hash &y) {
	// append y to end of x
	assert(x.l>=0&&y.l>=0&&x.l<=INT_MAX-y.l);init_hash(y.l);
	return {x.l + y.l,
			(x.h1 * pw1[y.l] + y.h1) % MOD1,
			(x.h2 * pw2[y.l] + y.h2) % MOD2};
}

inline Hash operator-(const Hash &x, const Hash &y) {
	// y is prefix of x, delete y from x
	assert(0<=y.l&&y.l<=x.l);init_hash(x.l-y.l);
	return {x.l - y.l,
			(x.h1 - y.h1 * pw1[x.l - y.l] % MOD1 + MOD1) % MOD1,
			(x.h2 - y.h2 * pw2[x.l - y.l] % MOD2 + MOD2) % MOD2};
}

template<typename T>
inline Hash operator+(const Hash &h, const T &x) {
	// append x to end of h
	assert(0<=h.l&&h.l<INT_MAX);
	return {h.l + 1,
			(h.h1 * B1 + x) % MOD1,
			(h.h2 * B2 + x) % MOD2};
}

template<typename T>
inline void operator+=(Hash &h, const T &x) {
	// append x to end of h
	h = h + x;
}

void init_hash(int n) {
	assert(0<=n&&n<INT_MAX);
	int m=pw1.size();if(n<m)return;
	n=max(n,(int)min(2ll*m,(ll)INT_MAX-1));
	pw1.resize(n+1);pw2.resize(n+1);
	for(int i=m;i<=n;i++){
		pw1[i]=pw1[i-1]*B1%MOD1;pw2[i]=pw2[i-1]*B2%MOD2;
	}
}

vector<Hash> Hash_of(const string &s) {
	int n = (int)s.length() - 1;
	assert(!s.empty()&&s[0]==' '&&s.size()<=INT_MAX);init_hash(n);
	vector<Hash> a(n+1);
	for (int i = 1; i <= n; ++i)
		a[i] = a[i-1] + (unsigned char)s[i];
	return a;
}
// end for string/rolling-hash.cpp
/////////////////////////
// !!!!! Call init_hash() first; hashing has collision risk; rnd may conflict with other templates. !!!!
