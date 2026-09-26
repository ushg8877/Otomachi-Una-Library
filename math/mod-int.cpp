////////////////////////////////////////////////////////////////
//
// template for mod-int
//
// usage: mint a,b; init(n); C(n,k); // fixed prime modulus
//
////////////////////////////////////////////////////////////////
template <unsigned _M>
struct ModInt{
static constexpr unsigned MOD=_M;
static_assert(1<MOD&&MOD<=INT_MAX);
unsigned x;
constexpr ModInt():x(0){}
constexpr ModInt(unsigned y):x(y%MOD){}
constexpr ModInt(int y):x((y%=static_cast<int>(MOD))<0?y+MOD:y){}
constexpr ModInt(unsigned long long y):x(y%MOD){}
constexpr ModInt(long long y):x((y%=static_cast<long long>(MOD))<0?y+MOD:y){}
ModInt& operator +=(const ModInt &a){x=(((x+=a.x)>=MOD)?x-MOD:x);return *this;}
ModInt& operator -=(const ModInt &a){x=(((x-=a.x)>=MOD)?x+MOD:x);return *this;}
ModInt& operator *=(const ModInt &a){
	x=static_cast<unsigned long long>(x)*a.x%MOD;
	return *this;
}
ModInt& operator /=(const ModInt &a){return (*this)*=a.inv();}
bool operator ==(const ModInt &a)const{return x==a.x;}
bool operator !=(const ModInt &a)const{return x!=a.x;}
explicit operator bool() const { return x != 0; }
bool operator!() const { return x == 0; }
ModInt inv() const {
	unsigned a=MOD,b=x;ll y=0,z=1;
	while(b){
		const unsigned q=a/b,c=a-q*b;
		a=b,b=c;
		const ll w=y-(ll)q*z;
		y=z,z=w;
	}
	assert(a==1);
	return ModInt(y);
}
ModInt pow(long long k) const {
	unsigned long long e=k;
	ModInt a=*this,b=1;
	if(k<0){a=a.inv();e=0-e;}
	for(;e;e>>=1){if(e&1)b*=a;a*=a;}
	return b;
}
ModInt& operator++() {*this+=1;return *this;}
ModInt& operator--() {*this-=1;return *this;}
ModInt operator--(int){ModInt t=*this;--*this;return t;}
ModInt operator++(int){ModInt t=*this;++*this;return t;}
ModInt operator +()const{return *this;}
ModInt operator -()const{ModInt a;a.x=(x?MOD-x:0u);return a;}
ModInt operator +(const ModInt &a)const{return ModInt(*this)+=a;}
ModInt operator -(const ModInt &a)const{return ModInt(*this)-=a;}
ModInt operator *(const ModInt &a)const{return ModInt(*this)*=a;}
ModInt operator /(const ModInt &a)const{return ModInt(*this)/=a;}

template<typename T>
ModInt friend operator +(T a,const ModInt &b){return ModInt(a)+=b;}
template<typename T>
ModInt friend operator -(T a,const ModInt &b){return ModInt(a)-=b;}
template<typename T>
ModInt friend operator *(T a,const ModInt &b){return ModInt(a)*=b;}
template<typename T>
ModInt friend operator /(T a,const ModInt &b){return ModInt(a)/=b;}
friend istream& operator >> (istream &i,ModInt &x){
	ll y;
	if(i>>y)x=ModInt(y);
	return i;
}
friend ostream& operator << (ostream &o,const ModInt &x){o<<x.x;return o;}
};
////////////////////////////////////////////////////////////////////
// Basic function, fac, ifac, binom
// init(n) prepares fac / ifac / inv through n
const int MOD=998244353;
using mint=ModInt<MOD>;
vector<mint> fac{1},ifac{1},inv{0,1};
void init(int n=0){
	assert(0<=n&&n<MOD);
	int m=fac.size();if(n<m)return;
	n=max(n,(int)min(2ll*m,(ll)MOD-1));
	fac.resize(n+1);
	ifac.resize(n+1);
	inv.resize(max(n+1,2));
	for(int i=m;i<=n;i++){
		fac[i]=fac[i-1]*i;
		if(i>1)inv[i]=-(MOD/i)*inv[MOD%i];
		ifac[i]=ifac[i-1]*inv[i];
	}
}
inline mint C(int x,int y){
	// choose x from y
	if(x<0||y<x)return 0;
	init(y);
	return fac[y]*ifac[x]*ifac[y-x];
}
inline mint binom(int y,int x){return C(x,y);}

inline mint Apple_in_Box(int n,int m){
	// n apples into m boxes, each box must be non-empty
	return C(m-1,n-1);
}
inline mint Apple_in_Box(int n,int m,int k){
	// n apples into m boxes, 
	// the first k boxes must be non-empty (others can be empty)
	assert(0<=k&&k<=m);
	return Apple_in_Box(n+(m-k),m);
}
// end for math/mod-int.cpp
/////////////////////////
// !!!!! Division needs an invertible divisor; factorials need a prime
// modulus. !!!!
