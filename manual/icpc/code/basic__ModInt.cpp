template<unsigned M>
struct ModInt{
static constexpr unsigned MOD=M;
static_assert(1<M&&M<=INT_MAX);
unsigned x=0;
ModInt()=default;
template<typename T> ModInt(T v){
	if constexpr(is_signed_v<T>){ll y=v%(ll)M;x=y<0?y+M:y;}
	else x=v%M;
}
ModInt &operator+=(ModInt a){if((x+=a.x)>=M)x-=M;return *this;}
ModInt &operator-=(ModInt a){if((x-=a.x)>=M)x+=M;return *this;}
ModInt &operator*=(ModInt a){x=(ull)x*a.x%M;return *this;}
ModInt &operator/=(ModInt a){return *this*=a.inv();}
friend ModInt operator+(ModInt a,ModInt b){return a+=b;}
friend ModInt operator-(ModInt a,ModInt b){return a-=b;}
friend ModInt operator*(ModInt a,ModInt b){return a*=b;}
friend ModInt operator/(ModInt a,ModInt b){return a/=b;}
bool operator==(ModInt a)const{return x==a.x;}
bool operator!=(ModInt a)const{return x!=a.x;}
explicit operator bool()const{return x!=0;}
bool operator!()const{return !x;}
ModInt operator+()const{return *this;}
ModInt operator-()const{ModInt a;a.x=x?M-x:0;return a;}
ModInt &operator++(){return *this+=1;}
ModInt &operator--(){return *this-=1;}
ModInt operator++(int){ModInt a=*this;++*this;return a;}
ModInt operator--(int){ModInt a=*this;--*this;return a;}
ModInt inv()const{
	ll a=M,b=x,y=0,z=1;
	while(b){ll q=a/b;a-=q*b;swap(a,b);y-=q*z;swap(y,z);}
	assert(a==1);return y;
}
ModInt pow(ll k)const{
	ull e=k;ModInt a=*this,b=1;
	if(k<0){a=a.inv();e=0-e;}
	for(;e;e>>=1,a*=a)if(e&1)b*=a;
	return b;
}
friend istream &operator>>(istream &o,ModInt &a){
	ll x;if(o>>x)a=x;return o;
}
friend ostream &operator<<(ostream &o,ModInt a){return o<<a.x;}
};

const int MOD=998244353;
using mint=ModInt<MOD>;
vector<mint> fac{1},ifac{1},inv{0,1};
void init(int n=0){
	assert(0<=n&&n<MOD);
	int m=fac.size();if(n<m)return;
	n=max(n,(int)min(2ll*m,(ll)MOD-1));
	fac.resize(n+1);ifac.resize(n+1);inv.resize(max(n+1,2));
	for(int i=m;i<=n;i++){
		fac[i]=fac[i-1]*i;
		if(i>1)inv[i]=-(MOD/i)*inv[MOD%i];
		ifac[i]=ifac[i-1]*inv[i];
	}
}
inline mint C(int x,int y){
	// choose x from y
	if(x<0||y<x)return 0;
	init(y);return fac[y]*ifac[x]*ifac[y-x];
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
