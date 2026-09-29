struct mint{
static inline unsigned M=998244353;
static inline ull R=-1ull/M;
static void setM(unsigned m){assert(1<m&&m<=INT_MAX);M=m;R=-1ull/M;}
unsigned x=0;
mint()=default;
template<typename T> mint(T v){
	if constexpr(is_signed_v<T>){ll y=v%(ll)M;x=y<0?y+M:y;}
	else x=v%M;
}
mint &operator+=(mint a){if((x+=a.x)>=M)x-=M;return *this;}
mint &operator-=(mint a){if((x-=a.x)>=M)x+=M;return *this;}
mint &operator*=(mint a){
	ull v=(ull)x*a.x,q=(__uint128_t)R*v>>64;
	v-=q*M;x=v>=M?v-M:v;return *this;
}
mint &operator/=(mint a){return *this*=a.inv();}
friend mint operator+(mint a,mint b){return a+=b;}
friend mint operator-(mint a,mint b){return a-=b;}
friend mint operator*(mint a,mint b){return a*=b;}
friend mint operator/(mint a,mint b){return a/=b;}
bool operator==(mint a)const{return x==a.x;}
bool operator!=(mint a)const{return x!=a.x;}
explicit operator bool()const{return x!=0;}
bool operator!()const{return !x;}
mint operator+()const{return *this;}
mint operator-()const{mint a;a.x=x?M-x:0;return a;}
mint &operator++(){return *this+=1;}
mint &operator--(){return *this-=1;}
mint operator++(int){mint a=*this;++*this;return a;}
mint operator--(int){mint a=*this;--*this;return a;}
mint inv()const{
	ll a=M,b=x,y=0,z=1;
	while(b){ll q=a/b;a-=q*b;swap(a,b);y-=q*z;swap(y,z);}
	assert(a==1);return y;
}
mint pow(ll k)const{
	ull e=k;mint a=*this,b=1;
	if(k<0){a=a.inv();e=0-e;}
	for(;e;e>>=1,a*=a)if(e&1)b*=a;
	return b;
}
friend istream &operator>>(istream &o,mint &a){
	ll x;if(o>>x)a=x;return o;
}
friend ostream &operator<<(ostream &o,mint a){return o<<a.x;}
};

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
