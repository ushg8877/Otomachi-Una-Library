////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/polynomial/mtt.cpp
//
// usage: mint::setM(mod); poly h=f*g;
//   Inv(f); Ln(f); Exp(f); init(n); C(n,k);
//
////////////////////////////////////////////////////////////////
template<unsigned _M>
struct ModInt{
static_assert(_M==0||(1<_M&&_M<=INT_MAX));
static inline conditional_t<_M==0,unsigned,const unsigned>MOD=_M?_M:1000000007;
static inline unsigned long long NEG_INV_M=-1ULL/MOD;
static void setM(unsigned m);
unsigned x;
constexpr ModInt():x(0){}
constexpr ModInt(unsigned y):x(y%MOD){}
constexpr ModInt(int y):x((y%=static_cast<int>(MOD))<0?y+MOD:y){}
constexpr ModInt(unsigned long long y):x(y%MOD){}
constexpr ModInt(long long y):x((y%=static_cast<long long>(MOD))<0?y+MOD:y){}
ModInt &operator+=(const ModInt &a){x=(((x+=a.x)>=MOD)?x-MOD:x);return *this;}
ModInt &operator-=(const ModInt &a){x=(((x-=a.x)>=MOD)?x+MOD:x);return *this;}
ModInt &operator*=(const ModInt &a){
	unsigned long long v=(unsigned long long)x*a.x;
	if constexpr(_M)x=v%MOD;
	else{
		unsigned long long q=(__uint128_t)NEG_INV_M*v>>64,r=v-q*MOD;
		x=r-MOD*(r>=MOD);
	}
	return *this;
}
ModInt &operator/=(const ModInt &a){return (*this)*=a.inv();}
bool operator==(const ModInt &a)const{return x==a.x;}
bool operator!=(const ModInt &a)const{return x!=a.x;}
explicit operator bool()const{return x!=0;}
bool operator!()const{return x==0;}
ModInt inv()const{
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
ModInt pow(long long k)const{
	unsigned long long e=k;
	ModInt a=*this,b=1;
	if(k<0){a=a.inv();e=0-e;}
	for(;e;e>>=1){if(e&1)b*=a;a*=a;}
	return b;
}
ModInt &operator++(){*this+=1;return *this;}
ModInt &operator--(){*this-=1;return *this;}
ModInt operator+()const{return *this;}
ModInt operator-()const{ModInt a;a.x=(x?MOD-x:0u);return a;}
ModInt operator+(const ModInt &a)const{return ModInt(*this)+=a;}
ModInt operator-(const ModInt &a)const{return ModInt(*this)-=a;}
ModInt operator*(const ModInt &a)const{return ModInt(*this)*=a;}
ModInt operator/(const ModInt &a)const{return ModInt(*this)/=a;}

template<typename T>
ModInt friend operator+(T a,const ModInt &b){return ModInt(a)+=b;}
template<typename T>
ModInt friend operator-(T a,const ModInt &b){return ModInt(a)-=b;}
template<typename T>
ModInt friend operator*(T a,const ModInt &b){return ModInt(a)*=b;}
template<typename T>
ModInt friend operator/(T a,const ModInt &b){return ModInt(a)/=b;}
friend istream&operator>>(istream&i,ModInt &x){
	ll y;
	if(i>>y)x=ModInt(y);
	return i;
}
friend ostream&operator<<(ostream&o,const ModInt &x){o<<x.x;return o;}
};
////////////////////////////////////////////////////////////////////
// Basic function, fac, ifac, binom
// init(n) prepares fac / ifac / inv through n
using mint=ModInt<0>;
const unsigned &MOD=mint::MOD;
using poly=vector<mint>;
using Poly=poly;
vector<mint>fac{1},ifac{1},inv{0,1};
// 2<=mod<=INT_MAX; clears fac / ifac / inv and invalidates old values.
template<unsigned _M>
void ModInt<_M>::setM(unsigned m){
	static_assert(_M==0,"only the target modulus can change");
	assert(1<m&&m<=INT_MAX);
	MOD=m;
	NEG_INV_M=-1ULL/m;
	fac.assign(1,1);
	ifac.assign(1,1);
	::inv={0,1};
}
// Prepare factorials through n; requires gcd(n!,MOD)=1.
void init(int n=0){
	assert(0<=n&&(unsigned)n<MOD);
	int m=fac.size();if(n<m)return;
	fac.resize(n+1);
	ifac.resize(n+1);
	inv.resize(max(n+1,2));
	for(int i=m;i<=n;i++)fac[i]=fac[i-1]*i;
	ifac[n]=fac[n].inv();
	for(int i=n;i>=m;i--){ifac[i-1]=ifac[i]*i;inv[i]=ifac[i]*fac[i-1];}
}
inline mint C(int x,int y){
	// choose x from y
	if(x<0||y<x)return 0;
	init(y);
	return fac[y]*ifac[x]*ifac[y-x];
}
inline mint binom(int y,int x){return C(x,y);}

////////////////////////////////////////////////////////////////
//
// template for MTT (three-prime CRT)
////////////////////////////////////////////////////////////////
const int FFT_MAX=23;
template<unsigned MO>
struct MTT_NTT{
using num=ModInt<MO>;
static constexpr unsigned MO2=MO*2;
vector<num>FFT_RATIOS,INV_FFT_RATIOS;
MTT_NTT():FFT_RATIOS(FFT_MAX,1),INV_FFT_RATIOS(FFT_MAX,1){
	static_assert(MO<(1u<<30)&&((MO-1)%(1<<FFT_MAX)==0));
	vector<num>root(FFT_MAX+1),iroot(FFT_MAX+1);
	root[FFT_MAX]=num(3).pow((MO-1)>>FFT_MAX);
	iroot[FFT_MAX]=root[FFT_MAX].inv();
	for(int i=FFT_MAX;i;i--)root[i-1]=root[i]*root[i],
		iroot[i-1]=iroot[i]*iroot[i];
	num x=1,y=1;
	for(int i=0;i<FFT_MAX-1;i++){
		FFT_RATIOS[i]=root[i+2]*x;INV_FFT_RATIOS[i]=iroot[i+2]*y;
		x*=iroot[i+2];y*=root[i+2];
	}
}
// as[rev(i)] <- \sum_j \zeta^(ij) as[j]
void fft(num *as,int n)const{
	assert(!(n&(n-1)));
	assert(1<=n);
	assert(n<=1<<FFT_MAX);
	int m=n;
	if(m>>=1){
		for(int i=0;i<m;++i){
			const unsigned x=as[i+m].x;
			as[i+m].x=as[i].x+MO-x;
			as[i].x+=x;
		}
	}
	if(m>>=1){
		num prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+m;++i){
				const unsigned x=(prod*as[i+m]).x;
				as[i+m].x=as[i].x+MO-x;
				as[i].x+=x;
			}
			prod*=FFT_RATIOS[__builtin_ctz(++h)];
		}
	}
	for(;m;){
		if(m>>=1){
			num prod=1U;
			for(int h=0,i0=0;i0<n;i0+=(m<<1)){
				for(int i=i0;i<i0+m;++i){
					const unsigned x=(prod*as[i+m]).x;
					as[i+m].x=as[i].x+MO-x;
					as[i].x+=x;
				}
				prod*=FFT_RATIOS[__builtin_ctz(++h)];
			}
		}
		if(m>>=1){
			num prod=1U;
			for(int h=0,i0=0;i0<n;i0+=(m<<1)){
				for(int i=i0;i<i0+m;++i){
					const unsigned x=(prod*as[i+m]).x;
					as[i].x=(as[i].x>=MO2)?(as[i].x-MO2):as[i].x;
					as[i+m].x=as[i].x+MO-x;
					as[i].x+=x;
				}
				prod*=FFT_RATIOS[__builtin_ctz(++h)];
			}
		}
	}
	for(int i=0;i<n;++i){
		as[i].x=(as[i].x>=MO2)?(as[i].x-MO2):as[i].x;
		as[i].x=(as[i].x>=MO)?(as[i].x-MO):as[i].x;
	}
}

// as[i] <- (1/n) \sum_j \zeta^(-ij) as[rev(j)]
void invFft(num *as,int n)const{
	assert(!(n&(n-1)));
	assert(1<=n);
	assert(n<=1<<FFT_MAX);
	int m=1;
	if(m<n>>1){
		num prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+m;++i){
				const unsigned long long y=as[i].x+MO-as[i+m].x;
				as[i].x+=as[i+m].x;
				as[i+m].x=(prod.x*y)%MO;
			}
			prod*=INV_FFT_RATIOS[__builtin_ctz(++h)];
		}
		m<<=1;
	}
	for(;m<n>>1;m<<=1){
		num prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+(m>>1);++i){
				const unsigned long long y=as[i].x+MO2-as[i+m].x;
				as[i].x+=as[i+m].x;
				as[i].x=(as[i].x>=MO2)?(as[i].x-MO2):as[i].x;
				as[i+m].x=(prod.x*y)%MO;
			}
			for(int i=i0+(m>>1);i<i0+m;++i){
				const unsigned long long y=as[i].x+MO-as[i+m].x;
				as[i].x+=as[i+m].x;
				as[i+m].x=(prod.x*y)%MO;
			}
			prod*=INV_FFT_RATIOS[__builtin_ctz(++h)];
		}
	}
	if(m<n){
		for(int i=0;i<m;++i){
			const unsigned y=as[i].x+MO2-as[i+m].x;
			as[i].x+=as[i+m].x;
			as[i+m].x=y;
		}
	}
	const num invN=num(n).inv();
	for(int i=0;i<n;++i){
		as[i]*=invN;
	}
}
};

template<unsigned MO>
vector<ModInt<MO>>convolution(const poly &f,const poly &g,int n){
	static const MTT_NTT<MO>ntt;
	vector<ModInt<MO>>a(n),b(n);
	for(int i=0;i<(int)f.size();i++)a[i]=f[i].x;
	for(int i=0;i<(int)g.size();i++)b[i]=g[i].x;
	ntt.fft(a.data(),n);ntt.fft(b.data(),n);
	for(int i=0;i<n;i++)a[i]*=b[i];
	ntt.invFft(a.data(),n);
	a.resize(f.size()+g.size()-1);
	return a;
}
// Composite moduli allowed; convolution length <=2^23.
poly operator*(const poly &f,const poly &g){
	if(f.empty()||g.empty())return{};
	assert(f.size()+g.size()-1<=(1<<FFT_MAX));
	int m=f.size()+g.size()-1;
	poly h(m);
	if(min(f.size(),g.size())<=32){
		for(int i=0;i<(int)f.size();i++)for(int j=0;j<(int)g.size();j++)
			h[i+j]+=f[i]*g[j];
		return h;
	}
	constexpr unsigned mod1=998244353,mod2=167772161,mod3=469762049;
	constexpr unsigned long long mod12=1ull*mod1*mod2;
	static const unsigned inv1=ModInt<mod2>(mod1).inv().x;
	static const unsigned inv2=ModInt<mod3>(mod12).inv().x;
	assert((__int128)min(f.size(),
		g.size())*(MOD-1)*(MOD-1)<(__int128)mod12*mod3);
	int n=1;
	while(n<m)n<<=1;
	auto a=convolution<mod1>(f,g,n);
	auto b=convolution<mod2>(f,g,n);
	auto c=convolution<mod3>(f,g,n);
	for(int i=0;i<m;i++){
		unsigned long long x=(b[i].x+mod2-a[i].x%mod2)*1ull*inv1%mod2;
		x=x*mod1+a[i].x;
		unsigned long long y=(c[i].x+mod3-x%mod3)*inv2%mod3;
		h[i]=x%MOD+(mod12%MOD)*y%MOD;
	}
	return h;
}
// Transposed multiplication: (f/g)[i] = sum_j f[i+j]*g[j].
poly operator/(const poly &f,poly g){
	if(f.empty()||g.empty())return poly(f.size());
	int m=g.size();
	reverse(g.begin(),g.end());
	g=f*g;
	return poly(g.begin()+m-1,g.end());
}
void operator+=(poly &f,const poly &g){
	if(f.size()<g.size())f.resize(g.size());
	for(int i=0;i<(int)g.size();i++)f[i]+=g[i];
}
poly operator+(poly f,const poly &g){f+=g;return f;}
void operator-=(poly &f,const poly &g){
	if(f.size()<g.size())f.resize(g.size());
	for(int i=0;i<(int)g.size();i++)f[i]-=g[i];
}
poly operator-(poly f,const poly &g){f-=g;return f;}
poly operator*(poly f,mint x){for(auto&y:f)y*=x;return f;}
poly operator*(mint x,poly f){for(auto&y:f)y*=x;return f;}
mint value(const poly &f,mint x){
	mint ans=0;
	for(int i=(int)f.size()-1;i>=0;i--)ans=ans*x+f[i];
	return ans;
}
// Formal inverse; constant term must be invertible.
// FPS length <=2^22; divisions require units modulo MOD.
poly Inv(const poly &f){
	assert(!f.empty()&&f[0]&&f.size()<=(1<<(FFT_MAX-1)));
	int n=f.size();poly g{f[0].inv()};
	while((int)g.size()<n){
		int m=min(n,(int)g.size()*2);
		poly h=poly(f.begin(),f.begin()+m)*g;h.resize(m);
		for(auto&x:h)x=-x;
		h[0]+=2;
		g=g*h;
		g.resize(m);
	}
	return g;
}
poly diff(poly f){
	if(f.empty())return{};
	for(int i=1;i<(int)f.size();i++)f[i-1]=f[i]*i;
	f.pop_back();
	return f;
}
poly integ(poly f){
	int n=f.size();
	init(n);
	f.resize(n+1);
	for(int i=n;i;i--)f[i]=f[i-1]*inv[i];
	f[0]=0;
	return f;
}
// Formal logarithm; f[0]=1. Integration denominators must be invertible.
poly Ln(const poly &f){
	assert(!f.empty()&&f[0].x==1&&f.size()<(size_t)MOD&&
		f.size()<=(1<<(FFT_MAX-1)));
	poly g=diff(f)*Inv(f);
	g.resize(f.size()-1);
	return integ(move(g));
}
// Formal exponential; f[0]=0. Integration denominators must be invertible.
poly Exp(const poly &f){
	assert(!f.empty()&&!f[0]&&f.size()<(size_t)MOD&&f.size()<=(1<<(FFT_MAX-1)));
	int n=f.size();poly g{1};
	while((int)g.size()<n){
		int m=min(n,(int)g.size()*2);
		poly h=g;
		h.resize(m);
		h=Ln(h);
		for(int i=0;i<m;i++)h[i]=f[i]-h[i];
		h[0]+=1;
		g=g*h;
		g.resize(m);
	}
	return g;
}
// Product tree shared by evaluation and interpolation; coefficients low first.
struct Poly_Tree{
private:
int n=0;
poly x,ig;
vector<poly> g;
void build(int u,int l,int r){
	if(l==r){g[u]={-x[l],1};return;}
	int m=(l+r)>>1,v=u+2*(m-l+1);
	build(u+1,l,m);build(v,m+1,r);
	g[u]=g[u+1]*g[v];
}
// First k coefficients of the correlation with reverse(b).
static poly down(const poly &a,const poly &b,int k){
	k=min(k,(int)a.size());
	int m=b.size()-1;
	if(k<=32){
		poly c(k);
		for(int i=0;i<k;i++)
			for(int j=0;j<=m&&i+j<(int)a.size();j++)c[i]+=a[i+j]*b[m-j];
		return c;
	}
	poly c=a*b;
	return poly(c.begin()+m,c.begin()+m+k);
}
// Keep the original transform limit when the root correlation is too large.
static poly rem(poly a,const poly &b){
	if(a.size()<b.size())return a;
	int k=a.size()-b.size()+1;
	poly q(a.rbegin(),a.rbegin()+k),v(b.rbegin(),b.rend());
	v.resize(k);q=q*Inv(v);q.resize(k);
	reverse(q.begin(),q.end());q=q*b;a.resize(b.size()-1);
	for(int i=0;i<(int)a.size();i++)a[i]-=q[i];
	return a;
}
void eval_rem(poly a,int u,int l,int r,poly &ans)const{
	if(r-l<32){
		for(int i=l;i<=r;i++)ans[i]=value(a,x[i]);
		return;
	}
	int m=(l+r)>>1,v=u+2*(m-l+1);
	eval_rem(rem(a,g[u+1]),u+1,l,m,ans);
	eval_rem(rem(move(a),g[v]),v,m+1,r,ans);
}
void eval(poly a,int u,int l,int r,poly &ans)const{
	if(a.empty())return;
	if(l==r){ans[l]=a[0];return;}
	int m=(l+r)>>1,v=u+2*(m-l+1);
	eval(down(a,g[v],m-l+1),u+1,l,m,ans);
	eval(down(a,g[u+1],r-m),v,m+1,r,ans);
}
poly join(const poly &a,int u,int l,int r)const{
	if(l==r)return {a[l]};
	int m=(l+r)>>1,v=u+2*(m-l+1);
	return join(a,u+1,l,m)*g[v]+join(a,v,m+1,r)*g[u+1];
}
public:
Poly_Tree()=default;
explicit Poly_Tree(poly a){set(move(a));}
// Replace the evaluation points; empty input clears the tree and cache.
void set(poly a={}){
	assert(a.size()<(1u<<23));
	x=move(a);n=x.size();ig.clear();g.clear();
	if(n){g.resize(2*n-1);build(0,0,n-1);}
}
// f(x[i]); repeated points allowed. O(M(m)+M(n)log(n)), m=f.size().
// One root inverse, then transposed products; cache reused for fixed points.
poly eval(const poly &f){
	poly ans(n);
	if(!n||f.empty())return ans;
	if(n<=32||f.size()<=32){
		for(int i=0;i<n;i++)ans[i]=value(f,x[i]);
		return ans;
	}
	int m=f.size();
	if(m>(1<<22)){eval_rem(f,0,0,n-1,ans);return ans;}
	if((int)ig.size()<m){
		poly a(g[0].rbegin(),g[0].rend());
		a.resize(m);ig=Inv(a);
	}
	poly a(ig.begin(),ig.begin()+m);reverse(a.begin(),a.end());
	eval(down(f,a,n),0,0,n-1,ans);
	return ans;
}
// n coefficients interpolating (x[i],y[i]); pairwise differences must be units.
// O(n log^2 n) time, O(n log n) space; empty input returns empty.
poly interpolate(const poly &y){
	assert(y.size()==x.size());
	if(!n)return{};
	poly d=eval(diff(g[0])),pre(n+1,1);
	for(int i=0;i<n;i++){assert(d[i]);pre[i+1]=pre[i]*d[i];}
	mint v=pre[n].inv();
	for(int i=n-1;i>=0;i--){
		mint t=v*pre[i];v*=d[i];d[i]=y[i]*t;
	}
	return join(d,0,0,n-1);
}
};
// One-shot evaluation; use Poly_Tree for repeated queries at the same points.
poly Eval(poly f,poly x){
	if(x.size()<=32||f.size()<=32){
		poly ans(x.size());
		for(int i=0;i<(int)x.size();i++)ans[i]=value(f,x[i]);
		return ans;
	}
	return Poly_Tree(move(x)).eval(f);
}

// Shortest order d=c.size()-1: c[0]=1, sum c[j]*a[i-j]=0 for i>=d.
// O(n*(d+1)) time, O(d+1) space; discrepancies must be invertible.
poly BM(const poly &a){
	poly c{1},b{1};
	int len=0,shift=1;
	mint inv_d=1;
	for(int i=0;i<(int)a.size();i++){
		mint d=a[i];
		for(int j=1;j<=len;j++)d+=c[j]*a[i-j];
		if(!d){shift++;continue;}
		bool grow=2*len<=i;
		poly old;
		if(grow)old=c;
		mint v=d*inv_d;
		if(c.size()<b.size()+shift)c.resize(b.size()+shift);
		for(int j=0;j<(int)b.size();j++)c[j+shift]-=v*b[j];
		if(grow){
			len=i+1-len;b=move(old);inv_d=d.inv();shift=1;
		}else shift++;
	}
	c.resize(len+1);
	return c;
}


mint FSPE(poly F,poly G,ll t){
	// find [x^t] F/G
	assert(t>=0&&!G.empty()&&G[0]);
	if(F.empty())return 0;
	while(t){
		poly G1=G;
		for(int i=1;i<(int)G1.size();i+=2)G1[i]=-G1[i];
		F=F*G1;G=G*G1;
		poly f,g;
		for(int i=t&1;i<(int)F.size();i+=2)f.push_back(F[i]);
		for(int i=0;i<(int)G.size();i+=2)g.push_back(G[i]);
		t>>=1;
		F=f;G=g;
		if(F.empty())return 0;
	}
	return F[0]/G[0];
}

mint RSPE(poly f,ll x){
	// find f[x] in O(n^2), note that f must be recursion
	poly g=BM(f);
	poly s(g.size()-1);
	for(int i=0;i<(int)s.size();i++)
		for(int j=0;j<=i&&j<(int)f.size();j++)
			s[i]+=f[j]*g[i-j];
	return FSPE(s,g,x);
}


// Interpolate arbitrary points; return n coefficients, degree < n.
// O(n log^2 n) time, O(n log n) space. Pairwise x differences must be units.
poly lagrange(const vector<mint>&x,const vector<mint>&y){
	return Poly_Tree(x).interpolate(y);
}
// end for polynomial/mtt.cpp
/////////////////////////
// !!!!! Call mint::setM(mod) before creating values; changing mod
// invalidates them. !!!!
