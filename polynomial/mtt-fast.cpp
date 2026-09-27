////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/polynomial/mtt-fast.cpp
//
// usage: mint::setM(mod); poly h=f*g;
//   Inv(f); Ln(f); Exp(f); init(n); C(k,n);
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
poly mtt_exact(const poly &f,const poly &g){
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
// Quadratic lattice lift follows the supplied cp-algorithms FFT.
#pragma GCC push_options
#pragma GCC target("avx2,fma")
#pragma GCC optimize("O3","no-unroll-loops")
struct MTT_FFT{
struct Point{
	vector<double>x,y;
	Point(int n):x(n),y(n){}
};
vector<double>w{1,1},v{0,0};
unsigned mod=0,root=0,d=0;
ll a=0,b=0,c=0,e=0;
double scale=1;
static unsigned power(unsigned x,unsigned k,unsigned p){
	unsigned y=1;
	for(;k;k>>=1,x=1ull*x*x%p)if(k&1)y=1ull*y*x%p;
	return y;
}
static unsigned sqrt(unsigned x,unsigned p){
	if(p%4==3)return power(x,(p+1)/4,p);
	unsigned q=p-1,z=2;
	int s=__builtin_ctz(q);q>>=s;
	while(power(z,(p-1)/2,p)!=p-1)z++;
	unsigned c=power(z,q,p),t=power(x,q,p),r=power(x,(q+1)/2,p);
	while(t!=1){
		int i=0;unsigned v=t;
		while(v!=1)v=1ull*v*v%p,i++;
		unsigned h=power(c,1u<<(s-i-1),p);
		r=1ull*r*h%p;c=1ull*h*h%p;t=1ull*t*c%p;s=i;
	}
	return r;
}
void setM(unsigned p){
	if(mod==p)return;
	mod=p;d=0;
	if(p<3||!(p&1))return;
	for(ll i=3;i*i<=p;i+=2)if(p%i==0)return;
	for(unsigned i=1;i<min(p,256u);i++)
		if(power(p-i,(p-1)/2,p)==1){d=i;break;}
	if(!d)return;
	root=sqrt(p-d,p);scale=std::sqrt(double(d));
	__int128 x=p,y=0,z=-(ll)root,t=1;
	auto norm=[&](__int128 x,__int128 y){return x*x+d*y*y;};
	while(true){
		if(norm(x,y)>norm(z,t))swap(x,z),swap(y,t);
		__int128 v=x*z+d*y*t,l=norm(x,y);
		__int128 k=(2*v+(v>=0?l:-l))/(2*l);
		if(!k||norm(z-k*x,t-k*y)>=norm(z,t))break;
		z-=k*x;t-=k*y;
	}
	a=x;b=y;c=z;e=t;
	assert(abs(a*e-b*c)==p);
}
void setN(int n){
	for(int k=w.size();k<n;k<<=1){
		w.resize(k*2);v.resize(k*2);
		long double step=acosl(-1)/k;
		for(int i=k/2;i<k;i++){
			w[i*2]=w[i];v[i*2]=v[i];
			long double x=step*(2*i-k+1);
			w[i*2+1]=cosl(x);v[i*2+1]=sinl(x);
		}
	}
}
// Fuse the last two forward / first two inverse levels in registers.
template<bool inv>void leaf(double *px,double *py,int n){
	for(int i=0;i<n;i+=4){
		auto x=_mm256_loadu_pd(px+i);
		auto y=_mm256_loadu_pd(py+i);
		if constexpr(inv){
			auto u=_mm256_permute4x64_pd(x,0xb1);
			auto v=_mm256_permute4x64_pd(y,0xb1);
			x=_mm256_blend_pd(x+u,u-x,0xa);
			y=_mm256_blend_pd(y+v,v-y,0xa);
			u=x;x=_mm256_blend_pd(x,y,8);y=_mm256_blend_pd(y,-u,8);
			u=_mm256_permute4x64_pd(x,0x4e);
			v=_mm256_permute4x64_pd(y,0x4e);
			x=_mm256_blend_pd(x+u,u-x,0xc);
			y=_mm256_blend_pd(y+v,v-y,0xc);
		}else{
			auto u=_mm256_permute4x64_pd(x,0x4e);
			auto v=_mm256_permute4x64_pd(y,0x4e);
			x=_mm256_blend_pd(x+u,u-x,0xc);
			y=_mm256_blend_pd(y+v,v-y,0xc);
			u=x;x=_mm256_blend_pd(x,-y,8);y=_mm256_blend_pd(y,u,8);
			u=_mm256_permute4x64_pd(x,0xb1);
			v=_mm256_permute4x64_pd(y,0xb1);
			x=_mm256_blend_pd(x+u,u-x,0xa);
			y=_mm256_blend_pd(y+v,v-y,0xa);
		}
		_mm256_storeu_pd(px+i,x);
		_mm256_storeu_pd(py+i,y);
	}
}
// DIF forward / DIT inverse. Split real/imaginary arrays avoid shuffles.
template<bool inv>void pass(double *px,double *py,int n,int k){
		for(int l=0;l<n;l+=k*2){
		int j=0;
		for(;j+3<k;j+=4){
			int i=l+j;
			auto x=_mm256_loadu_pd(px+i);
			auto y=_mm256_loadu_pd(py+i);
			auto z=_mm256_loadu_pd(px+i+k);
			auto t=_mm256_loadu_pd(py+i+k);
			auto u=_mm256_loadu_pd(w.data()+k+j);
			auto v0=_mm256_loadu_pd(v.data()+k+j);
			__m256d rx,ry,sx,sy;
			if(inv){
				auto p=_mm256_fmadd_pd(z,u,_mm256_mul_pd(t,v0));
				auto q=_mm256_fmsub_pd(t,u,_mm256_mul_pd(z,v0));
				rx=x+p;ry=y+q;sx=x-p;sy=y-q;
			}else{
				rx=x+z;ry=y+t;x-=z;y-=t;
				sx=_mm256_fmsub_pd(x,u,_mm256_mul_pd(y,v0));
				sy=_mm256_fmadd_pd(y,u,_mm256_mul_pd(x,v0));
			}
			_mm256_storeu_pd(px+i,rx);
			_mm256_storeu_pd(py+i,ry);
			_mm256_storeu_pd(px+i+k,sx);
			_mm256_storeu_pd(py+i+k,sy);
		}
	}
}
template<bool inv>void fft(Point &a){
	int n=a.x.size();setN(n);
	double *x=a.x.data(),*y=a.y.data();
	int block=min(n,8192);
	if constexpr(!inv)
		for(int k=n/2;k>=block;k>>=1)pass<false>(x,y,n,k);
	for(int i=0;i<n;i+=block){
		if constexpr(inv)leaf<true>(x+i,y+i,block);
		for(int k=inv?4:block/2;inv?k<block:k>=4;inv?k<<=1:k>>=1)
			pass<inv>(x+i,y+i,block,k);
		if constexpr(!inv)leaf<false>(x+i,y+i,block);
	}
	if constexpr(inv)
		for(int k=block;k<n;k<<=1)pass<true>(x,y,n,k);
}
static __m256d round(__m256d x){
	return _mm256_round_pd(x,_MM_FROUND_TO_NEAREST_INT|_MM_FROUND_NO_EXC);
}
// Independent stochastic rounding keeps the lift approximately centered.
void read(Point &f,const poly &g,ull seed){
	double det=double(a)*e-double(b)*c;
	auto qa=_mm256_set1_pd(e/det),qb=_mm256_set1_pd(-b/det);
	auto va=_mm256_set1_pd(a),vb=_mm256_set1_pd(b);
	auto vc=_mm256_set1_pd(c),ve=_mm256_set1_pd(e);
	auto vs=_mm256_set1_pd(scale),iv=_mm256_set1_pd(0x1p-32);
	__m256i state=_mm256_set_epi64x(seed+3,seed+2,seed+1,seed);
	int n=g.size(),i=0;
	for(;i+3<n;i+=4){
		state^=_mm256_slli_epi64(state,13);
		state^=_mm256_srli_epi64(state,7);
		state^=_mm256_slli_epi64(state,17);
		auto lo=_mm256_castsi256_si128(state);
		auto hi=_mm256_extracti128_si256(state,1);
		auto noise=_mm256_cvtepi32_pd(_mm_unpacklo_epi64(lo,hi))*iv;
		auto x=_mm256_cvtepi32_pd(_mm_loadu_si128(
			(const __m128i*)(g.data()+i)));
		auto q=round(x*qa+noise),t=round(x*qb+noise);
		auto re=_mm256_fnmadd_pd(t,vc,_mm256_fnmadd_pd(q,va,x));
		auto im=_mm256_fnmadd_pd(t,ve,-q*vb)*vs;
		_mm256_storeu_pd(f.x.data()+i,re);
		_mm256_storeu_pd(f.y.data()+i,im);
	}
	for(;i<n;i++){
		seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;
		double x=g[i].x,noise=int(seed)*0x1p-32;
		double q=nearbyint(x*e/det+noise),t=nearbyint(-x*b/det+noise);
		f.x[i]=x-q*a-t*c;f.y[i]=(-q*b-t*e)*scale;
	}
}
bool usable(int n,int m){
	setM(MOD);
	// Conservative empirical budget, not an exact-arithmetic guarantee.
	return d&&(__uint128_t)(n+m-1)*d*MOD*MOD<=((__uint128_t)1<<84);
}
poly mul(const poly &f,const poly &g){
	int m=f.size()+g.size()-1,n=1;
	while(n<m)n<<=1;
	Point a(n),b(n);
	read(a,f,0x123456789abcdefull);read(b,g,0xfedcba987654321ull);
	fft<false>(a);fft<false>(b);
	for(int i=0;i<n;i+=4){
		auto x=_mm256_loadu_pd(a.x.data()+i);
		auto y=_mm256_loadu_pd(a.y.data()+i);
		auto z=_mm256_loadu_pd(b.x.data()+i);
		auto t=_mm256_loadu_pd(b.y.data()+i);
		_mm256_storeu_pd(a.x.data()+i,_mm256_fmsub_pd(x,z,y*t));
		_mm256_storeu_pd(a.y.data()+i,_mm256_fmadd_pd(x,t,y*z));
	}
	fft<true>(a);
	poly h(m);
	for(int i=0;i<m;i++){
		ll x=llround(a.x[i]/n)%MOD,y=llround(a.y[i]/(n*scale))%MOD;
		h[i]=x+y*root;
	}
	return h;
}
};
#pragma GCC pop_options
MTT_FFT mtt;
// Odd primes use the lattice FFT; other moduli / large sizes use exact CRT.
poly operator*(const poly &f,const poly &g){
	if(f.empty()||g.empty())return{};
	assert(f.size()+g.size()-1<=(1<<FFT_MAX));
	if(min(f.size(),g.size())<=32||!mtt.usable(f.size(),g.size()))
		return mtt_exact(f,g);
	return mtt.mul(f,g);
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
poly BM(poly a){
	poly C{1},B{1};
	int L=0,m=1;
	mint b=1;
	for(int n=0;n<(int)a.size();n++){
		mint d=0;
		for(int i=0;i<=L;i++)d+=C[i]*a[n-i];
		if(!d){
			m++;
		}else{
			poly T=C;
			mint coef=d*b.inv();
			poly xmB(m,0);
			for(mint val:B)xmB.push_back(val);
			if(C.size()<xmB.size())C.resize(xmB.size(),0);
			for(int i=0;i<(int)xmB.size();i++)C[i]-=coef*xmB[i];
			if(2*L<=n){
				L=n+1-L;
				B=T;
				b=d;
				m=1;
			}else{
				m++;
			}
		}
	}
	return C;
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


poly lagrange(const vector<mint>&x,const vector<mint>&y){
	int n=x.size();assert(y.size()==x.size());
	poly g(n+1),f(n),q(n);g[0]=1;
	for(int i=0;i<n;i++){
		for(int j=i+1;j;j--)g[j]=g[j-1]-x[i]*g[j];
		g[0]*=-x[i];
	}
	for(int i=0;i<n;i++){
		q[n-1]=g[n];
		for(int j=n-2;j>=0;j--)q[j]=g[j+1]+x[i]*q[j+1];
		mint d=value(q,x[i]);assert(d);
		mint c=y[i]/d;
		for(int j=0;j<n;j++)f[j]+=q[j]*c;
	}
	return f;
}
// end for polynomial/mtt-fast.cpp
/////////////////////////
// !!!!! Call mint::setM(mod) first; requires AVX2/FMA and <immintrin.h>.
// FFT has rounding risk; use mtt_exact(f,g) for exact convolution. !!!!
