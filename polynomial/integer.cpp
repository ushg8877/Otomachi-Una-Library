////////////////////////////////////////////////////////////////
//
// template for integer polynomial / linear recurrence
// usage:
//   poly h=f*g; // exact signed ll coefficients, no target modulus
//   FSPE(F,G,i); // [x^i] F/G, G[0] must be 1 or -1
//   recurrence(a,c,i); // a[n]=sum c[j]*a[n-1-j], a.size()==c.size()
//   RSPE(a,i); // infer from samples; at least 2*d terms for order d
//   Indices start at 0; the true returned coefficients / term must fit ll.
//   Fast NTT core from ntt-fast.cpp; three-prime CRT; convolution O(n log n),
//   length <= 2^23.
//   Recurrence O(d log(d+1) log(i+1)); RSPE also needs O(a.size()*d).
//
////////////////////////////////////////////////////////////////
using poly=vector<ll>;
using Poly=poly;
template<unsigned MOD>
struct Int_Poly{
using V=vector<unsigned>;
static unsigned power(unsigned x,unsigned k){
	unsigned ans=1;
	for(;k;k>>=1,x=(ll)x*x%MOD)if(k&1)ans=(ll)ans*x%MOD;
	return ans;
}
static V read(const poly &f){
	V a(f.size());
	for(int i=0;i<(int)f.size();i++){
		ll x=f[i]%MOD;
		a[i]=x<0?x+MOD:x;
	}
	return a;
}
struct NTT{
static constexpr int FFT_MAX=23;
static constexpr unsigned MOD2=MOD*2;
V FFT_RATIOS,INV_FFT_RATIOS;
NTT():FFT_RATIOS(FFT_MAX,1),INV_FFT_RATIOS(FFT_MAX,1){
	static_assert(MOD<(1u<<30)&&((MOD-1)%(1<<FFT_MAX)==0));
	V root(FFT_MAX+1),iroot(FFT_MAX+1);
	root[FFT_MAX]=power(3,(MOD-1)>>FFT_MAX);
	iroot[FFT_MAX]=power(root[FFT_MAX],MOD-2);
	for(int i=FFT_MAX;i;i--)root[i-1]=1ull*root[i]*root[i]%MOD,
		iroot[i-1]=1ull*iroot[i]*iroot[i]%MOD;
	unsigned x=1,y=1;
	for(int i=0;i<FFT_MAX-1;i++){
		FFT_RATIOS[i]=1ull*root[i+2]*x%MOD;
		INV_FFT_RATIOS[i]=1ull*iroot[i+2]*y%MOD;
		x=1ull*x*iroot[i+2]%MOD;y=1ull*y*root[i+2]%MOD;
	}
}
// as[rev(i)] <- \sum_j \zeta^(ij) as[j]
void fft(unsigned *as,int n)const{
	assert(!(n&(n-1)));
	assert(1<=n);
	assert(n<=1<<FFT_MAX);
	int m=n;
	if(m>>=1){
		for(int i=0;i<m;++i){
			const unsigned x=as[i+m];
			as[i+m]=as[i]+MOD-x;
			as[i]+=x;
		}
	}
	if(m>>=1){
		unsigned prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+m;++i){
				const unsigned x=(1ull*prod*as[i+m])%MOD;
				as[i+m]=as[i]+MOD-x;
				as[i]+=x;
			}
			prod=1ull*prod*FFT_RATIOS[__builtin_ctz(++h)]%MOD;
		}
	}
	for(;m;){
		if(m>>=1){
			unsigned prod=1U;
			for(int h=0,i0=0;i0<n;i0+=(m<<1)){
				for(int i=i0;i<i0+m;++i){
					const unsigned x=(1ull*prod*as[i+m])%MOD;
					as[i+m]=as[i]+MOD-x;
					as[i]+=x;
				}
				prod=1ull*prod*FFT_RATIOS[__builtin_ctz(++h)]%MOD;
			}
		}
		if(m>>=1){
			unsigned prod=1U;
			for(int h=0,i0=0;i0<n;i0+=(m<<1)){
				for(int i=i0;i<i0+m;++i){
					const unsigned x=(1ull*prod*as[i+m])%MOD;
					as[i]=(as[i]>=MOD2)?(as[i]-MOD2):as[i];
					as[i+m]=as[i]+MOD-x;
					as[i]+=x;
				}
				prod=1ull*prod*FFT_RATIOS[__builtin_ctz(++h)]%MOD;
			}
		}
	}
	for(int i=0;i<n;++i){
		as[i]=(as[i]>=MOD2)?(as[i]-MOD2):as[i];
		as[i]=(as[i]>=MOD)?(as[i]-MOD):as[i];
	}
}

// as[i] <- (1/n) \sum_j \zeta^(-ij) as[rev(j)]
void invFft(unsigned *as,int n)const{
	assert(!(n&(n-1)));
	assert(1<=n);
	assert(n<=1<<FFT_MAX);
	int m=1;
	if(m<n>>1){
		unsigned prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+m;++i){
				const unsigned long long y=as[i]+MOD-as[i+m];
				as[i]+=as[i+m];
				as[i+m]=(prod*y)%MOD;
			}
			prod=1ull*prod*INV_FFT_RATIOS[__builtin_ctz(++h)]%MOD;
		}
		m<<=1;
	}
	for(;m<n>>1;m<<=1){
		unsigned prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+(m>>1);++i){
				const unsigned long long y=as[i]+MOD2-as[i+m];
				as[i]+=as[i+m];
				as[i]=(as[i]>=MOD2)?(as[i]-MOD2):as[i];
				as[i+m]=(prod*y)%MOD;
			}
			for(int i=i0+(m>>1);i<i0+m;++i){
				const unsigned long long y=as[i]+MOD-as[i+m];
				as[i]+=as[i+m];
				as[i+m]=(prod*y)%MOD;
			}
			prod=1ull*prod*INV_FFT_RATIOS[__builtin_ctz(++h)]%MOD;
		}
	}
	if(m<n){
		for(int i=0;i<m;++i){
			const unsigned y=as[i]+MOD2-as[i+m];
			as[i]+=as[i+m];
			as[i+m]=y;
		}
	}
	const unsigned invN=power(n,MOD-2);
	for(int i=0;i<n;++i){
		as[i]=1ull*as[i]*invN%MOD;
	}
}
};
static V mul(V a,V b){
	if(a.empty()||b.empty())return{};
	assert(a.size()+b.size()-1<=(1<<23));int m=a.size()+b.size()-1;
	if(min(a.size(),b.size())<=32){
		V c(m);
		for(int i=0;i<(int)a.size();i++)for(int j=0;j<(int)b.size();j++)
			c[i+j]=(c[i+j]+(ll)a[i]*b[j])%MOD;
		return c;
	}
	int n=1;
	while(n<m)n<<=1;
	a.resize(n);
	b.resize(n);
	static const NTT ntt;
	ntt.fft(a.data(),n);
	ntt.fft(b.data(),n);
	for(int i=0;i<n;i++)a[i]=(ll)a[i]*b[i]%MOD;
	ntt.invFft(a.data(),n);
	a.resize(m);
	return a;
}
static V bm(const V &a){
	V c{1},b{1};
	int len=0,m=1;
	unsigned last=1;
	for(int n=0;n<(int)a.size();n++){
		unsigned d=a[n];
		for(int i=1;i<=len;i++)d=(d+(ll)c[i]*a[n-i])%MOD;
		if(!d){m++;continue;}
		V old=c;unsigned v=(ll)d*power(last,MOD-2)%MOD;
		if(c.size()<b.size()+m)c.resize(b.size()+m);
		for(int i=0;i<(int)b.size();i++){
			unsigned x=(ll)v*b[i]%MOD;c[i+m]=c[i+m]>=x?c[i+m]-x:c[i+m]+MOD-x;
		}
		if(len<=n/2){len=n+1-len;b=move(old);last=d;m=1;}else m++;
	}
	c.resize(len+1);
	return c;
}
static unsigned nth(V f,V g,ll k){
	while(k&&!f.empty()){
		V h=g;
		for(int i=1;i<(int)h.size();i+=2)if(h[i])h[i]=MOD-h[i];
		f=mul(move(f),h);g=mul(move(g),move(h));
		int n=0;
		for(int i=k&1;i<(int)f.size();i+=2)f[n++]=f[i];
		f.resize(n);
		n=0;
		for(int i=0;i<(int)g.size();i+=2)g[n++]=g[i];
		g.resize(n);
		k>>=1;
	}
	return f.empty()?0:(ll)f[0]*power(g[0],MOD-2)%MOD;
}
static unsigned solve(V a,V g,ll k){
	int d=g.size()-1;if(!d)return 0;
	a.resize(d);
	a=mul(move(a),g);
	a.resize(d);
	return nth(move(a),move(g),k);
}
static unsigned fspe(const poly &f,const poly &g,ll k){
	return nth(read(f),read(g),k);
}
static unsigned recurrence(const poly &a,const poly &c,ll k){
	V g(c.size()+1);g[0]=1;
	for(int i=0;i<(int)c.size();i++){ll x=c[i]%MOD;g[i+1]=x>0?MOD-x:-x;}
	return solve(read(a),move(g),k);
}
static unsigned rspe(const poly &a,ll k){
	V f=read(a),g=bm(f);
	return solve(move(f),move(g),k);
}
};
ll int_crt(unsigned a,unsigned b,unsigned c){
	constexpr ll m1=998244353,m2=167772161,m3=469762049,m12=m1*m2;
	static const unsigned v1=Int_Poly<m2>::power(m1%m2,m2-2);
	static const unsigned v2=Int_Poly<m3>::power(m12%m3,m3-2);
	ll x=(b+m2-a%m2)*v1%m2*m1+a;
	ll y=(c+m3-x%m3)*v2%m3;
	__int128 ans=x+(__int128)m12*y,mod=(__int128)m12*m3;
	if(ans>mod/2)ans-=mod;
	assert(LLONG_MIN<=ans&&ans<=LLONG_MAX);
	return ans;
}
poly operator*(const poly &f,const poly &g){
	if(f.empty()||g.empty())return{};
	assert(f.size()+g.size()-1<=(1<<23));
	auto a=Int_Poly<998244353>::mul(Int_Poly<998244353>::read(f),
		Int_Poly<998244353>::read(g));
	auto b=Int_Poly<167772161>::mul(Int_Poly<167772161>::read(f),
		Int_Poly<167772161>::read(g));
	auto c=Int_Poly<469762049>::mul(Int_Poly<469762049>::read(f),
		Int_Poly<469762049>::read(g));
	poly h(a.size());
	for(int i=0;i<(int)h.size();i++)h[i]=int_crt(a[i],b[i],c[i]);
	return h;
}
ll FSPE(const poly &f,const poly &g,ll k){
	assert(k>=0&&!g.empty()&&(g[0]==1||g[0]==-1));
	if(f.empty())return 0;
	return int_crt(Int_Poly<998244353>::fspe(f,g,k),
		Int_Poly<167772161>::fspe(f,g,k),
		Int_Poly<469762049>::fspe(f,g,k));
}
ll recurrence(const poly &a,const poly &c,ll k){
	assert(k>=0&&a.size()==c.size());if(k<(ll)a.size())return a[k];
	if(a.empty())return 0;
	return int_crt(Int_Poly<998244353>::recurrence(a,c,k),
		Int_Poly<167772161>::recurrence(a,c,k),
		Int_Poly<469762049>::recurrence(a,c,k));
}
ll RSPE(const poly &a,ll k){
	assert(k>=0);if(k<(ll)a.size())return a[k];
	if(a.empty())return 0;
	return int_crt(Int_Poly<998244353>::rspe(a,k),
		Int_Poly<167772161>::rspe(a,k),
		Int_Poly<469762049>::rspe(a,k));
}
// end for polynomial/integer.cpp
/////////////////////////
// !!!!! Choose one poly implementation. True returned values must fit ll; CRT
// assertions cannot detect every overflow. RSPE needs enough samples of an
// integer-coefficient linear recurrence; internal intermediates need not fit
// ll. !!!!
