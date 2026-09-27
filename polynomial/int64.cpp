////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/polynomial/int64.cpp
//
// usage: poly h=f*g; ull x=recurrence(a,c,k);
//   FSPE(f,g,k); RSPE(a,k);
//
////////////////////////////////////////////////////////////////
using poly=vector<ull>;
using Poly=poly;
template<unsigned MOD>
struct Int_Poly{
using V=vector<unsigned>;
static unsigned power(unsigned x,unsigned k){
	unsigned ans=1;
	for(;k;k>>=1,x=1ull*x*x%MOD)if(k&1)ans=1ull*ans*x%MOD;
	return ans;
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
// Four 16-bit limbs; each integer coefficient is less than 2^56.
static vector<V> mul(const poly &f,const poly &g,int n){
	static const NTT ntt;
	vector<V>a(4,V(n)),b(4,V(n)),c(4);
	for(int k=0;k<4;k++){
		for(int i=0;i<(int)f.size();i++)a[k][i]=(f[i]>>(16*k))&65535;
		for(int i=0;i<(int)g.size();i++)b[k][i]=(g[i]>>(16*k))&65535;
		ntt.fft(a[k].data(),n);ntt.fft(b[k].data(),n);
	}
	V h(n);
	for(int k=0;k<4;k++){
		for(int i=0;i<n;i++){
			ull x=0;
			for(int j=0;j<=k;j++)x+=1ull*a[j][i]*b[k-j][i];
			h[i]=x%MOD;
		}
		ntt.invFft(h.data(),n);
		c[k].assign(h.begin(),h.begin()+f.size()+g.size()-1);
	}
	return c;
}
};
// Odd numbers are exactly the invertible elements modulo 2^64.
ull int_inv(ull x){
	assert(x&1);
	ull y=1;
	for(int i=0;i<6;i++)y*=2-x*y;
	return y;
}
// O(n log n); unsigned overflow is intentional throughout this template.
poly operator*(const poly &f,const poly &g){
	if(f.empty()||g.empty())return{};
	assert(f.size()+g.size()-1<=(1<<23));
	int m=f.size()+g.size()-1;
	poly h(m);
	if(min(f.size(),g.size())<=32){
		for(int i=0;i<(int)f.size();i++)
			for(int j=0;j<(int)g.size();j++)h[i+j]+=f[i]*g[j];
		return h;
	}
	int n=1;
	while(n<m)n<<=1;
	constexpr unsigned p=998244353,q=469762049;
	static const unsigned v=Int_Poly<q>::power(p%q,q-2);
	auto a=Int_Poly<p>::mul(f,g,n),b=Int_Poly<q>::mul(f,g,n);
	for(int k=0;k<4;k++)for(int i=0;i<m;i++){
		ull x=1ull*(b[k][i]+q-a[k][i]%q)*v%q;
		h[i]+=(x*p+a[k][i])<<(16*k);
	}
	return h;
}
// Return [x^k] f/g; g[0] must be odd. O(d log(d+1) log(k+1)).
ull FSPE(poly f,poly g,ll k){
	assert(k>=0&&!g.empty()&&(g[0]&1));
	while(k&&!f.empty()){
		poly h=g;
		for(int i=1;i<(int)h.size();i+=2)h[i]=0-h[i];
		f=f*h;g=g*h;
		int n=0;
		for(int i=k&1;i<(int)f.size();i+=2)f[n++]=f[i];
		f.resize(n);n=0;
		for(int i=0;i<(int)g.size();i+=2)g[n++]=g[i];
		g.resize(n);k>>=1;
	}
	return f.empty()?0:f[0]*int_inv(g[0]);
}
// a has d initial terms; a[n]=sum c[j]*a[n-1-j], 0-based k.
ull recurrence(const poly &a,const poly &c,ll k){
	assert(k>=0&&a.size()==c.size());
	if(k<(ll)a.size())return a[k];
	if(a.empty())return 0;
	poly g(c.size()+1);g[0]=1;
	for(int i=0;i<(int)c.size();i++)g[i+1]=0-c[i];
	poly f=a*g;f.resize(a.size());
	return FSPE(move(f),move(g),k);
}
// Solve for an order-d recurrence over Z/(2^64), including even pivots.
// Internal helper; O(a.size()*d^2) time, O(a.size()*d) space.
bool int_rec(const poly &a,int d,poly &c){
	int n=a.size()-d,r=0;
	vector<poly>b(n,poly(d+1));
	vector<int>p(d);
	iota(p.begin(),p.end(),0);
	for(int i=0;i<n;i++){
		for(int j=0;j<d;j++)b[i][j]=a[i+d-1-j];
		b[i][d]=a[i+d];
	}
	for(;r<min(n,d);r++){
		int x=-1,y=-1,v=64;
		for(int i=r;i<n;i++)for(int j=r;j<d;j++)if(b[i][j]){
			int t=__builtin_ctzll(b[i][j]);
			if(t<v)x=i,y=j,v=t;
		}
		if(x<0)break;
		swap(b[x],b[r]);swap(p[y],p[r]);
		for(poly &row:b)swap(row[y],row[r]);
		ull z=int_inv(b[r][r]>>v);
		for(int j=r;j<=d;j++)b[r][j]*=z;
		for(int i=r+1;i<n;i++){
			ull t=b[i][r]>>v;
			for(int j=r;j<=d;j++)b[i][j]-=t*b[r][j];
		}
	}
	for(int i=r;i<n;i++)if(b[i][d])return false;
	poly x(d);
	for(int i=r-1;i>=0;i--){
		ull z=b[i][d];
		for(int j=i+1;j<d;j++)z-=b[i][j]*x[j];
		if(z&(b[i][i]-1))return false;
		x[i]=z>>__builtin_ctzll(b[i][i]);
	}
	c.resize(d);
	for(int i=0;i<d;i++)c[p[i]]=x[i];
	return true;
}
// Infer a shortest recurrence; provide >=2*d samples of an order-d sequence.
// Ring elimination costs O(n^3 log(n+1)); prefer recurrence if c is known.
ull RSPE(const poly &a,ll k){
	assert(k>=0);
	if(k<(ll)a.size())return a[k];
	if(a.empty())return 0;
	int l=0,r=a.size()/2;
	poly c;
	bool ok=int_rec(a,r,c);assert(ok);
	while(l<r){
		int m=(l+r)>>1;
		if(int_rec(a,m,c))r=m;else l=m+1;
	}
	ok=int_rec(a,l,c);assert(ok);
	return recurrence(poly(a.begin(),a.begin()+l),c,k);
}
// end for polynomial/int64.cpp
/////////////////////////
// !!!!! ull arithmetic wraps modulo 2^64. FSPE requires odd g[0].
// Do not combine with another poly implementation. !!!!
