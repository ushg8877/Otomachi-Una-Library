////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/polynomial/ntt.cpp
//
// usage: poly h=f*g; Inv(f); Ln(f); Exp(f);
//   init(n); C(n,k); BM(a);
//
////////////////////////////////////////////////////////////////
template<unsigned _M>
struct ModInt{
static constexpr unsigned MOD=_M;
static_assert(1<MOD&&MOD<=INT_MAX);
unsigned x;
constexpr ModInt():x(0){}
constexpr ModInt(unsigned y):x(y%MOD){}
constexpr ModInt(int y):x((y%=static_cast<int>(MOD))<0?y+MOD:y){}
constexpr ModInt(unsigned long long y):x(y%MOD){}
constexpr ModInt(long long y):x((y%=static_cast<long long>(MOD))<0?y+MOD:y){}
ModInt &operator+=(const ModInt &a){x=(((x+=a.x)>=MOD)?x-MOD:x);return *this;}
ModInt &operator-=(const ModInt &a){x=(((x-=a.x)>=MOD)?x+MOD:x);return *this;}
ModInt &operator*=(const ModInt &a){
	x=static_cast<unsigned long long>(x)*a.x%MOD;
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
//////////////////////////////////////////////////////////////////

const int MOD=998244353;
using mint=ModInt<MOD>;
vector<mint>fac{1},ifac{1},inv{0,1};
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


////////////////////////////////////////////////////////////////
//
// template for polynomial, prepared by Otomachi Una
// version 6.1 (Last Update Jun 13th, 2026)
// NTT core replaced with the "tech" version
//
////////////////////////////////////////////////////////////////

using poly=vector<mint>;

// ---------- NTT constants and roots (tech version) ----------
constexpr unsigned MO=998244353U;
constexpr unsigned MO2=2U*MO;
constexpr int FFT_MAX=23;
constexpr array<mint,FFT_MAX+1> FFT_ROOTS={
	1U,998244352U,911660635U,372528824U,929031873U,452798380U,922799308U,
	781712469U,476477967U,166035806U,258648936U,584193783U,63912897U,
	350007156U,666702199U,968855178U,629671588U,24514907U,996173970U,
	363395222U,565042129U,733596141U,267099868U,15311432U
};
constexpr array<mint,FFT_MAX+1> INV_FFT_ROOTS={
	1U,998244352U,86583718U,509520358U,337190230U,87557064U,609441965U,
	135236158U,304459705U,685443576U,381598368U,335559352U,129292727U,
	358024708U,814576206U,708402881U,283043518U,3707709U,121392023U,
	704923114U,950391366U,428961804U,382752275U,469870224U
};
constexpr array<mint,FFT_MAX> FFT_RATIOS={
	911660635U,509520358U,369330050U,332049552U,983190778U,123842337U,
	238493703U,975955924U,603855026U,856644456U,131300601U,842657263U,
	730768835U,942482514U,806263778U,151565301U,510815449U,503497456U,
	743006876U,741047443U,56250497U,867605899U
};
constexpr array<mint,FFT_MAX> INV_FFT_RATIOS={
	86583718U,372528824U,373294451U,645684063U,112220581U,692852209U,155456985U,
	797128860U,90816748U,860285882U,927414960U,354738543U,109331171U,
	293255632U,535113200U,308540755U,121186627U,608385704U,438932459U,
	359477183U,824071951U,103369235U
};

// as[rev(i)] <- \sum_j \zeta^(ij) as[j]
void fft(mint *as,int n){
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
		mint prod=1U;
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
			mint prod=1U;
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
			mint prod=1U;
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
void invFft(mint *as,int n){
	assert(!(n&(n-1)));
	assert(1<=n);
	assert(n<=1<<FFT_MAX);
	int m=1;
	if(m<n>>1){
		mint prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+m;++i){
				// < 2 MO
				const unsigned long long y=as[i].x+MO-as[i+m].x;
				as[i].x+=as[i+m].x; // < 2 MO
				as[i+m].x=(prod.x*y)%MO; // < MO
			}
			prod*=INV_FFT_RATIOS[__builtin_ctz(++h)];
		}
		m<<=1;
	}
	for(;m<n>>1;m<<=1){
		mint prod=1U;
		for(int h=0,i0=0;i0<n;i0+=(m<<1)){
			for(int i=i0;i<i0+(m>>1);++i){
				// < 4 MO
				const unsigned long long y=as[i].x+MO2-as[i+m].x;
				as[i].x+=as[i+m].x; // < 4 MO
				// < 2 MO
				as[i].x=(as[i].x>=MO2)?(as[i].x-MO2):as[i].x;
				as[i+m].x=(prod.x*y)%MO; // < MO
			}
			for(int i=i0+(m>>1);i<i0+m;++i){
				// < 2 MO
				const unsigned long long y=as[i].x+MO-as[i+m].x;
				as[i].x+=as[i+m].x; // < 2 MO
				as[i+m].x=(prod.x*y)%MO; // < MO
			}
			prod*=INV_FFT_RATIOS[__builtin_ctz(++h)];
		}
	}
	if(m<n){
		for(int i=0;i<m;++i){
			const unsigned y=as[i].x+MO2-as[i+m].x; // < 4 MO
			as[i].x+=as[i+m].x; // < 4 MO
			as[i+m].x=y; // < 4 MO
		}
	}
	const mint invN=mint(n).inv();
	for(int i=0;i<n;++i){
		as[i]*=invN;
	}
}

void fft(poly &as){
	fft(as.data(),as.size());
}
void invFft(poly &as){
	invFft(as.data(),as.size());
}

// ----------------------------------------------------------------

poly operator*(poly f,poly g){
	if(f.empty()||g.empty())return{};
	assert(f.size()+g.size()-1<=(1<<FFT_MAX));
	if(min(f.size(),g.size())<=32){
		poly h(f.size()+g.size()-1);
		for(int i=0;i<(int)f.size();i++)
			for(int j=0;j<(int)g.size();j++)h[i+j]+=f[i]*g[j];
		return h;
	}
	int n=f.size()+g.size();
	int l=0;
	while((1<<l)<n-1)++l;
	f.resize(1<<l);g.resize(1<<l);
	fft(f);fft(g);
	for(int i=0;i<(1<<l);i++)f[i]*=g[i];
	invFft(f);
	f.resize(n-1);
	return f;
}

poly operator/(poly f,poly g){
	if(f.empty()||g.empty())return poly(f.size());
	int m=g.size();
	reverse(g.begin(),g.end());
	g=f*g;
	for(int i=0;i<(int)f.size();i++)f[i]=g[i+m-1];
	return f;
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

poly operator*(poly f,mint x){
	for(mint &i:f)i*=x;
	return f;
}

poly operator*(mint x,poly f){
	for(mint &i:f)i*=x;
	return f;
}

mint value(const poly &f,mint x){
	mint ans=0;
	for(int i=f.size()-1;i>=0;i--)ans=ans*x+f[i];
	return ans;
}

// Formal inverse; only compute the new upper half each round.
poly Inv(poly f){
	assert(!f.empty()&&f[0]&&f.size()<=(1<<FFT_MAX));
	int n=f.size();
	poly g{f[0].inv()},a,b;
	for(int m=1;m<n;m<<=1){
		int t=m<<1;
		a.assign(f.begin(),f.begin()+min(n,t));a.resize(t);
		b=g;b.resize(t);
		fft(a);fft(b);
		for(int i=0;i<t;i++)a[i]*=b[i];
		invFft(a);
		fill(a.begin(),a.begin()+m,0);
		fft(a);
		for(int i=0;i<t;i++)a[i]*=b[i];
		invFft(a);g.resize(min(n,t));
		for(int i=m;i<(int)g.size();i++)g[i]=-a[i];
	}
	return g;
}

poly integ(poly f){
	int n=f.size();
	::init(n);
	f.resize(n+1);
	for(int i=n;i>=1;i--)f[i]=f[i-1]*::inv[i];
	f[0]=0;
	return f;
}

poly diff(poly f){
	if(f.empty())return{};
	int n=f.size();
	for(int i=0;i<n-1;i++)
		f[i]=f[i+1]*(i+1);
	f.pop_back();
	return f;
}

// Formal logarithm; f[0]=1. Integration denominators must be invertible.
poly Ln(poly f){
	assert(!f.empty()&&f[0].x==1);
	int n=f.size();
	poly g=diff(f)*Inv(f);
	g.resize(n-1);
	return integ(move(g));
}

// Formal exponential; f[0]=0. Integration denominators must be invertible.
poly Exp(poly f){
	assert(!f.empty()&&!f[0]);
	int n=f.size();
	poly g{1};
	// exp(f)=g*(1+f-ln(g)) modulo x^(2*m).
	for(int m=1;m<n;m<<=1){
		int t=min(n,m<<1);
		g.resize(t);
		poly h=Ln(g),a(t-m);
		for(int i=m;i<t;i++)a[i-m]=f[i]-h[i];
		g.resize(m);a=g*a;g.resize(t);
		for(int i=m;i<t;i++)g[i]=a[i-m];
	}
	return g;
}

// f[0]=1; choose the square root with constant term 1.
poly Sqrt(poly f){
	assert(!f.empty()&&f[0].x==1);
	int n=f.size();
	poly g{1};
	mint v=mint(2).inv();
	for(int m=1;m<n;m<<=1){
		int t=min(n,m<<1);
		g.resize(t);
		poly h=poly(f.begin(),f.begin()+t)*Inv(g);
		for(int i=m;i<t;i++)g[i]=h[i]*v;
	}
	return g;
}

// f=q*g+r; g.back()!=0. Unlike operator/, this is division.
void Div(poly f,poly g,poly &q,poly &r){
	assert(!g.empty()&&g.back()&&(&q!=&r));
	if(f.size()<g.size()){q={0};r=move(f);return;}
	int n=f.size(),m=g.size(),k=n-m+1;
	poly a(f.rbegin(),f.rbegin()+k),b(g.rbegin(),g.rend());
	b.resize(k);q=a*Inv(b);q.resize(k);
	reverse(q.begin(),q.end());
	b=q*g;r.resize(m-1);
	for(int i=0;i<m-1;i++)r[i]=f[i]-b[i];
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



struct Poly_Comp{
	static void solve(poly &x,const poly &f,int m,int n,mint w){
		int o=m*n;
		if(n==1){
			if(!w){x=f;x.resize(o);return;}
			init(2*m-1);x.resize(m+1);
			mint v=ifac[m-1];
			for(int i=0;i<=m;i++){
				x[m-i]=v*fac[m-1+i]*ifac[i];v*=w;
			}
			x=x*f;x.erase(x.begin(),x.begin()+m);x.resize(o);
			return;
		}
		poly y(4*o),z(2*o),u(2*o);
		for(int i=0;i<m;i++)
			copy_n(x.begin()+i*n,n,y.begin()+2*i*n);
		y[2*o]=1;fft(y);
		for(int i=0;i<2*o;i++)z[i]=y[2*i]*y[2*i+1];
		invFft(z);z[0]-=1;
		for(int i=1;i<2*m;i++)
			copy_n(z.begin()+i*n,n/2,z.begin()+i*(n/2));
		z.resize(o);solve(z,f,2*m,n/2,w);
		for(int i=0;i<2*m;i++)
			copy_n(z.begin()+i*(n/2),n/2,u.begin()+i*n);
		fft(u);x.resize(4*o);
		for(int i=0;i<2*o;i++){
			x[2*i]=y[2*i+1]*u[i];x[2*i+1]=y[2*i]*u[i];
		}
		invFft(x);
		for(int i=0;i<m;i++)
			copy_n(x.begin()+2*(i+m)*n,n,x.begin()+i*n);
		x.resize(o);
	}
};

// f(g(x)) mod x^n; n defaults to f.size(). g[0] may be nonzero.
// O(k log^2 k), k=max(n,f.size()); at most 2^21 coefficients.
poly Compose(poly f,const poly &g,int n=-1){
	if(n<0)n=f.size();
	assert(0<=n&&n<=(1<<21));
	if(!n)return {};
	mint w=g.empty()?mint(0):g[0];
	if(!w&&f.size()>(size_t)n)f.resize(n);
	assert(f.size()<=(1<<21));
	if(f.empty())return poly(n);
	if(n<=32){
		poly x(n),y(n);
		for(int i=(int)f.size()-1;i>=0;i--){
			fill(y.begin(),y.end(),0);
			for(int j=0;j<n;j++)
				for(int k=0;k<min(n-j,(int)g.size());k++)
					y[j+k]+=x[j]*g[k];
			y[0]+=f[i];swap(x,y);
		}
		return x;
	}
	int k=1;
	while(k<max(n,(int)f.size()))k<<=1;
	poly x(k);
	for(int i=0;i<min(n,(int)g.size());i++)x[i]=-g[i];
	Poly_Comp::solve(x,f,1,k,w);x.resize(n);return x;
}

// g(0)=0, f(g(x))=x mod x^n; f[0]=0, f[1] and 1..n-1 invertible.
// O(n log^2 n) time, O(n) space; n defaults to f.size().
poly Compose_Inv(poly f,int n=-1){
	if(n<0)n=f.size();
	assert(0<=n&&n<=(1<<21));
	if(!n)return {};
	assert(f.empty()||!f[0]);
	if(n==1)return {0};
	assert(f.size()>1&&f[1]);
	if(n==2)return {0,f[1].inv()};
	f.resize(n);
	int len=1;
	while(len<n)len<<=1;
	poly x(2*len),y(2*len),u(4*len),v(4*len),rt(2*len);
	rt[0]=1;
	for(int i=1;i<2*len;i<<=1){
		mint w=mint(3).pow((MOD-1)/(4*i)).inv();
		for(int j=0;j<i;j++)rt[i+j]=rt[j]*w;
	}
	mint w=f[1].inv(),p=-1,c=mint(2).inv();
	for(int i=1;i<n;i++){p*=w;x[i]=f[i]*p;}
	y[0]=1;
	int j=1;
	for(int i=n;i>1;j<<=1){
		int o=(i+1)/2,m=1;
		while(m<4*i*j)m<<=1;
		u.assign(m,0);v.assign(m,0);
		for(int k=0;k<j;k++){
			copy_n(x.begin()+i*k,i,u.begin()+2*i*k);
			copy_n(y.begin()+i*k,i,v.begin()+2*i*k);
		}
		u[2*i*j]=1;fft(u);fft(v);x.resize(m/2);y.resize(m/2);
		for(int k=0;k<m/2;k++){
			x[k]=u[2*k]*u[2*k+1];
			mint a=u[2*k+1]*v[2*k],b=u[2*k]*v[2*k+1];
			y[k]=(i&1)?(a+b)*c:(a-b)*rt[k]*c;
		}
		invFft(x);invFft(y);
		if(4*i*j==m)x[0]-=1;
		for(int k=1;k<2*j;k++){
			copy_n(x.begin()+i*k,o,x.begin()+o*k);
			copy_n(y.begin()+i*k,o,y.begin()+o*k);
		}
		i=o;
	}
	y.resize(j);reverse(y.begin(),y.end());y.resize(n);init(n-1);
	for(int i=1;i<n;i++)y[i]*=(n-1)*inv[i];
	reverse(y.begin(),y.end());y.resize(n-1);
	y=Exp(Ln(move(y))*(-inv[n-1]))*w;
	y.insert(y.begin(),0);return y;
}


struct BM_Work{
	struct Mat{poly a{1},b,c,d{1};};
	static void trim(poly &a){while(!a.empty()&&!a.back())a.pop_back();}
	static poly high(const poly &a,int k){
		return poly(a.begin()+min(k,(int)a.size()),a.end());
	}
	static void apply(const Mat &m,poly &a,poly &b){
		int k=max({m.a.size(),m.b.size(),m.c.size(),m.d.size()});
		if(k>64&&min(a.size(),b.size())>64){
			int len=k+max(a.size(),b.size())-1,n=1;
			while(n<len)n<<=1;
			array<poly,4> f{m.a,m.b,m.c,m.d};
			for(poly &v:f){v.resize(n);fft(v);}
			a.resize(n);b.resize(n);fft(a);fft(b);
			for(int i=0;i<n;i++){
				mint x=f[0][i]*a[i]+f[1][i]*b[i];
				b[i]=f[2][i]*a[i]+f[3][i]*b[i];a[i]=x;
			}
			invFft(a);invFft(b);a.resize(len);b.resize(len);
			trim(a);trim(b);return;
		}
		poly x=m.a*a+m.b*b,y=m.c*a+m.d*b;
		trim(x);trim(y);a=move(x);b=move(y);
	}
	static Mat mul(const Mat &x,const Mat &y){
		int k=max({x.a.size(),x.b.size(),x.c.size(),x.d.size()});
		int l=max({y.a.size(),y.b.size(),y.c.size(),y.d.size()});
		if(min(k,l)>64){
			int len=k+l-1,n=1;
			while(n<len)n<<=1;
			array<poly,4> a{x.a,x.b,x.c,x.d},b{y.a,y.b,y.c,y.d},c;
			for(poly &v:a){v.resize(n);fft(v);}
			for(poly &v:b){v.resize(n);fft(v);}
			for(int i=0;i<4;i++){
				c[i].resize(n);
				for(int j=0;j<n;j++)
					c[i][j]=a[i/2*2][j]*b[i%2][j]
						+a[i/2*2+1][j]*b[i%2+2][j];
				invFft(c[i]);c[i].resize(len);trim(c[i]);
			}
			return {move(c[0]),move(c[1]),move(c[2]),move(c[3])};
		}
		Mat z{x.a*y.a+x.b*y.c,x.a*y.b+x.b*y.d,
			x.c*y.a+x.d*y.c,x.c*y.b+x.d*y.d};
		trim(z.a);trim(z.b);trim(z.c);trim(z.d);return z;
	}
	static poly div(poly &a,const poly &b){
		int n=a.size(),m=b.size(),k=n-m+1;
		poly q(k);
		if(min(m,k)<=32){
			mint v=b.back().inv();
			for(int i=k-1;i>=0;i--){
				q[i]=a[i+m-1]*v;
				for(int j=0;j<m;j++)a[i+j]-=q[i]*b[j];
			}
		}else{
			poly x(a.rbegin(),a.rbegin()+k),y(b.rbegin(),b.rend());
			y.resize(k);q=x*Inv(move(y));q.resize(k);
			reverse(q.begin(),q.end());a-=q*b;
		}
		a.resize(m-1);trim(a);return q;
	}
	static void step(Mat &m,poly &a,poly &b){
		poly q=div(a,b);swap(a,b);
		poly x=m.a-q*m.c,y=m.b-q*m.d;
		trim(x);trim(y);
		m.a=move(m.c);m.b=move(m.d);m.c=move(x);m.d=move(y);
	}
	static Mat solve(poly a,poly b){
		int k=(a.size()+1)/2;
		Mat m;
		if((int)b.size()<=k)return m;
		if(a.size()<=64){
			while((int)b.size()>k)step(m,a,b);
			return m;
		}
		m=solve(high(a,k),high(b,k));apply(m,a,b);
		if((int)b.size()<=k)return m;
		step(m,a,b);
		if((int)b.size()<=k)return m;
		int j=2*k-(int)a.size()+1;
		return mul(solve(high(a,j),high(b,j)),m);
	}
};

// Shortest order d=c.size()-1: c[0]=1, sum c[j]*a[i-j]=0 for i>=d.
// Half-GCD: O(n log^2 n). Work over a field; keep trailing zeroes.
poly BM(const poly &a){
	if(a.size()<=64){
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
	poly x(a.size()+1),y(a.rbegin(),a.rend());x.back()=1;
	BM_Work::trim(y);
	if(y.empty())return {1};
	auto m=BM_Work::solve(x,y);BM_Work::apply(m,x,y);
	while(y.size()>=m.d.size())BM_Work::step(m,x,y);
	poly c=move(m.d);
	assert(!c.empty()&&c.back());
	mint v=c.back().inv();reverse(c.begin(),c.end());
	for(mint &i:c)i*=v;
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
	// Infer a recurrence, then return f[x].
	assert(x>=0);
	if(x<(ll)f.size())return f[x];
	poly g=BM(f);
	f.resize(g.size()-1);
	poly s=f*g;s.resize(g.size()-1);
	return FSPE(move(s),move(g),x);
}

// Interpolate arbitrary points; return n coefficients, degree < n.
// O(n log^2 n) time, O(n log n) space. Pairwise x differences must be units.
poly lagrange(vector<mint>x,vector<mint>y){
	return Poly_Tree(move(x)).interpolate(y);
}

poly to_ex(poly f){ // f(e^x)
	int n=f.size();
	if(!n)return{};
	::init(n);
	function<pair<poly,poly>(int,int)>solve=[&](int l,int r){
		if(l==r)return make_pair(poly{f[l]},poly{1,MOD-l});
		int mid=l+r>>1;
		auto ls=solve(l,mid),rs=solve(mid+1,r);
		return make_pair(ls.first*rs.second+ls.second*rs.first,
			ls.second*rs.second);
	};
	auto ans=solve(0,n-1);
	poly g=ans.first*Inv(ans.second);
	g.resize(n);
	for(int i=0;i<n;i++)g[i]=g[i]*ifac[i];
	return g;
}

poly S2line(int n){
	// return S(n,i)
	assert(n>=0);::init(n);
	poly F(n+1),G(n+1);
	// S(n,i) * binom(j,i) -> F(j)
	for(int i=0;i<=n;i++){
		G[i]=ifac[i];F[i]=mint(i).pow(n)*ifac[i];
		if(i&1)F[i]*=-1;
	}
	F=F*G;
	F.resize(n+1);
	for(int i=0;i<=n;i++){
		if(i&1)F[i]*=-1;
	}
	return F;
}
// end for polynomial/ntt.cpp
/////////////////////////
// !!!!! Contains mint already; do not paste another mint or poly
// implementation. !!!!
