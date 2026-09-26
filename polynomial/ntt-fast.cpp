////////////////////////////////////////////////////////////////
//
// template for ModInt
// version 3.0 (Last Update Jun 15th, 2026)
//
// usage:
//   init(n); mint a,b; C(n,k); binom(n,k);
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
// template for polynomial
// version 7.0 (Last Update Jun 15th, 2026)
// NTT core replaced with the faster version
//
// usage:
//    init(n);
//   poly h=f*g; Inv, Ln, Exp, BM, FSPE, value;
//
////////////////////////////////////////////////////////////////

using ll=long long;
using poly=vector<mint>;

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
				const unsigned long long y=as[i].x+MO-as[i+m].x;
				as[i].x+=as[i+m].x;
				as[i+m].x=(prod.x*y)%MO;
			}
			prod*=INV_FFT_RATIOS[__builtin_ctz(++h)];
		}
		m<<=1;
	}
	for(;m<n>>1;m<<=1){
		mint prod=1U;
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

poly Inv(poly f){
	assert(!f.empty()&&f[0]);
	int n=f.size();
	int l=0;
	while((1<<l)<n)++l;
	f.resize(1<<l);
	poly g{f[0].inv()},_f;
	for(int i=1;i<=l;i++){
		_f=poly(begin(f),begin(f)+(1<<i));
		g.resize(1<<(i+1));_f.resize(1<<(i+1));
		fft(g);fft(_f);
		for(int j=0;j<(1<<(i+1));j++)
			g[j]=mint(2)*g[j]-g[j]*g[j]*_f[j];
		invFft(g);
		fill(begin(g)+(1<<i),end(g),0);
	}
	g.resize(n);
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

poly Ln(poly f){
	assert(!f.empty()&&f[0].x==1);
	poly f_=diff(f),_f=Inv(f);
	int n=f_.size(),m=_f.size();
	int l=0;
	while((1<<l)<n+m)++l;
	f_.resize(1<<l);_f.resize(1<<l);
	fft(f_);fft(_f);
	for(int i=0;i<(1<<l);i++)f_[i]*=_f[i];
	invFft(f_);
	f_=integ(f_);
	f_.resize(f.size());
	return f_;
}

poly Exp(poly f){
	assert(!f.empty()&&!f[0]);
	poly g{1},_f,_g;
	int n=f.size();
	int l=0;
	while((1<<l)<n)++l;
	f.resize(1<<l);
	for(int i=1;i<=l;i++){
		_f=poly(begin(f),begin(f)+(1<<i));
		_g=Ln(g);
		g.resize(1<<(i+1));
		_f.resize(1<<(i+1));
		_g.resize(1<<(i+1));
		fft(g);
		fft(_f);
		fft(_g);
		for(int j=0;j<(1<<(i+1));j++)
			g[j]*=mint(1)-_g[j]+_f[j];
		invFft(g);
		fill(begin(g)+(1<<i),end(g),0);
	}
	g.resize(n);
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

poly lagrange(vector<mint>x,vector<mint>y){
	int n=x.size();assert(y.size()==x.size());
	if(!n)return{};
	poly g(n+1),q(n),f(n);g[0]=1;
	for(int i=0;i<n;i++){
		for(int j=i+1;j>=1;j--)g[j]=g[j-1]-x[i]*g[j];
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
// use init(n) before accessing fac / ifac / inv directly
// end for polynomial/ntt-fast.cpp
/////////////////////////
// !!!!! Choose only one polynomial implementation; it already defines modular
// arithmetic. Check constant terms before inv/ln/exp and the transform size
// limit. !!!!
