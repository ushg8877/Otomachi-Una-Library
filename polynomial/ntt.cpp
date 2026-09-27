////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
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

// Return f(x[i]); repeated evaluation points are allowed.
poly Eval(poly f,poly x){
	int n=x.size();
	if(!n)return{};
	vector<poly>g(n*4);
	auto build=[&](auto &&self,int u,int l,int r)->void{
		if(l==r){g[u]={-x[l],1};return;}
		int m=(l+r)>>1;
		self(self,u*2,l,m);self(self,u*2+1,m+1,r);
		g[u]=g[u*2]*g[u*2+1];
	};
	build(build,1,0,n-1);
	poly ans(n);
	auto solve=[&](auto &&self,poly a,int u,int l,int r)->void{
		if(r-l<32){
			for(int i=l;i<=r;i++)ans[i]=value(a,x[i]);
			return;
		}
		int m=(l+r)>>1;
		poly q,b;
		Div(a,g[u*2],q,b);self(self,move(b),u*2,l,m);
		Div(move(a),g[u*2+1],q,b);self(self,move(b),u*2+1,m+1,r);
	};
	solve(solve,move(f),1,0,n-1);
	return ans;
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
