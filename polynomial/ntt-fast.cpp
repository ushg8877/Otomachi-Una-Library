////////////////////////////////////////////////////////////////
//
// template for polynomial modulo 998244353
//
// usage: poly h=f*g; Inv(f); Ln(f); Exp(f);
//   Sqrt(f); Div(f,g,q,r); Eval(f,x); FSPE(f,g,k);
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


// Newton transform reuse follows QedDust413 & Killer_joke.
using poly=vector<mint>;
using Poly=poly;

constexpr unsigned MO=998244353U;
constexpr unsigned MO2=2U*MO;
constexpr int FFT_MAX=23;
struct NTT{
struct Root{
	unsigned x,y;
	Root(unsigned v=1):x(v),y((1ull*v<<32)/MO){}
	unsigned mul(unsigned v)const{
		return v*x-unsigned((1ull*v*y)>>32)*MO;
	}
};
vector<Root>w{1,1},iw{1,1};
static unsigned norm(unsigned x){return x>=MO2?x-MO2:x;}
void setN(int n){
	for(int k=w.size();k<n;k<<=1){
		w.resize(k*2);iw.resize(k*2);
		mint x=mint(3).pow((MO-1)/(k*2)),y=x.inv();
		for(int i=k/2;i<k;i++){
			w[i*2]=w[i];w[i*2+1]=Root((mint(w[i].x)*x).x);
			iw[i*2]=iw[i];iw[i*2+1]=Root((mint(iw[i].x)*y).x);
		}
	}
}
// DIF/DIT, two stages per pass; Shoup products stay below 2*MOD.
void fft(mint *a,int n){
	setN(n);
	int len=n;
	for(;len>=4;len>>=2){
		int k=len>>2;
		for(int l=0;l<n;l+=len)for(int j=0;j<k;j++){
			unsigned x=a[l+j].x,y=a[l+j+k].x;
			unsigned z=a[l+j+2*k].x,t=a[l+j+3*k].x;
			unsigned u=norm(x+z),v=norm(y+t);
			x=w[2*k+j].mul(x+MO2-z);
			y=w[3*k+j].mul(y+MO2-t);
			a[l+j].x=norm(u+v);
			a[l+j+k].x=w[k+j].mul(u+MO2-v);
			a[l+j+2*k].x=norm(x+y);
			a[l+j+3*k].x=w[k+j].mul(x+MO2-y);
		}
	}
	if(len==2)for(int i=0;i<n;i+=2){
		unsigned x=a[i].x,y=a[i+1].x;
		a[i].x=norm(x+y);a[i+1].x=norm(x+MO2-y);
	}
	for(int i=0;i<n;i++)if(a[i].x>=MO)a[i].x-=MO;
}
void invFft(mint *a,int n){
	setN(n);
	int len=1;
	if(__builtin_ctz(unsigned(n))&1){
		for(int i=0;i<n;i+=2){
			unsigned x=a[i].x,y=a[i+1].x;
			a[i].x=x+y;a[i+1].x=x+MO-y;
		}
		len=2;
	}
	for(;len<n;len<<=2){
		int k=len;
		for(int l=0;l<n;l+=k*4)for(int j=0;j<k;j++){
			unsigned x=a[l+j].x,y=iw[k+j].mul(a[l+j+k].x);
			unsigned z=a[l+j+2*k].x,t=iw[k+j].mul(a[l+j+3*k].x);
			unsigned u=norm(x+y),v=norm(x+MO2-y);
			x=iw[2*k+j].mul(norm(z+t));
			y=iw[3*k+j].mul(norm(z+MO2-t));
			a[l+j].x=norm(u+x);a[l+j+2*k].x=norm(u+MO2-x);
			a[l+j+k].x=norm(v+y);a[l+j+3*k].x=norm(v+MO2-y);
		}
	}
	Root v(mint(n).inv().x);
	for(int i=0;i<n;i++){
		unsigned x=v.mul(a[i].x);
		a[i].x=x>=MO?x-MO:x;
	}
}
};
NTT ntt;
// Nonempty power-of-two length, <=2^23; spectrum uses bit-reversed order.
void fft(mint *a,int n){
	assert(n>0&&!(n&(n-1))&&n<=(1<<FFT_MAX));
	ntt.fft(a,n);
}
void invFft(mint *a,int n){
	assert(n>0&&!(n&(n-1))&&n<=(1<<FFT_MAX));
	ntt.invFft(a,n);
}
void fft(poly &a){fft(a.data(),a.size());}
void invFft(poly &a){invFft(a.data(),a.size());}

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

// Formal inverse; constant term must be invertible.
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

// f[0]=0; maintain g=exp(f), h=1/g and reuse their transforms.
poly Exp(poly f){
	assert(!f.empty()&&!f[0]&&f.size()<=(1<<FFT_MAX));
	int n=f.size(),l=1;
	while(l<n)l<<=1;
	init(l);f.resize(l);
	poly g(l),h(l),a(l),b(l),c(l);
	g[0]=h[0]=a[0]=1;
	if(l>1)a[1]=1;
	for(int m=1;m<l;m<<=1){
		int t=m<<1;
		for(int i=0;i<m;i++)c[i]=f[i]*i;
		fft(c.data(),m);
		for(int i=0;i<m;i++)c[i]*=a[i];
		invFft(c.data(),m);
		for(int i=0;i<m;i++)c[i+m]=g[i]*i-c[i],c[i]=0;
		copy(h.begin(),h.begin()+m,b.begin());
		fill(b.begin()+m,b.begin()+t,0);
		fft(c.data(),t);fft(b.data(),t);
		for(int i=0;i<t;i++)c[i]*=b[i];
		invFft(c.data(),t);
		fill(c.begin(),c.begin()+m,0);
		for(int i=m;i<t;i++)c[i]=c[i]*inv[i]-f[i];
		fft(c.data(),t);
		for(int i=0;i<t;i++)a[i]*=c[i];
		invFft(a.data(),t);
		for(int i=m;i<t;i++)g[i]=a[i]=-a[i];
		if(t==l)break;
		copy(g.begin(),g.begin()+m,a.begin());
		fill(a.begin()+t,a.begin()+2*t,0);
		fft(a.data(),2*t);
		for(int i=0;i<t;i++)c[i]=a[i]*b[i];
		invFft(c.data(),t);
		fill(c.begin(),c.begin()+m,0);
		fft(c.data(),t);
		for(int i=0;i<t;i++)c[i]*=b[i];
		invFft(c.data(),t);
		for(int i=m;i<t;i++)h[i]=-c[i];
	}
	g.resize(n);
	return g;
}

// f[0]=1; choose the square root with constant term 1.
poly Sqrt(poly f){
	assert(!f.empty()&&f[0].x==1&&f.size()<=(1<<FFT_MAX));
	int n=f.size(),l=1;
	while(l<n)l<<=1;
	f.resize(l);
	poly g(l),h(l),a(l),b(l),c(l);
	g[0]=h[0]=a[0]=1;
	mint v=499122177;
	for(int m=1;m<l;m<<=1){
		int t=m<<1;
		for(int i=0;i<m;i++)a[i]*=a[i];
		invFft(a.data(),m);
		for(int i=0;i<m;i++)a[i+m]=a[i]-f[i]-f[i+m],a[i]=0;
		copy(h.begin(),h.begin()+m,b.begin());
		fill(b.begin()+m,b.begin()+t,0);
		fft(a.data(),t);fft(b.data(),t);
		for(int i=0;i<t;i++)a[i]*=b[i];
		invFft(a.data(),t);
		for(int i=m;i<t;i++)g[i]=-a[i]*v;
		if(t==l)break;
		copy(g.begin(),g.begin()+t,a.begin());
		fft(a.data(),t);
		for(int i=0;i<t;i++)c[i]=a[i]*b[i];
		invFft(c.data(),t);
		fill(c.begin(),c.begin()+m,0);
		fft(c.data(),t);
		for(int i=0;i<t;i++)c[i]*=b[i];
		invFft(c.data(),t);
		for(int i=m;i<t;i++)h[i]=-c[i];
	}
	g.resize(n);
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
// use init(n) before accessing fac / ifac / inv directly
// end for polynomial/ntt-fast.cpp
/////////////////////////
// !!!!! Contains mint already; do not paste another mint or poly
// implementation. !!!!
