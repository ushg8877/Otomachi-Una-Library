#define LL __int128
#define MP make_pair
#define poly vector<long long>
////////////////////////////////////////////////////////////////
//
// template for polynomial (ll NTT)
// version 1.0 (last update: May 24th, 2026)
//
// usage:
//   init(n); poly f,g; f*g; ln(f);
//
////////////////////////////////////////////////////////////////
const ll MOD=998244353,inv2=(MOD+1)/2;
vector<ll> fac{1},ifac{1},Inv{0,1};
vector<int> bfly(1);vector<ll> inver{0,1};
int clogg(int x){assert(x>0);return x==1?0:32-__builtin_clz(x-1);} 
ll ksm(ll a,int b){ll res=1;while(b){if(b&1)res=res*a%MOD;a=a*a%MOD,b>>=1;}return res;}
void butterfly(int l){
	assert(0<=l&&l<=23);
	static int las=-1;
	if(las!=l){
		las=l;bfly.assign(1<<l,0); 
		for(int i=1;i<(1<<l);i++)
			bfly[i]=(bfly[i>>1]>>1)|((i&1)<<l-1);
	} 
}
void NTT(poly &f,int l,int typ){
	butterfly(l);f.resize(1<<l);
	for(int i=0;i<(1<<l);i++)
		if(bfly[i]<i) swap(f[i],f[bfly[i]]);
	for(int i=0;i<l;i++){
		ll step=ksm(3,MOD-1+(MOD-1>>i+1)*typ);
		for(int j=0;j<(1<<l);j+=(1<<i+1)){
			ll cur=1,l=j+(1<<i);
			for(int k=j;k<l;k++){
				ll u=f[k],v=f[k|(1<<i)]*cur%MOD;
				f[k]=(u+v>=MOD?u+v-MOD:u+v);
				f[k|(1<<i)]=(u>=v?u-v:u-v+MOD);
				cur=cur*step%MOD;
			}
		}
	}
	if(typ==-1){
		ll val=ksm(1<<l,MOD-2);
		for(int i=0;i<(1<<l);i++)
			f[i]=val*f[i]%MOD;
	}
	return;
}
poly operator *(poly f,poly g){
	// h[i+j] <- f[i]*g[j]
	if(f.empty()||g.empty())return {};
	assert(f.size()+g.size()-1<=(1<<23));
	int n=f.size()+g.size(),l=clogg(f.size()+g.size()-1);
	if(n<=64){
		poly h(n-1);
		for(int i=0;i<f.size();i++)
			for(int j=0;j<g.size();j++)
				h[i+j]+=f[i]*g[j]%MOD;
		for(ll &i:h) i%=MOD;
		return h;
	}
	NTT(f,l,1);NTT(g,l,1);
	for(int i=0;i<(1<<l);i++)
		f[i]=f[i]*g[i]%MOD;
	NTT(f,l,-1);f.resize(n-1);
	return f;
}
poly operator /(poly f,poly g){
	if(f.empty()||g.empty())return poly(f.size());
	// h[i-j] <- f[i]*g[j]
	int m=g.size();
	reverse(g.begin(),g.end());
	g=f*g;
	for(int i=0;i<f.size();i++) f[i]=g[i+m-1];
	return f;
}
void operator +=(poly &f,const poly &g){
	int n=max(f.size(),g.size());
	f.resize(n);
	for(int i=0;i<g.size();i++)
		f[i]=(f[i]+g[i]>=MOD?f[i]+g[i]-MOD:f[i]+g[i]);
	return;
}
poly operator +(poly f,poly g){
	f+=g;
	return f;
}
void operator -=(poly &f,const poly &g){
	int n=max(f.size(),g.size());
	f.resize(n);
	for(int i=0;i<g.size();i++)
		f[i]=(f[i]-g[i]<0?f[i]-g[i]+MOD:f[i]-g[i]);
	return;
}
poly operator -(poly f,poly g){
	f-=g;
	return f;
}
poly operator *(poly f,ll x){
	for(int i=0;i<f.size();i++) (f[i]*=x)%=MOD;
	return f;
}
poly operator *(ll x,poly f){
	for(int i=0;i<f.size();i++) (f[i]*=x)%=MOD;
	return f;
}
poly inv(poly f){
	assert(!f.empty()&&f[0]);
	int n=f.size(),l=clogg(n);
	f.resize(1<<l);
	poly g{ksm(f[0],MOD-2)},_f;
	for(int i=1;i<=l;i++){
		_f=poly(begin(f),begin(f)+(1<<i));
		NTT(g,i+1,1);NTT(_f,i+1,1);
		for(int j=0;j<(1<<i+1);j++)
			g[j]=(2*g[j]-g[j]*g[j]%MOD*_f[j]%MOD+MOD)%MOD;
		NTT(g,i+1,-1);
		fill(begin(g)+(1<<i),end(g),0);
	}
	g.resize(n);
	return g;
}
poly integ(poly f){
	int n=f.size();f.resize(n+1);
	assert(n<MOD);
	int m=inver.size();
	if(n>=m){
		int k=max(n,(int)min(2ll*m,MOD-1));inver.resize(k+1);
		for(int i=m;i<=k;i++)inver[i]=MOD-(MOD/i)*inver[MOD%i]%MOD;
	}
	for(int i=n;i>=1;i--) f[i]=f[i-1]*inver[i]%MOD;
	f[0]=0;
	return f;
}
poly diff(poly f){
	if(f.empty())return {};
	int n=f.size();
	for(int i=0;i<n-1;i++)
		f[i]=f[i+1]*(i+1)%MOD;
	f.pop_back();
	return f;
} 
poly ln(poly f){
	assert(!f.empty()&&f[0]==1);
	poly f_=diff(f),_f=inv(f);
	int n=f_.size(),m=_f.size();
	int l=clogg(n+m);
	NTT(f_,l,1);NTT(_f,l,1);
	for(int i=0;i<(1<<l);i++)
		f_[i]=f_[i]*_f[i]%MOD;
	NTT(f_,l,-1);
	f_=integ(f_);
	f_.resize(f.size());
	return f_;
}
poly exp(poly f){
	assert(!f.empty()&&f[0]==0);
	poly g{1},_f,_g;
	int n=f.size(),l=clogg(n);
	f.resize(1<<l);
	for(int i=1;i<=l;i++){
		_f=poly(begin(f),begin(f)+(1<<i));
		_g=ln(g);
		NTT(g,i+1,1);NTT(_f,i+1,1);NTT(_g,i+1,1);
		for(int j=0;j<(1<<i+1);j++)
			g[j]=g[j]*(1-_g[j]+_f[j]+MOD)%MOD;
		NTT(g,i+1,-1);
		fill(begin(g)+(1<<i),end(g),0);
	}
	g.resize(n);
	return g;
}
poly sqrt(poly f){
	assert(!f.empty()&&f[0]==1);
	poly g{1},_f,_g;
	int n=f.size(),l=clogg(n);
	f.resize(1<<l);
	for(int i=1;i<=l;i++){
		_f=poly(begin(f),begin(f)+(1<<i));
		_g=inv(g);
		NTT(_f,i+1,1);NTT(_g,i+1,1);NTT(g,i+1,1);
		for(int j=0;j<(1<<i+1);j++)
			g[j]=(_f[j]+g[j]*g[j]%MOD)*inv2%MOD*_g[j]%MOD;
		NTT(g,i+1,-1);
		fill(begin(g)+(1<<i),end(g),0);
	}
	g.resize(n);
	return g;
} 
void div(poly f,poly g,poly &q,poly &r){// ¹¹Ôì p,q Ê¹µÃ f=g*q+r; 
	assert(!g.empty()&&g.back());
	if(f.size()<g.size()){q={0};r=f;return;}
	int n=f.size()-1,m=g.size()-1;
	reverse(f.begin(),f.end());
	reverse(g.begin(),g.end());
	g.resize(n+1);
	q=(f*inv(g));q.resize(n-m+1);
	reverse(q.begin(),q.end());
	g.resize(m+1);
	reverse(g.begin(),g.end());
	reverse(f.begin(),f.end());
	poly h=q*g;
	r.resize(m);
	for(int i=0;i<m;i++) r[i]=(f[i]-h[i]+MOD)%MOD;
	return;
}
poly eval(poly f,poly x){
	if(x.empty())return {};
	vector<poly> func; map<pair<int,int>,int> mp;
	function<poly(int,int)> prod=[&](int l,int r){
		poly ans;
		if(l==r) ans=(poly){MOD-x[l],1};
		else {
			int mid=l+r>>1;
			ans=prod(l,mid)*prod(mid+1,r);
		} 
		func.push_back(ans);
		mp[MP(l,r)]=func.size()-1;
		return ans;
	};prod(0,x.size()-1);
	int t=-1;
	function<poly(poly,int,int)> solve=[&](poly f,int l,int r){
		if(l==r){
			ll val=0;
			for(int i=f.size()-1;i>=0;i--)
				val=(val*x[l]+f[i])%MOD;
			return (poly){val};
		}
		int mid=l+r>>1;
		poly p0=func[mp[MP(l,mid)]],p1=func[mp[MP(mid+1,r)]];
		poly d,q;
		div(f,p0,d,q);
		poly ans=solve(q,l,mid);
		div(f,p1,d,q);
		poly _ans=solve(q,mid+1,r);
		for(int i:_ans) ans.push_back(i);
		return ans;
	};
	return solve(f,0,x.size()-1);
}
poly BM(poly a){
	poly c{1},b{1};int l=0,m=1;ll last=1;
	for(int i=0;i<(int)a.size();i++){
		ll d=0;for(int j=0;j<=l;j++)d=(d+c[j]*a[i-j])%MOD;
		if(!d){m++;continue;}
		poly t=c;ll v=d*ksm(last,MOD-2)%MOD;
		if(c.size()<b.size()+m)c.resize(b.size()+m);
		for(int j=0;j<(int)b.size();j++)c[j+m]=(c[j+m]-v*b[j]%MOD+MOD)%MOD;
		if(2*l<=i){l=i+1-l;b=move(t);last=d;m=1;}
		else m++;
	}
	return c;
}
ll FSPE(poly F,poly G,ll t){// [x^t]F(x)/G(x)
	assert(t>=0&&!G.empty()&&G[0]);
	if(F.empty())return 0;
	while(t){
		poly G1=G;
		for(int i=1;i<G1.size();i+=2) G1[i]=MOD-G1[i];
		F=F*G1,G=G*G1;
		poly f,g;
		for(int i=t&1;i<F.size();i+=2) f.push_back(F[i]);
		for(int i=0;i<G.size();i+=2) g.push_back(G[i]);
		t>>=1;
		F=f,G=g;
		if(F.empty())return 0;
	}
	return F[0]*ksm(G[0],MOD-2)%MOD;

}
ll RSPE(poly f,ll x){
	poly g=BM(f),s(g.size()-1);
	for(int i=0;i<(int)s.size();i++)
		for(int j=0;j<=i&&j<(int)f.size();j++)s[i]=(s[i]+f[j]*g[i-j])%MOD;
	return FSPE(s,g,x);
}
mt19937 rnd(time(0));
void init(int n=0){
	assert(0<=n&&n<MOD);
	int m=fac.size();if(n<m)return;
	n=max(n,(int)min(2ll*m,MOD-1));
	fac.resize(n+1);ifac.resize(n+1);Inv.resize(max(n+1,2));
	for(int i=m;i<=n;i++){
		fac[i]=fac[i-1]*i%MOD;
		if(i>=2)Inv[i]=MOD-Inv[MOD%i]*(MOD/i)%MOD;
		ifac[i]=ifac[i-1]*Inv[i]%MOD;
	}
}
ll C(int x,int y){if(x<0||x>y)return 0;init(y);return fac[y]*ifac[x]%MOD*ifac[y-x]%MOD;}
// end for polynomial/ntt-int64.cpp
/////////////////////////
// !!!!! Choose only one polynomial implementation; it already defines modular arithmetic. Check constant terms before inv/ln/exp and the transform size limit. !!!!
