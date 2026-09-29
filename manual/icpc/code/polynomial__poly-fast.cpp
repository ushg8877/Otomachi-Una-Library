using poly=vector<mint>;
using Poly=poly;
const int FFT_MAX=23;
// Shared NTT; use in both the 998244353 and CRT versions.
template<unsigned P>
void ntt(vector<ModInt<P>> &a,bool rev=false){
	using num=ModInt<P>;
	int n=a.size();
	assert(n&&!(n&(n-1))&&(P-1)%n==0&&n<=(1<<FFT_MAX));
	for(int i=1,j=0;i<n;i++){
		int k=n>>1;
		for(;j&k;k>>=1)j^=k;
		j^=k;if(i<j)swap(a[i],a[j]);
	}
	static vector<num> w{0,1};
	for(int m=w.size();m<n;m<<=1){
		w.resize(2*m);num z=num(3).pow((P-1)/(2*m));
		for(int i=m/2;i<m;i++)w[2*i]=w[i],w[2*i+1]=w[i]*z;
	}
	for(int m=1;m<n;m<<=1)
		for(int i=0;i<n;i+=2*m)
			for(int j=0;j<m;j++){
				num x=a[i+j],y=a[i+j+m]*w[m+j];
				a[i+j]=x+y;a[i+j+m]=x-y;
			}
	if(rev){
		reverse(a.begin()+1,a.end());num v=num(n).inv();
		for(num &x:a)x*=v;
	}
}
// 998244353 convolution. Replace this function by the MTT core if needed.
poly operator*(poly f,poly g){
	if(f.empty()||g.empty())return {};
	int m=f.size()+g.size()-1,n=1;
	assert(m<=(1<<FFT_MAX));
	if(min(f.size(),g.size())<=32){
		poly h(m);
		for(int i=0;i<(int)f.size();i++)
			for(int j=0;j<(int)g.size();j++)h[i+j]+=f[i]*g[j];
		return h;
	}
	while(n<m)n<<=1;
	f.resize(n);g.resize(n);ntt(f);ntt(g);
	for(int i=0;i<n;i++)f[i]*=g[i];
	ntt(f,true);f.resize(m);return f;
}

// Shared FPS; also used with the MTT convolution.
// Transposed multiplication: (f/g)[i] = sum_j f[i+j]*g[j].
poly operator/(const poly &f,poly g){
	if(f.empty()||g.empty())return poly(f.size());
	int m=g.size();reverse(g.begin(),g.end());g=f*g;
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
poly operator*(poly f,mint x){for(auto &y:f)y*=x;return f;}
poly operator*(mint x,poly f){for(auto &y:f)y*=x;return f;}
mint value(const poly &f,mint x){
	mint ans=0;for(int i=(int)f.size()-1;i>=0;i--)ans=ans*x+f[i];return ans;
}
poly Inv(const poly &f){
	assert(!f.empty()&&f[0]&&f.size()<=(1<<(FFT_MAX-1)));
	int n=f.size();poly g{f[0].inv()};
	while((int)g.size()<n){
		int m=min(n,(int)g.size()*2);
		poly h=poly(f.begin(),f.begin()+m)*g;h.resize(m);
		for(auto &x:h)x=-x;
		h[0]+=2;g=g*h;g.resize(m);
	}
	return g;
}
poly diff(poly f){
	if(f.empty())return {};
	for(int i=1;i<(int)f.size();i++)f[i-1]=f[i]*i;
	f.pop_back();return f;
}
poly integ(poly f){
	int n=f.size();init(n);f.resize(n+1);
	for(int i=n;i;i--)f[i]=f[i-1]*inv[i];
	f[0]=0;return f;
}
poly Ln(const poly &f){
	assert(!f.empty()&&f[0].x==1&&f.size()<(size_t)MOD&&f.size()<=(1<<(FFT_MAX-1)));
	poly g=diff(f)*Inv(f);g.resize(f.size()-1);return integ(move(g));
}
poly Exp(const poly &f){
	assert(!f.empty()&&!f[0]&&f.size()<(size_t)MOD&&f.size()<=(1<<(FFT_MAX-1)));
	int n=f.size();poly g{1};
	while((int)g.size()<n){
		int m=min(n,(int)g.size()*2);poly h=g;h.resize(m);h=Ln(h);
		for(int i=0;i<m;i++)h[i]=f[i]-h[i];
		h[0]+=1;g=g*h;g.resize(m);
	}
	return g;
}
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

// Return [x^k] f/g; g[0] must be invertible.
mint FSPE(poly f,poly g,ll k){
	assert(k>=0&&!g.empty()&&g[0]);
	while(k&&!f.empty()){
		poly h=g;
		for(int i=1;i<(int)h.size();i+=2)h[i]=-h[i];
		f=f*h;g=g*h;
		int n=(f.size()+1-(k&1))/2,m=(g.size()+1)/2;
		for(int i=0;i<n;i++)f[i]=f[2*i+(k&1)];
		for(int i=0;i<m;i++)g[i]=g[2*i];
		f.resize(n);g.resize(m);k>>=1;
	}
	return f.empty()?mint(0):f[0]/g[0];
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

poly lagrange(const vector<mint> &x,const vector<mint> &y){
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
