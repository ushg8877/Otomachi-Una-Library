template<unsigned P>
vector<ModInt<P>> convolution(const poly &f,const poly &g,int n){
	vector<ModInt<P>> a(n),b(n);
	for(int i=0;i<(int)f.size();i++)a[i]=f[i].x;
	for(int i=0;i<(int)g.size();i++)b[i]=g[i].x;
	ntt(a);ntt(b);
	for(int i=0;i<n;i++)a[i]*=b[i];
	ntt(a,true);a.resize(f.size()+g.size()-1);return a;
}
poly operator*(const poly &f,const poly &g){
	if(f.empty()||g.empty())return {};
	int m=f.size()+g.size()-1,n=1;
	assert(m<=(1<<FFT_MAX));
	poly h(m);
	if(min(f.size(),g.size())<=32){
		for(int i=0;i<(int)f.size();i++)
			for(int j=0;j<(int)g.size();j++)h[i+j]+=f[i]*g[j];
		return h;
	}
	constexpr unsigned p1=998244353,p2=167772161,p3=469762049;
	constexpr ull p12=1ull*p1*p2;
	static const unsigned v1=ModInt<p2>(p1).inv().x;
	static const unsigned v2=ModInt<p3>(p12).inv().x;
	while(n<m)n<<=1;
	auto a=convolution<p1>(f,g,n);
	auto b=convolution<p2>(f,g,n);
	auto c=convolution<p3>(f,g,n);
	for(int i=0;i<m;i++){
		ull x=(b[i].x+p2-a[i].x%p2)*1ull*v1%p2;
		x=x*p1+a[i].x;
		ull y=(c[i].x+p3-x%p3)*v2%p3;
		h[i]=x%MOD+(p12%MOD)*y%MOD;
	}
	return h;
}
