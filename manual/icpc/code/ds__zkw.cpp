struct zkw_segt{
//////////////////////////////////////////////////////////////
// information
struct info{
	// modify here
	ll s,mx;
	info():s(0),mx(0){}
	info(ll _s,ll _mx):s(_s),mx(_mx){}
	inline info& operator +=(const info &y){
		ll os=s;
		s+=y.s;
		mx=max({mx,y.mx,os+y.mx});
		return *this;
	}
	inline info operator +(const info &y)const{
		return info(*this)+=y;
	}
};
/////////////////////////////////////////////////////////////////
// basic segment tree template, private operations
private:
int n=0,B=0;vector<info> a;
void modify(int x,info v){
	x+=B+1;
	a[x]=v;
	while(x>1) x>>=1,a[x]=a[x<<1]+a[x<<1|1]; 
}
info query(int l,int r){
	l+=B,r+=B+2;
	info ql,qr;
	while(l^r^1){
		if(~l&1) ql=ql+a[l^1];
		if(r&1) qr=a[r^1]+qr;
		l>>=1,r>>=1;
	}
	return ql+qr;
}
/////////////////////////////////////////////////////////////////
// basic segment tree template, public operations
// please transfer (int,int) to info here
public:
void setN(int _n){
	assert(1<=_n&&_n<(1<<29));n=_n;B=1;
	while(B<n+2)B*=2;
	a.assign(2*B,info());
}
void upd(int x,ll a,ll b){
	assert(0<=x&&x<n);
	modify(x,info(a,b));
}
ll ask(int l,int r){
	assert(0<=l&&l<=r&&r<n);
	return query(l,r).mx;
}
};
