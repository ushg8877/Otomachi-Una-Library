////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: BIT: setRange(l,r); add(i,v); ask(r); ask(l,r).
//   setN(n) = setRange(1,n); add(vr+1,v) is a difference sentinel.
//
////////////////////////////////////////////////////////////////
struct BIT{
private:
struct info{
	ll x=0;
	info(){}
	info(ll _x):x(_x){}
	info& operator +=(const info &y){x+=y.x;return *this;}
	info operator +(const info &y)const{return info(*this)+=y;}
};
int vl=1,vr=0,n=0;vector<info> a;
info query(int x)const{
	info ans;
	for(;x;x-=x&-x)ans+=a[x];
	return ans;
}
public:
void set(){vl=1;vr=n=0;a.clear();}
void setRange(int l,int r){
	ll m=(ll)r-l+1;assert(0<=m&&m<INT_MAX-1);
	vl=l;
	vr=r;
	n=m;
	a.assign(n+2,info());
}
void setN(int n){setRange(1,n);}
void add(ll x,ll v){
	assert(!a.empty()&&vl<=x&&x<=(ll)vr+1);
	for(ll i=x-vl+1;i<=n+1;i+=i&-i)a[i]+=info(v);
}
ll ask(ll x)const{assert((ll)vl-1<=x&&x<=vr);return query(x-vl+1).x;}
ll ask(int l,int r)const{assert(vl<=l&&l<=r&&r<=vr);return ask(r)-ask((ll)l-1);}
};
// Range add / range sum: setRange(l,r), add(l,r,v), ask(l,r).
struct rars{
int vl=1,vr=0;BIT T0,T1;
void set(){vl=1;vr=0;T0.set();T1.set();}
void setRange(int l,int r){vl=l;vr=r;T0.setRange(l,r);T1.setRange(l,r);}
void setN(int n){setRange(1,n);}
// Point add; vr+1 is allowed as a difference-array sentinel.
void add(int l,int r,ll x){
	assert(vl<=l&&l<=r&&r<=vr);
	ll a=(ll)l-vl+1,b=(ll)r-vl+2;
	T0.add(l,x);T0.add((ll)r+1,-x);
	T1.add(l,a*x);T1.add((ll)r+1,-b*x);
}
ll query(ll x)const{return (x-vl+2)*T0.ask(x)-T1.ask(x);}
ll ask(int l,int r)const{
	assert(vl<=l&&l<=r&&r<=vr);
	return query(r)-query((ll)l-1);
}
}T;
// end for data-structure/fenwick/fenwick.cpp
/////////////////////////
