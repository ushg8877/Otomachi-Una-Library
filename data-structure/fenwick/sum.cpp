////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: Sum_BIT t; t.setRange(l,r); // setN(n) = setRange(1,n)
//   add(i,v); ask(r); ask(l,r); initial values are zero.
//
////////////////////////////////////////////////////////////////
struct Sum_BIT{
private:
int vl=1,vr=0,n=0;vector<ll> a;
public:
void set(){vl=1;vr=n=0;a.clear();}
void setRange(int l,int r){
	ll m=(ll)r-l+1;assert(0<=m&&m<INT_MAX-1);
	vl=l;
	vr=r;
	n=m;
	a.assign(n+2,0);
}
void setN(int n){setRange(1,n);}
// Point add; vr+1 is allowed as a difference-array sentinel.
void add(ll x,ll v){
	assert(!a.empty()&&vl<=x&&x<=(ll)vr+1);
	for(ll i=x-vl+1;i<=n+1;i+=i&-i)a[i]+=v;
}
ll ask(ll x)const{
	assert((ll)vl-1<=x&&x<=vr);ll ans=0;
	for(int i=x-vl+1;i;i-=i&-i)ans+=a[i];
	return ans;
}
ll ask(int l,int r)const{assert(vl<=l&&l<=r&&r<=vr);return ask(r)-ask((ll)l-1);}
};
// end for data-structure/fenwick/sum.cpp
/////////////////////////
