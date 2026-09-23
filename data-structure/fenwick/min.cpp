////////////////////////////////////////////////////////////////
//
// template for BIT (min)
//
// usage:
//   Min_BIT t; t.setRange(l,r); // setN(n) = setRange(1,n)
//   upd(i,v) is chmin, not assignment; ask(r) queries [vl,r].
//   Initial / empty-prefix value: LLONG_MAX. No inverse for interval queries.
//
////////////////////////////////////////////////////////////////
struct Min_BIT{
private:
int vl=1,vr=0,n=0;vector<ll> a;
public:
void set(){vl=1;vr=n=0;a.clear();}
void setRange(int l,int r){
	ll m=(ll)r-l+1;assert(0<=m&&m<INT_MAX-1);
	vl=l;vr=r;n=m;a.assign(n+2,LLONG_MAX);
}
void setN(int n){setRange(1,n);}
void upd(ll x,ll v){
	assert(vl<=x&&x<=vr);
	for(ll i=x-vl+1;i<=n;i+=i&-i)a[i]=min(a[i],v);
}
ll ask(ll x)const{
	assert((ll)vl-1<=x&&x<=vr);ll ans=LLONG_MAX;
	for(int i=x-vl+1;i;i-=i&-i)ans=min(ans,a[i]);
	return ans;
}

};
// end for data-structure/fenwick/min.cpp
/////////////////////////
// !!!!! upd is chmin; arbitrary assignment and interval subtraction are not supported. !!!!
