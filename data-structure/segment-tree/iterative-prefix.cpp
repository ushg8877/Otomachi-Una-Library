////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/data-structure/segment-tree/iterative-prefix.cpp
//
// usage: zkw_segt z; z.setN(n); z.upd(i, s, mx); z.ask(l, r);
//   0-indexed; leaf mx=max(0ll,s); ask returns maximum prefix sum
//
////////////////////////////////////////////////////////////////
struct zkw_segt{
//////////////////////////////////////////////////////////////
// information
struct info{
	// sum and maximum prefix sum, including the empty prefix
	ll s,mx;
	info():s(0),mx(0){}
	info(ll _s,ll _mx):s(_s),mx(_mx){}
	inline info& operator +=(const info &y){
		ll os=s;
		s+=y.s;
		mx=max(mx,os+y.mx);
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
	assert(1<=_n&&_n<(1<<29));
	n=_n;
	B=1;
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
// Largest valid r in [l,r]; return {r,info}, empty r=l-1. O(log n).
// f(empty) must be true; extension may only change true to false.
template<class F>
pair<ll,info> max_right(ll l,F f)const{
	assert(0<=l&&l<=n);
	info s;bool ok=f(s);assert(ok);(void)ok;
	if(l==n) return {l-1,s};
	int p=l+B+1;
	do{
		while(!(p&1)) p>>=1;
		if(!f(s+a[p])){
			while(p<B){
				p<<=1;
				info t=s+a[p];
				if(f(t)){s=t;p++;}
			}
			return {p-B-2,s};
		}
		s=s+a[p];p++;
	}while((p&-p)!=p);
	return {n-1,s};
}
// Smallest valid l in [l,r]; return {l,info}, empty l=r+1.
// Same predicate requirements; the information keeps left-to-right order.
template<class F>
pair<ll,info> min_left(ll r,F f)const{
	assert(-1<=r&&r<n);
	info s;bool ok=f(s);assert(ok);(void)ok;
	if(r==-1) return {0,s};
	int p=r+B+2;
	do{
		p--;
		while(p>1&&(p&1)) p>>=1;
		if(!f(a[p]+s)){
			while(p<B){
				p=p<<1|1;
				info t=a[p]+s;
				if(f(t)){s=t;p--;}
			}
			return {p-B,s};
		}
		s=a[p]+s;
	}while((p&-p)!=p);
	return {0,s};
}

};
// end for data-structure/segment-tree/iterative-prefix.cpp
/////////////////////////
// !!!!! Search predicates must be monotone and accept empty info. !!!!
