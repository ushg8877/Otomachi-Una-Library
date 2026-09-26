////////////////////////////////////////////////////////////////
//
// template for static / compressed Li Chao tree
// usage:
//   struct F{
//       using value_type=ll;
//       ll k,b;int id;
//       ll operator()(ll x)const{return k*x+b;}
//   };
//   LiChao<F> lc; lc.set({1,4,9},F{0,LLONG_MAX,-1});
//   lc.add(F{2,3,0}); lc.add(F{-1,5,1},2,10); // [l,r)
//   auto [y,f]=lc.ask(4); // value and function; f.id is the source
//   LiChao<F,false,false> mx; mx.setRange(-100,101,F{0,LLONG_MIN,-1});
//   set() clears functions and keeps the coordinates / default function.
//   O(n) space; add / ask O(log n), interval add O(log^2 n).
//
////////////////////////////////////////////////////////////////
template<typename F,bool COMPRESS=true,bool MINIMIZE=true> struct LiChao{
using T=typename F::value_type;
int n=0,sz=1;
ll lo=0,hi=0;
vector<ll> X;
vector<F> a;
void set(const vector<ll> &x,const F &f){
	static_assert(COMPRESS);X=x;
	sort(X.begin(),X.end());X.erase(unique(X.begin(),X.end()),X.end());
	assert(!X.empty()&&X.size()<=(1<<29));n=X.size();build(f);
}
void setRange(ll l,ll r,const F &f){
	static_assert(!COMPRESS);
	assert(l<r&&(__int128)r-l<=(1<<29));lo=l;hi=r;n=r-l;build(f);
}
void set(){assert(n);fill(a.begin()+1,a.end(),a[0]);}
void add(F f){assert(n);add(1,move(f));}
void add(F f,ll l,ll r){
	assert(n&&l<=r);
	if constexpr(!COMPRESS){assert(lo<=l&&r<=hi);}
	int L=idx(l)+sz,R=idx(r)+sz;
	while(L<R){
		if(L&1)add(L++,f);
		if(R&1)add(--R,f);
		L>>=1;R>>=1;
	}
}
pair<T,F> ask(ll x)const{
	assert(n);int p=idx(x);
	if constexpr(COMPRESS){assert(p<n&&X[p]==x);}
	else{assert(lo<=x&&x<hi);}
	F f=a[0];T ans=f(x);
	for(int i=p+sz;i;i>>=1){
		T v=a[i](x);if(better(v,ans))ans=v,f=a[i];
	}
	return {ans,f};
}
private:
void build(const F &f){sz=1;while(sz<n)sz<<=1;a.assign(sz*2,f);}
static bool better(const T &x,const T &y){if constexpr(MINIMIZE)return x<y;else return y<x;}
int idx(ll x)const{
	if constexpr(COMPRESS)return lower_bound(X.begin(),X.end(),x)-X.begin();
	else{assert(lo<=x&&x<=hi);return x-lo;}
}
ll coord(int i)const{
	i=min(i,n-1);
	if constexpr(COMPRESS)return X[i];else return lo+i;
}
void add(int i,F f){
	int k=31-__builtin_clz((unsigned)i),len=sz>>k;
	int l=(i-(1<<k))*len,r=l+len;
	while(l<r){
		bool bl=better(f(coord(l)),a[i](coord(l)));
		bool br=better(f(coord(r-1)),a[i](coord(r-1)));
		if(bl&&br){a[i]=move(f);return;}
		if(!bl&&!br)return;
		int m=(l+r)>>1;
		bool bm=better(f(coord(m)),a[i](coord(m)));
		if(bm)swap(a[i],f);
		if(bl!=bm)i=i*2,r=m;else i=i*2+1,l=m;
	}
}
};
// end for data-structure/li-chao-tree-static.cpp
/////////////////////////
// !!!!! F needs value_type and const operator()(ll). Pairwise preference must
// switch at most once on the domain. Default F must be no better than any valid
// function; ties may return either. Compressed queries must belong to X;
// intervals and setRange use [l,r). All evaluations must fit value_type. !!!!
