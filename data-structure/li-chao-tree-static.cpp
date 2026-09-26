////////////////////////////////////////////////////////////////
//
// template for static / compressed Li Chao tree
// usage: LiChao<F> lc; lc.set(X,default_f); lc.add(f);
//   lc.add(f,l,r); auto [y,g]=lc.ask(x);
//
////////////////////////////////////////////////////////////////
// F needs value_type and const operator()(ll); for example:
// struct F{using value_type=ll;ll k,b;int id;
//     ll operator()(ll x)const{return k*x+b;}
// };
// COMPRESS=true uses discrete coordinates; MINIMIZE=true queries minima.
template<typename F,bool COMPRESS=true,bool MINIMIZE=true>
struct LiChao{
using T=typename F::value_type;
int n=0,sz=1;
ll lo=0,hi=0;
vector<ll> X;
vector<F> a;
// Initialize the query coordinates; queries must belong to X.
// f is the empty-node function, no better than any inserted function.
// For min use F{0,LLONG_MAX,-1}, for max F{0,LLONG_MIN,-1}. O(n log n).
void set(const vector<ll> &x,const F &f){
	static_assert(COMPRESS);X=x;
	sort(X.begin(),X.end());X.erase(unique(X.begin(),X.end()),X.end());
	assert(!X.empty()&&X.size()<=(1<<29));
	n=X.size();
	build(f);
}
// Initialize all integer coordinates in [l,r); requires COMPRESS=false.
// O(r-l) space; f is the empty-node function.
void setRange(ll l,ll r,const F &f){
	static_assert(!COMPRESS);
	assert(l<r&&(__int128)r-l<=(1<<29));
	lo=l;
	hi=r;
	n=r-l;
	build(f);
}
// Clear functions, retaining coordinates and the default function.
void set(){assert(n);fill(a.begin()+1,a.end(),a[0]);}
// Insert on the whole domain; O(log n).
void add(F f){assert(n);add(1,move(f));}
// Insert only on [l,r), left closed and right open; O(log^2 n).
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
// Return {best value, source function}; ties may choose either. O(log n).
// All evaluations must fit F::value_type.
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
static bool better(const T &x,const T &y){
	if constexpr(MINIMIZE)return x<y;
	else return y<x;
}
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
// !!!!! Pairwise preference must switch at most once on the domain. !!!!
