////////////////////////////////////////////////////////////////
//
// template for XOR equations (bitset Gaussian elimination)
// usage: Xor_Gauss<N> g; g.setN(n); g.add(a,b);
//   int res=g.solve(); auto x=g.x;
//
////////////////////////////////////////////////////////////////
template<size_t N>
struct Xor_Gauss{
static_assert(N<INT_MAX);
int n=N,rank=0;bool ok=true;
vector<bitset<N>> a=vector<bitset<N>>(N);
vector<ll> b=vector<ll>(N),x;
// Clear equations and keep n; setN(n) also changes the variable count.
void set(){
	rank=0;
	ok=true;
	x.clear();
	for(auto &v:a)v.reset();
	fill(b.begin(),b.end(),0);
}
void setN(int _n){
	assert(0<=_n&&(size_t)_n<=N);
	n=_n;
	rank=0;
	ok=true;
	a.assign(n,{});
	b.assign(n,0);
	x.clear();
}
// XOR of variables selected by v equals w; variables / RHS are ll.
// Return whether the whole system is consistent. W=ceil(N/64);
// O((n+1)*(W+1)) per insertion, including elimination.
bool add(bitset<N> v,ll w){
	assert((v>>n).none());x.clear();
	if(!ok)return false;
	for(size_t i=v._Find_first();i<(size_t)n;i=v._Find_next(i)){
		if(a[i][i])v^=a[i],w^=b[i];
		else{a[i]=v;b[i]=w;++rank;return true;}
	}
	return ok=(w==0);
}
// Return -1 / 0 / 1 for none / multiple / unique; x[0..n-1] is a solution.
// Free variables are zero; O(n*(n+W+1)), W=ceil(N/64).
int solve(){
	x.clear();if(!ok)return -1;
	x.assign(n,0);
	for(int i=n-1;i>=0;i--)if(a[i][i]){
		x[i]=b[i];
		for(size_t j=a[i]._Find_next(i);j<(size_t)n;j=a[i]._Find_next(j))
			x[i]^=x[j];
	}
	return rank==n?1:0;
}
};
// end for math/xor-equations.cpp
/////////////////////////
// !!!!! Requires GNU libstdc++ bitset extensions. !!!!
