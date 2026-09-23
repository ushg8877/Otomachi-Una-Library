////////////////////////////////////////////////////////////////
//
// template for XOR equations (bitset Gaussian elimination)
// usage:
//   Xor_Gauss<N> G; G.setN(n); // n<=N, default n=N
//   G.add(a,b); // xor of x[i] with a[i]=1 equals b
//   int res=G.solve(); // -1: no solution, 0: multiple, 1: unique
//   answer in G.x[0..n-1], free variables are set to zero
//   G.rank; G.set(); // clear equations, keep n
//   add returns whether the whole system is still consistent
//   variables and right-hand sides are ll, not individual bits
//   GNU bitset; W=ceil(N/64), O((n+1)*(W+1)) per add, O(n*(n+W+1)) to solve
//
////////////////////////////////////////////////////////////////
template<size_t N> struct Xor_Gauss{
static_assert(N<INT_MAX);
int n=N,rank=0;bool ok=true;
vector<bitset<N>> a=vector<bitset<N>>(N);
vector<ll> b=vector<ll>(N),x;
void set(){
	rank=0;ok=true;x.clear();
	for(auto &v:a)v.reset();
	fill(b.begin(),b.end(),0);
}
void setN(int _n){
	assert(0<=_n&&(size_t)_n<=N);n=_n;rank=0;ok=true;
	a.assign(n,{});b.assign(n,0);x.clear();
}
bool add(bitset<N> v,ll w){
	assert((v>>n).none());x.clear();
	if(!ok)return false;
	for(size_t i=v._Find_first();i<(size_t)n;i=v._Find_next(i)){
		if(a[i][i])v^=a[i],w^=b[i];
		else{a[i]=v;b[i]=w;++rank;return true;}
	}
	return ok=(w==0);
}
int solve(){
	x.clear();if(!ok)return -1;
	x.assign(n,0);
	for(int i=n-1;i>=0;i--)if(a[i][i]){
		x[i]=b[i];
		for(size_t j=a[i]._Find_next(i);j<(size_t)n;j=a[i]._Find_next(j))x[i]^=x[j];
	}
	return rank==n?1:0;
}
};
// end for math/xor-equations.cpp
/////////////////////////
// !!!!! Uses GNU libstdc++ bitset extensions; variables are ll and coordinates are 0..n-1. !!!!
