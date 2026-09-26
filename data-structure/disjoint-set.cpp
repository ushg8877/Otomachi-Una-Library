////////////////////////////////////////////////////////////////
//
// template for Disjoint Set Union (DSU)
//
// usage: DSU dsu; dsu.setN(n); dsu.merge(u, v);
//   dsu.same(u, v); int r = dsu.find(u);  // merge returns new root or 0
//
////////////////////////////////////////////////////////////////
struct DSU{
int n=0;
vector<int> fa,siz;
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);n=_n;
	fa.resize(n+1);
	iota(fa.begin(),fa.end(),0);
	siz.assign(n+1,1);
	siz[0]=0;
}
inline int find(int x){
	assert(1<=x&&x<=n);
	while(x^fa[x])x=fa[x]=fa[fa[x]];
	return x;
}
inline int merge(int x,int y){
	assert(1<=min(x,y)&&max(x,y)<=n);
	x=find(x),y=find(y);
	if(x==y) return 0;
	if(siz[x]<siz[y]) swap(x,y);
	siz[x]+=siz[y],fa[y]=x;
	return x;
}
inline bool same(int x,int y){
	assert(1<=min(x,y)&&max(x,y)<=n);
	return find(x)==find(y);
}
};
// end for data-structure/disjoint-set.cpp
/////////////////////////
