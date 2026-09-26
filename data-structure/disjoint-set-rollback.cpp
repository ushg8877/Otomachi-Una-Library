////////////////////////////////////////////////////////////////
//
// template for unionfind (accept undo, O(nlogn))
// saver.version 1.0 (Last Update Jun 28th, 2026)
//
// usage: unionfind uf; uf.setN(n); int t=version();
//   uf.merge(u,v); roll_back(t);
//
////////////////////////////////////////////////////////////////
struct unionfind{
int n=0;
vector<int> fa,siz;
inline void setN(int _n){
	// saved addresses must not survive setN()
	assert(version()==0);
	assert(0<=_n&&_n<INT_MAX);n=_n;
	fa.resize(n+1);
	iota(fa.begin(),fa.end(),0);
	siz.assign(n+1,1);
	siz[0]=0;
}
inline int find(int x){assert(1<=x&&x<=n);while(x^fa[x])x=fa[x];return x;}
inline int merge(int x,int y){
	assert(1<=min(x,y)&&max(x,y)<=n);
	x=find(x),y=find(y);
	if(x==y) return 0;
	if(siz[x]<siz[y]) swap(x,y);
	save(siz[x]);save(fa[y]);
	siz[x]+=siz[y];fa[y]=x;
	return x;
}
};
// end for data-structure/disjoint-set-rollback.cpp
/////////////////////////
// !!!!! Paste data-structure/rollback-int.cpp first. !!!!
