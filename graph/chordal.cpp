////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/graph/chordal.cpp
//
// usage: Chordal g; g.setN(n); g.add_edge(u,v); bool ok=g.solve();
//   If ok, peo[0..n-1] is the elimination order, rk[u] is 1-indexed.
//
////////////////////////////////////////////////////////////////
struct Chordal{
int n=0;
vector<vector<int>> edg;
vector<int> peo,rk;
void set(){n=0;edg.clear();peo.clear();rk.clear();}
void setN(int _n){assert(0<=_n&&_n<INT_MAX);set();n=_n;edg.resize(n+1);}
// Undirected, vertices 1..n; loops ignored, parallel edges merged.
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n);
	if(u!=v)edg[u].push_back(v),edg[v].push_back(u);
}
// O(n+m); success fills peo and 1-based rk; failure clears both.
bool solve(){
	vector<int> vis(n+1),cnt(n+1);vector<vector<int>> bucket(n+1),son(n+1);
	for(int u=1;u<=n;u++){
		int k=0;
		for(int v:edg[u])if(vis[v]!=u)vis[v]=u,edg[u][k++]=v;
		edg[u].resize(k);bucket[0].push_back(u);
	}
	peo.resize(n);
	rk.assign(n+1,0);
	int mx=0;
	for(int i=n;i>=1;i--){
		while(true){
			while(bucket[mx].empty())mx--;
			int u=bucket[mx].back();bucket[mx].pop_back();
			if(rk[u]||cnt[u]!=mx)continue;
			peo[i-1]=u;rk[u]=i;
			for(int v:edg[u])if(!rk[v]){
				bucket[++cnt[v]].push_back(v);mx=max(mx,cnt[v]);
			}
			break;
		}
	}
	for(int u=1;u<=n;u++){
		int p=0;
		for(int v:edg[u])if(rk[v]>rk[u]&&(!p||rk[v]<rk[p]))p=v;
		if(p)son[p].push_back(u);
	}
	fill(vis.begin(),vis.end(),0);
	for(int p=1;p<=n;p++){
		for(int v:edg[p])vis[v]=p;
		for(int u:son[p])for(int v:edg[u])if(v!=p&&rk[v]>rk[u]&&vis[v]!=p){
			peo.clear();
			rk.clear();
			return false;
		}
	}
	return true;
}
};
// end for graph/chordal.cpp
/////////////////////////
