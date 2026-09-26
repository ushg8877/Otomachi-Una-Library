////////////////////////////////////////////////////////////////
//
// template for Directed MST (Chu-Liu / Edmonds)
// usage:
//   Directed_MST g; g.setN(n); g.add_edge(u,v,w); auto a=g.solve(S);
//   a.ok; a.cost; a.id[v] is the incoming edge of v, a.id[S]=0.
//   Vertices / edge ids are 1-indexed; edges point away from the root S.
//   Negative weights / parallel edges / loops are allowed. No solution:
//   ok=false.
//   O((n+m)*log(n+1)) time, O(n+m) space; cost uses __int128.
//
////////////////////////////////////////////////////////////////
struct Directed_MST{
struct edge{int u,v;ll w;};
struct result{bool ok=false;__int128 cost=0;vector<int> id;};
int n=0;
vector<edge> e=vector<edge>(1);
void set(){n=0;e.assign(1,{});}
void setN(int _n){assert(0<=_n&&_n<INT_MAX);set();n=_n;}
int add_edge(int u,int v,ll w){
	assert(1<=min(u,v)&&max(u,v)<=n&&e.size()<INT_MAX);
	int id=e.size();e.push_back({u,v,w});return id;
}
private:
struct heap{
	struct node{int ls=0,rs=0,dis=0,id=0;__int128 w=0,lazy=0;};
	vector<node> a=vector<node>(1);
	void add(int u,__int128 w){if(u)a[u].w+=w,a[u].lazy+=w;}
	void down(int u){add(a[u].ls,a[u].lazy);add(a[u].rs,a[u].lazy);a[u].lazy=0;}
	int merge(int u,int v){
		if(!u||!v)return u|v;
		if(a[u].w>a[v].w||(a[u].w==a[v].w&&a[u].id>a[v].id))swap(u,v);
		down(u);a[u].rs=merge(a[u].rs,v);
		if(a[a[u].ls].dis<a[a[u].rs].dis)swap(a[u].ls,a[u].rs);
		a[u].dis=a[a[u].rs].dis+1;return u;
	}
	int push(int u,int id,ll w){
		a.push_back({0,0,1,id,w,0});return merge(u,(int)a.size()-1);
	}
	int pop(int u){down(u);return merge(a[u].ls,a[u].rs);}
};
struct DSU{
	vector<int> fa,siz;
	vector<pair<int,int>> stk;
	DSU(int n):fa(n+1),siz(n+1,1){iota(fa.begin(),fa.end(),0);}
	int find(int u)const{while(fa[u]!=u)u=fa[u];return u;}
	bool join(int u,int v){
		u=find(u);v=find(v);if(u==v)return false;
		if(siz[u]<siz[v])swap(u,v);
		stk.emplace_back(u,v);fa[v]=u;siz[u]+=siz[v];return true;
	}
	void rollback(int t){
		while((int)stk.size()>t){auto [u,v]=stk.back();stk.pop_back();fa[v]=v;siz[u]-=siz[v];}
	}
};
public:
result solve(int S)const{
	assert(1<=S&&S<=n);DSU dsu(n);heap h;
	vector<int> rt(n+1),seen(n+1),path(n),q(n),in(n+1);
	// Keep one cheapest edge per ordered pair, also bounding heap size by n*n.
	{
		vector<vector<int>> edg(n+1);vector<int> vis(n+1),best(n+1),a;
		for(int id=1;id<(int)e.size();id++)if(e[id].u!=e[id].v&&e[id].v!=S)edg[e[id].v].push_back(id);
		for(int v=1;v<=n;v++){
			a.clear();
			for(int id:edg[v]){
				int u=e[id].u;
				if(vis[u]!=v)vis[u]=v,best[u]=id,a.push_back(u);
				else if(e[id].w<e[best[u]].w)best[u]=id;
			}
			for(int u:a)rt[v]=h.push(rt[v],best[u],e[best[u]].w);
		}
	}
	struct cycle{int u,t;vector<int> e;};
	vector<cycle> cyc;__int128 cost=0;seen[S]=-1;
	for(int s=1;s<=n;s++){
		int u=dsu.find(s),top=0;
		while(!seen[u]){
			while(rt[u]&&dsu.find(e[h.a[rt[u]].id].u)==u)rt[u]=h.pop(rt[u]);
			if(!rt[u])return {};
			int id=h.a[rt[u]].id;__int128 w=h.a[rt[u]].w;
			cost+=w;h.add(rt[u],-w);rt[u]=h.pop(rt[u]);
			q[top]=id;path[top++]=u;seen[u]=s;u=dsu.find(e[id].u);
			if(seen[u]==s){
				int end=top,t=dsu.stk.size(),r=0,v;
				do{v=path[--top];r=h.merge(r,rt[v]);}while(dsu.join(u,v));
				u=dsu.find(u);rt[u]=r;seen[u]=0;
				cyc.push_back({u,t,vector<int>(q.begin()+top,q.begin()+end)});
			}
		}
		for(int i=0;i<top;i++)in[dsu.find(e[q[i]].v)]=q[i];
	}
	for(int i=(int)cyc.size()-1;i>=0;i--){
		auto &c=cyc[i];int id=in[c.u];assert(id);
		dsu.rollback(c.t);
		for(int x:c.e)in[dsu.find(e[x].v)]=x;
		in[dsu.find(e[id].v)]=id;
	}
	return {true,cost,move(in)};
}
};
// end for graph/directed-mst.cpp
/////////////////////////
// !!!!! The result cost is __int128; basic/fast-io.cpp provides decimal IO if
// needed. !!!!
