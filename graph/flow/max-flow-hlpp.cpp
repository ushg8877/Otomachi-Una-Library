//////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/graph/flow/max-flow-hlpp.cpp
//
// usage: Flow_Graph g; g.setN(n); g.setST(S,T); g.add_edge(u,v,cap);
//   ll f = g.max_flow(); g.edge_flow(id);
//
//////////////////////////////////////////////////////////////////
struct Flow_Graph{
private:

struct edge{int v,nxt;ll w;};
vector<edge> e=vector<edge>(2);bool is_flowed=false;
int n=0,m=1,S=0,T=0;ll ans=0;
vector<int> lst{0},level{0},cur{0},pos,adj,cnt,q;
vector<vector<int>> bucket;
vector<__int128> ex;
int hi=0;
ll work=0;
void add(int u,int v,ll w){
	assert(m<INT_MAX-1);
	e.push_back({v,lst[u],w});
	lst[u]=++m;
}
void active(int u){
	if(u==S||u==T||!ex[u])return;
	assert(level[u]<2*n);
	bucket[level[u]].push_back(u);hi=max(hi,level[u]);
}
void rebuild(){
	for(auto &v:bucket)v.clear();
	fill(cnt.begin(),cnt.end(),0);hi=0;
	for(int u=1;u<=n;u++){
		++cnt[level[u]];
		cur[u]=pos[u];
		active(u);
	}
}
void global(){
	fill(level.begin(),level.end(),2*n);
	level[T]=0;
	level[S]=n;
	int l=0,r=0;q[r++]=T;
	while(l<r){
		int u=q[l++];
		for(int i=pos[u];i<pos[u+1];i++){
			int id=adj[i],v=e[id].v;
			if(!e[id^1].w||level[v]!=2*n)continue;
			level[v]=level[u]+1;q[r++]=v;
		}
	}
	// Label the remaining vertices towards S, returning excess to the source.
	l=r=0;q[r++]=S;
	while(l<r){
		int u=q[l++];
		for(int i=pos[u];i<pos[u+1];i++){
			int id=adj[i],v=e[id].v;
			if(!e[id^1].w||level[v]!=2*n)continue;
			level[v]=level[u]+1;q[r++]=v;
		}
	}
	rebuild();work=0;
}
void push(int u,int id){
	int v=e[id].v;ll f=(ll)min(ex[u],(__int128)e[id].w);
	bool zero=!ex[v];
	e[id].w-=f;
	e[id^1].w+=f;
	ex[u]-=f;
	ex[v]+=f;
	if(zero&&f)active(v);
}
void discharge(int u){
	while(ex[u]){
		int &i=cur[u];
		while(i<pos[u+1]){
			int id=adj[i];++work;
			if(e[id].w&&level[u]==level[e[id].v]+1){
				push(u,id);if(!ex[u])return;
			}
			++i;
		}
		int old=level[u],d=2*n;
		for(int j=pos[u];j<pos[u+1];j++){
			int id=adj[j];++work;
			if(e[id].w)d=min(d,level[e[id].v]+1);
		}
		assert(d<2*n);
		level[u]=d;
		cur[u]=pos[u];
		++cnt[d];
		if(--cnt[old]==0&&old<n){
			for(int v=1;v<=n;v++)if(old<level[v]&&level[v]<n)level[v]=n+1;
			rebuild();
			return;
		}
		active(u);
		return;
	}
}
public:
void set(){
	n=S=T=0;
	m=1;
	ans=0;
	is_flowed=false;
	e.assign(2,{0,0,0});
	lst.assign(1,0);
	level.assign(1,0);
	cur.assign(1,0);
	pos.clear();
	adj.clear();
	cnt.clear();
	q.clear();
	bucket.clear();
	ex.clear();
}
// Clear the graph; use new_node() to keep existing edges.
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX/2);
	set();
	n=_n;
	lst.assign(n+1,0);
	level.assign(n+1,0);
	cur.assign(n+1,0);
}
void setST(int _S,int _T){
	assert(1<=min(_S,_T)&&max(_S,_T)<INT_MAX/2&&_S!=_T);
	if(S!=_S||T!=_T){
		for(int i=2;i<=m;i+=2)e[i].w+=e[i^1].w,e[i^1].w=0;
		ans=0;
	}
	S=_S;
	T=_T;
	n=max({n,S,T});
	lst.resize(n+1);
	level.resize(n+1);
	cur.resize(n+1);
	is_flowed=false;
}
int new_node(){
	assert(n<INT_MAX/2-1);
	++n;
	lst.push_back(0);
	level.push_back(0);
	cur.push_back(0);
	is_flowed=false;
	return n;
}
ll edge_flow(int id){
	assert(2<=id&&id<=m&&id%2==0);
	return e[id^1].w;
}
int add_edge(int u,int v,ll w){
	// return forward edge id
	assert(1<=min(u,v)&&max(u,v)<=n&&0<=w);
	assert(m&1);
	is_flowed=false;
	add(u,v,w);
	add(v,u,0);
	return m-1;
}
// Return total flow; adding edges continues on the residual graph.
// Excess is returned to S, so edge_flow(id) gives a feasible flow.
// Capacities / total flow must fit ll; excess uses __int128.
ll max_flow(){
	assert(1<=min(S,T)&&max(S,T)<=n&&S!=T);
	if(is_flowed) return ans;
	pos.assign(n+2,0);adj.resize(m-1);
	for(int u=1;u<=n;u++){
		pos[u+1]=pos[u];
		for(int id=lst[u];id;id=e[id].nxt)adj[pos[u+1]++]=id;
	}
	ex.assign(n+1,0);
	cnt.resize(2*n+1);
	q.resize(n+1);
	bucket.resize(2*n+1);
	for(int i=pos[S];i<pos[S+1];i++){
		int id=adj[i],v=e[id].v;if(v==S)continue;
		ll f=e[id].w;
		e[id].w=0;
		e[id^1].w+=f;
		ex[S]-=f;
		ex[v]+=f;
	}
	global();
	while(true){
		while(hi>=0&&bucket[hi].empty())--hi;
		if(hi<0)break;
		int u=bucket[hi].back();bucket[hi].pop_back();
		discharge(u);
		if(work>4LL*((ll)n+m))global();
	}
	assert(ex[T]<=LLONG_MAX-ans);
	ans+=(ll)ex[T];
	is_flowed=true;
	return ans;
}
};
// end for graph/flow/max-flow-hlpp.cpp
/////////////////////////
