//////////////////////////////////////////////////////////////////
//
// template for min cost max flow (potential + dijkstra + ISAP)
//
// usage:
//   Cost_Flow_Graph g; g.setN(n); g.setST(S,T); g.add_edge(u,v,cap,cost);
//   auto [flow,cost]=g.min_cost_flow(); g.edge_flow(id);
//   Alternative to min-cost-flow-simplex.cpp; use only one of the two
//   templates.
//   Negative costs allowed, input graph must have no negative cost cycle.
//   Editing the graph makes the next solve start over.
//   Flow and intermediate total cost must fit ll; distances use __int128.
//   Initial SPFA: O(n*m); each cost phase:
//   O(n*n*m+(n+m)*log(n+m+2)), O(n+m) space.
//
//////////////////////////////////////////////////////////////////
struct Cost_Flow_Graph{
private:
static constexpr __int128 inf=(__int128)1<<120;
struct edge{int v,nxt;ll w,c;};
vector<edge> e=vector<edge>(2);bool is_flowed=false;
int n=0,m=1,S=0,T=0;ll ans=0,cost=0;
vector<int> lst{0},pre{0},dep,cur,q,cnt,pos,adj;
vector<char> ok;
vector<ll> lim,rem;
vector<__int128> h{0},dis{0};
void add(int u,int v,ll w,ll c){
	assert(m<INT_MAX-1);e.push_back({v,lst[u],w,c});lst[u]=++m;
}
void SPFA(){
	fill(h.begin(),h.end(),0);
	bool neg=false;
	for(int i=2;i<=m;i++)if(e[i].w&&e[i].c<0){neg=true;break;}
	if(!neg)return;
	vector<int> cnt(n+1);vector<char> inq(n+1,true);
	queue<int> q;for(int i=1;i<=n;i++)q.push(i);
	while(!q.empty()){
		int u=q.front();q.pop();inq[u]=false;
		for(int i=pos[u];i<pos[u+1];i++){
			int id=adj[i];if(!e[id].w)continue;
			int v=e[id].v;
			if(h[v]<=h[u]+e[id].c)continue;
			h[v]=h[u]+e[id].c;cnt[v]=cnt[u]+1;
			assert(cnt[v]<n);
			if(!inq[v])inq[v]=true,q.push(v);
		}
	}
}
bool Dijkstra(){
	fill(dis.begin(),dis.end(),inf);dis[S]=0;
	priority_queue<pair<__int128,int>,vector<pair<__int128,int>>,greater<pair<__int128,int>>> q;
	q.push({0,S});
	while(!q.empty()){
		auto [d,u]=q.top();q.pop();
		if(d!=dis[u])continue;
		if(u==T)break;
		for(int i=pos[u];i<pos[u+1];i++){
			int id=adj[i];if(!e[id].w)continue;
			int v=e[id].v;__int128 w=e[id].c+h[u]-h[v];
			assert(w>=0);
			if(dis[v]<=d+w)continue;
			dis[v]=d+w;q.push({dis[v],v});
		}
	}
	if(dis[T]==inf)return false;
	for(int i=1;i<=n;i++)h[i]+=min(dis[i],dis[T]);
	return true;
}
void build(){
	pos.assign((size_t)n+2,0);adj.resize(m-1);
	for(int u=1;u<=n;u++){
		pos[u+1]=pos[u];
		for(int id=lst[u];id;id=e[id].nxt)adj[pos[u+1]++]=id;
	}
}
ll augment(){
	// Reverse BFS and gap only use edges with zero reduced cost.
	for(int id=2;id<=m;id+=2)
		ok[id]=ok[id^1]=(h[e[id].v]==h[e[id^1].v]+e[id].c);
	fill(dep.begin(),dep.end(),n);dep[T]=0;
	int l=0,r=0;q[r++]=T;
	while(l<r){
		int u=q[l++];
		for(int i=pos[u];i<pos[u+1];i++){
			int id=adj[i],v=e[id].v;
			if(dep[v]!=n||!ok[id]||!e[id^1].w)continue;
			dep[v]=dep[u]+1;q[r++]=v;
		}
	}
	if(dep[S]==n)return 0;
	fill(cnt.begin(),cnt.end(),0);
	for(int u=1;u<=n;u++)++cnt[dep[u]];
	cur=pos;int u=S;bool stop=false;
	lim[u]=rem[u]=LLONG_MAX;
	while(true){
		if(u==T)rem[u]=0;
		int &i=cur[u];
		while(!stop&&rem[u]&&i<pos[u+1]){
			int id=adj[i],v=e[id].v;
			if(e[id].w&&ok[id]&&dep[u]==dep[v]+1)break;
			++i;
		}
		if(!stop&&rem[u]&&i==pos[u+1]){
			int d=n;
			// The incoming reverse edge also gains the flow sent below.
			if(u!=S&&lim[u]!=rem[u])d=min(d,dep[e[pre[u]^1].v]+1);
			for(int j=pos[u];j<pos[u+1];j++){
				int id=adj[j];
				if(e[id].w&&ok[id])d=min(d,dep[e[id].v]+1);
			}
			stop=(--cnt[dep[u]]==0);
			dep[u]=d;++cnt[d];cur[u]=pos[u];
			if(u==S){if(d==n)stop=true;if(!stop)continue;}
		}else if(!stop&&rem[u]){
			int id=adj[i],v=e[id].v;pre[v]=id;
			lim[v]=rem[v]=min(rem[u],e[id].w);u=v;continue;
		}
		// Unwind partial pushes as well when gap terminates the search.
		ll f=lim[u]-rem[u];
		if(u==S)return f;
		int id=pre[u];u=e[id^1].v;
		e[id].w-=f;e[id^1].w+=f;rem[u]-=f;
		if(rem[u])++cur[u];
	}
}
public:
void set(){
	n=S=T=0;m=1;ans=cost=0;is_flowed=false;e.assign(2,{0,0,0,0});
	dep.clear();cur.clear();q.clear();cnt.clear();pos.clear();adj.clear();ok.clear();
	lim.clear();rem.clear();
	lst.assign(1,0);pre.assign(1,0);h.assign(1,0);dis.assign(1,0);
}
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);set();n=_n;
	lst.assign(n+1,0);pre.assign(n+1,0);h.assign(n+1,0);dis.assign(n+1,0);
}
void setST(int _S,int _T){
	assert(1<=min(_S,_T)&&max(_S,_T)<INT_MAX&&_S!=_T);
	S=_S;T=_T;n=max({n,S,T});
	lst.resize(n+1);pre.resize(n+1);h.resize(n+1);dis.resize(n+1);is_flowed=false;
}
int new_node(){
	assert(n<INT_MAX-1);++n;lst.push_back(0);pre.push_back(0);h.push_back(0);dis.push_back(0);
	is_flowed=false;return n;
}
ll edge_flow(int id){
	assert(2<=id&&id<=m&&id%2==0);return e[id^1].w;
}
int add_edge(int u,int v,ll w,ll c){
	assert(1<=min(u,v)&&max(u,v)<=n&&w>=0&&c!=LLONG_MIN);
	assert(m&1);is_flowed=false;add(u,v,w,c);add(v,u,0,-c);return m-1;
}
pair<ll,ll> min_cost_flow(){
	assert(1<=min(S,T)&&max(S,T)<=n&&S!=T);
	if(is_flowed)return {ans,cost};
	for(int i=2;i<=m;i+=2)e[i].w+=e[i^1].w,e[i^1].w=0;
	ans=cost=0;build();SPFA();
	dep.resize(n+1);cur.resize(n+1);q.resize(n+1);
	lim.resize(n+1);rem.resize(n+1);cnt.resize(n+1);ok.resize(m+1);
	while(Dijkstra()){
		ll f=augment();assert(f>0);
		__int128 c=h[T]-h[S];
		assert(ans<=LLONG_MAX-f);
		assert(c>=((__int128)LLONG_MIN-cost)/f&&c<=((__int128)LLONG_MAX-cost)/f);
		ans+=f;cost=(ll)(cost+c*f);
	}
	is_flowed=true;return {ans,cost};
}
};
// end for graph/flow/min-cost-flow-isap.cpp
/////////////////////////
// !!!!! Only paste one Cost_Flow_Graph variant; input must have no negative
// cost cycle; flow and total cost must fit ll. !!!!
