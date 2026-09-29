struct Flow_Graph{
private:
static const ll inf=2.02e18;
struct edge{int v,nxt;ll w;};
vector<edge> e=vector<edge>(2);bool is_flowed=false;
int n=0,m=1,S=0,T=0;ll ans=0;
vector<int> lst{0},level{0},cur{0};
void add(int u,int v,ll w){
	assert(m<INT_MAX-1);e.push_back({v,lst[u],w});lst[u]=++m;
}
bool BFS(){
	for(int i=0;i<=n;i++) level[i]=-1;
	queue<int> q;q.push(S);level[S]=0;
	while(!q.empty()){
		int u=q.front();q.pop();
		if(u==T) return true;
		for(int id=lst[u];id;id=e[id].nxt) if(e[id].w&&level[e[id].v]==-1){
			level[e[id].v]=level[u]+1;
			q.push(e[id].v);
		}
	}
	return false;
}
ll DFS(int u,ll flow){
	if(u==S) return flow;
	ll now=0;
	for(int &i=cur[u];i;i=e[i].nxt) if(level[e[i].v]==level[u]-1&&e[i^1].w){
		ll fo=DFS(e[i].v,min(e[i^1].w,flow-now));
		e[i^1].w-=fo;e[i].w+=fo;now+=fo;
		if(now==flow) break;
	}
	return now;
}
public:
void set(){
	n=S=T=0;m=1;ans=0;is_flowed=false;e.assign(2,{0,0,0});
	lst.assign(1,0);level.assign(1,0);cur.assign(1,0);
}
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);set();n=_n;
	lst.assign(n+1,0);level.assign(n+1,0);cur.assign(n+1,0);
}
void setST(int _S,int _T){
	assert(1<=min(_S,_T)&&max(_S,_T)<INT_MAX&&_S!=_T);
	if(S!=_S||T!=_T){
		for(int i=2;i<=m;i+=2)e[i].w+=e[i^1].w,e[i^1].w=0;
		ans=0;
	}
	S=_S;T=_T;n=max({n,S,T});
	lst.resize(n+1);level.resize(n+1);cur.resize(n+1);is_flowed=false;
}
int new_node(){
	assert(n<INT_MAX-1);++n;lst.push_back(0);level.push_back(0);
		cur.push_back(0);
	is_flowed=false;return n;
}
ll edge_flow(int id){
	assert(2<=id&&id<=m&&id%2==0);
	return e[id^1].w;
}
int add_edge(int u,int v,ll w){
	// return forward edge id
	assert(1<=min(u,v)&&max(u,v)<=n&&0<=w);
	assert(m&1);
	is_flowed=false;add(u,v,w);add(v,u,0);return m-1;
}
ll max_flow(){
	assert(1<=min(S,T)&&max(S,T)<=n&&S!=T);
	if(is_flowed) return ans;
	is_flowed=true;
	while(BFS()){
		for(int i=1;i<=n;i++) cur[i]=lst[i];
		ll f=DFS(T,inf);assert(ans<=LLONG_MAX-f);ans+=f;
	}
	return ans;
}
};
