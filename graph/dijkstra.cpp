//////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/graph/dijkstra.cpp
//
// usage: graph_SSSP G; G.setN(n); G.setS(s); G.add_edge(u,v,w);
//   auto d = G.SSSP();  // unreachable 1e18, non-negative weights
//
//////////////////////////////////////////////////////////////////
struct graph_SSSP{
private:
int n=0,s=0;
vector<vector<pair<int,ll>>> edg;
public:
void setN(int N){assert(0<=N&&N<INT_MAX);n=N;s=0;edg.assign(n+1,{});}
void setS(int S){s=S;assert(1<=s&&s<=n);}
int new_node(){edg.push_back({});return ++n;}
void add_edge(int u,int v,ll w){
	assert(1<=min(u,v)&&max(u,v)<=n&&0<=w);
	edg[u].push_back(make_pair(v,w));
}
// Unreachable distances are 1e18; finite distances must stay below it.
vector<ll> SSSP(){
	// return a vector of length n+1 --- ans[i] shortest path between s and i
	// if no path, ans[i] will be 1e18
	assert(1<=s&&s<=n);
	priority_queue<pair<ll,int>,vector<pair<ll,int>>,greater<pair<ll,int>>> q;
	vector<ll> d(n+1,1e18);
	d[s]=0;q.push(make_pair(0,s));
	while(!q.empty()){
		auto [dd,u]=q.top();q.pop();
		if(dd!=d[u]) continue;
		for(auto [v,w]:edg[u]) if(w<d[v]-d[u]){
			d[v]=d[u]+w;
			q.push(make_pair(d[v],v));
		}
	}
	return d;
}
};
// end for graph/dijkstra.cpp
/////////////////////////
// !!!!! Edge weights must be nonnegative. !!!!
