////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/graph/shortest-path-label.cpp
//
// usage: Shortest_Path<string> g; g.add_edge("a","b",2); g.setS("a");
//   g.print("b"); // setS(vector<string>{"a","c"}) for multiple sources
//
////////////////////////////////////////////////////////////////
template<typename T=string>
struct Shortest_Path{
private:
map<T,int> id;
vector<T> name;
vector<vector<pair<int,ll>>> edg;
vector<int> s,pre;
vector<ll> d;
bool ready=false;
int node(const T &x){
	auto it=id.find(x);
	if(it!=id.end()) return it->second;
	assert(name.size()<INT_MAX);
	int u=name.size();
	id.emplace(x,u);name.push_back(x);edg.emplace_back();
	ready=false;
	return u;
}
void build(){
	if(ready) return;
	d.assign(name.size(),LLONG_MAX);pre.assign(name.size(),-1);
	priority_queue<pair<ll,int>,vector<pair<ll,int>>,greater<pair<ll,int>>> q;
	for(int u:s) if(d[u]){d[u]=0;q.push({0,u});}
	while(!q.empty()){
		auto [dd,u]=q.top();q.pop();
		if(dd!=d[u]) continue;
		for(auto [v,w]:edg[u]) if(w<d[v]-dd){
			d[v]=dd+w;pre[v]=u;
			q.push({d[v],v});
		}
	}
	ready=true;
}
public:
// Clear the graph and sources.
void set(){
	id.clear();name.clear();edg.clear();s.clear();pre.clear();d.clear();
	ready=false;
}
// Add a directed edge; labels are created automatically, default cost 1.
// For an undirected edge, also call add_edge(v,u,w).
void add_edge(const T &u,const T &v,ll w=1){
	assert(0<=w&&w<LLONG_MAX);
	int a=node(u),b=node(v);
	edg[a].push_back({b,w});ready=false;
}
// Start a new round with one source; the graph is retained.
void setS(const T &x){s.clear();s.push_back(node(x));ready=false;}
// Start a new round with multiple sources, all at distance 0.
void setS(const vector<T> &a){
	s.clear();
	for(const T &x:a) s.push_back(node(x));
	ready=false;
}
// Return {distance, full path}; unreachable gives {-1,{}}.
// First query after changes: O((n+m)log(n+m+1)); later: O(log n+path).
pair<ll,vector<T>> solve(const T &t){
	auto it=id.find(t);
	if(it==id.end()) return {-1,{}};
	build();
	int u=it->second;
	if(d[u]==LLONG_MAX) return {-1,{}};
	vector<T> path;
	for(int v=u;v!=-1;v=pre[v]) path.push_back(name[v]);
	reverse(path.begin(),path.end());
	return {d[u],move(path)};
}
// Print distance, then the labels on one line; unreachable prints only -1.
void print(const T &t,ostream &out=cout){
	auto [dis,path]=solve(t);
	out<<dis<<'\n';
	if(dis<0) return;
	for(int i=0;i<(int)path.size();i++){
		if(i) out<<" -> ";
		out<<path[i];
	}
	out<<'\n';
}
};
// end for graph/shortest-path-label.cpp
/////////////////////////
// !!!!! Nonnegative weights; finite distances must be < LLONG_MAX. !!!!
