////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/math/sprague-grundy.cpp
//
// usage: SG_Game g; g.add_edge("a","b"); g.add("c",{"a","b"});
//   g.print(); int sg=g.get("a"); // terminal nodes have SG 0
//
////////////////////////////////////////////////////////////////
struct SG_Game{
private:
map<string,int> id;
vector<string> name;
vector<vector<int>> edg;
vector<int> sg;
bool ready=false;
int node(const string &s){
	auto it=id.find(s);
	if(it!=id.end()) return it->second;
	assert(name.size()<INT_MAX);
	int u=name.size();
	id[s]=u;name.push_back(s);edg.emplace_back();ready=false;
	return u;
}
public:
// Clear the graph and computed values.
void set(){id.clear();name.clear();edg.clear();sg.clear();ready=false;}
// Add a directed move u -> v; labels are created automatically.
void add_edge(const string &u,const string &v){
	int x=node(u),y=node(v);
	edg[x].push_back(y);ready=false;
}
// Set/replace all moves from s. add(s,{}) creates a terminal node.
// Referenced labels without their own moves are also terminal nodes.
void add(const string &s,const vector<string> &to){
	int u=node(s);
	vector<int> e;
	for(auto &t:to) e.push_back(node(t));
	edg[u]=move(e);ready=false;
}
// Compute all SG values in O(n+m) time and space; repeated calls are cached.
void solve(){
	if(ready) return;
	int n=name.size();
	vector<size_t> deg(n);
	vector<int> q,vis(n+1,-1);
	vector<vector<int>> rev(n);
	for(int u=0;u<n;u++){
		deg[u]=edg[u].size();
		for(int v:edg[u]) rev[v].push_back(u);
		if(!deg[u]) q.push_back(u);
	}
	sg.assign(n,0);
	for(int i=0;i<(int)q.size();i++){
		int u=q[i];
		for(int v:edg[u]) vis[sg[v]]=u;
		while(vis[sg[u]]==u) sg[u]++;
		for(int v:rev[u]) if(!--deg[v]) q.push_back(v);
	}
	assert((int)q.size()==n);ready=true;
}
// Query an existing label; automatically recompute after graph changes.
int get(const string &s){
	auto it=id.find(s);assert(it!=id.end());
	solve();return sg[it->second];
}
// Print each label and its SG value, in the order labels were created.
void print(ostream &o=cout){
	solve();
	for(int i=0;i<(int)name.size();i++) o<<name[i]<<" = "<<sg[i]<<'\n';
}
};
// end for math/sprague-grundy.cpp
/////////////////////////
// !!!!! DAG only; both players share the moves, and no move means losing. !!!!
