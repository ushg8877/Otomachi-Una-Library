//////////////////////////////////////////////////////////////////
//
// template for minimum Steiner tree
//
// usage: Steiner_Tree g; g.setN(n); g.add_edge(u,v,w);
//   ll ans=g.solve({s1,s2,s3}); // minimum cost, -1 if disconnected
//
//////////////////////////////////////////////////////////////////
struct Steiner_Tree{
private:
static constexpr ll inf=LLONG_MAX;
int n=0;
vector<vector<pair<int,ll>>> edg;
public:
void set(){n=0;edg.clear();}
void setN(int _n){assert(0<=_n&&_n<INT_MAX);n=_n;edg.assign(n+1,{});}
// Undirected, vertices 1..n; weights must be nonnegative and < LLONG_MAX.
void add_edge(int u,int v,ll w){
	assert(1<=min(u,v)&&max(u,v)<=n&&0<=w&&w<inf);
	edg[u].push_back({v,w});edg[v].push_back({u,w});
}
// Minimum cost, -1 if disconnected; repeated terminals ignored.
// Empty / single-terminal set returns 0; finite answers < LLONG_MAX.
// For k terminals: O(n*3^k+m*log(m+2)*2^k) time, O(n*2^k+n+m) space.
ll solve(vector<int> s)const{
	for(int u:s)assert(1<=u&&u<=n);
	sort(s.begin(),s.end());s.erase(unique(s.begin(),s.end()),s.end());
	if(s.size()<=1)return 0;
	// Fix one terminal as the root, saving half of the subset states.
	int rt=s.back();
	s.pop_back();
	assert(s.size()<31);
	int k=s.size(),N=1<<k;size_t len=(size_t)n+1;
	vector<ll> f;assert((size_t)N<=f.max_size()/len);
	f.assign((size_t)N*len,inf);
	for(int i=0;i<k;i++)f[(size_t)(1<<i)*len+s[i]]=0;
	vector<pair<ll,int>> q;greater<pair<ll,int>> cmp;
	for(int mask=1;mask<N;mask++){
		ll *d=f.data()+(size_t)mask*len;
		for(int sub=(mask-1)&mask;sub;sub=(sub-1)&mask){
			int other=mask^sub;if(sub>other)continue;
			const ll *a=f.data()+(size_t)sub*len,*b=f.data()+(size_t)other*len;
			for(int u=1;u<=n;u++)if(a[u]<d[u]&&b[u]<d[u]-a[u])d[u]=a[u]+b[u];
		}
		q.clear();
		for(int u=1;u<=n;u++)if(d[u]!=inf)q.push_back({d[u],u});
		make_heap(q.begin(),q.end(),cmp);
		while(!q.empty()){
			auto [dd,u]=q.front();
			pop_heap(q.begin(),q.end(),cmp);
			q.pop_back();
			if(dd!=d[u])continue;
			if(mask==N-1&&u==rt)return dd;
			for(auto [v,w]:edg[u])if(w<d[v]-dd){
				d[v]=dd+w;
				q.push_back({d[v],v});
				push_heap(q.begin(),q.end(),cmp);
			}
		}
	}
	return -1;
}
};
// end for graph/steiner-tree.cpp
/////////////////////////
