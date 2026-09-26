////////////////////////////////////////////////////////////////
//
// template for Regular Bipartite Matching
// usage:
//   Regular_Bipartite g; g.setN(n); g.add_edge(u,v);
//   auto a=g.solve(); // d perfect matchings, a[c][u] is the right partner
//   Each row has n+1 entries, a[c][0]=0; c is 0-indexed, u is 1-indexed.
//   Both sides are 1..n; every vertex must have the same degree d.
//   Parallel edges are allowed and counted with multiplicity.
//   Deterministic O(n+nd*log(nd+1)) time, O(n+nd) space.
//
////////////////////////////////////////////////////////////////
struct Regular_Bipartite{
struct edge{int u,v;};
int n=0;
vector<edge> e=vector<edge>(1);
void set(){n=0;e.assign(1,{});}
void setN(int _n){assert(0<=_n&&_n<INT_MAX/4);set();n=_n;}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n&&e.size()+n<INT_MAX/2);
	e.push_back({u,v});
}
private:
// Split an even-degree bipartite graph into two equal-degree subgraphs.
// Negative id -u denotes a temporary edge (u,u).
vector<char> split(const vector<int> &a)const{
	assert(a.size()<INT_MAX/2);int m=a.size();
	vector<int> head(2*n+1,-1),nxt(2*m);
	vector<char> vis(m),col(m);
	for(int i=0;i<m;i++){
		int id=a[i],u=id>0?e[id].u:-id,v=(id>0?e[id].v:-id)+n;
		nxt[i*2]=head[u];head[u]=i*2;
		nxt[i*2+1]=head[v];head[v]=i*2+1;
	}
	for(int s=1;s<=2*n;s++){
		int u=s,c=0;
		while(head[u]!=-1){
			int x=head[u];
			head[u]=nxt[x];
			int i=x/2;
			if(vis[i])continue;
			vis[i]=true;
			col[i]=c;
			c^=1;
			int id=a[i];u=x&1?(id>0?e[id].u:-id):(id>0?e[id].v:-id)+n;
		}
		assert(u==s&&!c);
	}
	return col;
}
vector<int> perfect(const vector<int> &ids,int d)const{
	ll q=1;
	while(q<(ll)ids.size())q*=2;
	vector<pair<int,ll>> a; a.reserve(ids.size()+n);
	for(int id:ids)a.emplace_back(id,q/d);
	if(q%d)for(int u=1;u<=n;u++)a.emplace_back(-u,q%d);
	while(q>1){
		vector<int> odd;
		for(auto [id,w]:a)if(w&1)odd.push_back(id);
		auto col=split(odd);int c0=0,c1=0;
		for(int i=0;i<(int)odd.size();i++)if(odd[i]<0)(col[i]?c1:c0)++;
		int keep=c1<c0,j=0,k=0;
		for(auto [id,w]:a){
			ll v=w/2;if(w&1)v+=col[k++]==keep;
			if(v)a[j++]={id,v};
		}
		a.resize(j);q/=2;
	}
	vector<int> ans(n);
	for(auto [id,w]:a){assert(id>0&&w==1);ans[e[id].u-1]=id;}
	assert((int)a.size()==n);
	return ans;
}
void divide(vector<int> a,int d,vector<vector<int>> &ans)const{
	if(!d)return;
	if(d==1){
		vector<int> b(n);
		for(int id:a)b[e[id].u-1]=id;
		ans.push_back(move(b));
		return;
	}
	if(d&1){
		ans.push_back(perfect(a,d));const auto &b=ans.back();
		a.erase(remove_if(a.begin(),a.end(),
			[&](int id){return b[e[id].u-1]==id;}),a.end());
		d--;
	}
	auto col=split(a);
	vector<int> l,r;
	l.reserve(a.size()/2);
	r.reserve(a.size()/2);
	for(int i=0;i<(int)a.size();i++)(col[i]?r:l).push_back(a[i]);
	divide(move(l),d/2,ans);
	// Borrow some completed matchings, so the other degree is a power of two.
	int k=1;
	while(k<=((d-1)/2))k*=2;
	for(int i=0;i<k-d/2;i++){
		r.insert(r.end(),ans.back().begin(),ans.back().end());ans.pop_back();
	}
	divide(move(r),k,ans);
}
public:
vector<vector<int>> solve()const{
	if(!n)return {};
	vector<int> l(n+1),r(n+1),ids(e.size()-1);
	for(int id=1;id<(int)e.size();id++)l[e[id].u]++,r[e[id].v]++,ids[id-1]=id;
	int d=l[1];
	for(int u=1;u<=n;u++)assert(l[u]==d&&r[u]==d);
	vector<vector<int>> ans;
	ans.reserve(d);
	divide(move(ids),d,ans);
	for(auto &a:ans){
		vector<int> b(n+1);
		for(int u=1;u<=n;u++)b[u]=e[a[u-1]].v;
		a=move(b);
	}
	return ans;
}
};
// end for graph/matching/regular-bipartite.cpp
/////////////////////////
// !!!!! Both sides must have the same size and every vertex the same degree;
// parallel edges count separately. !!!!
