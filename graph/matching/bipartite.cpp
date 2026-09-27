////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: Bipartite_Matching g; g.setN(n,m); g.add_edge(u,v);
//   auto a=g.solve(); // pairs {u,v}, ordered by the left endpoint
//
////////////////////////////////////////////////////////////////
struct Bipartite_Matching{
struct edge{int u,v;};
int n=0,m=0;
vector<edge> e;
vector<int> matchL,matchR;
void set(){n=m=0;e.clear();matchL.clear();matchR.clear();}
void setN(int _n,int _m){
	assert(0<=_n&&_n<INT_MAX-1&&0<=_m&&_m<INT_MAX);
	set();
	n=_n;
	m=_m;
	matchL.assign(n+1,0);matchR.assign(m+1,0);
}
void add_edge(int u,int v){
	assert(1<=u&&u<=n&&1<=v&&v<=m&&e.size()<INT_MAX);
	e.push_back({u,v});
}
// Rebuild matching; return pairs ordered by left endpoint.
// matchL / matchR store partners, 0 if unmatched; left 1..n, right 1..m.
// Greedy + multi-source BFS: O(n+E) per round, O(n+m+E) space.
vector<pair<int,int>> solve(){
	matchL.assign(n+1,0);matchR.assign(m+1,0);
	vector<int> head(n+2),to(e.size()),rt(n+1),pre(n+1),q;
	for(auto [u,v]:e)head[u+1]++;
	partial_sum(head.begin(),head.end(),head.begin());
	copy(head.begin(),head.begin()+n+1,rt.begin());
	for(auto [u,v]:e)to[rt[u]++]=v;
	q.reserve(n);
	for(int u=1;u<=n;u++)for(int i=head[u];i<head[u+1];i++)if(!matchR[to[i]]){
		matchL[u]=to[i];
		matchR[to[i]]=u;
		break;
	}
	while(true){
		fill(rt.begin(),rt.end(),0);
		q.clear();
		bool found=false;
		for(int u=1;u<=n;u++)if(!matchL[u]&&head[u]!=head[u+1]){
			rt[u]=u;
			pre[u]=0;
			q.push_back(u);
		}
		for(int h=0;h<(int)q.size();h++){
			int u=q[h];if(matchL[rt[u]])continue;
			for(int i=head[u];i<head[u+1];i++){
				int v=to[i],w=matchR[v];
				if(!w){
					found=true;
					for(int x=u;v;x=pre[x]){
						int t=matchL[x];
						matchL[x]=v;
						matchR[v]=x;
						v=t;
					}
					break;
				}
				if(!rt[w])rt[w]=rt[u],pre[w]=u,q.push_back(w);
			}
		}
		if(!found)break;
	}
	vector<pair<int,int>> ans;ans.reserve(min(n,m));
	for(int u=1;u<=n;u++)if(matchL[u])ans.emplace_back(u,matchL[u]);
	return ans;
}
};
// end for graph/matching/bipartite.cpp
/////////////////////////
