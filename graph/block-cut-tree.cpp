////////////////////////////////////////////////////////////////
//
// template for Block Cut Tree
// use add_edge to add (u,v), tr[n+i] for vertex in i-th component
// ec[i] for edge in i-th component
// version 1.0 (Last Update Jun 17th, 2026)
//
// usage: block_cut bc; bc.setN(n); bc.add_edge(u,v); bc.build();
//   vertex u -> tree u; BCC i -> tree n+i; bc.ec[i] edges
//
////////////////////////////////////////////////////////////////
struct block_cut{
int n=0,tot=0,scc=0,top=0;
vector<char> inq;
vector<int> low,dfn,stk,dep,ff;
vector<vector<int>> edg,tr;
vector<vector<array<int,2>>> ec;
void setN(int _n){
	assert(0<=_n&&_n<(INT_MAX-2)/2);n=_n;
	low.assign(n+2,0);
	dfn.assign(n+2,0);
	stk.assign(n+2,0);
	inq.assign(n+2,0);
	edg.assign(n+2,{});
	tr.assign(2*n+2,{});
	ec.assign(n+2,{});
	dep.assign(2*n+2,0);
	ff.assign(2*n+2,0);
	tot=scc=top=0;
}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n&&u!=v);
	edg[u].push_back(v);
	edg[v].push_back(u);
}
void tarjan(int u,int fa=0){
	// frame: vertex, parent, next edge, skipped parent edge
	vector<array<int,4>> q{{u,fa,0,0}};
	low[u]=dfn[u]=++tot;
	stk[++top]=u;
	inq[u]=true;
	while(!q.empty()){
		auto &[x,p,i,skip]=q.back();
		if(i<(int)edg[x].size()){
			int v=edg[x][i++];
			if(v==p&&!skip){skip=1;continue;}
			if(!dfn[v]){
				low[v]=dfn[v]=++tot;
				stk[++top]=v;
				inq[v]=true;
				q.push_back({v,x,0,0});
			}else if(inq[v])low[x]=min(low[x],dfn[v]);
		}else{
			int v=x,parent=p;q.pop_back();
			if(q.empty())continue;
			low[parent]=min(low[parent],low[v]);
			if(low[v]>=dfn[parent]){
				++scc;
				tr[scc].push_back(parent);
				tr[parent].push_back(scc);
				int y;
				do{
					y=stk[top--];
					inq[y]=false;
					tr[y].push_back(scc);
					tr[scc].push_back(y);
				}
				while(y!=v);
			}
		}
	}
}
void dfs(int u,int fa){
	vector<int> q{u};
	ff[u]=fa;
	dep[u]=dep[fa]+1;
	for(int i=0;i<(int)q.size();i++){
		int x=q[i];
		for(int v:tr[x])if(v!=ff[x]){ff[v]=x;dep[v]=dep[x]+1;q.push_back(v);}
	}
}
void build(){
	fill(dfn.begin(),dfn.end(),0);fill(inq.begin(),inq.end(),0);
	fill(dep.begin(),dep.end(),0);fill(ff.begin(),ff.end(),0);
	for(auto &v:tr)v.clear();
	for(auto &v:ec)v.clear();
	tot=top=0;scc=n;
	for(int i=1;i<=n;i++)if(!dfn[i]){
		tarjan(i);
		dfs(i,0);
		inq[i]=false;
		top=0;
	}
	for(int i=1;i<=n;i++)for(int j:edg[i])if(i<j){
		int x=ff[i]==ff[j]?ff[i]:(dep[i]>dep[j]?ff[i]:ff[j]);
		assert(n<x&&x<=scc);ec[x-n].push_back({i,j});
	}
}
};
// end for graph/block-cut-tree.cpp
/////////////////////////
