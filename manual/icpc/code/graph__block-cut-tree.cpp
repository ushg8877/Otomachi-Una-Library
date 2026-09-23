struct block_cut{
int n=0,tot=0,scc=0,top=0;
vector<char> inq;
vector<int> low,dfn,stk,dep,ff;
vector<vector<int>> edg,tr;
vector<vector<array<int,2>>> ec;
void init(int _n){
	assert(0<=_n&&_n<(INT_MAX-2)/2);n=_n;
	low.assign(n+2,0);dfn.assign(n+2,0);stk.assign(n+2,0);inq.assign(n+2,0);
	edg.assign(n+2,{});tr.assign(2*n+2,{});ec.assign(n+2,{});
	dep.assign(2*n+2,0);ff.assign(2*n+2,0);tot=scc=top=0;
}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n&&u!=v);
	edg[u].push_back(v);
	edg[v].push_back(u);
}
void tarjan(int u,int fa=0){
	low[u]=dfn[u]=++tot;stk[++top]=u;inq[u]=true;
	bool vis_fa=false;
	for(int v:edg[u]){
		if(v==fa&&!vis_fa){vis_fa=true;continue;}
		if(!dfn[v]){
			tarjan(v,u);
			if(low[v]>=dfn[u]){
				++scc;
				tr[scc].push_back(u);
				tr[u].push_back(scc);
				while(1){
					int w=stk[top--];
					tr[w].push_back(scc);
					tr[scc].push_back(w);
					inq[w]=false;
					if(w==v)break;
				}
			}
			low[u]=min(low[u],low[v]);
		}else if(inq[v]) low[u]=min(low[u],dfn[v]);
	}
}
void dfs(int u,int fa){
	ff[u]=fa;dep[u]=dep[fa]+1;
	for(int v:tr[u]) if(v!=fa) dfs(v,u);
}
void build(){
	fill(dfn.begin(),dfn.end(),0);fill(inq.begin(),inq.end(),0);
	fill(dep.begin(),dep.end(),0);fill(ff.begin(),ff.end(),0);
	for(auto &v:tr)v.clear();for(auto &v:ec)v.clear();
	tot=top=0;scc=n;
	for(int i=1;i<=n;i++)if(!dfn[i]){
		tarjan(i);dfs(i,0);inq[i]=false;top=0;
	}
	for(int i=1;i<=n;i++)for(int j:edg[i])if(i<j){
		int x=ff[i]==ff[j]?ff[i]:(dep[i]>dep[j]?ff[i]:ff[j]);
		assert(n<x&&x<=scc);ec[x-n].push_back({i,j});
	}
}
};
