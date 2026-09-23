struct bi_edge{
int n=0,top=0,tot=0,scnt=0;
vector<vector<int>> edg;
vector<int> low,dfn,stk,bel;
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);n=_n;edg.assign(n+1,{});
	low.assign(n+1,0);dfn.assign(n+1,0);stk.assign(n+1,0);bel.assign(n+1,0);
	top=tot=scnt=0;
}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n);
	edg[u].push_back(v);
	edg[v].push_back(u);
}
vector<array<int,2>> build(){
	fill(dfn.begin(),dfn.end(),0);top=tot=scnt=0;
	auto tarjan=[&](auto&& self,int u,int fa)->void{
		low[u]=dfn[u]=++tot;stk[++top]=u;
		bool vis_fa=false;
		for(int v:edg[u]){
			if(v==fa){if(!vis_fa){vis_fa=true;continue;}}
			if(!dfn[v]) self(self,v,u),chkmin(low[u],low[v]);
			else chkmin(low[u],dfn[v]);
		}
		if(low[u]==dfn[u]){
			++scnt;
			int x;
			do{
				x=stk[top--];
				bel[x]=scnt;
			}while(x!=u);
		}
	};
	for(int i=1;i<=n;i++) if(!dfn[i]) tarjan(tarjan,i,0);
	vector<array<int,2>> E;
	for(int i=1;i<=n;i++) for(int j:edg[i]) if(bel[i]<bel[j])
		E.push_back({bel[i],bel[j]});
	return E;
}
};
