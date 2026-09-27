////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/graph/edge-biconnected.cpp
//
// usage: bi_edge b; b.setN(n); b.add_edge(u,v); auto E = b.build();
//   bel[u] is EBCC id; E lists edges between different blocks
//
////////////////////////////////////////////////////////////////
struct bi_edge{
int n=0,top=0,tot=0,scnt=0;
vector<vector<int>> edg;
vector<int> low,dfn,stk,bel;
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);
	n=_n;
	edg.assign(n+1,{});
	low.assign(n+1,0);
	dfn.assign(n+1,0);
	stk.assign(n+1,0);
	bel.assign(n+1,0);
	top=tot=scnt=0;
}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n);
	edg[u].push_back(v);
	edg[v].push_back(u);
}
vector<array<int,2>> build(){
	fill(dfn.begin(),dfn.end(),0);top=tot=scnt=0;
	for(int u=1;u<=n;u++)if(!dfn[u]){
		vector<array<int,4>> q{{u,0,0,0}};
		low[u]=dfn[u]=++tot;stk[++top]=u;
		while(!q.empty()){
			auto &[x,p,i,skip]=q.back();
			if(i<(int)edg[x].size()){
				int v=edg[x][i++];
				if(v==p&&!skip){skip=1;continue;}
				if(!dfn[v]){
					low[v]=dfn[v]=++tot;
					stk[++top]=v;
					q.push_back({v,x,0,0});
				}
				else low[x]=min(low[x],dfn[v]);
			}else{
				int v=x,parent=p;q.pop_back();
				if(parent)low[parent]=min(low[parent],low[v]);
				if(low[v]==dfn[v]){
					scnt++;
					int y;
					do{y=stk[top--];bel[y]=scnt;}while(y!=v);
				}
			}
		}
	}
	vector<array<int,2>> E;
	for(int i=1;i<=n;i++) for(int j:edg[i]) if(bel[i]<bel[j])
		E.push_back({bel[i],bel[j]});
	return E;
}
};
// end for graph/edge-biconnected.cpp
/////////////////////////
