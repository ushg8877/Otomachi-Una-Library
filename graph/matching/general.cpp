////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/graph/matching/general.cpp
//
// usage: blossom g; g.setN(n); g.add_edge(u,v); int cnt=g.solve();
//   cnt*2==n if perfect; match[u]=0 if unmatched; 1-index
//
////////////////////////////////////////////////////////////////
struct blossom{
int n=0;
vector<vector<int>> edg;
vector<int> match,pre,bel,q;
vector<char> inq,vis,flower;
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);n=_n;
	edg.assign(n+1,{});match.assign(n+1,0);
	pre.resize(n+1);
	bel.resize(n+1);
	q.resize(n+1);
	inq.resize(n+1);
	vis.resize(n+1);
	flower.resize(n+1);
}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n);
	if(u==v)return;
	edg[u].push_back(v);edg[v].push_back(u);
}
int lca(int u,int v){
	fill(vis.begin(),vis.end(),false);
	while(1){
		u=bel[u];vis[u]=true;
		if(!match[u])break;
		u=pre[match[u]];
	}
	while(!vis[bel[v]])v=pre[match[bel[v]]];
	return bel[v];
}
void shrink(int u,int v,int w){
	while(bel[u]!=w){
		flower[bel[u]]=flower[bel[match[u]]]=true;
		pre[u]=v;
		v=match[u];
		u=pre[v];
	}
}
bool BFS(int s){
	fill(inq.begin(),inq.end(),false);
	fill(pre.begin(),pre.end(),0);
	iota(bel.begin(),bel.end(),0);
	int h=0,t=0;
	q[t++]=s;
	inq[s]=true;
	while(h<t){
		int u=q[h++];
		for(int v:edg[u]){
			if(bel[u]==bel[v]||match[u]==v)continue;
			if(v==s||(match[v]&&pre[match[v]])){
				int w=lca(u,v);
				fill(flower.begin(),flower.end(),false);
				shrink(u,v,w);shrink(v,u,w);
				for(int i=1;i<=n;i++)if(flower[bel[i]]){
					bel[i]=w;
					if(!inq[i])inq[i]=true,q[t++]=i;
				}
			}else if(!pre[v]){
				pre[v]=u;
				if(!match[v]){
					while(v){
						int x=pre[v],y=match[x];
						match[v]=x;
						match[x]=v;
						v=y;
					}
					return true;
				}
				int x=match[v];
				inq[x]=true;
				q[t++]=x;
			}
		}
	}
	return false;
}
int solve(){
	fill(match.begin(),match.end(),0);
	int ans=0;
	for(int i=1;i<=n;i++)if(!match[i])ans+=BFS(i);
	return ans;
}
};
// end for graph/matching/general.cpp
/////////////////////////
