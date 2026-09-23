using namespace std;
int n,m,cnt,ans;vector<int> h,pre,fa,c,mak;
vector<vector<int>> edg;queue<int> q;
inline int find(int x){return x==fa[x]?x:fa[x]=find(fa[x]);}
int lca(int x,int y){
	for(cnt++;mak[x]!=cnt;swap(x,y)) if(x) mak[x]=cnt,x=find(pre[h[x]]);
	return x;
}
void blo(int u,int v,int rt){
	for(;find(u)!=rt;u=pre[v]){
		pre[u]=v;fa[u]=fa[v=h[u]]=rt;
		if(c[v]==1)q.push(v),c[v]=2;
	}
}
bool bfs(int rt){
	for(int i=1;i<=n;i++) fa[i]=i,pre[i]=c[i]=0;
	while(!q.empty()) q.pop();
	c[rt]=2;q.push(rt);
	while(!q.empty()){
		int u=q.front();q.pop();
		for(int v:edg[u]) if(!c[v]){
			pre[v]=u;c[v]=1;
			if(!h[v]){
				while(u) u=h[pre[v]],h[h[v]=pre[v]]=v,v=u;
				return true;
			}else c[h[v]]=2,q.push(h[v]);
		}else if(c[v]==2){
			int l=lca(u,v);
			blo(u,v,l);blo(v,u,l);
		}
	}
	return false;
}
