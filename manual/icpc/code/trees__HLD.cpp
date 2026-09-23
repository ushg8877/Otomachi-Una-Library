struct Tree{
////////////////////////////////// basic tree variable
int n=0,rt=1,tot=0;
vector<vector<int>> edg,st;
vector<int> in,out,siz,ff,dep,dfn,hson,top,ed;
////////////////////////////////// self variable
void setN(int _n){
	assert(1<=_n&&_n<INT_MAX-1);n=_n;tot=0;rt=1;
	edg.assign(n+1,{});st.assign(__lg(n)+1,vector<int>(n+1));
	in.assign(n+1,0);
	out.assign(n+1,0);
	siz.assign(n+1,0);
	ff.assign(n+1,0);
	dep.assign(n+1,0);
	dfn.assign(n+1,0);
	hson.assign(n+1,0);
	top.assign(n+1,0);
	ed.assign(n+1,0);
}
void setRt(int _rt){assert(1<=_rt&&_rt<=n);rt=_rt;}
void add_edge(int u,int v){
	assert(1<=min(u,v)&&max(u,v)<=n&&u!=v);
	edg[u].push_back(v);
	edg[v].push_back(u);
}
void HLD(int u,int fa){
	assert(!siz[u]);siz[u]=1;ff[u]=fa;
	for(int v:edg[u]) if(v!=fa){
		dep[v]=dep[u]+1;HLD(v,u);siz[u]+=siz[v];
		if(siz[hson[u]]<siz[v]) hson[u]=v;
	}
} 
void HLD1(int u,int fa){
	++tot;in[u]=tot;dfn[tot]=u;st[0][tot]=fa;
	if(hson[u]){top[hson[u]]=top[u];HLD1(hson[u],u);ed[u]=ed[hson[u]];}
	else ed[u]=u;
	for(int v:edg[u]) if(v!=fa&&v!=hson[u]){
		top[v]=v;HLD1(v,u);
	}
	out[u]=tot;
}
inline int higher(int x,int y){return dep[x]<dep[y]?x:y;}
void build(){
	// build all basic information in O(nlogn) (without lca O(n))
	assert(1<=rt&&rt<=n);tot=0;
	fill(hson.begin(),hson.end(),0);fill(siz.begin(),siz.end(),0);
	dep[rt]=1;top[rt]=rt;HLD(rt,0);HLD1(rt,0);assert(tot==n);
	for(int i=1;i<(int)st.size();i++) for(int j=1;j<=n-(1<<i)+1;j++)
		st[i][j]=higher(st[i-1][j],st[i-1][j+(1<<(i-1))]);
}
inline int lca(int x,int y){
	// find lca ancestor in O(1)
	assert(1<=min(x,y)&&max(x,y)<=n);
	if(x==y) return x;
	x=in[x],y=in[y];
	if(x>y) swap(x,y);
	int k=__lg(y-(x++));
	return higher(st[k][x],st[k][y-(1<<k)+1]);
}
inline int dist(int u,int v){
	return dep[u]+dep[v]-2*dep[lca(u,v)];
}
inline int kth_anc(int u,int k){
	// find k-th ancestor in O(logn)
	assert(1<=u&&u<=n&&0<=k&&k<dep[u]);
	while(1){
		if(dep[u]-dep[top[u]]>=k) return dfn[in[u]-k];
		k-=dep[u]-dep[top[u]]+1;
		u=ff[top[u]];
		if(!u) return rt;
	}
}
inline vector<array<int,2>> index_of_path(int x,int y){
	// return all in ids in path (x,y)
	assert(1<=min(x,y)&&max(x,y)<=n);
	vector<array<int,2>> id;
	while(top[x]!=top[y]){
		if(dep[top[x]]<dep[top[y]]) swap(x,y);
		id.push_back({in[top[x]],in[x]});
		x=ff[top[x]];
	}
	if(dep[x]>dep[y]) swap(x,y);
	id.push_back({in[x],in[y]});
	return id;
}
void output(){
	cerr<<"***********************************\n";
	cerr<<"structure:"<<endl;
	cerr<<"root = "<<rt<<endl;
	for(int i=1;i<=n;i++) for(int j:edg[i]) if(i<j)
		cerr<<i<<' '<<j<<'\n';
	cerr<<"***********************************\n";
	cerr<<"dfn-order: ";
	for(int i=1;i<=n;i++) cerr<<dfn[i]<<" \n"[i==n];
	cerr<<"top: ";
	for(int i=1;i<=n;i++) cerr<<top[i]<<" \n"[i==n];
}
};
//!!!!!!!!!! do not forget setN() and build() !!!!!!!!!!
