struct Tree{
////////////////////////////////// basic tree variable
int n=0,rt=1,tot=0;
vector<vector<int>> edg,st;
vector<int> in,out,siz,ff,dep,dfn,hson,top,ed,msiz,csiz;
vector<char> vis;int tsiz=0,cen=0;
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
	msiz.assign(n+1,0);
	csiz.assign(n+1,0);
	vis.assign(n+1,0);tsiz=cen=0;
}
void init(int _n){setN(_n);}
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
void HLD_build(){
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
void get_size(int u,int fa){
	csiz[u]=msiz[u]=1;
	for(int v:edg[u]) if(v!=fa&&!vis[v]){
		get_size(v,u);
		csiz[u]+=csiz[v];
		msiz[u]=max(msiz[u],csiz[v]);
	}
	msiz[u]=max(msiz[u],tsiz-csiz[u]);
	if(msiz[u]<msiz[cen]) cen=u;
}
void get(int u,int fa){
	for(int v:edg[u]) if(v!=fa&&!vis[v]){
		get(v,u);
	}
}
void DC(int u){
	vis[u]=true;
	for(int v:edg[u]) if(!vis[v]){
		// get some info?
		get(v,0);
	}
	// solve whole
	for(int v:edg[u]) if(!vis[v]){
		tsiz=csiz[v];cen=0;
		get_size(v,0);get_size(cen,0);DC(cen);
	}
	// cerr<<"finish"<<endl;
}
void DC_build(){
	assert(1<=rt&&rt<=n);fill(vis.begin(),vis.end(),0);
	tsiz=n;msiz[cen=0]=1e9;get_size(rt,0);get_size(cen,0);DC(cen);
}
void debug(){
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
//////////////////////////////////////////////////////
// modify from here
void solve(){
	// do not forget init
	// use add_edge(u,v) to add a edge
	// do not forget HLD_build()
	int N;cin>>N;setN(N);
	for(int i=1,f;i<=n;i++){
		cin>>f;
		if(f) add_edge(i,f);
		else rt=i;
	}
	HLD_build();
	DC_build();
}
};
