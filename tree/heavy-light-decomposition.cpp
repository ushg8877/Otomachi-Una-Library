////////////////////////////////////////////////////////////////
//
// template for Heavy-Light Decomposition, prepared by Otomachi Una
// version 1.0 (Last Update Jun 13th, 2026)
//
// usage: Tree tr; tr.setN(n); tr.add_edge(u,v); tr.build();
//   tr.lca(u,v); tr.index_of_path(u,v);
//
////////////////////////////////////////////////////////////////
struct Tree{
////////////////////////////////// basic tree variable
int n=0,rt=1,tot=0;
vector<vector<int>> edg,st;
vector<int> in,out,siz,ff,dep,dfn,hson,top,ed;
////////////////////////////////// self variable
void setN(int _n){
	assert(1<=_n&&_n<INT_MAX-1);
	n=_n;
	tot=0;
	rt=1;
	edg.assign(n+1,{});st.resize(__lg(n)+1);
	for(int i=0;i<(int)st.size();i++)st[i].assign(n-(1<<i)+2,0);
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
	vector<int> q{u};
	assert(!siz[u]);
	siz[u]=1;
	ff[u]=fa;
	for(int i=0;i<(int)q.size();i++){
		int x=q[i];
		for(int v:edg[x])if(v!=ff[x]){
			assert(!siz[v]);
			siz[v]=1;
			ff[v]=x;
			dep[v]=dep[x]+1;
			q.push_back(v);
		}
	}
	for(int i=(int)q.size()-1;i>0;i--){
		int x=q[i],p=ff[x];siz[p]+=siz[x];
		if(siz[x]>siz[hson[p]])hson[p]=x;
	}
} 
void HLD1(int u,int fa){
	int start=tot;
	vector<int> q{u};
	ff[u]=fa;
	while(!q.empty()){
		int x=q.back();
		q.pop_back();
		in[x]=++tot;
		dfn[tot]=x;
		st[0][tot]=ff[x];
		for(int i=(int)edg[x].size()-1;i>=0;i--){
			int v=edg[x][i];if(v==ff[x]||v==hson[x])continue;
			top[v]=v;q.push_back(v);
		}
		if(hson[x]){top[hson[x]]=top[x];q.push_back(hson[x]);}
		out[x]=in[x]+siz[x]-1;
	}
	for(int i=tot;i>start;i--){int x=dfn[i];ed[x]=hson[x]?ed[hson[x]]:x;}
}
inline int higher(int x,int y){return dep[x]<dep[y]?x:y;}
// Call setN and add edges before build; queries require build first.
void build(){
	// build all basic information in O(nlogn) (without lca O(n))
	assert(1<=rt&&rt<=n);tot=0;
	fill(hson.begin(),hson.end(),0);fill(siz.begin(),siz.end(),0);
	dep[rt]=1;
	top[rt]=rt;
	HLD(rt,0);
	HLD1(rt,0);
	assert(tot==n);
	for(int i=1;i<(int)st.size();i++) for(int j=1;j<=n-(1<<i)+1;j++)
		st[i][j]=higher(st[i-1][j],st[i-1][j+(1<<(i-1))]);
}
inline int lca(int x,int y){
	// find lca ancestor in O(1)
	assert(1<=min(x,y)&&max(x,y)<=n);
	if(x==y) return x;
	x=in[x],y=in[y];
	if(x>y) swap(x,y);
	x++;int k=__lg(y-x+1);
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
// end for tree/heavy-light-decomposition.cpp
/////////////////////////
