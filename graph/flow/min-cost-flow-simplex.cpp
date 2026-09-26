//////////////////////////////////////////////////////////////////
//
// template for min cost max flow (network simplex)
//
// usage:
//   Cost_Flow_Graph g; g.setN(n); g.setST(S,T); g.add_edge(u,v,cap,cost);
//   auto [flow,cost]=g.min_cost_flow(); g.edge_flow(id);
//   negative costs allowed, input graph must have no negative cost cycle
//   setN clears the graph; editing the graph makes the next solve start over
//   flow and cost must fit ll; potentials use ll when safe, otherwise __int128
//
//////////////////////////////////////////////////////////////////
struct Cost_Flow_Graph{
private:
struct edge{int u,v;ll w,c;};
vector<edge> e=vector<edge>(2);
vector<vector<int>> edg=vector<vector<int>>(1);
int n=0,S=0,T=0;ll ans=0,cost=0;bool is_flowed=false;
template<class C>
void simplex(C big,int m){
	auto weight=[&](int id)->C{return id<m?e[id].c:id==m?-big:big;};
	vector<int> fa(n+1),fe(n+1),cur(n+1),path,cir;
	vector<ll> tag(n+1);vector<C> h(n+1);
	path.reserve(n);cir.reserve(n+1);
	// Start with a residual spanning forest, joined to virtual root 0.
	auto build=[&](int rt){
		tag[rt]=1;path.push_back(rt);
		while(!path.empty()){
			int u=path.back(),&i=cur[u];
			if(i==(int)edg[u].size()){path.pop_back();continue;}
			int id=edg[u][i++],v=e[id].v;
			if(!e[id].w||tag[v])continue;
			fa[v]=u;fe[v]=id;h[v]=h[u]+weight(id);tag[v]=1;path.push_back(v);
		}
	};
	build(T);for(int u=1;u<=n;u++)if(!tag[u])build(u);
	ll stamp=1;tag[0]=stamp;
	auto sum=[&](int u){
		int v=u;
		while(tag[v]!=stamp){path.push_back(v);v=fa[v];}
		while(!path.empty()){
			v=path.back();path.pop_back();h[v]=h[fa[v]]+weight(fe[v]);tag[v]=stamp;
		}
		return h[u];
	};
	bool changed=true;
	while(changed){
		changed=false;
		for(int id=2;id<m+2;id++){
			if(!e[id].w||weight(id)+sum(e[id].u)-sum(e[id].v)>=0)continue;
			changed=true;++stamp;tag[0]=stamp;
			for(int u=e[id].u;u;u=fa[u])tag[u]=stamp;
			int lca=e[id].v;
			while(tag[lca]!=stamp)tag[lca]=stamp,lca=fa[lca];
			ll f=e[id].w;int del=0,side=2;cir.clear();
			for(int u=e[id].u;u!=lca;u=fa[u]){
				int x=fe[u];cir.push_back(x);
				if(e[x].w<f)f=e[x].w,del=u,side=0;
			}
			for(int u=e[id].v;u!=lca;u=fa[u]){
				int x=fe[u]^1;cir.push_back(x);
				if(e[x].w<f)f=e[x].w,del=u,side=1;
			}
			cir.push_back(id);
			if(f)for(int x:cir)e[x].w-=f,e[x^1].w+=f;
			if(side==2)continue;
			int u=e[id].u,v=e[id].v,x=id^side;
			if(side)swap(u,v);
			// Only invalidate the part of the tree whose parents change.
			while(v!=del){
				x^=1;tag[u]=0;swap(fe[u],x);int p=fa[u];fa[u]=v;v=u;u=p;
			}
		}
	}
}
public:
void set(){n=S=T=0;ans=cost=0;is_flowed=false;e.assign(2,{});edg.assign(1,{});}
void setN(int _n){
	assert(0<=_n&&_n<INT_MAX);set();n=_n;edg.resize(n+1);
}
void setST(int _S,int _T){
	assert(1<=min(_S,_T)&&max(_S,_T)<INT_MAX&&_S!=_T);
	S=_S;T=_T;n=max({n,S,T});edg.resize(n+1);is_flowed=false;
}
int new_node(){
	assert(n<INT_MAX-1);++n;edg.push_back({});is_flowed=false;return n;
}
ll edge_flow(int id){
	assert(2<=id&&id<(int)e.size()&&id%2==0);return e[id^1].w;
}
int add_edge(int u,int v,ll w,ll c){
	assert(1<=min(u,v)&&max(u,v)<=n&&w>=0&&c!=LLONG_MIN);
	assert(e.size()<INT_MAX-3);
	int id=e.size();e.push_back({u,v,w,c});e.push_back({v,u,0,-c});
	edg[u].push_back(id);edg[v].push_back(id^1);is_flowed=false;return id;
}
pair<ll,ll> min_cost_flow(){
	assert(1<=min(S,T)&&max(S,T)<=n&&S!=T);
	if(is_flowed)return {ans,cost};
	int m=e.size();__int128 big=1;
	for(int id=2;id<m;id+=2){
		e[id].w+=e[id^1].w;e[id^1].w=0;
		big+=e[id].c<0?-(__int128)e[id].c:e[id].c;
	}
	add_edge(T,S,LLONG_MAX,0);
	// Every tree path has at most n edges; bound the reduced-cost arithmetic.
	if(big<=LLONG_MAX/(2LL*n+1))simplex<ll>((ll)big,m);
	else simplex<__int128>(big,m);
	ans=e[m^1].w;edg[T].pop_back();edg[S].pop_back();e.resize(m);
	// Sum opposite signs alternately; never multiply the artificial edge cost.
	constexpr __int128 mx=(((__int128)1<<126)-1)*2+1;
	__int128 c=0;int pos=2,neg=2;
	while(true){
		while(pos<m&&e[pos].c<=0)pos+=2;
		while(neg<m&&e[neg].c>=0)neg+=2;
		if(pos==m&&neg==m)break;
		int &id=(pos==m||(c>0&&neg<m))?neg:pos;
		__int128 x=(__int128)e[id].c*e[id^1].w;
		if(x>0)assert(c<=mx-x);else assert(c>=-mx-1-x);
		c+=x;id+=2;
	}
	assert(LLONG_MIN<=c&&c<=LLONG_MAX);
	cost=(ll)c;is_flowed=true;return {ans,cost};
}
};
// end for graph/flow/min-cost-flow-simplex.cpp
/////////////////////////
// !!!!! Only paste one Cost_Flow_Graph variant; input must have no negative
// cost cycle; flow and total cost must fit ll. !!!!
