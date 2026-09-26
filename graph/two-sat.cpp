////////////////////////////////////////////////////////////////
//
// template for 2-SAT
//
// usage:
//   two_sat sat; sat.set(); sat.setN(n); sat.add(x1,o1,x2,o2);
//   auto ans = sat.solve();  // empty if unsat; var i: 2i-1=T, 2i=F
//   setN clears old clauses
//
////////////////////////////////////////////////////////////////
struct two_sat{
int tot=0,scnt=0,top=0,n=0;
vector<int> low,dfn,stk,bel;
vector<vector<int>> edg;vector<char> inq;
void set(){setN(0);}
void setN(int _n){
	// 1-index, 2i-1 true 2i false
	assert(0<=_n&&_n<(INT_MAX-2)/2);n=_n;
	low.assign(2*n+2,0);
	dfn.assign(2*n+2,0);
	stk.assign(2*n+2,0);
	bel.assign(2*n+2,0);
	inq.assign(2*n+2,0);
	edg.assign(2*n+2,{});
	tot=scnt=top=0;
}
void tarjan(int u){
	vector<pair<int,int>> q{{u,0}};
	low[u]=dfn[u]=++tot;
	stk[++top]=u;
	inq[u]=true;
	while(!q.empty()){
		auto &[x,i]=q.back();
		if(i<(int)edg[x].size()){
			int v=edg[x][i++];
			if(!dfn[v]){
				low[v]=dfn[v]=++tot;
				stk[++top]=v;
				inq[v]=true;
				q.push_back({v,0});
			}else if(inq[v])low[x]=min(low[x],dfn[v]);
		}else{
			int v=x;q.pop_back();
			if(!q.empty()){int p=q.back().first;low[p]=min(low[p],low[v]);}
			if(low[v]==dfn[v]){
				scnt++;int y;
				do{y=stk[top--];inq[y]=false;bel[y]=scnt;}while(y!=v);
			}
		}
	}
}
vector<int> solve(){
	// if impossible return empty
	// else return a vector with length n+1, the a[i] - any valid answer
	tot=scnt=top=0;
	for(int i=1;i<=2*n;i++) low[i]=dfn[i]=bel[i]=inq[i]=0;
	for(int i=1;i<=2*n;i++) if(!dfn[i]) tarjan(i); 
	for(int i=1;i<=n;i++) if(bel[2*i-1]==bel[2*i]) return {};
	vector<int> ans(n+1);
	for(int i=1;i<=n;i++) ans[i]=bel[2*i-1]<bel[2*i];
	return ans;
}
void add(int x1,int o1,int x2,int o2){
	// x1 = o1 or x2 = o2
	assert(1<=min(x1,x2)&&max(x1,x2)<=n&&0<=min(o1,o2)&&max(o1,o2)<=1);
	edg[2*x1-(!o1)].push_back(2*x2-o2);
	edg[2*x2-(!o2)].push_back(2*x1-o1);
}
};
// end for graph/two-sat.cpp
/////////////////////////
