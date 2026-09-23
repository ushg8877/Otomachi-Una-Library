using namespace std;
#define ll long long
#define MP make_pair
const ll inf=2020000000000000000ll;
struct edge{int v,w,c,nxt;};
struct Dinic{
vector<edge> e=vector<edge>(2);
vector<int> lst{0},cur{0};vector<ll> dis{0};vector<char> inq{0};
int tot=1,cnt=0,S=0,T=0;
void set_node_number(int x){
	assert(0<=x&&x<INT_MAX);cnt=x;tot=1;e.assign(2,{0,0,0,0});
	lst.assign(cnt+1,0);dis.assign(cnt+1,0);cur.assign(cnt+1,0);inq.assign(cnt+1,0);
}
void init(){set_node_number(0);S=T=0;}
int newnode(){
	assert(cnt<INT_MAX-1);lst.push_back(0);dis.push_back(0);cur.push_back(0);inq.push_back(0);
	return ++cnt;
}
void add(int u,int v,int w,int c){
	assert(1<=min(u,v)&&max(u,v)<=cnt&&w>=0&&tot<INT_MAX-1);
	e.push_back({v,w,c,lst[u]});lst[u]=++tot;
}
int Add(int u,int v,int w,int c){
	assert(c!=INT_MIN);add(u,v,w,c);add(v,u,0,-c);
	return tot-1;
}
bool flowed(int i){assert(2<=i&&i<=tot&&i%2==0);return e[i^1].w;}
bool SPFA(){
	for(int i=1;i<=cnt;i++) dis[i]=inf,inq[i]=false;
	dis[S]=0;queue<int> q;q.push(S);
	while(!q.empty()){
		int u=q.front();q.pop();inq[u]=false;
		for(int id=lst[u];id;id=e[id].nxt) if(e[id].w&&dis[e[id].v]>dis[u]+e[id].c){
			dis[e[id].v]=dis[u]+e[id].c;
			if(!inq[e[id].v]) inq[e[id].v]=true,q.push(e[id].v);
		}
	}
	return dis[T]<inf;
}
int flow(int u,int fl){
	if(u==T) return fl;
	int rmn=fl;
	inq[u]=true;
	for(int id=cur[u];id&&rmn;id=e[id].nxt){
		cur[u]=id;
		if(!inq[e[id].v]&&e[id].w&&dis[e[id].v]==dis[u]+e[id].c){
			int r=flow(e[id].v,min(e[id].w,rmn));
			rmn-=r;e[id].w-=r;e[id^1].w+=r;
		}
	}
	inq[u]=false;
	return fl-rmn;
}
void MCMF(ll &fl,ll &co){
	assert(1<=min(S,T)&&max(S,T)<=cnt&&S!=T);
	fl=co=0;
	while(SPFA()){
		for(int i=1;i<=cnt;i++) cur[i]=lst[i];
		ll r=flow(S,INT_MAX);
		assert(fl<=LLONG_MAX-r);fl+=r;
		__int128 c=(__int128)co+r*(__int128)dis[T];
		assert(LLONG_MIN<=c&&c<=LLONG_MAX);co=(ll)c;
	}
	return;
}
}G;
