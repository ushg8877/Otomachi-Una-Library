struct du_jiao{
vector<int> mu;
vector<ll> phi;
unordered_map<ll,ll> mem_mu;
unordered_map<ll,__int128> mem_phi;
void init_mu(int m){
	assert(m>=1&&m<INT_MAX);
	vector<int> prime;vector<bool> vis(m+1,false);
	mu.assign(m+1,0);mem_mu.clear();mu[1]=1;
	for(int i=2;i<=m;i++){
		if(!vis[i])prime.push_back(i),mu[i]=-1;
		for(int p:prime){
			if(p>m/i)break;
			int x=i*p;vis[x]=true;
			if(i%p==0)break;
			mu[x]=-mu[i];
		}
	}
	for(int i=1;i<=m;i++)mu[i]+=mu[i-1];
}
void init_phi(int m){
	assert(m>=1&&m<INT_MAX);
	vector<int> prime;vector<bool> vis(m+1,false);
	phi.assign(m+1,0);mem_phi.clear();phi[1]=1;
	for(int i=2;i<=m;i++){
		if(!vis[i])prime.push_back(i),phi[i]=i-1;
		for(int p:prime){
			if(p>m/i)break;
			int x=i*p;vis[x]=true;
			if(i%p==0){phi[x]=phi[i]*p;break;}
			phi[x]=phi[i]*(p-1);
		}
	}
	for(int i=1;i<=m;i++)phi[i]+=phi[i-1];
}
ll sum_mu(ll n){
	assert(n>=0&&!mu.empty());
	if(n<(ll)mu.size())return mu[n];
	auto it=mem_mu.find(n);
	if(it!=mem_mu.end())return it->second;
	__int128 ans=1;
	for(ll l=2,r;l<=n;){
		r=n/(n/l);ans-=(__int128)(r-l+1)*sum_mu(n/l);
		if(r==n)break;
		l=r+1;
	}
	return mem_mu[n]=(ll)ans;
}
__int128 sum_phi(ll n){
	assert(n>=0&&!phi.empty());
	if(n<(ll)phi.size())return phi[n];
	auto it=mem_phi.find(n);
	if(it!=mem_phi.end())return it->second;
	__int128 ans=(__int128)n*(n+(__int128)1)/2;
	for(ll l=2,r;l<=n;){
		r=n/(n/l);ans-=(__int128)(r-l+1)*sum_phi(n/l);
		if(r==n)break;
		l=r+1;
	}
	return mem_phi[n]=ans;
}
};
