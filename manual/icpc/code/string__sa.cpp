struct SA{
int n=0,m=0;
vector<int> sa,rk;
vector<unsigned char> s;
void set(){n=m=0;sa.clear();rk.clear();s.clear();}
void build(string S,int M=128){
	assert(!S.empty()&&S[0]==' '&&1<=M&&M<=256);
	assert(S.size()<INT_MAX/2);
	n=S.size()-1;m=M;
	sa.assign(n+1,0);rk.assign(n+1,0);s.assign(n+1,0);
	vector<int> y(n+1),old(n+1),cnt(max(n+1,M),0);
	for(int i=1;i<=n;i++){
		s[i-1]=(unsigned char)S[i];
		assert(0<s[i-1]&&s[i-1]<M);
		rk[i]=s[i-1];cnt[rk[i]]++;
	}
	for(int i=1;i<m;i++)cnt[i]+=cnt[i-1];
	for(int i=n;i;i--)sa[cnt[rk[i]]--]=i;
	for(int k=1;k<n;k<<=1){
		int p=0;
		for(int i=n-k+1;i<=n;i++)y[++p]=i;
		for(int i=1;i<=n;i++)if(sa[i]>k)y[++p]=sa[i]-k;
		fill(cnt.begin(),cnt.end(),0);
		for(int i=1;i<=n;i++)cnt[rk[i]]++;
		for(int i=1;i<m;i++)cnt[i]+=cnt[i-1];
		for(int i=n;i;i--)sa[cnt[rk[y[i]]]--]=y[i];
		old=rk;p=0;
		for(int i=1;i<=n;i++){
			int x=sa[i],z=sa[i-1];
			if(i==1||old[x]!=old[z]||
				(x+k<=n?old[x+k]:0)!=(z+k<=n?old[z+k]:0))p++;
			rk[x]=p;
		}
		m=p+1;if(p==n)break;
	}
	for(int i=1;i<=n;i++)rk[sa[i]]=i;
}
void output(){
	cerr<<"string: ";for(int i=0;i<n;i++)cerr<<s[i];cerr<<'\n';
	cerr<<"sa: ";for(int i=1;i<=n;i++)cerr<<sa[i]<<" \n"[i==n];
	cerr<<"rk: ";for(int i=1;i<=n;i++)cerr<<rk[i]<<" \n"[i==n];
}
void solve(){for(int i=1;i<=n;i++)cout<<sa[i]<<" \n"[i==n];}
};
// SA_LCP L; L.build(sa); L.lcp(i,j); // 1-indexed
// O(n) LCP build, O(1) query; paste Linear_RMQ first.
struct SA_LCP{
int n=0;
vector<int> rk,height;
Linear_RMQ<int> st;
void build(const SA &S){
	n=S.n;rk=S.rk;height.assign(n+1,0);st.set();
	assert((int)S.sa.size()==n+1&&(int)rk.size()==n+1);
	for(int i=1,k=0;i<=n;i++){
		int r=rk[i];
		if(r==1){k=0;continue;}
		int j=S.sa[r-1];
		while(i+k<=n&&j+k<=n&&S.s[i+k-1]==S.s[j+k-1])k++;
		height[r]=k;if(k)k--;
	}
	st.build(height);
}
void build(const string &s,int M=128){SA S;S.build(s,M);build(S);}
int lcp(int i,int j)const{
	assert(1<=min(i,j)&&max(i,j)<=n);
	if(i==j)return n-i+1;
	int l=rk[i],r=rk[j];if(l>r)swap(l,r);l++;
	return st.ask(l,r);
}
int ask(int i,int j)const{return lcp(i,j);}
};
