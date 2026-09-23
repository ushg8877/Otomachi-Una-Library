struct SA{
// this is template for Suffix Array, build sa in O(n) by SAIS
int n=0,m=0;
vector<int> lms,cnt,sa,rk;
vector<unsigned char> s;
vector<char> typ;
void init(){
	n=m=0;lms.clear();cnt.clear();sa.clear();rk.clear();s.clear();typ.clear();
}

template<typename T>
void sais(int n,int m,T s[],char typ[],int lms[],int cnt[]){
	int n1=0,ch=-1;
	typ[n-1]=1;
	for(int i=n-2;i>=0;i--) typ[i]=s[i]==s[i+1]?typ[i+1]:s[i]<s[i+1];
	memset(cnt,0,sizeof(int)*m);
	for(int i=0;i<n;i++) cnt[(int)s[i]]++;
	partial_sum(cnt,cnt+m,cnt);
	for(int i=1;i<n;i++) if(typ[i-1]==0&&typ[i]==1)
		typ[i]=2,lms[n1++]=i;
	auto induced_sort=[&](const int v[]){
		fill(sa.begin(),sa.begin()+n,0);
		int *cur=cnt+m;
		auto pL=[&](const int x){sa[cur[(int)s[x]]++]=x;};
		auto pR=[&](const int x){sa[--cur[(int)s[x]]]=x;};
		copy(cnt,cnt+m,cur);
		for(int i=n1-1;i>=0;i--) pR(v[i]);
		copy(cnt,cnt+m-1,cur+1);
		for(int i=0;i<n;i++) if(sa[i]>0&&typ[sa[i]-1]==0) pL(sa[i]-1);
		copy(cnt,cnt+m,cur);
		for(int i=n-1;i>=0;i--) if(sa[i]>0&&typ[sa[i]-1]) pR(sa[i]-1);
	};
	auto lms_equal=[&](int x,int y){
		if(s[x]==s[y]) while(s[++x]==s[++y]&&typ[x]==typ[y])
			if(typ[x]==2||typ[y]==2) return true;
		return false;
	};
	induced_sort(lms);
	int *s1=remove_if(sa.data(),sa.data()+n,[&](int x){return typ[x]!=2;});
	for(int i=0;i<n1;i++) s1[sa[i]>>1]=ch+=ch<=0||!lms_equal(sa[i],sa[i-1]);
	for(int i=0;i<n1;i++) s1[i]=s1[lms[i]>>1];
	if(ch+1<n1) sais(n1,ch+1,s1,typ+n,lms+n1,cnt+m);
	else for(int i=0;i<n1;i++) sa[s1[i]]=i;
	for(int i=0;i<n1;i++) lms[n1+i]=lms[sa[i]];
	induced_sort(lms+n1);
}
void output(){
	// this allow you to see S, sa, rk
	cerr<<"string: ";for(int i=1;i<=n;i++) cerr<<s[i-1];cerr<<'\n';
	cerr<<"sa: ";for(int i=1;i<=n;i++) cerr<<sa[i]<<" \n"[i==n];
	cerr<<"rk: ";for(int i=1;i<=n;i++) cerr<<rk[i]<<" \n"[i==n];
}
void build(string S,int M=128){
	// M for the size of alphabet; 0 is reserved for the sentinel
	assert(!S.empty()&&S[0]==' '&&1<=M&&M<=256&&S.size()<(INT_MAX-1024)/2);
	n=(int)S.size()-1;m=M;
	lms.assign(2*n+2,0);cnt.assign(2*(n+m)+4,0);sa.assign(2*n+2,0);
	rk.assign(n+1,0);s.assign(n+1,0);typ.assign(2*n+2,0);
	for(int i=0;i<n;i++){
		s[i]=(unsigned char)S[i+1];assert(0<s[i]&&s[i]<m);
	}
	sais(n+1,m,s.data(),typ.data(),lms.data(),cnt.data());
	for(int i=1;i<=n;i++)sa[i]++,rk[sa[i]]=i;
	sa.resize(n+1);
}
void solve(){
	// please, write your code from here
	for(int i=1;i<=n;i++) cout<<sa[i]<<" \n"[i==n];
}
};

// SA_LCP L; L.build(sa); L.lcp(i,j); // 1-indexed
// O(n log n) build, O(1) query
struct SA_LCP{
int n=0;
vector<int> rk,height;
vector<vector<int>> st;
void build(const SA &S){
	n=S.n;rk=S.rk;height.assign(n+1,0);st.clear();
	assert((int)S.sa.size()==n+1&&(int)rk.size()==n+1);
	for(int i=1,k=0;i<=n;i++){
		int r=rk[i];
		if(r==1){k=0;continue;}
		int j=S.sa[r-1];
		while(i+k<=n&&j+k<=n&&S.s[i+k-1]==S.s[j+k-1])k++;
		height[r]=k;if(k)k--;
	}
	if(!n)return;
	st.push_back(height);
	for(int k=1;(1<<k)<=n;k++){
		st.emplace_back(n-(1<<k)+2);
		for(int i=1;i+(1<<k)-1<=n;i++)st[k][i]=min(st[k-1][i],st[k-1][i+(1<<(k-1))]);
	}
}
void build(const string &s,int M=128){SA S;S.build(s,M);build(S);}
int lcp(int i,int j)const{
	assert(1<=min(i,j)&&max(i,j)<=n);
	if(i==j)return n-i+1;
	int l=rk[i],r=rk[j];if(l>r)swap(l,r);l++;
	int k=__lg(r-l+1);return min(st[k][l],st[k][r-(1<<k)+1]);
}
int ask(int i,int j)const{return lcp(i,j);}
};
