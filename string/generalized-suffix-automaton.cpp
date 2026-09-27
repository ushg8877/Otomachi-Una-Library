////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/string/generalized-suffix-automaton.cpp
//
// usage: GSAM sam; sam.insert(" aba"); sam.insert(" bab"); sam.build();
//   sam.count(" ab"); sam.distinct();
//
////////////////////////////////////////////////////////////////
struct GSAM{
int n=0,tot=0;
ll strings=0;
bool built=false;
vector<array<int,26>> ch;
vector<int> fail,len;
vector<ll> cnt,siz;
GSAM(){set();}
void set(int m=0){
	assert(0<=m&&m<(INT_MAX-1)/2);
	n=tot=0;
	strings=0;
	built=false;
	ch.assign(1,{});
	fail.assign(1,-1);
	len.assign(1,0);
	cnt.assign(1,0);
	siz.clear();
	ch.reserve(2*m+1);
	fail.reserve(2*m+1);
	len.reserve(2*m+1);
	cnt.reserve(2*m+1);
}
int new_node(){
	assert(tot<INT_MAX-1);
	ch.push_back({});
	fail.push_back(0);
	len.push_back(0);
	cnt.push_back(0);
	return ++tot;
}
int clone(int x,int l){
	int y=new_node();
	ch[y]=ch[x];
	fail[y]=fail[x];
	len[y]=l;
	return y;
}
int extend(int p,int c){
	if(ch[p][c]){
		int q=ch[p][c];if(len[q]==len[p]+1)return q;
		int r=clone(q,len[p]+1);
		while(p!=-1&&ch[p][c]==q)ch[p][c]=r,p=fail[p];
		fail[q]=r;
		return r;
	}
	int cur=new_node();len[cur]=len[p]+1;
	while(p!=-1&&!ch[p][c])ch[p][c]=cur,p=fail[p];
	if(p==-1)fail[cur]=0;
	else{
		int q=ch[p][c];
		if(len[q]==len[p]+1)fail[cur]=q;
		else{
			int r=clone(q,len[p]+1);
			while(p!=-1&&ch[p][c]==q)ch[p][c]=r,p=fail[p];
			fail[q]=fail[cur]=r;
		}
	}
	return cur;
}
// String has a leading space and lowercase letters.
int insert(const string &s){
	assert(!s.empty()&&s[0]==' '&&s.size()<INT_MAX);
	int u=0;
	n=max(n,(int)s.size()-1);
	strings++;
	built=false;
	for(int i=1;i<(int)s.size();i++){
		assert('a'<=s[i]&&s[i]<='z');
		u=extend(u,s[i]-'a');
		cnt[u]++;
	}
	return u;
}
void build(){
	vector<int> c(n+1),q(tot+1);siz=cnt;
	for(int i=0;i<=tot;i++)c[len[i]]++;
	for(int i=1;i<=n;i++)c[i]+=c[i-1];
	for(int i=0;i<=tot;i++)q[--c[len[i]]]=i;
	for(int i=tot;i;i--)siz[fail[q[i]]]+=siz[q[i]];
	siz[0]+=strings;built=true;
}
int find(const string &s)const{
	assert(!s.empty()&&s[0]==' '&&s.size()<INT_MAX);int u=0;
	for(int i=1;i<(int)s.size();i++){
		assert('a'<=s[i]&&s[i]<='z');
		u=ch[u][s[i]-'a'];
		if(!u)return -1;
	}
	return u;
}
// Count occurrences over all inserted strings, including repetitions.
ll count(const string &s)const{
	assert(built);
	int u=find(s);
	return u==-1?0:siz[u];
}
ll distinct()const{
	ll ans=0;
	for(int i=1;i<=tot;i++)ans+=len[i]-len[fail[i]];
	return ans;
}
};
// end for string/generalized-suffix-automaton.cpp
/////////////////////////
// !!!!! Insert every string before build(). !!!!
