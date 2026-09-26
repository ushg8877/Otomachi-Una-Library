////////////////////////////////////////////////////////////////
// string algorithm template by Otomachi Una (Junlin Ye)
// version 1.0 (last update May 24th 2026)
//
// usage: SAM sam; sam.build(" abaab");  // leading space, lowercase
//   sam.len[u], sam.fail[u], sam.siz[u]
//
////////////////////////////////////////////////////////////////
struct SAM{
// this SAM rooted at 0
string s;
int n=0,tot=0,lst=0;
vector<array<int,26>> ch;
vector<int> fail,len,siz;
vector<vector<int>> son;
SAM(){set();}
void set(int m=0){
	assert(0<=m&&m<(INT_MAX-1)/2);n=tot=lst=0;
	ch.assign(1,{});
	fail.assign(1,-1);
	len.assign(1,0);
	siz.assign(1,0);
	son.assign(1,{});
	ch.reserve(2*m+1);
	fail.reserve(2*m+1);
	len.reserve(2*m+1);
	siz.reserve(2*m+1);son.reserve(2*m+1);
}
int new_node(){
	assert(tot<INT_MAX-1);
	ch.push_back({});
	fail.push_back(0);
	len.push_back(0);
	siz.push_back(0);
	son.emplace_back();
	return ++tot;
}
inline int clone(int x){
	new_node();ch[tot]=ch[x];
	fail[tot]=fail[x];
	return tot;
}
void extend(char c){
	assert('a'<=c&&c<='z');
	c-='a';
	// attention the alphabet
	int cur=new_node();len[cur]=len[lst]+1;
	while(~lst&&!ch[lst][c]) ch[lst][c]=tot,lst=fail[lst];
	if(lst==-1) fail[cur]=0;
	else{
		int p=ch[lst][c];
		if(len[p]==len[lst]+1) fail[cur]=p;
		else{
			int q=clone(p);len[q]=len[lst]+1;
			while(~lst&&ch[lst][c]==p) ch[lst][c]=q,lst=fail[lst];
			fail[p]=fail[cur]=q;
		}
	}
	lst=cur;
}
void output(){
	// this allow you to see the structure of SAM
	// output fail, ed, len etc.
	cerr<<"string:"<<s<<endl;
	cerr<<"# node = "<<tot<<endl;
	cerr<<"fail: ";
	for(int i=0;i<=tot;i++) cerr<<fail[i]<<" \n"[i==tot];
	cerr<<"len: ";
	for(int i=0;i<=tot;i++) cerr<<len[i]<<" \n"[i==tot];
	// other else?
	cerr<<"siz: ";
	for(int i=0;i<=tot;i++) cerr<<siz[i]<<" \n"[i==tot];
}
void dfs(int u){
	vector<int> q{u};
	for(int i=0;i<(int)q.size();i++)for(int v:son[q[i]])q.push_back(v);
	for(int i=(int)q.size()-1;i>0;i--)siz[fail[q[i]]]+=siz[q[i]];
}
void build(string S){
	assert(!S.empty()&&S[0]==' '&&S.size()<(INT_MAX-1)/2);
	set((int)S.size()-1);
	s=move(S);
	n=(int)s.size()-1;
	for(int i=1;i<=n;i++){
		extend(s[i]);
		// something record here
		siz[lst]=1;
	}
	for(int i=1;i<=tot;i++) son[fail[i]].push_back(i);
	dfs(0);
	// output();
}
void solve(){
	// this is the main work
	ll ans=0;
	for(int i=1;i<=tot;i++) if(siz[i]>1) ans=max(ans,1ll*siz[i]*len[i]);
	cout<<ans<<'\n';
}
};
// end for string/suffix-automaton.cpp
/////////////////////////
