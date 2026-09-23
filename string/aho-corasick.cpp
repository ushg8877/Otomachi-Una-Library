////////////////////////////////////////////////////////////////
//
// template for Aho-Corasick
// usage:
//   AC ac; int id=ac.insert(" aba"); ac.build();
//   auto cnt=ac.count(" ababa"); cnt[id];
//   leading space, lowercase; insert before build
//
////////////////////////////////////////////////////////////////
struct AC{
int tot=0;
bool built=false;
vector<array<int,26>> ch;
vector<int> fail,q;
AC(){set();}
void set(int m=0){
	assert(0<=m&&m<INT_MAX-1);tot=0;built=false;
	ch.assign(1,{});fail.assign(1,0);q.clear();
	ch.reserve(m+1);fail.reserve(m+1);
}
int new_node(){
	assert(tot<INT_MAX-1);ch.push_back({});fail.push_back(0);return ++tot;
}
int insert(const string &s){
	assert(!built&&!s.empty()&&s[0]==' '&&s.size()<INT_MAX);
	int u=0;
	for(int i=1;i<(int)s.size();i++){
		assert('a'<=s[i]&&s[i]<='z');int c=s[i]-'a';
		if(!ch[u][c]){int v=new_node();ch[u][c]=v;}
		u=ch[u][c];
	}
	return u;
}
void build(){
	if(built)return;
	q.clear();q.reserve(tot+1);q.push_back(0);
	for(int i=0;i<(int)q.size();i++){
		int u=q[i];
		for(int c=0;c<26;c++){
			int v=ch[u][c];
			if(v){fail[v]=u?ch[fail[u]][c]:0;q.push_back(v);}
			else if(u)ch[u][c]=ch[fail[u]][c];
		}
	}
	built=true;
}
vector<ll> count(const string &s)const{
	assert(built&&!s.empty()&&s[0]==' '&&s.size()<INT_MAX);
	vector<ll> cnt(tot+1);cnt[0]=1;int u=0;
	for(int i=1;i<(int)s.size();i++){
		assert('a'<=s[i]&&s[i]<='z');u=ch[u][s[i]-'a'];cnt[u]++;
	}
	for(int i=(int)q.size()-1;i;i--)cnt[fail[q[i]]]+=cnt[q[i]];
	return cnt;
}
};
// end for string/aho-corasick.cpp
/////////////////////////
// !!!!! Leading space and lowercase letters; insert every pattern before build(). !!!!
