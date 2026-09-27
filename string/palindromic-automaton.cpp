////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: PAM pam; pam.build(" ababa"); auto cnt=pam.count();
//   Or pam.set(n), then pam.extend(c); returns the longest suffix node.
//
////////////////////////////////////////////////////////////////
struct PAM{
int n=0,tot=1,lst=0;
string s;
vector<array<int,26>> ch;
vector<int> len,fail,num,pos,siz;
PAM(){set();}
void set(int m=0){
	assert(0<=m&&m<INT_MAX-2);
	n=0;
	tot=1;
	lst=0;
	s=" ";
	ch.assign(2,{});
	len={0,-1};
	fail={1,1};
	num.assign(2,0);
	pos.assign(2,0);
	siz.assign(2,0);
	s.reserve(m+1);
	ch.reserve(m+2);
	len.reserve(m+2);
	fail.reserve(m+2);
	num.reserve(m+2);
	pos.reserve(m+2);
	siz.reserve(m+2);
}
int get_fail(int u)const{
	while(s[n-len[u]-1]!=s[n])u=fail[u];
	return u;
}
// Append one lowercase letter; return the longest palindromic suffix node.
int extend(char c){
	assert('a'<=c&&c<='z'&&n<INT_MAX-3);
	s+=c;
	n++;
	c-='a';
	int u=get_fail(lst);
	if(!ch[u][c]){
		int v=++tot;
		ch.push_back({});len.push_back(len[u]+2);
		fail.push_back(len[v]==1?0:ch[get_fail(fail[u])][c]);
		num.push_back(num[fail[v]]+1);
		pos.push_back(n);
		siz.push_back(0);
		ch[u][c]=v;
	}
	lst=ch[u][c];
	siz[lst]++;
	return lst;
}
// Leading space and lowercase letters; O(n) time / space.
// Nodes 0,1 are virtual; real nodes 2..tot. len/fail/num/pos store
// length, suffix link, suffix count and first occurrence end.
void build(string S){
	assert(!S.empty()&&S[0]==' '&&S.size()<INT_MAX-2);
	set((int)S.size()-1);
	for(int i=1;i<(int)S.size();i++)extend(S[i]);
}
// Occurrence counts; repeatable, O(n).
vector<int> count()const{
	vector<int> cnt=siz;
	for(int u=tot;u>=2;u--)cnt[fail[u]]+=cnt[u];
	cnt[0]=cnt[1]=0;
	return cnt;
}
};
// end for string/palindromic-automaton.cpp
/////////////////////////
