struct unionfind{
int n=0;
vector<int> fa,siz,h;
void setN(int m){
	assert(0<=m&&m<INT_MAX);n=m;h.clear();
	fa.resize(n+1);iota(fa.begin(),fa.end(),0);
	siz.assign(n+1,1);siz[0]=0;
}
int find(int x)const{
	assert(1<=x&&x<=n);
	while(x!=fa[x])x=fa[x];
	return x;
}
int merge(int x,int y){
	x=find(x);y=find(y);if(x==y)return 0;
	if(siz[x]<siz[y])swap(x,y);
	h.push_back(y);fa[y]=x;siz[x]+=siz[y];return x;
}
int version()const{return h.size();}
void roll_back(int t){
	assert(0<=t&&t<=(int)h.size());
	while((int)h.size()>t){
		int y=h.back();h.pop_back();siz[fa[y]]-=siz[y];fa[y]=y;
	}
}
};
