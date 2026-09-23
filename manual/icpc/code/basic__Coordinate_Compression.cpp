template<typename T>
struct compress{
vector<T> I;
vector<T*> II;
// make compress values to 1 ~ size
void set(){I.clear();II.clear();}
void pb(T &x){I.push_back(x);II.push_back(&x);}
int id(T x)const{
	auto it=lower_bound(I.begin(),I.end(),x);
	assert(it!=I.end()&&*it==x);
	return it-I.begin()+1;
}
void build(){
	sort(I.begin(),I.end());
	I.erase(unique(I.begin(),I.end()),I.end());
	for(auto i:II) *i=id(*i);
	II.clear();
}
int size()const{return I.size();}
array<int,2> range(int l,int r){
	l=lower_bound(I.begin(),I.end(),l)-I.begin();
	r=upper_bound(I.begin(),I.end(),r)-I.begin()-1;
	return {l+1,r+1};
}
T operator [](int x)const{
	assert(1<=x&&x<=I.size());
	return I[x-1];
}
};
// !!!!!!!!!! before each use, use cc.set() to clear !!!!!!!
// !!!!!!!!!! 如果元素被插入，cc.build() 会自动重编号 !!!!!!!!!!!!!
