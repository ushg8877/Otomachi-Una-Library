////////////////////////////////////////////////////////////////
//
// template for xor prefix linear basis
//
// usage:
// Prefix_Xor_Basis B; int removed=B.insert(x,id);
// id must increase; 0: rank increased, otherwise removed id
//
////////////////////////////////////////////////////////////////
struct Prefix_Xor_Basis{
int rank=0,last=0;
vector<ull> a=vector<ull>(64);
vector<int> pos=vector<int>(64);
void set(){rank=last=0;fill(a.begin(),a.end(),0);fill(pos.begin(),pos.end(),0);}
int insert(ull x,int id){
	assert(id>last);last=id;
	for(int i=63;i>=0;i--)if(x>>i&1){
		if(!a[i]){a[i]=x;pos[i]=id;rank++;return 0;}
		if(pos[i]<id){swap(a[i],x);swap(pos[i],id);}
		x^=a[i];
	}
	return id;
}
int get_rank(int l=1)const{
	assert(l>=1);int ans=0;
	for(int i=0;i<64;i++)ans+=pos[i]>=l;
	return ans;
}
bool contains(ull x,int l=1)const{
	assert(l>=1);
	for(int i=63;i>=0;i--)if((x>>i&1)&&pos[i]>=l)x^=a[i];
	return x==0;
}
ull ask(ull x=0,int l=1)const{
	assert(l>=1);
	for(int i=63;i>=0;i--)if(pos[i]>=l)x=max(x,x^a[i]);
	return x;
}
};
// end for math/linear-basis/xor-prefix.cpp
/////////////////////////
// !!!!! Ids must be positive and strictly increasing; a dependent insertion may remove its own id. !!!!
