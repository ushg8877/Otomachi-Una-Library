////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/math/linear-basis/xor-prefix.cpp
//
// usage: Prefix_Xor_Basis B; B.setM(m); int removed=B.insert(x,id);
//   m bits, 0<=m<=64; default 64. setM clears; set keeps m.
//
////////////////////////////////////////////////////////////////
struct Prefix_Xor_Basis{
int m=64,rank=0,last=0;
vector<ull> a=vector<ull>(64);
vector<int> pos=vector<int>(64);
// Set bit width 0..64 and clear; default 64. set() retains the width.
void setM(int _m){
	assert(0<=_m&&_m<=64);
	m=_m;
	rank=last=0;
	a.assign(m,0);
	pos.assign(m,0);
}
void set(){rank=last=0;fill(a.begin(),a.end(),0);fill(pos.begin(),pos.end(),0);}
// Ids must be positive and strictly increasing. Return 0 if rank grows;
// otherwise return the removed id, possibly the new element itself.
int insert(ull x,int id){
	assert(m==64||!(x>>m));
	assert(id>last);last=id;
	for(int i=m-1;i>=0;i--)if(x>>i&1){
		if(!a[i]){a[i]=x;pos[i]=id;rank++;return 0;}
		if(pos[i]<id){swap(a[i],x);swap(pos[i],id);}
		x^=a[i];
	}
	return id;
}
int get_rank(int l=1)const{
	assert(l>=1);int ans=0;
	for(int i=0;i<m;i++)ans+=pos[i]>=l;
	return ans;
}
bool contains(ull x,int l=1)const{
	assert(l>=1);
	for(int i=m-1;i>=0;i--)if((x>>i&1)&&pos[i]>=l)x^=a[i];
	return x==0;
}
ull ask(ull x=0,int l=1)const{
	assert(l>=1);
	for(int i=m-1;i>=0;i--)if(pos[i]>=l)x=max(x,x^a[i]);
	return x;
}
};
// end for math/linear-basis/xor-prefix.cpp
/////////////////////////
// !!!!! Inserted values must fit m bits; ids must be positive and
// increasing. !!!!
