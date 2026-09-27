////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/math/linear-basis/xor.cpp
//
// usage: Xor_Basis B; B.setM(m); B.insert(x); B.contains(x); B.ask(x);
//   m bits, 0<=m<=64; default 64. setM clears; set keeps m.
//
////////////////////////////////////////////////////////////////
struct Xor_Basis{
int m=64,rank=0;
vector<ull> a=vector<ull>(64);
// Set bit width 0..64 and clear; default 64. set() retains the width.
void setM(int _m){assert(0<=_m&&_m<=64);m=_m;rank=0;a.assign(m,0);}
void set(){rank=0;fill(a.begin(),a.end(),0);}
bool insert(ull x){
	assert(m==64||!(x>>m));
	for(int i=m-1;i>=0;i--)if(x>>i&1){
		if(a[i])x^=a[i];
		else{a[i]=x;rank++;return true;}
	}
	return false;
}
bool contains(ull x)const{
	for(int i=m-1;i>=0;i--)if(x>>i&1)x^=a[i];
	return x==0;
}
// Maximize x xor span; bits above m in x remain unchanged.
ull ask(ull x=0)const{
	for(int i=m-1;i>=0;i--)x=max(x,x^a[i]);
	return x;
}
};
// end for math/linear-basis/xor.cpp
/////////////////////////
// !!!!! Inserted values must fit m bits. !!!!
