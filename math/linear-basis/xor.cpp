////////////////////////////////////////////////////////////////
//
// template for xor linear basis
//
// usage:
// Xor_Basis B; B.insert(x); B.contains(x); B.ask(x);
//
////////////////////////////////////////////////////////////////
struct Xor_Basis{
int rank=0;
vector<ull> a=vector<ull>(64);
void set(){rank=0;fill(a.begin(),a.end(),0);}
bool insert(ull x){
	for(int i=63;i>=0;i--)if(x>>i&1){
		if(a[i])x^=a[i];
		else{a[i]=x;rank++;return true;}
	}
	return false;
}
bool contains(ull x)const{
	for(int i=63;i>=0;i--)if(x>>i&1)x^=a[i];
	return x==0;
}
ull ask(ull x=0)const{
	for(int i=63;i>=0;i--)x=max(x,x^a[i]);
	return x;
}
};
// end for math/linear-basis/xor.cpp
/////////////////////////
