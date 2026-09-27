////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library/blob/main/data-structure/rollback-int.cpp
//
// usage: int t=version(); save(x); x=1; roll_back(t);
//
////////////////////////////////////////////////////////////////
vector<pair<int*,int>> buf;
void set(){buf.clear();}
void save(int &x){buf.push_back(make_pair(&x,x));}
int version(){return buf.size();}
void roll_back(int t){
	assert(0<=t&&t<=version());
	while(buf.size()>t){
		auto [a,b]=buf.back();buf.pop_back();
		*a=b;
	}
}
// end for data-structure/rollback-int.cpp
/////////////////////////
// !!!!! Saved elements must keep their addresses until roll_back(). !!!!
