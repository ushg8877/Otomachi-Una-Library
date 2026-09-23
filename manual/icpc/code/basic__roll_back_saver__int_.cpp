vector<pair<int*,int>> buf;
void saver_clear(){buf.clear();}
void save(int &x){buf.push_back(make_pair(&x,x));}
int version(){return buf.size();}
void roll_back(int t){
	assert(0<=t&&t<=version());
	while(buf.size()>t){
		auto [a,b]=buf.back();buf.pop_back();
		*a=b;
	}
}
