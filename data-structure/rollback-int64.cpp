////////////////////////////////////////////////////////////////
//
// template for roll back saver
//
// usage:
//   ::set(); save(x); save(y); roll_back(t);
//   supports int and ll; do not paste with (int) saver
//   saved elements must stay at the same address until rollback
//
////////////////////////////////////////////////////////////////
vector<pair<int*,int>> buf;
vector<pair<ll*,ll>> buf1;
int saved_ele=0;
vector<char> R=vector<char>(1);
void set(){saved_ele=0;buf.clear();buf1.clear();R.assign(1,false);}
void save(int &x){buf.push_back(make_pair(&x,x));R.push_back(false);++saved_ele;}
void save(ll &x){buf1.push_back(make_pair(&x,x));R.push_back(true);++saved_ele;}
int version(){return saved_ele;}
void roll_back(int t){
	assert(0<=t&&t<=version());
	while(saved_ele>t){
		if(R[saved_ele]){
			auto [a,b]=buf1.back();buf1.pop_back();
			*a=b;
		}else{
			auto [a,b]=buf.back();buf.pop_back();
			*a=b;
		}
		saved_ele--;R.pop_back();
	}
}
// end for data-structure/rollback-int64.cpp
/////////////////////////
// !!!!! Do not paste both rollback variants; saved elements must keep their
// addresses until roll_back(). !!!!
