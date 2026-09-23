struct seg{
	int l,r,c;
	seg():l(0),r(0),c(0){}
	seg(int _l,int _r,int _c):l(_l),r(_r),c(_c){}
	inline bool operator < (const seg &x) const {return l<x.l;}
};
struct ODT{
private:
set<seg> S;
void cut(int y){
	auto it=S.upper_bound(seg(y,0,0));
	if(it==S.begin()) return;
	auto [l,r,c]=*(--it);
	if(l<y&&y<=r) S.erase(it),S.insert(seg(l,y-1,c)),S.insert(seg(y,r,c));
}
public:
void init(){S.clear();}
vector<seg> erase(int l,int r){
	// erase segments from [l,r]
	assert(l<=r);
	cut(l);if(r<INT_MAX)cut(r+1);
	auto it=S.lower_bound(seg(l,0,0));
	vector<seg> I;
	while(1){
		if(it==S.end()||it->l>r) break;
		I.push_back(*it);
		it=S.erase(it);
	}
	return I;
}
vector<seg> extract(int l,int r){
	// extract segments from [l,r]
	// before you use this function, pay attention is the TC correct?
	assert(l<=r);
	cut(l);if(r<INT_MAX)cut(r+1);
	auto it=S.lower_bound(seg(l,0,0));
	vector<seg> I;
	while(1){
		if(it==S.end()||it->l>r) break;
		I.push_back(*it);
		it++;
	}
	return I;
}
void add(int l,int r,int c){
	assert(l<=r);
	S.insert(seg(l,r,c));
}
};
// !!!!!!!!!! before add segments, erase the range first!!!!!!!!!!
