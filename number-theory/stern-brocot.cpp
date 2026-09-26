////////////////////////////////////////////////////////////////
//
// template for Stern-Brocot search
// usage:
//   SBT t; t.setM(M);
//   while(!t.done()){Rat x=t.ask();t.tell(cmp(x));}
//   cmp: -1 if x is too small, 0 if equal, +1 if too large.
//   Search 0<=p<=M, 1<=q<=M; feedback must be monotone.
//   found(): l==r is the answer; otherwise l,r bracket the target.
//   l,r are {p,q}; {-1,0} / {1,0} mean no lower / upper candidate.
//   O(log(M+1)) queries / local time, O(1) space; setM clears the search.
//
////////////////////////////////////////////////////////////////
struct SBT{
pair<ll,ll> l,r;
private:
ll M=0,lo=0,hi=0,k=0,cap=0;
int phase=-1,dir=0;
bool stop=true;
Rat cur;
pair<ll,ll> at(ll v)const{
	auto a=dir<0?l:r,b=dir<0?r:l;
	return {a.first+v*b.first,a.second+v*b.second};
}
void next(){
	// The bounds have determinant 1, so every query is already reduced.
	if(phase==0){
		if(l.first>M-r.first||l.second>M-r.second){stop=true;return;}
		cur.p=l.first+r.first;cur.q=l.second+r.second;
	}else{auto [p,q]=at(k);cur.p=p;cur.q=q;}
}
void split(){
	if(hi-lo>1){k=lo+(hi-lo)/2;next();return;}
	auto a=at(lo),b=at(hi);
	if(dir<0)l=a,r=b;else r=a,l=b;
	phase=0;next();
}
public:
void setM(ll m){
	assert(m>=1);M=m;l={0,1};r={1,0};phase=-1;stop=false;cur=Rat();
}
bool done()const{return stop;}
bool found()const{return stop&&M&&l==r;}
Rat ask()const{assert(!stop);return cur;}
void tell(int res){
	assert(!stop&&-1<=res&&res<=1);
	if(!res){l=r={cur.p,cur.q};stop=true;return;}
	if(phase==-1){
		if(res>0){l={-1,0};r={0,1};stop=true;return;}
		phase=0;next();return;
	}
	if(phase==0){
		dir=res;auto a=dir<0?l:r,b=dir<0?r:l;cap=M;
		if(b.first)cap=min(cap,(M-a.first)/b.first);
		if(b.second)cap=min(cap,(M-a.second)/b.second);
		lo=k=1;phase=1;
	}else if(res==dir)lo=k;
	else{hi=k;phase=2;}
	if(phase==2){split();return;}
	if(lo==cap){if(dir<0)l=at(lo);else r=at(lo);stop=true;return;}
	k=lo>cap/2?cap:lo*2;next();
}
};
// end for number-theory/stern-brocot.cpp
/////////////////////////
// !!!!! Paste math/rational.cpp first. tell(-1/0/1) describes the last ask()
// result, not the target; l/r with q=0 are sentinels. !!!!
