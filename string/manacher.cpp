////////////////////////////////////////////////////////////////
//
// template for Manacher
// usage:
//   Manacher t; t.build(" abba"); t.ask(1,4);
//   Leading space, 1-indexed; ask(l,r) checks a nonempty closed interval.
//   Odd: [i-d1[i]+1,i+d1[i]-1]; even: [i-d2[i],i+d2[i]-1].
//   O(n) build / space, O(1) query.
//
////////////////////////////////////////////////////////////////
struct Manacher{
int n=0;
vector<int> d1,d2;
void set(){n=0;d1.clear();d2.clear();}
void build(const string &s){
	assert(!s.empty()&&s[0]==' '&&s.size()<INT_MAX/2);
	n=(int)s.size()-1;d1.assign(n+1,0);d2.assign(n+1,0);
	for(int i=1,l=1,r=0;i<=n;i++){
		int k=i>r?1:min(d1[l+r-i],r-i+1);
		while(i-k>=1&&i+k<=n&&s[i-k]==s[i+k])k++;
		d1[i]=k--;
		if(i+k>r)l=i-k,r=i+k;
	}
	for(int i=1,l=1,r=0;i<=n;i++){
		int k=i>r?0:min(d2[l+r-i+1],r-i+1);
		while(i-k-1>=1&&i+k<=n&&s[i-k-1]==s[i+k])k++;
		d2[i]=k;
		if(i+k-1>r)l=i-k,r=i+k-1;
	}
}
bool ask(int l,int r)const{
	assert(1<=l&&l<=r&&r<=n);int k=r-l+1;
	return k&1?d1[(l+r)/2]>=k/2+1:d2[(l+r+1)/2]>=k/2;
}
};
// end for string/manacher.cpp
/////////////////////////
