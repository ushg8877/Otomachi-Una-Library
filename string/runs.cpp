////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/string/runs.cpp
//
// usage: Runs t; auto a=t.build(" ababab"); // {1,6,2}
//
////////////////////////////////////////////////////////////////
struct Runs{
public:
// Leading space, bytes 1..255. Return sorted unique {l,r,p}:
// [l,r] is maximal with least period p, r-l+1>=2*p; 1-based.
// O(n log n) time and O(n) space.
vector<array<int,3>> build(const string &s)const{
	assert(!s.empty()&&s[0]==' '&&s.size()<(INT_MAX-1024)/2);
	int n=(int)s.size()-1;
	for(int i=1;i<=n;i++)assert((unsigned char)s[i]>0);
	vector<array<int,3>> ans;if(n<2)return ans;
	SA_LCP a;
	SA_LCS b;
	a.build(s,256);
	b.build(s,256);
	vector<int> st;st.reserve(n);
	for(int op=0;op<2;op++){
		st.clear();
		for(int i=n;i>=1;i--){
			while(!st.empty()){
				int j=st.back(),k=a.ask(i,j);
				// Reversing the alphabet does not reverse the prefix rule.
				if(j+k>n)break;
				bool less=(unsigned char)s[i+k]<(unsigned char)s[j+k];
				if(less==bool(op))break;
				st.pop_back();
			}
			if(!st.empty()){
				int j=st.back(),p=j-i;
				int l=b.ask(i-1,j-1),r=a.ask(i,j);
				if(l+r>=p)ans.push_back({i-l,j+r-1,p});
			}
			st.push_back(i);
		}
	}
	sort(ans.begin(),ans.end());
	ans.erase(unique(ans.begin(),ans.end()),ans.end());
	return ans;
}
};
// end for string/runs.cpp
/////////////////////////
// !!!!! Paste data-structure/linear-rmq.cpp, then string/suffix-array.cpp. !!!!
