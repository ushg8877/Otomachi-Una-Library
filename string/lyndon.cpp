////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
// blob/main/string/lyndon.cpp
//
// usage: auto a=lyndon(" ababab"); // [1,2], [3,4], [5,6]
//   Leading space, 1-indexed closed intervals, ordered from left to right.
//
////////////////////////////////////////////////////////////////
// Leading space; return 1-based closed intervals from left to right.
// Factors are non-increasing lexicographically; O(n) time and O(1) extra space.
vector<pair<int,int>> lyndon(const string &s){
	assert(!s.empty()&&s[0]==' '&&s.size()<INT_MAX);
	int n=(int)s.size()-1;vector<pair<int,int>> ans;
	for(int i=1;i<=n;){
		int j=i+1,k=i;
		while(j<=n&&(unsigned char)s[k]<=(unsigned char)s[j]){
			k=s[k]==s[j]?k+1:i;j++;
		}
		int p=j-k;
		while(i<=k){ans.emplace_back(i,i+p-1);i+=p;}
	}
	return ans;
}
// end for string/lyndon.cpp
/////////////////////////
