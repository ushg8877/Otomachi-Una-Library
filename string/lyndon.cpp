////////////////////////////////////////////////////////////////
//
// template for Lyndon Factorization
// usage:
//   auto a=lyndon(" ababab"); // [1,2], [3,4], [5,6]
//   Leading space, 1-indexed closed intervals, ordered from left to right.
//   The factors are Lyndon words in non-increasing lexicographic order.
//   O(n) time, O(1) extra space excluding the returned intervals.
//
////////////////////////////////////////////////////////////////
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
