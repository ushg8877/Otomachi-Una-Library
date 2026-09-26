////////////////////////////////////////////////////////////////
//
// template for Cartesian Tree (min heap)
//
// usage: auto [ls,rs]=cartesian(a); // 0-indexed, absent child = -1
//   O(n); equal values keep the earlier index above the later one.
//
////////////////////////////////////////////////////////////////
// O(n); 0-based child indices, -1 if absent. Equal values keep the earlier
// index above the later one; root is the first minimum. Empty input is valid.
template<typename T>
pair<vector<int>,vector<int>> cartesian(const vector<T> &a){
	assert(a.size()<(size_t)INT_MAX);int n=a.size();
	vector<int> ls(n,-1),rs(n,-1),st;st.reserve(n);
	for(int i=0;i<n;i++){
		int last=-1;
		while(!st.empty()&&a[i]<a[st.back()]){last=st.back();st.pop_back();}
		ls[i]=last;
		if(!st.empty())rs[st.back()]=i;
		st.push_back(i);
	}
	return {move(ls),move(rs)};
}
// end for tree/cartesian-tree.cpp
/////////////////////////
