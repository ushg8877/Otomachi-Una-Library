////////////////////////////////////////////////////////////////
//
// https://github.com/ushg8877/Otomachi-Una-Library
//
// usage: init(n); is_prime[x]; prime[i]  // i-th prime, 1-index
//
////////////////////////////////////////////////////////////////
vector<char> is_prime;
vector<int> prime;
int m=0;
void init(int n){
	assert(1<=n&&n<INT_MAX);
	is_prime.assign(n+1,true);
	prime.assign(1,0);
	m=0;
	is_prime[0]=is_prime[1]=false;
	for(int i=2;i<=n;i++){
		if(is_prime[i])prime.push_back(i),++m;
		for(int j=1;j<=m&&prime[j]<=n/i;j++){
			int x=i*prime[j];is_prime[x]=false;
			if(i%prime[j]==0)break;
		}
	}
}
// end for number-theory/linear-sieve.cpp
/////////////////////////
