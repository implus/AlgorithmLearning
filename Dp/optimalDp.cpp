/*
key: Dp, unordered_map, n^3 -> n^2
http://codeforces.com/contest/489/problem/F
n x n 01矩阵，每行每列和都是2的方案数目。
可以发现，对于某一行，已知的列和0的个数、1的个数；
就可以推出2的个数，同时可以推出行数。

所以两个参量就可以表示所有的有效状态。
dp[n][n], 记忆化即可;
很容易有n^3的方法，可以发现n^2中的有效状态其实很稀疏，所以可以用一个有效状态数组优化成n的。
*/

#include <tr1/unordered_map>
using namespace std::tr1;

int n;
typedef pair<int, int> pii;
struct hash_func{
	size_t operator()(const pii & a)const{
		return a.first * n + a.second;
	}
};
/*
struct cmp_func{
	bool operator()(const pii &a, const pii &b)const{
		return a.first == b.first && a.second == b.second;
	}
}*/
unordered_map<pii, int, hash_func/*, cmp_func*/> mp;


char a[555][555];
void add(int & x, int a){
	if( (x += a) >= mod ) x -= mod;
}
int gao(int j, int k){
	pii p = pii(j, k);
	if(mp.find(p) != mp.end()) return mp[p];
	if(j > n || k > n || j < 0 || k < 0) return 0;
	int ans = 0;
	add(ans, 1LL * (j + 2) * (j + 1) / 2 * gao(j + 2, k - 2) % mod);
	add(ans, 1LL * (j + 1) * k * gao(j + 1, k) % mod);
	add(ans, 1LL * (k + 2) * (k + 1) / 2 * gao(j, k + 2) % mod);
	return mp[p] = ans;
}

int main() {
	cin>>n>>m>>mod;
	for(int i = 0; i < m; i++){
		cin>>a[i];
	}
	int c[3] = {0};
	for(int j = 0; j < n; j++){
		int cnt = 0;
		for(int i = 0; i < m; i++)
			cnt += a[i][j] == '1';
		if(cnt > 2) {
			puts("0"); return 0;
		}
		c[cnt]++;
	}
	mp[pii(c[0], c[1])] = 1;
	cout<<gao(0, 0)<<endl;
	return 0;
}
