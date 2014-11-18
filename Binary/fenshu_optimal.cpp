/*
key: BinarySearch, Dp

http://codeforces.com/contest/489/problem/E
从n个点中选取一些点最优化
f = sigma(sqrt(x[i] - x[j] - l)) / sigma(b[k])
i, j 为相邻选择点, n 必须选。

简化问题：
如果去掉分母，dp很容易解决；
如何去掉分母？二分答案。
*/

const int MAX_N = 1010;
const double eps = 1e-8;
double dp[MAX_N];
int n, l;
int x[MAX_N], b[MAX_N], nxt[MAX_N];
vector<int> ans;
double gao(double v){
	ans.clear();
	fill(dp, dp + n + 1, 1e99);
	fill(nxt, nxt + n + 1, 0);
	dp[0] = 0;
	for(int i = 1; i <= n; i++){
		for(int j = 0; j < i; j++){
			double ndp = dp[j] + sqrt(1.0 * abs(x[i] - x[j] - l)) - v * b[i];
			if(ndp < dp[i]){
				dp[i] = ndp;
				nxt[i] = j;
			}
		}
	}
	int p = n;
	while(p){
		ans.push_back(p); p = nxt[p];
	}
	reverse(ans.begin(), ans.end());
	return dp[n];
}

int main() {
	cin>>n>>l;
	for(int i = 1; i <= n; i++){
		cin>>x[i]>>b[i];
	}

	double l = 0, r = 1e10, m;
	while(r - l > eps){
		m = (l + r) / 2.0;
		if(gao(m) < 0) r = m;
		else l = m;
	}
	gao(m);
	cerr<<"m = "<<m<<endl;
	for(int i =0 ; i< ans.size(); i++)
		cout<<ans[i]<<" ";
	return 0;
}
