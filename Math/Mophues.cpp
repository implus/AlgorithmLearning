/*
mophues base:

sigma[d|n] (u(d)) = [n == 1]

sigma[d|n] (u(d) * n / d) = phi(n) 

sigma[d|n] (phi(d)) = n



0.	http://acm.hdu.edu.cn/showproblem.php?pid=1695
sigma sigma (gcd(a, b) == 1)


1.	http://acm.hdu.edu.cn/showproblem.php?pid=4746
expand use of sigma sigma (gcd(a, b) == 1) 
sigma sigma ( gcd(a, b) == k ) || k must satisfy decomposed by prime numbers' amount <= p

ans = sigma u(d) * a/(k*d) * b/(k*d)
for each i == k * d, we precalculate sum of u(i / k), k meets some conditions

2.	http://acm.hust.edu.cn/vjudge/problem/viewproblem.action?id=37193
solve sigma^3 gcd(a, b, c) == 1

ans = 3 + solvetriple(a, b, c) + solvedouble(a, b) + solvedouble(a, c) + solvedouble(b, c)
	  3 means (0, 0, 1), (0, 1, 0) (1, 0, 0)

3.	http://acm.hust.edu.cn/vjudge/problem/viewproblem.action?id=10581

4.	http://acm.seu.edu.cn/oj/problem.php?id=123
find xishu needs a multi-function calculate, very good problem
when induce, you should be very careful, for example i = k*i', do not ignore k in every place.

*/


//0.	http://acm.hdu.edu.cn/showproblem.php?pid=1695
/*
#include <stdio.h> 
#include <string.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <set>
#include <map>
using namespace std;
typedef long long lint;
const int max_n = 1e5 + 10;

int u[max_n], prime[max_n], tot, notp[max_n];
void init(){
	u[1] = 1, tot = 0;
	for(int i = 2; i < max_n; i++){
		if(!notp[i]){
			u[i] = -1; prime[tot++] = i;
		}
		for(int j = 0; j < tot; j++){
			if(i * prime[j] > max_n) break;
			notp[i * prime[j]] = 1;
			if(i % prime[j] == 0) {
				u[i * prime[j]] = 0;
				break;
			}else {
				u[i * prime[j]] = -u[i];
			}
		}
	}
	for(int i = 1; i < max_n; i++) u[i] += u[i - 1];
}

int main(){
	init();
	int t, a, b, c, d, k, icase = 1;
	scanf("%d", &t);
	while(t--){
		scanf("%d%d%d%d%d", &a, &b, &c, &d, &k);
		lint ans = 0;
		if(k == 0) {} else{
			b /= k, d /= k;
			if( b > d ) swap(b, d);

			for(int i = 1, j = 1; i <= b; i = j + 1){
				j = min(b / (b / i) , d / (d / i));
				lint x = b / i, y = d / i;
				ans += 1ll * (u[j] - u[i - 1]) * (x * y - x * (x - 1) / 2);
			}
		}
		printf("case %d: %lld\n", icase++, ans);
	}
	return 0;
}
*/


//1.	http://acm.hdu.edu.cn/showproblem.php?pid=4746
/*
#include <stdio.h> 
#include <string.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <set>
#include <map>
using namespace std;
typedef long long lint;

const int max_n = 5e5 + 10;
int u[max_n], prime[max_n], tot, cnt[max_n], notp[max_n];
lint dp[max_n][20];
void init(){
	u[1] = 1; tot = 0; cnt[1] = 0;
	for(int i = 2; i < max_n; i++){
		if(!notp[i]){
			prime[tot++] = i;
			cnt[i] = 1;
			u[i] = -1;
		}
		for(int j = 0; j < tot; j++){
			if(i * prime[j] > max_n) break;
			notp[i * prime[j]] = 1;
			cnt[i * prime[j]] = cnt[i] + 1; 
			if(i % prime[j] == 0){
				u[i * prime[j]] = 0;
				break;
			}else{
				u[i * prime[j]] = -u[i];
			}
		}
	}

	memset(dp, 0, sizeof(dp));
	for(int i = 1; i < max_n; i++){
		for(int j = i; j < max_n; j += i){
			dp[j][cnt[i]] += u[j / i];
		}
	}

	for(int i = 0; i < max_n; i++){
		for(int j = 1; j < 19; j++){
			dp[i][j] += dp[i][j - 1];
		}
	}

	for(int i = 1; i < max_n; i++){
		for(int j = 0; j < 19; j++){
			dp[i][j] += dp[i - 1][j];
		}
	}
}

int main() {
	init();
	int t, m, n, p;
	while(scanf("%d", &t) == 1){
		while(t--){
			scanf("%d%d%d", &n, &m, &p);
			lint ans = 0;
			if(n > m) swap(n, m);
			if(p >= 19){
				ans = 1ll * n * m;
			}else
			for(int i = 1, j = 1; i <= n; i = j + 1){
				j = min(n / (n / i), m / (m / i));
				ans += (dp[j][p] - dp[i - 1][p]) * (n / i) * (m / i);
			}
			printf("%lld\n", ans);
		}
	}
	return 0;
}
*/

//2.	http://acm.hust.edu.cn/vjudge/problem/viewproblem.action?id=37193
/*
#include <stdio.h> 
#include <string.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <set>
#include <map>
using namespace std;
typedef long long lint;

const int max_n = 1e6 + 10;
int prime[max_n], tot, notp[max_n], u[max_n];

void init(){
	u[1] = 1, tot = 0;
	for(int i = 2; i < max_n; i++){
		if(!notp[i]) prime[tot++] = i, u[i] = -1;
		for(int j = 0; j < tot; j++){
			int v = i * prime[j];
			if(v >= max_n) break;
			notp[v] = 1;
			if(i % prime[j] == 0){
				u[v] = 0;
				break;
			}else u[v] = -u[i];
		}
	}
	for(int i = 1; i < max_n; i++) u[i] += u[i - 1];
}

lint solvetriple(int a, int b, int c){
	lint ans = 0;
	for(int i = 1, j = 1; i <= a; i = j + 1){
		j = (a / (a/i));
		ans += 1ll* (u[j] - u[i - 1]) * (a/i) * (b/i) * (c/i);
	}
	return ans;
}

lint solvedouble(int a, int b){
	lint ans = 0;
	for(int i = 1, j = 1; i <= a; i = j + 1){
		j = a / (a/i);
		ans += 1ll * (u[j] - u[i - 1]) * (a/i) * (b/i);
	}
	return ans;
}

int main(){
	init();
	int t, n;
	scanf("%d", &t);
	while(t--){
		scanf("%d", &n);
		lint ans = 3 + solvetriple(n, n, n) + 3 * solvedouble(n, n);
		printf("%lld\n", ans);
	}
	return 0;
}
*/



//3.	http://acm.hust.edu.cn/vjudge/problem/viewproblem.action?id=10581
/*
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <string>
#include <set>
#include <map>
using namespace std;
typedef long long lint;

const int max_n = 1e7 + 5;
int prime[max_n], tot, notp[max_n],u[max_n];
int dp[max_n];
// dp[i] means sum of u[d] d == i / k, k is a prime
// use Eular Shai method

void init(){
	u[1] = 1, tot = 0;
	for(int i = 2; i < max_n; i++){
		if(!notp[i]) prime[tot++] = i, u[i] = -1, dp[i] = u[1];
		for(int j = 0; j < tot; j++){
			int v = i * prime[j];
			if(v >= max_n) break;
			notp[v] = 1;
			if(i % prime[j] == 0){
				u[v] = 0;
				dp[v] = u[i]; // i / k (k != prime[j]) then i must have two prime[j], then all u = 0, except k = prime[j], that's just u[i];
				break;
			}else{
				u[v] = -u[i];
				dp[v] = -dp[i] + u[i]; // seems like add 1 factor, in dp[i], all u[d]'s d must add a new factor, they become -u[d], so dp[i] become -dp[i], also consider erase prime[j], we get u[i];
			}
		}
	}
	for(int i = 1; i < max_n; i++){ dp[i] += dp[i - 1]; }
}

lint gao(int a, int b){
	lint ans = 0;
	if( a > b ) swap(a, b);
	for(int i = 1, j = 1; i <= a; i = j + 1){
		j = min(a /(a/i), b / (b/i));
		ans += 1ll * (dp[j] - dp[i - 1]) * (a/i) * (b/i);
	}
	return ans;
}

int main() {
	init();
	int t, a, b;
	scanf("%d", &t);
	while(t--){
		scanf("%d%d", &a, &b);
		printf("%lld\n", gao(a, b));
	}
	return 0;
}
*/

//4.	http://acm.seu.edu.cn/oj/problem.php?id=123
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <string>
#include <set>
#include <map>
using namespace std;
typedef unsigned int lint;

const int MAX_N = 1e7 + 20;
int prime[MAX_N], tot;
bool notp[MAX_N];
lint dp[MAX_N], f[MAX_N];

inline lint po(lint p, int v){
	lint ans = 1;
	while(v){
		if(v & 1) ans = ans * p;
		p = p * p;
		v >>= 1;
	}
	return ans;
}

void init(){
	tot = 0; f[1] = 1;
	for(int i = 2; i < MAX_N; i++){
		if(!notp[i]){
			prime[tot++] = i;
			f[i] = i - 1;
		}
		for(int j = 0; j < tot; j++){
			int t = i * prime[j];
			if(t >= MAX_N) break;
			notp[t] = 1;
			if( i % prime[j] == 0 ){
				f[t] = f[i] * prime[j];
				break;
			}else
				f[t] = f[i] * (prime[j] - 1);
		}
	}
}

void initdp(int v, int nm){
	for(int i = 0; i < nm; i++){
		dp[i] = f[i] * po(i, v);
	}
	for(int i = 1; i < nm; i++){
		dp[i] += dp[i - 1];
	}
}

lint gao(int n, int m, int a, int b){
	lint ans = 0, da = 0, db = 0;
	for(int i = 1; i <= n; i++) da += po(i, a);
	for(int i = 1; i <= m; i++) db += po(i, b);
	lint pa = n, pb = m;
	for(int i = 1, j = 1; i <= n; i = j + 1){
		j = min(n / (n/i), m / (m/i));
		while(pa > n/i) da -= po(pa, a), pa--;
		while(pb > m/i) db -= po(pb, b), pb--;
		ans += (dp[j] - dp[i - 1]) * da * db;
	}
	return ans;
}

int main() {
	init();
	int t, a, b, n, m, icase = 1;
	scanf("%d", &t);
	while(t--){
		scanf("%d%d%d%d", &n, &m, &a, &b);
		// remember swap n,m a,b the same time
		if(n > m) swap(n, m), swap(a, b);

		initdp(a + b - 2, max(n,m) + 2);
		printf("Case #%d: %u\n",icase++, gao(n, m, a - 1, b - 1));
	}
	return 0;
}
