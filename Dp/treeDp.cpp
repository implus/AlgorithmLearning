/*
http://acm.hdu.edu.cn/showproblem.php?pid=5103
how many rooted tree form ? 
limit: every point has subtree points number from l[i] to r[i]

状压dp，先考虑二叉树的情况，设f[opt]为以opt（二进制状态）中所有节点构成一棵有根树的方案。则我们每次需要枚举以哪个节点为根和哪些节点放在左子树，然后剩下的节点放在右子树，最后将两棵子树的方案相乘累加到f[opt]。但是这样做可能会重复，我们可以随便固定某个节点一定在左子树，这样就可以解决重复的问题。多叉树可以通过左儿子右兄弟的方法转化成二叉树的问题，这时只需要加一维状态f[opt][0/1], 0的意义和二叉树一样，1的意义是以opt构成一个有根树森林的方法，转移方式和二叉树类似。dp的状态总数为O(2^n)，转移先枚举根节点再枚举子集，所以总复杂度为O(n3^n)。
*/
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


const int MAX_S = 14;
int cnt[1<<MAX_S];
int l[MAX_S], r[MAX_S];
const int MOD = 1e9 + 7;
int dp[1<<MAX_S][2];

int main() {
	cnt[0] = 0;
	for(int i = 1; i < (1<<MAX_S); i++){
		cnt[i] = cnt[i >> 1] + (i & 1);
	}

	int T, n;
	scanf("%d", &T);
	while(T--){
		scanf("%d", &n);
		for(int i = 0; i < n; i++) {
			scanf("%d%d", &l[i], &r[i]);
		}
		memset(dp, 0, sizeof(dp));

		int MS = (1<<n);
		for(int s = 0; s < MS; s++){
			// s find rt
			int p = cnt[s];
			if(p == 0){
				dp[s][0] = dp[s][1] = 1;
			}else if(p == 1){
				int i;
				for(i = 0; i < n; i++) if(s & (1<<i)) break;
				if(l[i] <= 1 && 1 <= r[i]) dp[s][0] = dp[s][1] = 1;
				else dp[s][0] = dp[s][1] = 0;
			}else {
				for(int i = 0; i < n; i++) if(s&(1<<i)){ // rt i
					if(cnt[s] < l[i] || cnt[s] > r[i]) continue;
					// find left tree
					int j;
					for(j = 0; j < n; j++) if(j - i && (s & (1<<j))) break;
					int ls = (s ^ (1<<i));

					for(int zs = ls; ; zs = (zs - 1) & ls) {
						if(zs & (1<<j)){ // make j must be in left tree
							dp[s][0] += 1LL * dp[zs][0] * dp[ls ^ zs][1] % MOD;
							dp[s][0] %= MOD;
						}
						if(zs == 0) break;
					}
				}

				int i;
				for(i = 0; i < n; i++) if(s & (1<<i)) break; // left forest rt

				for(int zs = s; ;zs = (zs - 1) & s) {
					if(zs & (1<<i)){ // make i must be in left tree
						dp[s][1] += 1LL * dp[zs][0] * dp[s ^ zs][1] % MOD;
						dp[s][1] %= MOD;
					}
					if(zs == 0) break;
				}
			}
		}

		printf("%d\n", dp[MS - 1][0]);
	}
	return 0;
}

