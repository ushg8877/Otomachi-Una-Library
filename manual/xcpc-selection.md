# XCPC 手册保留清单

按当前 `icpc-handbook.tex` 的 72 个条目整理，编号固定。
这是考场手册现有目录，路径沿用手册旧名称，不代表源码库的最新目录。

## 怎么选

- 把需要保留的条目前 `[ ]` 改成 `[x]`，保存后告诉我“按清单改”。
- 也可以直接在对话里回复“保留 01、03、14–18……”或“删除 02、05，其余保留”。
- 只留部分功能时，在该条目后追加说明，例如“只留凸包、旋转卡壳”。
- 尚未反馈前，未勾选只表示待选；不会据此删改手册，也不会删除源码。


## 基础与数值类型

- [x] **01 公共底稿** — `basic/src.cpp`
- [ ] **02 静态区间最值** — `basic/ST table.cpp`
- [x] **03 四毛子 RMQ：线性预处理与空间** — `basic/ST table (linear).cpp`
- [x] **04 普通并查集** — `basic/DSU.cpp`
- [ ] **05 128 位整数输入输出** — `basic/int128_IO.cpp`
- [x] **06 高精度有符号整数** — `basic/BigInt.cpp`
- [ ] **07 普通文件上的 mmap 快读快写** — `basic/FastIO.cpp`
- [x] **08 离散化并原地回写** — `basic/Coordinate Compression.cpp`
- [x] **09 固定模数及组合数** — `basic/ModInt.cpp`
- [x] **10 运行时模数及组合数** — `basic/ModInt-nonstaticMod.cpp`
- [ ] **11 记录并撤销变量修改** — `basic/roll back saver (int).cpp`
- [ ] **12 记录并撤销变量修改** — `basic/roll back saver (int or ll).cpp`
- [x] **13 可撤销并查集** — `basic/DSU (with undo).cpp`
## 数据结构

- [ ] **14 自定义下标的树状数组** — `ds/BIT.cpp`
- [ ] **15 树状数组：和特化** — `ds/BIT/Sum.cpp`
- [ ] **16 树状数组：最大值特化** — `ds/BIT/Max.cpp`
- [ ] **17 树状数组：最小值特化** — `ds/BIT/Min.cpp`
- [ ] **18 区间加、区间和** — `ds/Segment Tree/AddSum.cpp`
- [ ] **19 区间加、区间最大值** — `ds/Segment Tree/MaxAdd.cpp`
- [ ] **20 区间加、区间最小值** — `ds/Segment Tree/MinAdd.cpp`
- [ ] **21 通用线段树（带懒标记）** — `ds/Segment Tree(with lazytag).cpp`
- [ ] **22 通用线段树（不带懒标记）** — `ds/Segment Tree(without lazytag).cpp`
- [x] **23 非递归线段树框架** — `ds/zkw.cpp`
- [x] **24 大值域稀疏单点统计** — `ds/Dynamic Segment Tree.cpp`
- [ ] **25 可合并的多棵线段树** — `ds/Segment Tree Merge.cpp`
- [x] **26 线段李超树** — `ds/lichao.cpp`
- [ ] **27 有序区间集合** — `ds/ODT.cpp`
- [ ] **28 有界整数集合** — `ds/FastSet.cpp`
- [ ] **29 位集加速的区间覆盖** — `ds/FastODT.cpp`
- [ ] **30 按键或按序列位置分裂的 Treap** — `ds/FHQ.cpp`
- [x] **31 把优先队列删除转为栈撤销** — `ds/pq2stack.cpp`
- [x] **32 析合树：返回完整结构** — `ds/xihe_tree.cpp`
## 图论

- [x] **33 Dinic 最大流** — `graph/flow.cpp`
- [x] **34 最小费用最大流** — `graph/mcmf.cpp`
- [ ] **35 非负权单源最短路** — `graph/sssp.cpp`
- [x] **36 布尔约束求解** — `graph/2-sat.cpp`
- [x] **37 点双连通分量与圆方树** — `graph/block-cut-tree.cpp`
- [ ] **38 边双连通分量及桥森林** — `graph/Edge BCC.cpp`
- [x] **39 无权一般图最大匹配** — `graph/blossom.cpp`
- [x] **40 二分图最大匹配** — `graph/bipartite-matching.cpp`
- [ ] **41 正则二分图的完美匹配分解** — `graph/regular-bipartite.cpp`
- [x] **42 指定根的最小树形图** — `graph/directed-mst.cpp`
- [x] **43 弦图判定与完美消除序列** — `graph/chordal.cpp`
## 树上算法

- [ ] **44 LCA、祖先和路径分段** — `trees/HLD.cpp`
- [ ] **45 点分治框架** — `trees/DC.cpp`
- [x] **46 小根堆笛卡尔树** — `trees/cartesian.cpp`
## 数论

- [x] **47 模逆元和扩展 CRT** — `number theory/basic.cpp`
- [ ] **48 线性筛素数** — `number theory/linear sieve.cpp`
- [x] **49 大整数判素与质因数分解** — `number theory/pollard-rho.cpp`
- [x] **50 类欧几里得整除求和** — `number theory/floor sum.cpp`
- [x] **51 \ensuremath{\mu}、\ensuremath{\varphi} 前缀和** — `number theory/du jiao sieve.cpp`
## 多项式

- [x] **52 998244353 下的常用多项式** — `polynomial/poly-fast.cpp`
- [ ] **53 带长除、多点求值的多项式** — `polynomial/poly-new.cpp`
- [ ] **54 ll 系数的旧版多项式** — `polynomial/poly-old.cpp`
- [x] **55 非 NTT 模数下的多项式** — `polynomial/MTT.cpp`
- [x] **56 按位 AND / OR / XOR 卷积** — `polynomial/FWT.cpp`
## 字符串

- [x] **57 后缀数组和 O(1) LCP** — `string/sa.cpp`
- [x] **58 单串后缀自动机** — `string/sam.cpp`
- [x] **59 多模式串 AC 自动机** — `string/ac.cpp`
- [x] **60 广义后缀自动机** — `string/gsam.cpp`
- [x] **61 前缀函数和 Z 函数** — `string/kmp.cpp`
- [x] **62 双模前缀哈希** — `string/hash.cpp`
- [x] **63 回文自动机** — `string/pam.cpp`
- [x] **64 Manacher 回文半径** — `string/manacher.cpp`
- [x] **65 Runs 极大周期子串** — `string/runs.cpp`
- [x] **66 Lyndon 分解** — `string/lyndon.cpp`
## 线性代数

- [ ] **67 线性方程组与行列式** — `math/gauss.cpp`
- [ ] **68 异或线性基** — `math/linear basis/xor.cpp`
- [ ] **69 前缀异或线性基** — `math/linear basis/xor-prefix.cpp`
- [ ] **70 模意义线性基** — `math/linear basis/linear.cpp`
- [ ] **71 前缀模意义线性基** — `math/linear basis/linear-prefix.cpp`
## 计算几何

- [x] **72 整数二维几何基础** — `geometry/geometry-main.cpp`

## 额外要求

- 希望补进手册的模板：
- 希望压缩或合并的重复模板：
- 页数上限（没有就留空）：
