# Otomachi Una Library

Otomachi Una 的 XCPC 算法模板库。使用 GNU C++17 或更新标准。

以 `basic/template.cpp` 为底稿，把所需模板及其依赖粘贴到 `main` 前。除底稿外，模板不带头文件和 `main`。只选需要的实现，同名结构和全局变量不能同时粘贴。

文件顶部说明接口，末尾 `// !!!!! ... !!!!` 标出依赖和重要前提。下标、初值和复杂度以对应模板说明为准。

- [使用手册](manual/manual.pdf)：接口、调用顺序、例子和限制；源码为 [manual.tex](manual/manual.tex)。根目录 `manual.pdf` 是同一份文件。
- [XCPC 考场手册](manual/icpc-handbook.pdf)：保留的独立旧版，部分实现更短；本次没有更新，其中的旧路径不代表当前库路径。
- [路径迁移记录](PATHS.md)：旧路径与新路径对应关系，结构体及函数名尽量保留。

`basic/fast-io.cpp` 使用 fread/fwrite，适用于 Windows/Linux 的 GNU C++，同时提供 `__int128` 的普通流输入输出。SA/LCP/LCS 和 Runs 的粘贴顺序是 `linear-rmq.cpp` → `suffix-array.cpp` → `runs.cpp`。

`manual/icpc/code` 是旧手册的历史代码快照。`manual/tools/build_icpc_handbook.py` 可从当前模板生成另一本 `generated-handbook.tex`，不会覆盖考场手册，也不会自动替换成简化算法。使用手册在 manual 目录运行 XeLaTeX 三次即可重建。

## 模板目录

### 基础

- [basic/coordinate-compression.cpp](basic/coordinate-compression.cpp)
- [basic/fast-io.cpp](basic/fast-io.cpp)
- [basic/template.cpp](basic/template.cpp)
- [basic/vector.cpp](basic/vector.cpp)

### 数据结构

- [data-structure/disjoint-set-rollback.cpp](data-structure/disjoint-set-rollback.cpp)
- [data-structure/disjoint-set.cpp](data-structure/disjoint-set.cpp)
- [data-structure/fast-set.cpp](data-structure/fast-set.cpp)
- [data-structure/fenwick/fenwick.cpp](data-structure/fenwick/fenwick.cpp)
- [data-structure/fenwick/max.cpp](data-structure/fenwick/max.cpp)
- [data-structure/fenwick/min.cpp](data-structure/fenwick/min.cpp)
- [data-structure/fenwick/sum.cpp](data-structure/fenwick/sum.cpp)
- [data-structure/fhq-treap.cpp](data-structure/fhq-treap.cpp)
- [data-structure/li-chao-tree.cpp](data-structure/li-chao-tree.cpp)
- [data-structure/li-chao-tree-static.cpp](data-structure/li-chao-tree-static.cpp)
- [data-structure/linear-rmq.cpp](data-structure/linear-rmq.cpp)
- [data-structure/ordered-disjoint-interval-tree-fast.cpp](data-structure/ordered-disjoint-interval-tree-fast.cpp)
- [data-structure/ordered-disjoint-interval-tree.cpp](data-structure/ordered-disjoint-interval-tree.cpp)
- [data-structure/priority-queue-two-stacks.cpp](data-structure/priority-queue-two-stacks.cpp)
- [data-structure/rollback-int.cpp](data-structure/rollback-int.cpp)
- [data-structure/rollback-int64.cpp](data-structure/rollback-int64.cpp)
- [data-structure/segment-tree/dynamic.cpp](data-structure/segment-tree/dynamic.cpp)
- [data-structure/segment-tree/iterative-prefix.cpp](data-structure/segment-tree/iterative-prefix.cpp)
- [data-structure/segment-tree/lazy.cpp](data-structure/segment-tree/lazy.cpp)
- [data-structure/segment-tree/merge.cpp](data-structure/segment-tree/merge.cpp)
- [data-structure/segment-tree/point.cpp](data-structure/segment-tree/point.cpp)
- [data-structure/segment-tree/range-add-max.cpp](data-structure/segment-tree/range-add-max.cpp)
- [data-structure/segment-tree/range-add-min.cpp](data-structure/segment-tree/range-add-min.cpp)
- [data-structure/segment-tree/range-add-sum.cpp](data-structure/segment-tree/range-add-sum.cpp)
- [data-structure/sparse-table.cpp](data-structure/sparse-table.cpp)

### 图论

- [graph/block-cut-tree.cpp](graph/block-cut-tree.cpp)
- [graph/chordal.cpp](graph/chordal.cpp)
- [graph/dijkstra.cpp](graph/dijkstra.cpp)
- [graph/directed-mst.cpp](graph/directed-mst.cpp)
- [graph/edge-biconnected.cpp](graph/edge-biconnected.cpp)
- [graph/flow/max-flow-hlpp.cpp](graph/flow/max-flow-hlpp.cpp)
- [graph/flow/min-cost-flow-isap.cpp](graph/flow/min-cost-flow-isap.cpp)
- [graph/flow/min-cost-flow-simplex.cpp](graph/flow/min-cost-flow-simplex.cpp)
- [graph/matching/bipartite.cpp](graph/matching/bipartite.cpp)
- [graph/matching/general.cpp](graph/matching/general.cpp)
- [graph/matching/regular-bipartite.cpp](graph/matching/regular-bipartite.cpp)
- [graph/steiner-tree.cpp](graph/steiner-tree.cpp)
- [graph/two-sat.cpp](graph/two-sat.cpp)

### 树

- [tree/cartesian-tree.cpp](tree/cartesian-tree.cpp)
- [tree/centroid-decomposition.cpp](tree/centroid-decomposition.cpp)
- [tree/divide-combine-tree.cpp](tree/divide-combine-tree.cpp)
- [tree/heavy-light-decomposition.cpp](tree/heavy-light-decomposition.cpp)

### 字符串

- [string/aho-corasick.cpp](string/aho-corasick.cpp)
- [string/generalized-suffix-automaton.cpp](string/generalized-suffix-automaton.cpp)
- [string/kmp-z.cpp](string/kmp-z.cpp)
- [string/lyndon.cpp](string/lyndon.cpp)
- [string/manacher.cpp](string/manacher.cpp)
- [string/palindromic-automaton.cpp](string/palindromic-automaton.cpp)
- [string/rolling-hash.cpp](string/rolling-hash.cpp)
- [string/runs.cpp](string/runs.cpp)
- [string/suffix-array.cpp](string/suffix-array.cpp)
- [string/suffix-automaton.cpp](string/suffix-automaton.cpp)

### 数学

- [math/big-int.cpp](math/big-int.cpp)
- [math/dynamic-mod-int.cpp](math/dynamic-mod-int.cpp)
- [math/gaussian-elimination.cpp](math/gaussian-elimination.cpp)
- [math/linear-basis/field-prefix.cpp](math/linear-basis/field-prefix.cpp)
- [math/linear-basis/field.cpp](math/linear-basis/field.cpp)
- [math/linear-basis/xor-prefix.cpp](math/linear-basis/xor-prefix.cpp)
- [math/linear-basis/xor.cpp](math/linear-basis/xor.cpp)
- [math/mod-int.cpp](math/mod-int.cpp)
- [math/rational.cpp](math/rational.cpp)
- [math/xor-equations.cpp](math/xor-equations.cpp)

### 数论

- [number-theory/dujiao-sieve.cpp](number-theory/dujiao-sieve.cpp)
- [number-theory/floor-sum.cpp](number-theory/floor-sum.cpp)
- [number-theory/linear-sieve.cpp](number-theory/linear-sieve.cpp)
- [number-theory/modular-arithmetic.cpp](number-theory/modular-arithmetic.cpp)
- [number-theory/pollard-rho.cpp](number-theory/pollard-rho.cpp)
- [number-theory/stern-brocot.cpp](number-theory/stern-brocot.cpp)

### 多项式

- [polynomial/fwt.cpp](polynomial/fwt.cpp)
- [polynomial/integer.cpp](polynomial/integer.cpp)
- [polynomial/mtt.cpp](polynomial/mtt.cpp)
- [polynomial/ntt-fast.cpp](polynomial/ntt-fast.cpp)
- [polynomial/ntt-int64.cpp](polynomial/ntt-int64.cpp)
- [polynomial/ntt.cpp](polynomial/ntt.cpp)

### 计算几何

- [geometry/geometry-2d.cpp](geometry/geometry-2d.cpp)

辅助工具：[交互测试器](basic/interact.py)。
