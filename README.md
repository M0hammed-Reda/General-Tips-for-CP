# Competitive Programming Notes

## General Tips for Competitive Programming
1. **Read the Problem Statement Carefully**: Misunderstanding the problem can lead to incorrect solutions.
2. **Start with Simple Examples**: Use small test cases to understand the problem.
3. **Optimize for Time and Space**: Be aware of time and space complexity.
4. **Edge Cases**: Consider boundary conditions and edge cases (e.g., empty input, large input).
5. **Use Fast I/O**: Use `ios_base::sync_with_stdio(false); cin.tie(NULL);` for faster input/output.

## C++ Specific Tips
1. **Template Usage**:
    ```cpp
    #define fast_io ios_base::sync_with_stdio(false); cin.tie(NULL);
    #define pb push_back
    #define mp make_pair
    #define fi first
    #define se second
    typedef long long ll;
    typedef pair<int, int> pii;
    typedef vector<int> vi;
    typedef vector<pii> vpii;
    ```

2. **Common Libraries**:
    ```cpp
    #include <bits/stdc++.h>
    using namespace std;
    ```

3. **STL Containers and Algorithms**:
    - **Vector**: Dynamic array, random access.
    - **Set/Multiset**: Unique elements, automatic sorting.
    - **Map/Multimap**: Key-value pairs, automatic sorting by keys.
    - **Priority Queue**: Implements heap.
    - **Deque**: Double-ended queue.

4. **Common STL Functions**:
    - `sort(v.begin(), v.end())`
    - `reverse(v.begin(), v.end())`
    - `next_permutation(v.begin(), v.end())`
    - `lower_bound(v.begin(), v.end(), value)`
    - `upper_bound(v.begin(), v.end(), value)`
    - `binary_search(v.begin(), v.end(), value)`
    - `accumulate(v.begin(), v.end(), init_value)`

## Data Structures
1. **Array**: Fixed-size sequence of elements.
    - Use for simple, static collections of data.
    - O(1) access time.

2. **Linked List**: Sequence of elements, each containing a reference to the next element.
    - Use for dynamic data structures where insertion and deletion are frequent.
    - O(1) insertion/deletion time, O(n) access time.

3. **Stack**: LIFO data structure.
    - Use for backtracking, parsing expressions.
    - O(1) push/pop time.

4. **Queue**: FIFO data structure.
    - Use for scheduling, buffering.
    - O(1) enqueue/dequeue time.

5. **Binary Tree**: Tree data structure with at most two children.
    - Use for hierarchical data.
    - Traversals: Inorder, Preorder, Postorder.

6. **Binary Search Tree (BST)**: Binary tree with ordered elements.
    - Use for dynamic sets, maintains sorted order.
    - Average O(log n) time for insert/search/delete.

7. **Heap**: Complete binary tree, min-heap or max-heap.
    - Use for priority queues.
    - O(log n) insertion/deletion time.

8. **Graph**: Collection of nodes (vertices) and edges.
    - Representations: Adjacency matrix, adjacency list.
    - Traversals: BFS (Breadth-First Search), DFS (Depth-First Search).
    - Algorithms: Dijkstra's (shortest path), Kruskal's (MST), Prim's (MST).

## Algorithms
1. **Sorting Algorithms**:
    - **Quick Sort**: Average O(n log n), worst O(n^2).
    - **Merge Sort**: O(n log n), stable.
    - **Heap Sort**: O(n log n), in-place.

2. **Search Algorithms**:
    - **Binary Search**: O(log n), requires sorted array.
    - **Breadth-First Search (BFS)**: O(V + E), shortest path in unweighted graph.
    - **Depth-First Search (DFS)**: O(V + E), path finding, topological sort.

3. **Dynamic Programming**:
    - **Memoization**: Top-down approach, store intermediate results.
    - **Tabulation**: Bottom-up approach, build solution iteratively.
    - **Common Problems**: Fibonacci, knapsack, longest common subsequence, matrix chain multiplication.

4. **Greedy Algorithms**:
    - **Characteristics**: Make the locally optimal choice at each step.
    - **Common Problems**: Activity selection, Huffman coding, Prim's and Kruskal's for MST.

5. **Graph Algorithms**:
    - **Dijkstra's Algorithm**: Shortest path in weighted graph, O(V^2) with adjacency matrix, O((V + E) log V) with adjacency list + priority queue.
    - **Floyd-Warshall Algorithm**: All-pairs shortest path, O(V^3).
    - **Bellman-Ford Algorithm**: Shortest path with negative weights, O(VE).
    - **Kruskal's and Prim's Algorithms**: Minimum Spanning Tree.

## Problem-Solving Strategies
1. **Brute Force**: Try all possibilities, often leads to inefficiency.
2. **Divide and Conquer**: Break problem into smaller subproblems, solve recursively.
3. **Backtracking**: Build solution incrementally, discard solutions that fail to satisfy constraints.
4. **Greedy**: Make a series of choices, each of which is locally optimal.
5. **Dynamic Programming**: Solve overlapping subproblems, store results to avoid recomputation.

## Contest Preparation
1. **Practice**: Solve problems on platforms like Codeforces, LeetCode, HackerRank, and AtCoder.
2. **Time Management**: Work on improving speed and accuracy under timed conditions.
3. **Debugging**: Learn to debug quickly, use `assert` statements to catch errors early.
4. **Team Coordination**: If in a team contest, practice dividing and coordinating work efficiently.
