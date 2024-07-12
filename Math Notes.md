
# Competitive Programming Mathematics Notes

This repository contains essential mathematical concepts and functions for competitive programming. These notes will help you prepare for contests, especially those focusing on the C++ programming language and data structures and algorithms.

## Table of Contents
- [Greatest Common Divisor (GCD) and Least Common Multiple (LCM)](#gcd-and-lcm)
- [Prime Numbers and Prime Factorization](#prime-numbers-and-prime-factorization)
- [Modular Arithmetic](#modular-arithmetic)
- [Combinatorics](#combinatorics)
- [Geometry](#geometry)
- [Bit Manipulation](#bit-manipulation)
- [Matrix Exponentiation](#matrix-exponentiation)
- [Fibonacci Numbers](#fibonacci-numbers)
- [Powers and Logarithms](#powers-and-logarithms)
- [Fast Fourier Transform (FFT)](#fast-fourier-transform-fft)
- [Graph Theory](#graph-theory)
- [Number Theory](#number-theory)
- [Dynamic Programming](#dynamic-programming)
- [Sorting Algorithms](#sorting-algorithms)

## GCD and LCM

### Greatest Common Divisor (GCD)
The GCD is the largest number that divides two or more numbers without leaving a remainder. It can be calculated using the Euclidean algorithm.
```cpp
int gcd(int a, int b) {
    if (b == 0) return a;
    return gcd(b, a % b);
}
```

### Least Common Multiple (LCM)
The LCM is the smallest number that is a multiple of two or more numbers. It can be found using the relation between GCD and LCM:
```cpp
int lcm(int a, int b) {
    return (a / gcd(a, b)) * b;
}
```

## Prime Numbers and Prime Factorization

### Prime Check
To check if a number is prime, iterate from 2 to the square root of the number.
```cpp
bool isPrime(int n) {
    if (n <= 1) return false;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}
```

### Sieve of Eratosthenes
An efficient algorithm to find all prime numbers up to a given limit.
```cpp
void sieve(int n) {
    vector<bool> isPrime(n+1, true);
    isPrime[0] = isPrime[1] = false;
    for (int i = 2; i * i <= n; i++) {
        if (isPrime[i]) {
            for (int j = i * i; j <= n; j += i) {
                isPrime[j] = false;
            }
        }
    }
}
```

## Modular Arithmetic

### Modular Exponentiation
Efficiently computes `(base^exp) % mod`.
```cpp
int modExp(int base, int exp, int mod) {
    int result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        exp = exp >> 1;
        base = (base * base) % mod;
    }
    return result;
}
```

### Modular Inverse
For a number `a` and modulus `m`, the modular inverse is a number `x` such that `(a * x) % m = 1`. It can be found using the Extended Euclidean Algorithm.
```cpp
int modInverse(int a, int m) {
    int m0 = m, t, q;
    int x0 = 0, x1 = 1;
    if (m == 1) return 0;
    while (a > 1) {
        q = a / m;
        t = m;
        m = a % m, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
}
```

## Combinatorics

### Factorials and Binomial Coefficients
Factorials can be used to compute combinations and permutations.
```cpp
int factorial(int n) {
    if (n == 0) return 1;
    return n * factorial(n - 1);
}

int binomialCoeff(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    return binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k);
}
```

## Geometry

### Distance between two points
```cpp
double distance(double x1, double y1, double x2, double y2) {
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
```

### Area of a triangle using coordinates
```cpp
double areaOfTriangle(double x1, double y1, double x2, double y2, double x3, double y3) {
    return abs((x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)) / 2.0);
}
```

## Bit Manipulation

### Checking if a number is a power of 2
```cpp
bool isPowerOfTwo(int n) {
    return (n && !(n & (n - 1)));
}
```

### Counting the number of 1s in binary representation (Hamming Weight)
```cpp
int hammingWeight(int n) {
    int count = 0;
    while (n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
}
```

## Matrix Exponentiation

Useful for solving linear recurrence relations efficiently.
```cpp
const int MOD = 1e9 + 7;
vector<vector<int>> multiply(vector<vector<int>>& A, vector<vector<int>>& B) {
    int n = A.size();
    vector<vector<int>> C(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                C[i][j] = (C[i][j] + 1LL * A[i][k] * B[k][j] % MOD) % MOD;
            }
        }
    }
    return C;
}

vector<vector<int>> matrixExpo(vector<vector<int>>& A, int exp) {
    int n = A.size();
    vector<vector<int>> result(n, vector<int>(n));
    for (int i = 0; i < n; i++) result[i][i] = 1;
    while (exp > 0) {
        if (exp % 2 == 1) result = multiply(result, A);
        A = multiply(A, A);
        exp /= 2;
    }
    return result;
}
```

## Fibonacci Numbers

### Recursive Method
```cpp
int fibonacci(int n) {
    if (n <= 1) return n;
    return fibonacci(n - 1) + fibonacci(n - 2);
}
```

### Dynamic Programming
```cpp
int fibonacci(int n) {
    if (n <= 1) return n;
    int fib[n + 1];
    fib[0] = 0;
    fib[1] = 1;
    for (int i = 2; i <= n; i++) {
        fib[i] = fib[i - 1] + fib[i - 2];
    }
    return fib[n];
}
```

### Matrix Exponentiation
```cpp
void multiply(int F[2][2], int M[2][2]) {
    int x = F[0][0] * M[0][0] + F[0][1] * M[1][0];
    int y = F[0][0] * M[0][1] + F[0][1] * M[1][1];
    int z = F[1][0] * M[0][0] + F[1][1] * M[1][0];
    int w = F[1][0] * M[0][1] + F[1][1] * M[1][1];

    F[0][0] = x;
    F[0][1] = y;
    F[1][0] = z;
    F[1][1] = w;
}

void power(int F[2][2], int n) {
    if (n == 0 || n == 1) return;
    int M[2][2] = {{1, 1}, {1, 0}};

    power(F, n / 2);
    multiply(F, F);

    if (n % 2 != 0) multiply(F, M);
}

int fibonacci(int n) {
    int F[2][2] = {{1, 1}, {1, 0}};
    if (n == 0) return 0;
    power(F, n - 1);
    return F[0][0];
}
```

## Powers and Logarithms

### Power Function
```cpp
int power(int base, int exp) {
    int result

 = 1;
    while (exp > 0) {
        if (exp % 2 == 1) result *= base;
        base *= base;
        exp /= 2;
    }
    return result;
}
```

### Logarithm Base 2
```cpp
int log2(int n) {
    int res = 0;
    while (n >>= 1) res++;
    return res;
}
```

## Fast Fourier Transform (FFT)

Used for polynomial multiplication.
```cpp
#include <complex>
#include <vector>
#include <cmath>

using namespace std;

typedef complex<double> cd;
const double PI = acos(-1);

void fft(vector<cd>& a, bool invert) {
    int n = a.size();
    if (n == 1) return;

    vector<cd> a0(n / 2), a1(n / 2);
    for (int i = 0; 2 * i < n; i++) {
        a0[i] = a[i * 2];
        a1[i] = a[i * 2 + 1];
    }
    fft(a0, invert);
    fft(a1, invert);

    double angle = 2 * PI / n * (invert ? -1 : 1);
    cd w(1), wn(cos(angle), sin(angle));
    for (int i = 0; 2 * i < n; i++) {
        a[i] = a0[i] + w * a1[i];
        a[i + n / 2] = a0[i] - w * a1[i];
        if (invert) {
            a[i] /= 2;
            a[i + n / 2] /= 2;
        }
        w *= wn;
    }
}
```

## Graph Theory

### Breadth-First Search (BFS)
```cpp
void bfs(int start, vector<vector<int>>& adj) {
    vector<bool> visited(adj.size(), false);
    queue<int> q;
    q.push(start);
    visited[start] = true;
    while (!q.empty()) {
        int node = q.front();
        q.pop();
        for (int neighbor : adj[node]) {
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }
}
```

### Depth-First Search (DFS)
```cpp
void dfs(int node, vector<vector<int>>& adj, vector<bool>& visited) {
    visited[node] = true;
    for (int neighbor : adj[node]) {
        if (!visited[neighbor]) {
            dfs(neighbor, adj, visited);
        }
    }
}
```

## Number Theory

### Sieve of Eratosthenes
```cpp
void sieve(int n) {
    vector<bool> isPrime(n+1, true);
    isPrime[0] = isPrime[1] = false;
    for (int i = 2; i * i <= n; i++) {
        if (isPrime[i]) {
            for (int j = i * i; j <= n; j += i) {
                isPrime[j] = false;
            }
        }
    }
}
```

### Euler's Totient Function
```cpp
int eulerTotient(int n) {
    int result = n;
    for (int p = 2; p * p <= n; p++) {
        if (n % p == 0) {
            while (n % p == 0) n /= p;
            result -= result / p;
        }
    }
    if (n > 1) result -= result / n;
    return result;
}
```

## Dynamic Programming

### Longest Common Subsequence (LCS)
```cpp
int lcs(string s1, string s2) {
    int m = s1.size(), n = s2.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1));
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (s1[i - 1] == s2[j - 1]) dp[i][j] = dp[i - 1][j - 1] + 1;
            else dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
        }
    }
    return dp[m][n];
}
```

### Knapsack Problem
```cpp
int knapsack(int W, vector<int>& weights, vector<int>& values) {
    int n = weights.size();
    vector<vector<int>> dp(n + 1, vector<int>(W + 1));
    for (int i = 1; i <= n; i++) {
        for (int w = 0; w <= W; w++) {
            if (weights[i - 1] <= w)
                dp[i][w] = max(dp[i - 1][w], dp[i - 1][w - weights[i - 1]] + values[i - 1]);
            else dp[i][w] = dp[i - 1][w];
        }
    }
    return dp[n][W];
}
```

## Sorting Algorithms

### Quick Sort
```cpp
int partition(vector<int>& arr, int low, int high) {
    int pivot = arr[high];
    int i = low - 1;
    for (int j = low; j < high; j++) {
        if (arr[j] < pivot) {
            i++;
            swap(arr[i], arr[j]);
        }
    }
    swap(arr[i + 1], arr[high]);
    return i + 1;
}

void quickSort(vector<int>& arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
```

### Merge Sort
```cpp
void merge(vector<int>& arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;
    vector<int> L(n1), R(n2);
    for (int i = 0; i < n1; i++) L[i] = arr[l + i];
    for (int j = 0; j < n2; j++) R[j] = arr[m + 1 + j];
    int i = 0, j = 0, k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) arr[k++] = L[i++];
        else arr[k++] = R[j++];
    }
    while (i < n1) arr[k++] = L[i++];
    while (j < n2) arr[k++] = R[j++];
}

void mergeSort(vector<int>& arr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}
```
```

This is a comprehensive overview of essential mathematics topics for competitive programming. These notes should help you quickly reference important algorithms and concepts during your contest preparation.
