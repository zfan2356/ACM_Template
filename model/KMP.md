```c++
template<class T>
std::vector<int> prefixFunction(T S) {
    int n = S.size();
    std::vector<int> P(n);
    for (int i = 2, j = 0; i <= n; i++) {
        while (j && S[i - 1] != S[j]) {
            j = P[j - 1];
        }
        if (S[i - 1] == S[j]) {
            j++;
        }
        P[i - 1] = j;
    }
    return P;
}

template<class T>
std::vector<int> KMP(T Text, T Pattern) {
    int n = Pattern.size(), m = Text.size();
    std::vector<int> P = prefixFunction(Pattern);
    std::vector<int> Pos;
    for (int i = 1, j = 0; i <= m; i++) {
        while (j && Text[i - 1] != Pattern[j]) {
            j = P[j - 1];
        }
        if (Text[i - 1] == Pattern[j]) {
            j++;
        }
        if (j == n) {
            Pos.push_back(i);
            j = P[j - 1];
        }
    }
    return Pos;
}
```