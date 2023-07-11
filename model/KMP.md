```c++
template<class T>
std::vector<int> prefixFunction(T S) {
    int n = S.size();
    std::vector<int> P(n);
    for (int i = 2, j = 0; i < n; i++) {
        while (j && S[i] != S[j + 1]) {
            j = P[j];
        }
        if (S[i] == S[j + 1]) {
            j++;
        }
        P[i] = j;
    }
    return P;
}

template<class T>
std::vector<int> KMP(T Text, T Pattern) {
    int n = Pattern.size(), m = Text.size();
    std::vector<int> P = prefixFunction(Pattern);
    std::vector<int> Pos;
    for (int i = 1, j = 0; i < m; i++) {
        while (j && Text[i] != Pattern[j + 1]) {
            j = P[j];
        }
        if (Text[i] == Pattern[j + 1]) {
            j++;
        }
        if (j == n - 1) {
            Pos.push_back(i);
            j = P[j];
        }
    }
    return Pos;
}
```