#include "testlib.h"
#include <bits/stdc++.h>

#define ll long long
using namespace std;

int main() {
    for (int i = 1; i <= 10000; i++) {
        system("Generator.exe>data.txt");
        system("myCode.exe<data.txt>myANS.txt");
        double st = clock();
        system("std.exe<data.txt>stdANS.txt");
        double ed = clock();
        if (system("fc myANS.txt stdANS.txt")) {
            printf("WA\n");
            break;
        } else printf("AC #%d Time:%.3lfms\n", i, ed - st);
    }
    return 0;
}

