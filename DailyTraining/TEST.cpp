#include <fstream>
#include <iostream>

int main() {
    std::ifstream file1("D:\\mypile\\acm\\ICPC\\DailyTraining\\data1.txt");
    std::ifstream file2("D:\\mypile\\acm\\ICPC\\DailyTraining\\data2.txt");

    if (!file1.is_open() || !file2.is_open()) {
        std::cout << "error" << std::endl;
        return 1;
    }

    std::string line1, line2;
    int lineCount = 1; // 行号计数器
    bool hasDifference = false; // 记录是否存在差异

    while (std::getline(file1, line1) && std::getline(file2, line2)) {
        if (line1 != line2) {
            std::cout << "difference at " << lineCount << " line" << std::endl;
            std::cout << line1 << std::endl;
            std::cout << line2 << std::endl;
            hasDifference = true;
            return 0;
        }
        lineCount++;
    }

    // 如果两个文件的行数不同，也会被认为是有差异的
    if (std::getline(file1, line1) || std::getline(file2, line2)) {
        std::cout << "line number deference" << std::endl;
        hasDifference = true;
    }

    file1.close();
    file2.close();

    if (!hasDifference) {
        std::cout << "Yes" << std::endl;
    }

    return 0;
}
