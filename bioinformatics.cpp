// bioinformatics.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>

int patternCount(std::string genom, std::string pattern) {
    int count = 0;

    for (auto i = 0; i <= genom.length() - pattern.length();) {
        auto position = genom.find(pattern, i);
        if (position == std::string::npos) break;
        count++;
        i = position + 1;
    }
    return count;
}

std::vector<std::string> frequentWords(std::string text, int k) {
    std::vector<std::string> frequenPatterns;
    std::vector<size_t> count;
    for (auto i = 0; i <= text.length() - k; i++) {
        std::string pattern = text.substr(i, k);
        count.push_back(patternCount(text, pattern));
    }
    int maxCount = *std::max_element(count.begin(), count.end());
    for (auto i = 0; i <= text.length() - k; i++) {
        if (count[i] == maxCount)
            frequenPatterns.push_back(text.substr(i, k));
    }

    std::sort(std::begin(frequenPatterns), std::end(frequenPatterns));
    frequenPatterns.erase(std::unique(frequenPatterns.begin(), frequenPatterns.end()), frequenPatterns.end());
    return frequenPatterns;
}

std::string reverseComplement(std::string pattern) {
    std::string answer;
    for (auto i = 0; i < pattern.length(); i++) {
        switch (pattern[pattern.length()-1-i])
        {
        case 'A':
            answer.push_back('T');
            break;
        case 'T':
            answer.push_back('A');
            break;
        case 'C':
            answer.push_back('G');
            break;
        case 'G':
            answer.push_back('C');
            break;
        default:
            break;
        }
    }

    return answer;
}
std::vector<size_t> patternMatching(std::string genom, std::string pattern) {
    std::vector<size_t> answer;
    for (auto i = 0; i <= genom.length() - pattern.length();) {
        auto position = genom.find(pattern, i);
        if (position == std::string::npos) break;
        answer.push_back(position);
        i = position + 1;
    }

    return answer;
}
std::map<std::string, size_t> frequencyTable(std::string text, int k) {
    std::map<std::string, size_t> freqMap;
    auto n = text.length();
    for (size_t i = 0; i < n - k; i++) {
        std::string pattern = text.substr(i, k);
        if (freqMap.count(pattern) == 0)
            freqMap[pattern] = 1;
        else
            freqMap[pattern] = freqMap[pattern] + 1;
    }
    return freqMap;
}
std::vector<std::string> findClumps(std::string text, int k, int L, int t) {
    std::vector<std::string> patterns;
    size_t n = text.length();
    for (size_t i = 0; i < n - L; i++) {
        std::string window = text.substr(i, L);
        std::map<std::string, size_t> freqMap = frequencyTable(window, k);
        for (const auto& element : freqMap)
            if (freqMap[element.first] >= t)
                patterns.push_back(element.first);
    }

    std::sort(std::begin(patterns), std::end(patterns));
    patterns.erase(std::unique(patterns.begin(), patterns.end()), patterns.end());
    return patterns;
}
size_t hammingDistance(std::string s1, std::string s2) {
    size_t count = 0;
    if (s1.length() != s2.length())
        return -1;
    for (size_t i = 0; i < s1.length(); i++) {
        if (s1[i] != s2[i])
            count++;
    }

    return count;
}
std::vector<size_t> approximatePatternMatching(std::string pattern, std::string genom, int d) {
    std::vector<size_t> answer;
    for (auto i = 0; i <= genom.length() - pattern.length();) {
        std::string pattern_ = genom.substr(i, pattern.length());
        auto position = genom.find(pattern_, i);
        if (position == std::string::npos) break;
        size_t hammingDistance_ = hammingDistance(pattern, pattern_);
        if (hammingDistance_ <= d) {
            answer.push_back(position);
        }
        i = position + 1;
    }

    return answer;
}
size_t approximatePatternCount(std::string pattern, std::string genom, int d) {
    size_t count = 0;
    for (auto i = 0; i <= genom.length() - pattern.length();) {
        std::string pattern_ = genom.substr(i, pattern.length());
        auto position = genom.find(pattern_, i);
        if (position == std::string::npos) break;
        size_t hammingDistance_ = hammingDistance(pattern, pattern_);
        if (hammingDistance_ <= d) {
            count++;
        }
        i = position + 1;
    }

    return count;
}
std::unordered_map<std::string, int> getNeighbours(std::string text, int k, int d) {
    std::unordered_map<std::string, int> neighbours;
    if (d == 0) {
        for (size_t i = 0; i < text.length(); i++)
            neighbours[text.substr(i, k)] = 0;
        return neighbours;
    }
    if (k == 1)
        return{ {"A,",0},{"C",0},{"G",0},{"T",0} };

    for (size_t i = 0; i < text.length(); i++) {
        std::string pattern = text.substr(i, k);
        for (size_t j = 0; j < pattern.length(); j++) {
            for (char nucliotide : {'A', 'C', 'G', 'T'}) {
                std::string neighbour = pattern;
                neighbour[j] = nucliotide;
                neighbours[neighbour] = 0;
            }
        }
    }
    return neighbours;
}
std::vector<std::string> frequentWordsWithMismatches(std::string text, int k, int d) {
    std::unordered_map<std::string, int> neighbours = getNeighbours(text, k, d);

    int max_ = 0;
    for (const auto& element : neighbours) {
        for (size_t i = 0; i < text.length() - k; i++) {
            if (hammingDistance(element.first, text.substr(i, k)) <= d) {
                neighbours[element.first]++;
                if (neighbours[element.first] > max_)
                    max_ = neighbours[element.first];
            }
        }
    }

    std::vector<std::string> answer;
    for (const auto& element : neighbours)
        if (element.second == max_)
            answer.push_back(element.first);

    return answer;
}

int main()
{
    //1.1 Pattern Count
    /*
    std::string genom;
    std::string pattern;

    std::cin >> pattern;
    std::cin >> genom;

    
    int count = patternCount(genom, pattern);
    std::cout << count;
    */

    //1.2 Frequent Words
   /* std::string text;
    int k;

    std::cin >> text;
    std::cin >> k;
    std::vector<std::string> frequenPatterns = frequentWords(text, k);
    
    for (auto i = 0; i < frequenPatterns.size(); i++)
        std::cout << frequenPatterns[i] << " ";*/

    //1.3 Reverse Complement
    /*
    std::string pattern;
    std::cin >> pattern;

    std::cout << reverseComplement(pattern);
    */

    //1.4 Pattern Matching
    /*
    std::string genom;
    std::string pattern;

    std::cin >> pattern;
    std::cin >> genom;

    std::vector<size_t> answer = patternMatching(genom, pattern);

    for (auto i = 0; i < answer.size(); i++)
        std::cout << answer[i] << " ";
    */

    //1.5 Clump Finding 
    /*
    std::string text;
    int k, L, t;
    std::cin >> text;
    std::cin >> k >> L >> t;

    std::vector<std::string> answer = findClumps(text, k, L, t);
    for (auto i = 0; i < answer.size(); i++)
        std::cout << answer[i] << " ";
    */  

    //1.6 Hamming Distance
    /*
    std::string s1;
    std::string s2;
    std::cin >> s1;
    std::cin >> s2;

    size_t count = hammingDistance(s1, s2);
    std::cout << count;
    */

    //1.7 Approximate Pattern Matching
    /*
    std::string pattern;
    std::string genom;
    int d;
    std::cin >> pattern;
    std::cin >> genom;
    std::cin >> d;

    std::vector<size_t> answer = approximatePatternMatching(pattern, genom, d);
    for (size_t i = 0; i < answer.size(); i++)
        std::cout << answer[i] << " ";
    */

    //1.8 Approximate Pattern Count
    /*
    std::string pattern;
    std::string genom;
    int d;
    std::cin >> pattern;
    std::cin >> genom;
    std::cin >> d;

    size_t count = approximatePatternCount(pattern, genom, d);
    std::cout << count;
    */

    //1.9 Frequent Words with Mismatches
    std::string text;
    int k, d;
    std::cin >> text;
    std::cin >> k >> d;

    std::vector<std::string> answer = frequentWordsWithMismatches(text, k, d);
    for (size_t i = 0; i < answer.size(); i++)
        std::cout << answer[i] << " ";

    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
