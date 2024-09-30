#pragma once
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
        switch (pattern[pattern.length() - 1 - i])
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

