// bioinformatics.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Introdaction.h"
#include "PeptideSequencing.h"
#include "CompareBiologicalSequence.h"


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
    /*
    std::string text;
    int k, d;
    std::cin >> text;
    std::cin >> k >> d;

    std::vector<std::string> answer = frequentWordsWithMismatches(text, k, d);
    for (size_t i = 0; i < answer.size(); i++)
        std::cout << answer[i] << " ";
    */

    /*======= SECTION TWO =======*/

    //2.1 Protein Translation Problem
    /*
    std::string RNA_pattern;
    std::cin >> RNA_pattern;
    std::string protein = proteinTranslationProblem(RNA_pattern);
    std::cout << protein;
    */

    //2.2 Peptide Encoding Problem
    /*
    std::string DNA;
    std::cin >> DNA;
    std::string geneticCode;
    std::cin >> geneticCode;

    std::vector<std::string> answer = peptideEncodingProblem(DNA, geneticCode);
    for (size_t i = 0; i < answer.size(); i++)
        std::cout << answer[i] << std::endl;
    */

    //2.3 Subpeptides Count Problem
    /*
    int n;
    std::cin >> n;
    long long int ans = subpeptidesCountProblem(n);
    std::cout << ans;
    */

    //2.4 Generating Theoretical Spectrum Problem
    /*
    std::string peptide;
    std::cin >> peptide;
    std::vector<int> cyclospectrum = generatingTheoreticalSpectrumProblem(peptide);
    for (int el : cyclospectrum)
        std::cout << el << " ";
    */

    //2.5 Counting Peptides with Given Mass Problem
    
    int n;
    std::cin >> n;
    int ans = countingPeptidesWithGivenMassProblem(n);
    std::cout << ans;


	    /*======= SECTION FORE =======*/

    //4.1 Change Problem
    /*
    int money;
    std::vector<int> coins;
    std::cin >> money;
    //std::copy(std::istream_iterator<int>(std::cin), std::istream_iterator<int>(), std::back_inserter(coins));
    //std::copy(coins.begin(), coins.end(), std::ostream_iterator<size_t>(std::cout, " "));
    std::cin.ignore();
    coins = ReadVectorFromConsol();

    int MinNumCoins = MoneyChangeProblem(money, coins);
    std::cout << MinNumCoins;
    */

    //4.2 Manhattan Tourist Problem
    /*
    size_t n, m;
    std::cin >> n >> m;

    std::vector<std::vector<size_t>> Down(n, std::vector<size_t>(m + 1));
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m + 1; ++j)
            std::cin >> Down[i][j];

    std::cin.ignore();
    std::cin.ignore();

    std::vector < std::vector<size_t>> Right(n + 1, std::vector<size_t>(m));
    for (size_t i = 0; i < n + 1; ++i)
        for (size_t j = 0; j < m; ++j)
            std::cin >> Right[i][j];
    size_t answer = ManhattanTouristProblem(n, m, Down, Right);
    std::cout << answer;
    */

    //4.3 Longest Common Subsequence Problem
    /*
    std::string v, w;
    std::cin >> v >> w;
    v = '-' + v;
    w = '-' + w;
    std::string answer;
    std::vector<std::vector<ARROW>> BackTrack = LCSBackTrack(v, w);
    OutPutLSC(BackTrack, v, v.length() - 1, w.length() - 1, answer);
    std::reverse(answer.begin(), answer.end());
    //answer = OutPutLSC(BackTrack, v);
    std::cout << answer;
    */

    //4.4 Topological Ordering
    /*
    std::unordered_map<int, std::vector<int>> graph;
    std::unordered_map<int, int> degreeEachNode;
    
    readGraph(graph, degreeEachNode);
    std::vector<int> result = topologicalOrdering(graph, degreeEachNode);
    
    if (!result.empty()) {
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << result[i];
            if (i != result.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << std::endl;
    }
    */

    //4.5 Longest Path in a DAG
    
    int source, sink;
    std::cin >> source >> sink;
    std::unordered_map<int, std::vector<Edge>> graph;
    std::unordered_map<int, int> in_degree;

    readWeightGraph(graph, in_degree);

    std::vector<int> topo_order = topologicalOrderingWeight(graph, in_degree);

    auto result = findLongestPath(graph, source, sink, topo_order);
    int longestPathLength = result.first;
    std::vector<int> longestPath = result.second;

    std::cout << longestPathLength << std::endl;
    if (!longestPath.empty()) {
        printPath(longestPath);
    }

    //4.6 Multiple Longest Common Subsequence
    /*
    std::string v, w, u;
    std::cin >> v >> w >> u;
    
    AlignmentResult result = LCS3D(v, w, u);
    
    std::cout << result.score << std::endl;
    std::cout << result.alignedV << std::endl;
    std::cout << result.alignedW << std::endl;
    std::cout << result.alignedU << std::endl;
    */

    //4.7 Edit Distance
    /*
    std::string v, w;
    std::cin >> v >> w;
    
    int result = editDistance(v, w);
    std::cout << result << std::endl;*/
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
