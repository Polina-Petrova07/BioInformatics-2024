// bioinformatics.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <sstream>
#include "Introdaction.h"
#include "PeptideSequencing.h"
#include "CellularClocks.h"
#include "CompareBiologicalSequence.h"
#include "GraphAlgorithms.h"


int main()
{
    /*======= SECTION ONE =======*/

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
    /*
    int n;
    std::cin >> n;
    int ans = countingPeptidesWithGivenMassProblem(n);
    std::cout << ans;
    */

    //2.6 Counting Spectrum of the Linear Peptide Problem
    /*
    long long int n;
    std::cin >> n;
    long long int result = countingSpectrumOfTheLinearPeptideProblem(n);
    std::cout << result;
    */

    //2.7 Cyclopeptide Sequencing Problem
    /*
    std::string spectrum;
    std::getline(std::cin, spectrum);
    
    
    std::istringstream iss(spectrum);
    std::vector<int> spectrum_int((std::istream_iterator<int>(iss)), std::istream_iterator<int>());
    
    
    //spectrum_int.emplace(spectrum_int.begin(), 0);
    
    auto answerPeptides = cyclopeptideSequencingProblem(spectrum_int);
    
    for (size_t i = 0; i < answerPeptides.size(); i++) {
        for (size_t j = 0; j < answerPeptides[i].length(); j++) {
            if (answerPeptides[i][j] == 'Q')
                answerPeptides[i][j] = 'K';
            if (answerPeptides[i][j] == 'L')
                answerPeptides[i][j] = 'I';
        }
    }
    
    std::sort(answerPeptides.begin(), answerPeptides.end());
    std::vector<std::string> uniqueAnswerPeptides;
    for (const auto& element : answerPeptides) {
        if (uniqueAnswerPeptides.empty() || uniqueAnswerPeptides.back() != element) {
            uniqueAnswerPeptides.push_back(element);
        }
    }
    for (const auto& peptide : uniqueAnswerPeptides) {
        for (size_t currentPeptide = 0; currentPeptide < peptide.size(); ++currentPeptide) {
            if (currentPeptide == peptide.size() - 1) {
                std::cout << TheoreticalSpectrum[peptide[currentPeptide]];
            }
            else {
                std::cout << TheoreticalSpectrum[peptide[currentPeptide]] << "-";
            }
        }
        std::cout << " ";
    }
    std::cout << std::endl;
    */

    //2.8 Cyclopeptide Scoring Problem
    /*std::string peptide;
    std::cin >> peptide;

    std::string spectrum;
    std::cin >> spectrum;
    std::getline(std::cin, spectrum);

    std::istringstream iss(spectrum);
    std::vector<int> sprctrum_Int((std::istream_iterator<int>(iss)), std::istream_iterator<int>());
    sprctrum_Int.emplace(sprctrum_Int.begin(), 0);

    int score = cyclopeptideScoringProblem(peptide, sprctrum_Int);
    std::cout << score << std::endl;
    */

    //2.9 Leaderboard Cyclopeptide Sequencing Problem
    /*
    int N;
    std::cin >> N;

    std::string spectrum;
    std::cin >> spectrum;
    std::getline(std::cin, spectrum);

    std::istringstream iss(spectrum);
    std::vector<int> sprctrum_Int((std::istream_iterator<int>(iss)), std::istream_iterator<int>());

    std::cout << LeaderboardCyclopeptideSequencingProblem(sprctrum_Int, N) << std::endl;
    */

    /*======= SECTION THREE =======*/

    //3.1 Motif Enumeration Problem
    /*
    unsigned int k;
    unsigned int d;
    std::vector<std::string> DNA;
    
    std::cin >> k >> d;
    for (size_t i = 0; i < 4; ++i) {
        std::string tmp;
        std::cin >> tmp;
        DNA.push_back(tmp);
    }

    std::vector <std::string> answer = MotifEnumeration(DNA, k, d);
    for (std::string a :answer)
        std::cout << a << " ";
    */

    //3.2 Median String Problem
    /*
    std::vector<std::string> DNA;
    unsigned int k;
    std::cin >> k;

    std::string line;

    while (std::cin >> line) {
        DNA.push_back(line);
    }

    std::string answer = MedianString(DNA, k);
    std::cout << answer;
    */

    //3.3 Profile-most Probable k-mer Problem
    /*
    std::string text;
    int k;
    std::cin >> text >> k;

    Matrix mat;
    mat = ReadMatrix(k);

    std::string answer = ProfileMostProbableKMer(text, k, mat);
    std::cout << answer;
    
    return 0;
    */

    //3.4 Greedy Motif Search
    /*
    size_t k, t;
    std::cin >> k >> t;
    
    std::vector<std::string> DNA(t);
    for (int i = 0; i < t; ++i) {
        std::cin >> DNA[i];
    }
    
    std::vector<std::string> BestMotif = GreedyMotifSearch(DNA, k, t);
    
    for (std::string motif : BestMotif) {
        std::cout << motif << std::endl;
    }
    */

    //3.5 Greedy Motif Search with pseudocounts   
    /*size_t k, t;
    std::cin >> k >> t;
    
    std::vector<std::string> DNA(t);
    for (int i = 0; i < t; ++i) {
        std::cin >> DNA[i];
    }
    
    std::vector<std::string> BestMotif = GreedyMotifSearchLaplass(DNA, k, t);
    
    for (std::string motif : BestMotif) {
        std::cout << motif << std::endl;
    }*/

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
    /*std::unordered_map<int, std::vector<int>> graph;
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
    }*/
    
    //4.5 Longest Path in a DAG
    /*int source, sink;
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
    }*/

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
    
    /*======= SECTION FIVE =======*/

    //5.1 String Composition Problem
    /*int k;
    std::string text;
    std::cin >> k >> text;
    std::vector<std::string> answer = Composition(text, k);

    for (std::string str : answer)
        std::cout << str << std::endl;*/

    //5.2 String Spelled by a Genome Path Problem
    /*std::vector<std::string> sequencePatterns;
    std::string line;
    while (std::cin >> line)
        sequencePatterns.push_back(line);

    std::string answer = StringSpelled(sequencePatterns);
    std::cout << answer;*/

    //5.3 Overlap Graph Problem
    /*
    std::vector<std::string> nodes;
    std::string line;

    while (std::cin >> line)
        nodes.push_back(line);

    std::map<std::string, std::string> answer = OverlapGraph(nodes);
    printVectorStringPairs(answer);*/

    //5.4 DeBruijn Graph from k-mers Problem
    /*
    std::vector<std::string> nodes;
    std::string line;
    
    while (std::cin >> line)
        nodes.push_back(line);

    DeBruijnGraph(nodes);*/

    //5.5 Eulerian Cycle Problem
    /*
    std::unordered_map<int, std::deque<int>> graph;
    readGraphQ(graph);

    std::vector<int> answer = findEulerianCycle(graph);
    printEulerianCycle(answer);*/

    //5.6 Eulerian Path Problem
    /*
    std::unordered_map<int, std::deque<int>> graph;
    readGraphQ(graph);

    std::vector<int> answer = findEulerianPath(graph);
    printEulerianCycle(answer);*/

    
    //5.7 String Reconstruction Problem
    int k;
    std::cin >> k;
    std::vector<std::string> patterns;
    std::string pattern;
    while (std::cin >> pattern)
        patterns.push_back(pattern);

    std::unordered_map<std::string, std::deque<std::string>> graph;
    DeBruijnGraphNew(patterns, graph);
    std::vector<std::string> eulerianPath = findEulerianPathString(graph);
    std::string answer = StringSpelled(eulerianPath);

    std::cout << answer;
    
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
