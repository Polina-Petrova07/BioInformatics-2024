// bioinformatics.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include "Introdaction.h"
#include "PeptideSequencing.h"


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
    
    size_t k, t;
    std::cin >> k >> t;
    
    std::vector<std::string> DNA(t);
    for (int i = 0; i < t; ++i) {
        std::cin >> DNA[i];
    }
    
    std::vector<std::string> BestMotif = GreedyMotifSearchLaplass(DNA, k, t);
    
    for (std::string motif : BestMotif) {
        std::cout << motif << std::endl;
    }
    
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
