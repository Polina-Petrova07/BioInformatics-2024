#pragma once
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include "Introdaction.h"

std::map<std::string, char> RNA_codon_table{
	{"AAA", 'K'}, {"AAC", 'N'}, {"AAG", 'K'}, {"AAU", 'N'}, {"ACA", 'T'},
	{"ACC", 'T'}, {"ACG", 'T'}, {"ACU", 'T'}, {"AGA", 'R'}, {"AGC", 'S'},
	{"AGG", 'R'}, {"AGU", 'S'}, {"AUA", 'I'}, {"AUC", 'I'}, {"AUG", 'M'},
	{"AUU", 'I'}, {"CAA", 'Q'}, {"CAC", 'H'}, {"CAG", 'Q'}, {"CAU", 'H'},
	{"CCA", 'P'}, {"CCC", 'P'}, {"CCG", 'P'}, {"CCU", 'P'}, {"CGA", 'R'},
	{"CGC", 'R'}, {"CGG", 'R'}, {"CGU", 'R'}, {"CUA", 'L'}, {"CUC", 'L'},
	{"CUG", 'L'}, {"CUU", 'L'}, {"GAA", 'E'}, {"GAC", 'D'}, {"GAG", 'E'},
	{"GAU", 'D'}, {"GCA", 'A'}, {"GCC", 'A'}, {"GCG", 'A'}, {"GCU", 'A'},
	{"GGA", 'G'}, {"GGC", 'G'}, {"GGG", 'G'}, {"GGU", 'G'}, {"GUA", 'V'},
	{"GUC", 'V'}, {"GUG", 'V'}, {"GUU", 'V'}, {"UAA", ' '}, {"UAC", 'Y'},
	{"UAG", ' '}, {"UAU", 'Y'}, {"UCA", 'S'}, {"UCC", 'S'}, {"UCG", 'S'},
	{"UCU", 'S'}, {"UGA", ' '}, {"UGC", 'C'}, {"UGG", 'W'}, {"UGU", 'C'},
	{"UUA", 'L'}, {"UUC", 'F'}, {"UUG", 'L'}, {"UUU", 'F'}
};

std::map<char, int> TheoreticalSpectrum{
	{'G', 57},{'A', 71},{'S', 87},{'P', 97},{'V', 99},
	{'T', 101},{'C', 103},{'I', 113},{'L', 113},{'N', 114},
	{'D', 115},{'K',128},{'Q', 128},{'E', 129},{'M', 131},
	{'H', 137},{'F', 147},{'R', 156},{'Y', 163},{'W', 186}
};

std::string proteinTranslationProblem(std::string RNA_pattern) {
	std::string protein;
	// here was genetic code table
	for (size_t i = 0; i < RNA_pattern.length();) {
		std::string codon = RNA_pattern.substr(i, 3);
		if (RNA_codon_table.contains(codon))
			protein.push_back(RNA_codon_table[codon]);
		i += 3;
	}
	return protein;
}
std::string dnaToRna(std::string DNA) {
	std::string RNA = DNA;
	for (size_t i = 0; i < RNA.length(); ++i)
		if (RNA[i] == 'T')
			RNA[i] = 'U';

	return RNA;
}
std::vector<std::string> peptideEncodingProblem(std::string DNA, std::string geneticCode) {
	std::vector<std::string> answer;
	size_t lenghtSubStr = geneticCode.length() * 3;
	for (size_t i = 0; i < DNA.length(); i++) {
		std::string currentSubString_ = DNA.substr(i, lenghtSubStr);
		std::string currentSubString = dnaToRna(currentSubString_);
		std::string currentProtein = proteinTranslationProblem(currentSubString);
		if (currentProtein == geneticCode)
			answer.push_back(currentSubString_);
		else {
			std::string reversSubStrDNA = reverseComplement(currentSubString_);
			std::string reversSubStrRNA = dnaToRna(reversSubStrDNA);
			std::string currentProteinForRevers = proteinTranslationProblem(reversSubStrRNA);
			if (currentProteinForRevers == geneticCode)
				answer.push_back(currentSubString_);
		}
	}
	return answer;
}
long long int subpeptidesCountProblem(long long int n) {
	return (n * (n - 1));
}
std::vector<std::string> findAllSubStrings(const std::string& s) {
	std::string ss = s + s;
	std::vector<std::string> allSubStrings;
	for (size_t i = 0; i < s.length(); i++) {
		for (size_t len = 1; len < s.length(); len++) {
			std::string substr = ss.substr(i, len);
			allSubStrings.push_back(substr);
		}
	}
	return allSubStrings;
}
std::vector<int> generatingTheoreticalSpectrumProblem(std::string peptide) {
	std::vector<int> cyclospectrum;
	std::vector<std::string> allSubStrings = findAllSubStrings(peptide);
	for (std::string s : allSubStrings) {
		int currentValue = 0;
		for (size_t j = 0; j < s.length(); j++) {
			currentValue += TheoreticalSpectrum[s[j]];
		}
		cyclospectrum.push_back(currentValue);
	}
	cyclospectrum.push_back(0);
	int currentValue = 0;
	for (size_t i = 0; i < peptide.length(); i++) {
		currentValue += TheoreticalSpectrum[peptide[i]];
	}
	cyclospectrum.push_back(currentValue);
	std::sort(cyclospectrum.begin(), cyclospectrum.end());
	return cyclospectrum;
}
int countingPeptidesWithGivenMassProblem(int n) { // wrong answer
	std::vector<int> masses = { 57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186 };
	std::vector<int> count(static_cast<std::vector<int, std::allocator<int>>::size_type>(n + 1), 0);
	count[0] = 1;
	for (int i = 57; i <= n; ++i)
		for (int m : masses)
			if ((i - m) >= 0)
				count[i] += count[static_cast<std::vector<int, std::allocator<int>>::size_type>(i) - m];
	return count[n];
}
long long int countingSpectrumOfTheLinearPeptideProblem(long long int n) {
	return (((1 + n) * n) / 2 + 1);
}