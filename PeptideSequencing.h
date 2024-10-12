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
std::vector<int> CyclicSpectrum(const std::string& peptide) {
	std::vector<int> prefix_mass(1, 0);
	for (size_t i = 0; i < peptide.size(); ++i) {
		prefix_mass.push_back(prefix_mass[i] + TheoreticalSpectrum[peptide[i]]);
	}
	int peptide_mass = prefix_mass.back();
	std::vector<int> cyclo_spectrum = { 0 };
	for (size_t i = 0; i < peptide.size(); ++i) {
		for (size_t j = i + 1; j <= peptide.size(); ++j) {
			cyclo_spectrum.push_back(prefix_mass[j] - prefix_mass[i]);
			if (i > 0 && j < peptide.size()) {
				int cur = peptide_mass - (prefix_mass[j] - prefix_mass[i]);
				cyclo_spectrum.push_back(cur);
			}
		}
	}
	std::sort(cyclo_spectrum.begin(), cyclo_spectrum.end());
	return cyclo_spectrum;
}

std::vector<int> LinearSpectrum_1(const std::string& peptide) {
	std::vector<int> prefix_mass(1, 0);
	for (size_t i = 0; i < peptide.size(); ++i) {
		prefix_mass.push_back(prefix_mass[i] + TheoreticalSpectrum[peptide[i]]);
	}
	std::vector<int> linear_spectrum = { 0 };
	for (size_t i = 0; i < peptide.size(); ++i) {
		for (size_t j = i + 1; j <= peptide.size(); ++j) {
			linear_spectrum.push_back(prefix_mass[j] - prefix_mass[i]);
		}
	}
	std::sort(linear_spectrum.begin(), linear_spectrum.end());
	return linear_spectrum;
}

std::vector<std::string> Expand_1(const std::vector<std::string>& peptides) {
	std::vector<std::string> expanded_peptides;
	for (const auto& peptide : peptides) {
		for (const auto& key : TheoreticalSpectrum) {
			expanded_peptides.push_back(peptide + key.first);
		}
	}
	return expanded_peptides;
}

int onePeptideMass(const std::string& peptide) {
	int m = 0;
	for (char pep : peptide) {
		m += TheoreticalSpectrum[pep];
	}
	return m;
}

int ParentMass(const std::vector<int>& spectrum) {
	return spectrum.back();
}

bool Inconsistent(const std::string& peptide, const std::vector<int>& spectrum) {
	auto a = LinearSpectrum_1(peptide);
	std::multiset<int> s(spectrum.begin(), spectrum.end());
	for (int pep : a) {
		auto it = s.find(pep);
		if (it == s.end()) {
			return true;
		}
		s.erase(it);
	}
	return false;
}
std::vector<std::string> cyclopeptideSequencingProblem(const std::vector<int>& spectrum) {
	std::vector<std::string> peptides = { "" };
	std::vector<std::string> ans;

	auto es = Expand_1(peptides);

	for (const auto& p : es) {
		if (!Inconsistent(p, spectrum)) {
			peptides.push_back(p);
		}
	}

	while (!peptides.empty()) {
		peptides = Expand_1(peptides);
		for (size_t i = 0; i < peptides.size(); ++i) {
			const auto& p = peptides[i];
			if (onePeptideMass(p) == ParentMass(spectrum)) {
				if (CyclicSpectrum(p) == spectrum) {
					ans.push_back(p);
				}
				peptides.erase(peptides.begin() + i);
				--i;
			}
			else if (Inconsistent(p, spectrum)) {
				peptides.erase(peptides.begin() + i);
				--i;
			}
		}
	}

	std::set<std::string> unique_ans(ans.begin(), ans.end());
	return std::vector<std::string>(unique_ans.begin(), unique_ans.end());
}
std::map<int, int> cyclicSpectrumDict(std::string peptide) {
	int n = peptide.length();
	std::vector<int> PrefixMass(n + 1, 0);

	for (int i = 0; i < n; ++i) {
		PrefixMass[i + 1] = PrefixMass[i] + TheoreticalSpectrum[peptide[i]];
	}

	int peptideMass = PrefixMass[n];
	std::map<int, int> cSpectrumDict;
	cSpectrumDict[0] = 1;

	for (int i = 0; i < n; ++i) {
		for (int j = i + 1; j <= n; ++j) {
			int s = PrefixMass[j] - PrefixMass[i];
			cSpectrumDict[s]++;
			if (i > 0 && j < n) {
				s = peptideMass - (PrefixMass[j] - PrefixMass[i]);
				cSpectrumDict[s]++;
			}
		}
	}

	return cSpectrumDict;
}

int cyclopeptideScoringProblem(const std::string& peptide, const std::vector<int>& spectrum) {
	auto theoSpectrumDict = cyclicSpectrumDict(peptide);
	int score = 0;
	std::map<int, int> spectrumDict;

	for (int s : spectrum) {
		spectrumDict[s]++;
	}

	for (const auto& [s, v] : theoSpectrumDict) {
		int v0 = spectrumDict[s];
		score += std::min(v0, v);
	}

	return score;
}
std::vector<std::string> allMasses = { "57", "71", "87", "97", "99", "101", "103", "113", "114", "115", "128", "129", "131", "137", "147", "156", "163", "186" };

std::vector<int> SplitStringToList(const std::string& peptide, char delimiter = '-') {
	std::vector<int> result;
	std::stringstream ss(peptide);
	std::string item;
	while (getline(ss, item, delimiter)) {
		result.push_back(stoi(item));
	}
	return result;
}

int GetMass(const std::string& peptide) {
	std::vector<int> masses = SplitStringToList(peptide);
	int total = 0;
	for (int mass : masses) {
		total += mass;
	}
	return total;
}

std::map<int, int> Spectrum(const std::vector<int>& spectrum) {
	std::map<int, int> result;
	for (int mass : spectrum) {
		result[mass]++;
	}
	return result;
}

std::map<int, int> CycloSpectrumProblem(const std::string& peptide) {
	std::vector<int> peptide_list = SplitStringToList(peptide);
	int rank = peptide_list.size();
	std::vector<int> prefix(rank + 1, 0);
	for (int i = 0; i < rank; ++i) {
		prefix[i + 1] = prefix[i] + peptide_list[i];
	}
	int peptide_mass = prefix[rank];
	std::vector<int> cspectrum = { 0 };

	for (int i = 0; i < rank; ++i) {
		for (int j = i + 1; j <= rank; ++j) {
			int diff = prefix[j] - prefix[i];
			cspectrum.push_back(diff);
			if (i > 0 && j < rank) {
				cspectrum.push_back(peptide_mass - diff);
			}
		}
	}
	return Spectrum(cspectrum);
}

int CyclopeptideScoringProblem(const std::string& peptide, const std::map<int, int>& spectrum) {
	if (peptide.empty()) {
		return 0;
	}
	std::map<int, int> cyclo_spectrum = CycloSpectrumProblem(peptide);
	int score = 0;
	for (const auto& [s, v] : cyclo_spectrum) {
		int v0 = spectrum.find(s) != spectrum.end() ? spectrum.at(s) : 0;
		score += std::min(v0, v);
	}
	return score;
}

std::map<int, int> LinearSpectrum(const std::string& peptide) {
	std::vector<int> peptide_list = SplitStringToList(peptide);
	int rank = peptide_list.size();
	std::vector<int> prefix(rank + 1, 0);
	for (int i = 0; i < rank; ++i) {
		prefix[i + 1] = prefix[i] + peptide_list[i];
	}
	std::vector<int> lspectrum = { 0 };
	for (int i = 0; i < rank; ++i) {
		for (int j = i + 1; j <= rank; ++j) {
			int diff = prefix[j] - prefix[i];
			lspectrum.push_back(diff);
		}
	}
	return Spectrum(lspectrum);
}

int LinearScore(const std::string& peptide, const std::map<int, int>& spectrum) {
	if (peptide.empty()) {
		return 0;
	}
	std::map<int, int> linear_spectrum = LinearSpectrum(peptide);
	int score = 0;
	for (const auto& [mass, count] : linear_spectrum) {
		int value = spectrum.find(mass) != spectrum.end() ? spectrum.at(mass) : 0;
		score += std::min(count, value);
	}
	return score;
}

std::set<std::string> Expand(const std::set<std::string>& peptides) {
	std::set<std::string> expanded;
	for (const auto& peptide : peptides) {
		for (const auto& mass : allMasses) {
			std::string new_peptide = peptide.empty() ? mass : peptide + "-" + mass;
			expanded.insert(new_peptide);
		}
	}
	return expanded;
}

std::set<std::string> Trim(const std::set<std::string>& leaderboard, const std::map<int, int>& spectrum, int count_leaders) {
	std::vector<std::pair<std::string, int>> scores;
	for (const auto& peptide : leaderboard) {
		scores.emplace_back(peptide, LinearScore(peptide, spectrum));
	}
	sort(scores.begin(), scores.end(), [](const auto& a, const auto& b) {
		return a.second > b.second;
		});
	std::set<std::string> result;
	for (int i = 0; i < std::min((int)scores.size(), count_leaders); ++i) {
		result.insert(scores[i].first);
	}
	return result;
}

std::string LeaderboardCyclopeptideSequencingProblem(const std::vector<int>& spectrum, int count_leaders) {
	int parentMass = *max_element(spectrum.begin(), spectrum.end());
	std::map<int, int> spectrum_dict = Spectrum(spectrum);
	std::set<std::string> peptides = { "" };
	std::string max_peptide = "";
	int max_score = 0;

	while (!peptides.empty()) {
		peptides = Expand(peptides);
		std::set<std::string> removes;
		for (const auto& peptide : peptides) {
			if (GetMass(peptide) == parentMass) {
				int score = CyclopeptideScoringProblem(peptide, spectrum_dict);
				if (score > max_score) {
					max_peptide = peptide;
					max_score = score;
				}
			}
			else if (GetMass(peptide) > parentMass) {
				removes.insert(peptide);
			}
		}
		for (const auto& peptide : removes) {
			peptides.erase(peptide);
		}
		peptides = Trim(peptides, spectrum_dict, count_leaders);
	}

	return max_peptide;
}