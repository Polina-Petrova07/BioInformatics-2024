#pragma once

#include <algorithm>
#include <vector>
#include <string>

#include "Introdaction.h"

std::vector<char> Nucleotide{ 'A', 'T', 'G', 'C' };
std::unordered_map<char, int> map = { {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3} };

using Matrix = std::vector<std::vector<double>>;

std::vector<std::string> ImmediateNeighbors(std::string Pattern) {
	std::vector<std::string> Neighborhood;
	Neighborhood.push_back(Pattern);
	for (size_t i = 0; i < Pattern.length(); ++i) {
		auto symbol = Pattern[i];
		for (char sym : Nucleotide) {
			if (sym != symbol) {
				std::string neighbor = Pattern;
				neighbor[i] = sym;
				Neighborhood.push_back(neighbor);
			}
		}
	}
	return Neighborhood;
}

std::vector<std::string> IterativeNeighbors(std::string Pattern, size_t d) {
    std::vector<std::string> Neighborhood;
    Neighborhood.push_back(Pattern);

    for (size_t j = 0; j < d; ++j) {
        std::vector<std::string> newNeighbors;

        for (const std::string& ss : Neighborhood) {
            std::vector<std::string> tmp = ImmediateNeighbors(ss);
            newNeighbors.insert(newNeighbors.end(), tmp.begin(), tmp.end());
        }

        Neighborhood.insert(Neighborhood.end(), newNeighbors.begin(), newNeighbors.end());
        std::sort(Neighborhood.begin(), Neighborhood.end());
        Neighborhood.erase(std::unique(Neighborhood.begin(), Neighborhood.end()), Neighborhood.end());
    }
    return Neighborhood;
}
 
std::vector<std::string> MotifEnumeration(std::vector<std::string> Dna, unsigned int k, unsigned int d) {
	std::vector<std::string> Patterns;
	for (size_t i = 0; i < Dna[0].length(); i++) {
		std::string pattern = Dna[0].substr(i, k);
		std::vector<std::string> Neighborhood = IterativeNeighbors(pattern, d);
		for (std::string pattern_ : Neighborhood) {
			bool is_motif = true;
			for (std::string sub_Dna : Dna) {
				bool found = false;
				for (size_t j = 0; j <= sub_Dna.length() - k; j++) {
					std::string sub_kmer = sub_Dna.substr(j, k);
					if (hammingDistance(pattern_, sub_kmer) <= d) {
						found = true;
						break;
					}
				}
				if (!found) {
					is_motif = false;
					break;
				}
			}
			if (is_motif)
				Patterns.push_back(pattern_);
		}
	}
	std::sort(Patterns.begin(), Patterns.end());
	Patterns.erase(std::unique(Patterns.begin(), Patterns.end()), Patterns.end());

	return Patterns;
}

void GenerateWords(std::string current_word, size_t k, std::vector<std::string>& answer) {
	if (current_word.length() == k) {
		answer.push_back(current_word);
		return;
	}
	for (const char& letter : Nucleotide) {
		GenerateWords(current_word + letter, k, answer);
	}
}

//Find best position in DNA by minimaze Hamming distance
unsigned int FindBestPosition(std::vector<std::string> DNA, std::string pattern) {
	unsigned int bestHammingDistance = 0;
	size_t patternSize = pattern.length();
	for (std::string string_ : DNA) {
		unsigned int min = INT_MAX;
		for (size_t i = 0; i < string_.length(); ++i) {
			std::string sub_strDNA = string_.substr(i, patternSize);
			if (hammingDistance(sub_strDNA, pattern) < min)
				min = hammingDistance(sub_strDNA, pattern);
		}
		bestHammingDistance += min;
	}
	return bestHammingDistance;
}

std::string MedianString(std::vector<std::string> DNA, unsigned int k) {
	std::string median;
	unsigned int distance = INT_MAX;
	std::vector<std::string> allKMers;
	GenerateWords("", k, allKMers);
	std::sort(allKMers.begin(), allKMers.end());
	for (std::string kMer : allKMers) {
		if (distance > FindBestPosition(DNA, kMer)) {
			distance = FindBestPosition(DNA, kMer);
			median = kMer;
		}
	}
	return median;
}

Matrix ReadMatrix(int k) {
	Matrix mat(4, std::vector<double>(k));

	for (size_t i = 0; i < 4; ++i)
		for (size_t j = 0; j < k; ++j)
			std::cin >> mat[i][j];

	return mat;
}

std::string ProfileMostProbableKMer(const std::string &text, int k, const Matrix &mat) {
	size_t lenText = text.length();
	std::string answer;
	double max_prob = 0.0;
	for (size_t i = 0; i <= lenText - k; ++i) {
		std::string subStr = text.substr(i, k);
		double probability = 1.0;
		for (size_t j = 0; j < k; ++j) {
			auto sym = subStr[j];
			switch (sym)
			{
			case 'A':
				probability *= mat[0][j];
				break;
			case 'C':
				probability *= mat[1][j];
				break;
			case 'G':
				probability *= mat[2][j];
				break;
			case 'T':
				probability *= mat[3][j];
				break;
			default:
				break;
			}
		}
		if (probability > max_prob) {
			max_prob = probability;
			answer = subStr;
		}
	}
	return answer;
}
int Score(const std::vector<std::string>& MotifMatrix) {
	int score = 0;
	int len = MotifMatrix[0].size();

	for (int col = 0; col < len; ++col) {
		std::unordered_map<char, int> dic = { {'A', 0}, {'C', 0}, {'G', 0}, {'T', 0} };
		char key_max = 'A';
		int max = 0;

		for (std::string row : MotifMatrix) {
			dic[row[col]]++;
			if (dic[row[col]] > max) {
				max = dic[row[col]];
				key_max = row[col];
			}
		}

		for (std::string row : MotifMatrix) {
			if (row[col] != key_max) score++;
		}
	}

	return score;
}

std::string GetPopularMotifFromDNA(const std::vector<std::string>& DNA, int k, const Matrix& Profil, int i) {
	double p_max = 0.0;
	std::string res = DNA[i].substr(0, k);

	for (int j = 0; j <= DNA[0].size() - k; ++j) {
		std::string str = DNA[i].substr(j, k);
		double p = 1.0;

		for (int ind_elem = 0; ind_elem < k; ++ind_elem) {
			int ind_in_map = map[str[ind_elem]];
			p *= Profil[ind_in_map][ind_elem];
		}

		if (p > p_max) {
			p_max = p;
			res = str;
		}
	}
	return res;
}

std::vector<std::string> GreedyMotifSearch(const std::vector<std::string>& DNAs, int k, int t) {
	std::vector<std::string> BestMotif(t);
	Matrix Profil(4, std::vector<double>(k, 0.0));

	for (int i = 0; i < t; ++i) {
		BestMotif[i] = DNAs[i].substr(0, k);
	}

	for (int el = 0; el <= DNAs[0].size() - k; ++el) {
		std::string dna1 = DNAs[0].substr(el, k);
		std::fill(Profil.begin(), Profil.end(), std::vector<double>(k, 0.0));

		for (int i = 0; i < k; ++i) {
			Profil[map[dna1[i]]][i] = 1;
		}

		std::vector<std::string> MotifMatrix = { dna1 };

		for (int i = 1; i < t; ++i) {
			std::string MostPopularMotif = GetPopularMotifFromDNA(DNAs, k, Profil, i);
			MotifMatrix.push_back(MostPopularMotif);

			for (int c = 0; c < k; ++c) {
				int ind_in_map = map[MostPopularMotif[c]];

				for (int nuc = 0; nuc < 4; ++nuc) {
					if (nuc == ind_in_map) {
						Profil[nuc][c] = (Profil[nuc][c] * i + 1) / (i + 1);
					}
					else {
						Profil[nuc][c] = (Profil[nuc][c] * i) / (i + 1);
					}
				}
			}
		}

		if (Score(MotifMatrix) < Score(BestMotif)) {
			BestMotif = MotifMatrix;
		}
	}

	return BestMotif;
}

std::vector<std::string> GreedyMotifSearchLaplass(const std::vector<std::string>& DNA, int k, int t) {
	std::vector<std::string> BestMotif(t);
	Matrix Profil(4, std::vector<double>(k, 0.0));

	for (int i = 0; i < t; ++i) {
		BestMotif[i] = DNA[i].substr(0, k);
	}

	for (int el = 0; el <= DNA[0].size() - k; ++el) {
		std::string dna1 = DNA[0].substr(el, k);

		for (int i = 0; i < k; ++i)
			for (int j = 0; j < 4; ++j)
				Profil[j][i] = 1;

		for (int i = 0; i < k; ++i) {
			int ind = map[dna1[i]];
			Profil[ind][i] += 1;
		}

		for (int i = 0; i < k; ++i)
			for (int j = 0; j < 4; ++j)
				Profil[j][i] /= 5.0;

		std::vector<std::string> MotifMatrix = { dna1 };
		for (int i = 1; i < t; ++i) {
			std::string MostPopularMotif = GetPopularMotifFromDNA(DNA, k, Profil, i);
			MotifMatrix.push_back(MostPopularMotif);

			for (int c = 0; c < k; ++c) {
				int ind_in_map = map[MostPopularMotif[c]];
				for (int nuc = 0; nuc < 4; ++nuc) {
					if (nuc == ind_in_map) {
						Profil[nuc][c] = (Profil[nuc][c] * (i + 4) + 1) / (i + 5);
					}
					else {
						Profil[nuc][c] = (Profil[nuc][c] * (i + 4)) / (i + 5);
					}
				}
			}
		}

		if (Score(MotifMatrix) < Score(BestMotif)) {
			BestMotif = MotifMatrix;
		}
	}

	return BestMotif;
}