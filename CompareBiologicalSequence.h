#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <string>

enum ARROW {GREEN, BLUE, RED};
std::vector<int> ReadVectorFromConsol() {
	std::string input;
	std::vector<int> numbers;
	std::getline(std::cin, input);

	std::stringstream ss(input);
	std::string item;

	while (std::getline(ss, item, ',')) {
        int num = std::stoi(item);
        numbers.push_back(num);
	}

	return numbers;
}
int MoneyChangeProblem(int money, std::vector<int> coins) {
	std::vector<int> MinNumCoins(money + 1, INT_MAX);
	MinNumCoins[0] = 0;
	for (int m = 1; m <= money; ++m) {
		for (int i = 0; i < coins.size(); ++i) {
			if (m >= coins[i]) {
				if ((MinNumCoins[m - coins[i]] + 1) < MinNumCoins[m]) {
					MinNumCoins[m] = MinNumCoins[m - coins[i]] + 1;
				}
			}
		}
	}
	return MinNumCoins[money];
}

size_t ManhattanTouristProblem(size_t n, size_t m, const std::vector<std::vector<size_t>> &Down, const std::vector<std::vector<size_t>>& Right) {
	std::vector<std::vector<size_t>> S(n + 1);
	for (size_t i = 0; i < n + 1; ++i)
		S[i].resize(m + 1);
	S[0][0] = 0;

	for (size_t i = 1; i <= n; ++i)
		S[i][0] = S[i - 1][0] + Down[i - 1][0];

	for (size_t j = 1; j <= m; ++j)
		S[0][j] = S[0][j - 1] + Right[0][j - 1];

	for (size_t i = 1; i <= n; ++i)
		for (size_t j = 1; j <= m; ++j)
			S[i][j] = std::max(S[i - 1][j] + Down[i - 1][j], S[i][j - 1] + Right[i][j - 1]);

	return S[n][m];
}

size_t MaxFromThree(size_t first, size_t second, size_t third) {
	return std::max(first, std::max(second, third));
}

std::vector<std::vector<ARROW>> LCSBackTrack(const std::string &v, const std::string &w) {
	size_t n = v.length();
	size_t m = w.length();

	std::vector<std::vector<ARROW>> BackTrack(n , std::vector<ARROW>(m));

	std::vector<std::vector<size_t>> S(n, std::vector<size_t>(m, 0));

	for (size_t i = 1; i < n; ++i) {
		for (size_t j = 1; j < m; ++j) {
			S[i][j] = std::max({ S[i - 1][j], S[i][j - 1], S[i - 1][j - 1] + (v[i] == w[j] ? 1 : 0) });
			if (S[i][j] == S[i - 1][j]) {
				BackTrack[i][j] = GREEN;
			}
			else if (S[i][j] == S[i][j - 1]) {
				BackTrack[i][j] = BLUE;
			}
			else if (S[i][j] == S[i - 1][j - 1] + 1 && v[i] == w[j]) {
				BackTrack[i][j] = RED;
			}
		}
	}
	return BackTrack;
}

void OutPutLSC(std::vector<std::vector<ARROW>> BackTrack, std::string v, size_t lenV, size_t lenW, std::string &answer) {
	if ((lenV == 0) || (lenW == 0))
		return;
	if (BackTrack[lenV][lenW] == GREEN)
		return OutPutLSC(BackTrack, v, lenV - 1, lenW, answer);
	else if (BackTrack[lenV][lenW] == BLUE)
		return OutPutLSC(BackTrack, v, lenV, lenW - 1, answer);
	else {
		//std::cout << v[lenV];
		answer.push_back(v[lenV]);
		return OutPutLSC(BackTrack, v, lenV - 1, lenW - 1, answer);
	}

}

//std::string OutPutLSC(const std::vector<std::vector<ARROW>>& BackTrack, const std::string &v) {
//	std::string answer;
//	size_t i = BackTrack.size() - 1;
//	size_t j = BackTrack[0].size() - 1;
//
//	while (i > 0 && j > 0) {
//		if (BackTrack[i][j] == GREEN) {
//			--i;
//		}
//		else if (BackTrack[i][j] == BLUE) {
//			--j;
//		}
//		else if (BackTrack[i][j] == RED) {
//			answer += v[i];
//			--i;
//			--j;
//		}
//	}
//
//	std::reverse(answer.begin(), answer.end());
//	return answer;
//}

void readGraph(std::unordered_map<int, std::vector<int>>& graph, std::unordered_map<int, int>& degreeEachNode) {
	std::string line;
	while (std::getline(std::cin, line) && !line.empty()) {
		std::istringstream iss(line);
		int node;
		char arrow;
		iss >> node >> arrow >> arrow;

		std::vector<int> neighbors;
		int neighbor;
		while (iss >> neighbor) {
			neighbors.push_back(neighbor);
			if (iss.peek() == ',') iss.ignore();
		}

		graph[node] = neighbors;
		for (int currentNeighbor : neighbors) {
			degreeEachNode[currentNeighbor]++;
			if (degreeEachNode.find(node) == degreeEachNode.end()) {
				degreeEachNode[node] = 0;
			}
		}
	}
}

std::vector<int> topologicalOrdering(const std::unordered_map<int, std::vector<int>>& graph, std::unordered_map<int, int>& degreeEachNode) {
	std::queue<int> candidates;
	std::vector<int> list;

	for (const auto& pair : degreeEachNode) {
		if (pair.second == 0) {
			candidates.push(pair.first);
		}
	}

	while (!candidates.empty()) {
		int node = candidates.front();
		candidates.pop();
		list.push_back(node);

		if (graph.find(node) != graph.end()) {
			for (int neighbor : graph.at(node)) {
				degreeEachNode[neighbor]--;
				if (degreeEachNode[neighbor] == 0) {
					candidates.push(neighbor);
				}
			}
		}
	}

	if (list.size() != degreeEachNode.size()) {
		std::cerr << "The input graph is not a DAG graph" << std::endl;
		return {};
	}

	return list;
}

struct Edge {
	int to;
	int weight;
};

void readWeightGraph(std::unordered_map<int, std::vector<Edge>>& graph, std::unordered_map<int, int>& degreeEachNode) {
	std::string line;
	std::getline(std::cin, line);
	while (std::getline(std::cin, line) && !line.empty()) {
		std::istringstream iss(line);
		int from, to, weight;
		char arrow, colon;

		iss >> from >> arrow >> arrow >> to >> colon >> weight;
		graph[from].push_back({ to, weight });
		degreeEachNode[to]++;
		if (degreeEachNode.find(from) == degreeEachNode.end()) {
			degreeEachNode[from] = 0;
		}
	}
}

std::vector<int> topologicalOrderingWeight(const std::unordered_map<int, std::vector<Edge>>& graph, std::unordered_map<int, int>& degreeEachNode) {
	std::queue<int> candidates;
	std::vector<int> list;

	for (const auto& pair : degreeEachNode) {
		if (pair.second == 0) {
			candidates.push(pair.first);
		}
	}

	while (!candidates.empty()) {
		int node = candidates.front();
		candidates.pop();
		list.push_back(node);

		if (graph.find(node) != graph.end()) {
			for (const auto& edge : graph.at(node)) {
				degreeEachNode[edge.to]--;
				if (degreeEachNode[edge.to] == 0) {
					candidates.push(edge.to);
				}
			}
		}
	}

	return list;
}

std::pair<int, std::vector<int>> findLongestPath(const std::unordered_map<int, std::vector<Edge>>& graph, int source, int sink, const std::vector<int>& topo_order) {
	std::unordered_map<int, int> S;
	std::unordered_map<int, int> previous;

	for (int node : topo_order) {
		S[node] = std::numeric_limits<int>::min();
	}
	S[source] = 0;

	for (int node : topo_order) {
		if (S[node] != std::numeric_limits<int>::min()) {
			if (graph.find(node) != graph.end()) {
				for (const auto& edge : graph.at(node)) {
					int newDist = S[node] + edge.weight;
					if (newDist > S[edge.to]) {
						S[edge.to] = newDist;
						previous[edge.to] = node;
					}
				}
			}
		}
	}

	int lenLongestPath = S[sink];
	std::vector<int> path;
	if (lenLongestPath == std::numeric_limits<int>::min()) {
		return { 0, path };
	}

	for (int at = sink; at != source; at = previous[at]) {
		path.push_back(at);
	}
	path.push_back(source);
	std::reverse(path.begin(), path.end());

	return { lenLongestPath, path };
}

void printPath(const std::vector<int>& path) {
	for (size_t i = 0; i < path.size(); ++i) {
		std::cout << path[i];
		if (i != path.size() - 1) {
			std::cout << "->";
		}
	}
	std::cout << std::endl;
}

struct AlignmentResult {
	int score;
	std::string alignedV;
	std::string alignedW;
	std::string alignedU;
};

AlignmentResult LCS3D(const std::string& v, const std::string& w, const std::string& u) {
	size_t n = v.length();
	size_t m = w.length();
	size_t o = u.length();

	std::vector<std::vector<std::vector<int>>> S(n + 1, std::vector<std::vector<int>>(m + 1, std::vector<int>(o + 1, 0)));

	for (size_t i = 1; i <= n; ++i) {
		for (size_t j = 1; j <= m; ++j) {
			for (size_t k = 1; k <= o; ++k) {
				int matchScore = (v[i - 1] == w[j - 1] && w[j - 1] == u[k - 1]) ? 1 : 0;

				S[i][j][k] = std::max({
					S[i - 1][j][k],
					S[i][j - 1][k],
					S[i][j][k - 1],
					S[i - 1][j - 1][k],
					S[i - 1][j][k - 1],
					S[i][j - 1][k - 1],
					S[i - 1][j - 1][k - 1] + matchScore
					});
			}
		}
	}

	std::string alignedV, alignedW, alignedU;
	size_t i = n, j = m, k = o;
	while (i > 0 || j > 0 || k > 0) {
		if (i > 0 && S[i][j][k] == S[i - 1][j][k]) {
			alignedV.push_back(v[i - 1]);
			alignedW.push_back('-');
			alignedU.push_back('-');
			--i;
		}
		else if (j > 0 && S[i][j][k] == S[i][j - 1][k]) {
			alignedV.push_back('-');
			alignedW.push_back(w[j - 1]);
			alignedU.push_back('-');
			--j;
		}
		else if (k > 0 && S[i][j][k] == S[i][j][k - 1]) {
			alignedV.push_back('-');
			alignedW.push_back('-');
			alignedU.push_back(u[k - 1]);
			--k;
		}
		else if (i > 0 && j > 0 && S[i][j][k] == S[i - 1][j - 1][k]) {
			alignedV.push_back(v[i - 1]);
			alignedW.push_back(w[j - 1]);
			alignedU.push_back('-');
			--i; --j;
		}
		else if (i > 0 && k > 0 && S[i][j][k] == S[i - 1][j][k - 1]) {
			alignedV.push_back(v[i - 1]);
			alignedW.push_back('-');
			alignedU.push_back(u[k - 1]);
			--i; --k;
		}
		else if (j > 0 && k > 0 && S[i][j][k] == S[i][j - 1][k - 1]) {
			alignedV.push_back('-');
			alignedW.push_back(w[j - 1]);
			alignedU.push_back(u[k - 1]);
			--j; --k;
		}
		else if (i > 0 && j > 0 && k > 0 && S[i][j][k] == S[i - 1][j - 1][k - 1] + 1 && v[i - 1] == w[j - 1] && w[j - 1] == u[k - 1]) {
			alignedV.push_back(v[i - 1]);
			alignedW.push_back(w[j - 1]);
			alignedU.push_back(u[k - 1]);
			--i; --j; --k;
		}
	}

	std::reverse(alignedV.begin(), alignedV.end());
	std::reverse(alignedW.begin(), alignedW.end());
	std::reverse(alignedU.begin(), alignedU.end());

	return { S[n][m][o], alignedV, alignedW, alignedU };
}

int editDistance(const std::string& v, const std::string& w) {
	size_t n = v.length();
	size_t m = w.length();

	std::vector<std::vector<int>> S(n + 1, std::vector<int>(m + 1, 0));

	for (size_t i = 0; i <= n; ++i) S[i][0] = i;
	for (size_t j = 0; j <= m; ++j) S[0][j] = j;

	for (size_t i = 1; i <= n; ++i) {
		for (size_t j = 1; j <= m; ++j) {
			if (v[i - 1] == w[j - 1]) {
				S[i][j] = S[i - 1][j - 1];
			}
			else {
				S[i][j] = std::min({
					S[i - 1][j] + 1,
					S[i][j - 1] + 1,
					S[i - 1][j - 1] + 1
					});
			}
		}
	}

	return S[n][m];
}
