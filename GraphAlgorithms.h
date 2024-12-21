#pragma once
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <iterator>
#include <deque>
#include <unordered_map>

std::vector<std::string>Composition(const std::string& text, int k) {
	std::vector<std::string> answer;
	for (size_t i = 0; i <= text.length() - k; ++i) {
		std::string tmp = text.substr(i, k);
		answer.push_back(tmp);
	}
	std::sort(answer.begin(), answer.end());
	return answer;
}

std::string StringSpelled(const std::vector<std::string>& sequensePatterns) {
	size_t k = sequensePatterns[0].length();
	std::string answer = sequensePatterns[0];
	for (size_t i = 1; i < sequensePatterns.size(); ++i) {
		answer += sequensePatterns[i][k - 1];
	}
	return answer;
}

std::string GetPrefix(const std::string& str) {
	size_t k = str.length();
	return (str.substr(0, k - 1));
}

std::string GetSuffix(const std::string& str) {
	size_t k = str.length();
	return (str.substr(1, k - 1));
}

void printMapGraph(std::map<std::string, std::string> graph) {
	for (const auto& [node, neighbour] : graph) {
		if (neighbour != "NULL") {
			std::cout << node << " -> " << neighbour << std::endl;
		} else {
			continue;
		}
	}
}
std::map<std::string, std::string> OverlapGraph(const std::vector<std::string>& nodes) {
	std::map<std::string, std::string> answer;

	for (size_t i = 0; i < nodes.size();++i) {
		std::string currentNode = nodes[i];
		for (size_t j = 0; j < nodes.size(); ++j) {
			if (GetSuffix(currentNode) == GetPrefix(nodes[j])) {
				answer[currentNode] = nodes[j];
			} else {
				answer.insert({ currentNode, "NULL"});
			}
		}
	}
	return answer;
}

std::vector<std::string> generateUniqueKmers(const std::vector<std::string>& nodes) {
	std::vector<std::string> uniqueKmers;
	for (const auto& str : nodes) {
		std::string prefix = GetPrefix(str);
		std::string suffix = GetSuffix(str);
		if (std::find(uniqueKmers.begin(), uniqueKmers.end(), prefix) == uniqueKmers.end()) {
			uniqueKmers.push_back(prefix);
		}
		if (std::find(uniqueKmers.begin(), uniqueKmers.end(), suffix) == uniqueKmers.end()) {
			uniqueKmers.push_back(suffix);
		}
	}
	return uniqueKmers;
}

int countOccurrences(const std::vector<std::string>& nodes, const std::string& subStr) {
	return std::count(nodes.begin(), nodes.end(), subStr);
}

void DeBruijnGraph (const std::vector<std::string>& nodes) {
	std::vector<std::string> uniqueKmers  = generateUniqueKmers(nodes);
	for (const auto& kmer : uniqueKmers) {
		std::string tmp = kmer.substr(1);
		std::vector<std::string> answer;

		for (const auto& uniqKmer : uniqueKmers) {
			std::string kmerOgir = kmer + uniqKmer.back();
			int count = countOccurrences(nodes, kmerOgir);
			if (tmp == uniqKmer.substr(0, uniqKmer.size() - 1) && count != 0) {
				for (int j = 0; j < count; ++j) {
					answer.push_back(uniqKmer);
				}
			}
		}

		if (!answer.empty()) {
			std::cout << kmer << " -> ";
			for (size_t i = 0; i < answer.size(); ++i) {
				std::cout << answer[i];
				if (i < answer.size() - 1) std::cout << ",";
			}
			std::cout << std::endl;
		}
	}
}

void readGraphQ(std::unordered_map<int, std::deque<int>>& graph) {
	std::string line;
	while (std::getline(std::cin, line) && !line.empty()) {
		std::istringstream iss(line);
		int node;
		char arrow;
		iss >> node >> arrow >> arrow;

		std::deque<int> neighbors;
		int neighbor;
		while (iss >> neighbor) {
			neighbors.push_back(neighbor);
			if (iss.peek() == ',') iss.ignore();
		}

		graph[node] = neighbors;
	}
}
int getRandomNode(const std::unordered_map<int, std::deque<int>>& graph) {
	auto it = graph.begin();
	std::advance(it, rand() % graph.size());
	return it->first;
}

std::vector<int> findEulerianCycle(std::unordered_map<int, std::deque<int>>& graph) {
	std::vector<int> tmpArray;
	tmpArray.push_back(getRandomNode(graph));
	std::vector<int> answer;

	while (!tmpArray.empty()) {
		int curr = tmpArray.back();
		if (!graph[curr].empty()) {
			int next = graph[curr].front();
			graph[curr].pop_front();
			tmpArray.push_back(next);
		}
		else {
			answer.push_back(tmpArray.back());
			tmpArray.pop_back();
		}
	}

	std::reverse(answer.begin(), answer.end());

	return answer;
}

void printEulerianCycle (std::vector<int> answer) {
    for (size_t i = 0; i < answer.size(); ++i) {
    	std::cout << answer[i];
    	if (i < answer.size() - 1) {
    		std::cout << "->";
    	}
    }
}

std::vector<int> findEulerianPath(std::unordered_map<int, std::deque<int>>& graph) {
	std::unordered_map<int, int> degree;
	for (const auto& [node, neighbors] : graph) {
		degree[node] += neighbors.size();
		for (int neighbor : neighbors) {
			degree[neighbor]--;
		}
	}

	int unbalancedStart = -1, unbalancedEnd = -1;
	for (const auto& [node, d] : degree) {
		if (d == 1) {
			unbalancedStart = node;
		}
		else if (d == -1) {
			unbalancedEnd = node;
		}
	}

	if (unbalancedEnd != -1 && unbalancedStart != -1) {
		graph[unbalancedEnd].push_back(unbalancedStart);
	}

	std::vector<int> path = findEulerianCycle(graph);

	if (unbalancedEnd != -1 && unbalancedStart != -1) {
		for (size_t i = 0; i < path.size() - 1; ++i) {
			if (path[i] == unbalancedEnd && path[i + 1] == unbalancedStart) {
				std::vector<int> eulerianPath(path.begin() + i + 1, path.end());
				eulerianPath.insert(eulerianPath.end(), path.begin() + 1, path.begin() + i + 1);
				return eulerianPath;
			}
		}
	}

	return path;
}

void DeBruijnGraphNew(const std::vector<std::string>& nodes, std::unordered_map<std::string, std::deque<std::string>>& graph) {
	for (const auto& node : nodes) {
		std::string prefix = GetPrefix(node);
		std::string suffix = GetSuffix(node);
		graph[prefix].push_back(suffix);
	}
}

std::string getRandomNodeString(const std::unordered_map<std::string, std::deque<std::string>>& graph) {
	auto it = graph.begin();
	std::advance(it, rand() % graph.size());
	return it->first;
}

std::vector<std::string> findEulerianCycleString(std::unordered_map<std::string, std::deque<std::string>>& graph) {
	std::vector<std::string> tmpArray;
	tmpArray.push_back(getRandomNodeString(graph));
	std::vector<std::string> answer;

	while (!tmpArray.empty()) {
		std::string curr = tmpArray.back();
		if (!graph[curr].empty()) {
			std::string next = graph[curr].front();
			graph[curr].pop_front();
			tmpArray.push_back(next);
		}
		else {
			answer.push_back(curr);
			tmpArray.pop_back();
		}
	}

	std::reverse(answer.begin(), answer.end());

	return answer;
}

std::vector<std::string> findEulerianPathString(std::unordered_map<std::string, std::deque<std::string>>& graph) {
	std::unordered_map<std::string, int> degree;
	for (const auto& [node, neighbors] : graph) {
		degree[node] += neighbors.size();
		for (std::string neighbor : neighbors) {
			degree[neighbor]--;
		}
	}

	std::string unbalancedStart, unbalancedEnd;
	for (const auto& [node, d] : degree) {
		if (d == 1) {
			unbalancedStart = node;
		}
		else if (d == -1) {
			unbalancedEnd = node;
		}
	}

	if (!unbalancedEnd.empty() && !unbalancedStart.empty()) {
		graph[unbalancedEnd].push_back(unbalancedStart);
	}

	std::vector<std::string> path = findEulerianCycleString(graph);

	if (!unbalancedEnd.empty() && !unbalancedStart.empty()) {
		for (size_t i = 0; i < path.size() - 1; ++i) {
			if (path[i] == unbalancedEnd && path[i + 1] == unbalancedStart) {
				std::vector<std::string> eulerianPath(path.begin() + i + 1, path.end());
				eulerianPath.insert(eulerianPath.end(), path.begin() + 1, path.begin() + i + 1);
				return eulerianPath;
			}
		}
	}

	return path;
}
