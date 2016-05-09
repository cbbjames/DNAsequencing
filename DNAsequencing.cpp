#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <fstream>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cstdio>
#include <vector>
#include <cmath>
#include <queue>
#include <set>
using namespace std;
class DNASequencing{
public:
	int out[1000000];
	int f[1000000];
	int g[1000000][5];
	bool first = true;
	int ref_cnt = 0;
	int ref_map[26];//id>ref
	string main_ref[25];
	int main_has[25];
	class REF {
	public:
		char** ref;
		vector<int> table[26][4];
		void ini_table(int id) {
			for (int i = 0; i < 4; i++)
				table[id][i].reserve(500000);
		}
	};
	class cREF {
	public:
		char** cref;
		vector<int> table[26][4];
		void ini_table(int id) {
			for (int i = 0; i < 4; i++)
				table[id][i].reserve(500000);
		}
	};
	vector<pair<int, int> > keypos[33];
	vector<pair<int, int> > keypos2[33];
public:
	REF Ref;
	cREF cRef;
	char cp(char c) {
		if (c == 'A') return 'T';
		else if (c == 'T') return 'A';
		else if (c == 'G')return 'C';
		else if (c == 'C') return 'G';
	}
	void ini_pass(int chromatidSequenceId) {
		if (first) { 
			first = false; 
			Ref.ref = new char*[26]; 
			cRef.cref = new char*[24];
			memset(main_has, 0, sizeof(main_has));
		}
		ref_map[chromatidSequenceId] = ref_cnt;
		Ref.ref[ref_cnt] = new char[650000000];
		cRef.cref[ref_cnt] = new char[650000000];
		Ref.ini_table(ref_cnt);
		cRef.ini_table(ref_cnt);
		ref_cnt++;
	}
	int jump(string& s,int i) {
		if (s[i] == '>') {
			while (s[i] != '\n') { i++; }
			i++;
		}
		if (s[i] == 'N') {
			while (s[i] != 'N') { i++; }
			i++;
		};
		return i;
	}
	int passReferenceGenome(int chromatidSequenceId, const vector<string>& chromatidSequence) {
		ini_pass(chromatidSequenceId);
		string s = chromatidSequence[chromatidSequenceId];
		main_ref[chromatidSequenceId]= chromatidSequence[chromatidSequenceId];
		main_has[chromatidSequenceId] = 1;
		for (int i = 0; i < s.size(); i++) {
			if (s[i] == '>' || s[i] == 'N') i=jump(s, i);
			char cps = cp(s[i]);
			Ref.ref[ref_map[chromatidSequenceId]][i] = s[i];
			cRef.cref[ref_map[chromatidSequenceId]][i] = cps;
			Ref.table[s[i] - 'A'][ref_map[chromatidSequenceId]].push_back(i);
			cRef.table[cps - 'A'][ref_map[chromatidSequenceId]].push_back(i);
		}
		return 0;
	}
	int initTest(int testDifficulty) {
		return 0;
	}
	int preProcessing() {

	}
	int makec(char c) {
		if (c == 'A') return 0;
		else if (c == 'T')return 1;
		else if (c == 'G')return 2;
		else if (c == 'C')return 3;
	}
	int buildMatchingMachine(string *words) {
		memset(out, 0, sizeof out);
		memset(f, -1, sizeof f);
		memset(g, -1, sizeof g);

		int states = 1; // Initially, we just have the 0 state

		for (int i = 0; i < 33; ++i) {
			const string &keyword = words[i];
			int currentState = 0;
			for (int j = 0; j < keyword.size(); ++j) {
				int c = makec(keyword[j]);
				if (g[currentState][c] == -1) { // Allocate a new node
					g[currentState][c] = states++;
				}
				currentState = g[currentState][c];
			}
			out[currentState] |= (1 << i); // There's a match of keywords[i] at node currentState.
		}

		// State 0 should have an outgoing edge for all characters.
		for (int c = 0; c < 5; ++c) {
			if (g[0][c] == -1) {
				g[0][c] = 0;
			}
		}

		// Now, let's build the failure function
		queue<int> q;
		for (int c = 0; c <= 3; ++c) {  // Iterate over every possible input
										// All nodes s of depth 1 have f[s] = 0
			if (g[0][c] != -1 && g[0][c] != 0) {
				f[g[0][c]] = 0;
				q.push(g[0][c]);
			}
		}
		while (q.size()) {
			int state = q.front();
			q.pop();
			for (int c = 0; c <= 3; ++c) {
				if (g[state][c] != -1) {
					int failure = f[state];
					while (g[failure][c] == -1) {
						failure = f[failure];
					}
					failure = g[failure][c];
					f[g[state][c]] = failure;
					out[g[state][c]] |= out[failure]; // Merge out values
					q.push(g[state][c]);
				}
			}
		}

		return states;
	}
	int findNextState(int currentState, char nextInput) {
		int answer = currentState;
		int c = makec(nextInput);
		while (g[answer][c] == -1) answer = f[answer];
		return g[answer][c];
	}
	void matching(string* keywords) {
		buildMatchingMachine(keywords);
		int currentState = 0;
		for (int main = 1; main < 25; main++) {
			if (!main_has[main]) continue;
			for (int i = 0; i < main_ref[main].size(); ++i) {
				currentState = findNextState(currentState, main_ref[main][i]);
				if (out[currentState] == 0) continue; // Nothing new, let's move on to the next character.
				for (int j = 0; j < 33; ++j) {
					if (out[currentState] & (1 << j)) { 
						// Matched keywords[j]
						//cout << "Keyword " << keywords[j] << " appears from "
						//<< i - keywords[j].size() + 1 << " to " << i << endl;
						pair<int, int>temp;
						temp.first = i - keywords[j].size() + 1 + 1;
						temp.second = i + 1;
						keypos[j].push_back(temp);
					}
				}
			}
		}
	}
	void matching2(string* keywords) {
		buildMatchingMachine(keywords);
		int currentState = 0;
		for (int main = 1; main < 25; main++) {
			if (!main_has[main]) continue;
			for (int i = 0; i < main_ref[main].size(); ++i) {
				currentState = findNextState(currentState, main_ref[main][i]);
				if (out[currentState] == 0) continue; // Nothing new, let's move on to the next character.
				for (int j = 0; j < 33; ++j) {
					if (out[currentState] & (1 << j)) {
						// Matched keywords[j]
						//cout << "Keyword " << keywords[j] << " appears from "
						//<< i - keywords[j].size() + 1 << " to " << i << endl;
						pair<int, int>temp;
						temp.first = i - keywords[j].size() + 1 + 1;
						temp.second = i + 1;
						keypos[j].push_back(temp);
					}
				}
			}
		}
	}
	vector<string> getAlignment(int N, double NormA, double NormS, const vector<string>& readName, const vector<string>& readSequence) {
		vector<string> ret(N, "");
		string keywords[33];
		string keywords2[33];
		string s1, s2;
		int keycnt = 0;
		for (int i = 0; i < 2*N; i += 4) {
			s1 = readSequence[i + 1];
			s2 = readSequence[i + 3];
			if (keycnt >= 32) { 
				matching(keywords);
				//implement nice 
				for (int j = 0; j < 33; j++) {
					keypos[j].clear();
					keypos2[j].clear();
				}
				keycnt = 0; 
			}
			keywords[keycnt] = s1;
			keywords2[keycnt] = s2;
			keycnt ++;
		}
		for (int i = 0; i<N; ++i) {
			string qname = "sim" + to_string(1 + i / 2) + '/' + ((i % 2) ? '2' : '1');
			ret[i] = qname + ",20,1,150,+,1";
		}
		return ret;
	}
};