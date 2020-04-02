
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "bio.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;

bool IsValidDNASequence(const std::string& input) {

	for (auto c : input) {
		if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
			return false;
		}
	}
	return true;	
}

void GetReverseComplementSequence(const std::string& input, std::string* const output) {

	string copy_input = input;
	std::reverse(copy_input.begin(), copy_input.end());
	string temp = "";
	for (auto c : input) {
		if (c == 'A'){
			temp +='T';
		}
		else if (c == 'T') {
			temp += 'A';
		}
		else if (c == 'C') {
			temp += 'G';
		}
		else if (c == 'G') {
			temp += 'C';
		}
	}
	std::reverse(temp.begin(), temp.end());
	*output = temp;
}

std::string GetRNATranscript(const std::string& input) {

	string copy_input = input;
	string string_output = "";
	GetReverseComplementSequence(input, &string_output);
	for (auto &c : string_output) {
		if (c == 'T'){
			c = 'U';
		}
	}
	return string_output;
}

std::vector<std::vector<std::string>> GetReadingFramesAsCodons(const std::string& input) {

	string original1 = input;
	string original2;
	string original3;
	string antill1;
	string antill2;
	string antill3;
	
	GetReverseComplementSequence(original1, &antill1);
	original1 = GetRNATranscript(original1);
	antill1 = GetRNATranscript(antill1);
	

	/* 
		these next for loops are being used to offset the RNA sqeunce
	*/
	for (int i = 1; i < static_cast<int>(input.size()); ++i) {
		antill2 += antill1[i];
	}
	for (int i = 2; i < static_cast<int>(input.size()); ++i) {
		antill3 += antill1[i];
	}
	for (int i = 1; i < static_cast<int>(input.size()); ++i) {
		original2 += original1[i];
	}
	for (int i = 2; i < static_cast<int>(input.size()); ++i) {
		original3 += original1[i];
	}
	
	vector<string> v1 = { original1, original2, original3, antill1, antill2, antill3 };
	vector<vector<string>> return_vector;
	vector<string> current;

	/*
		these next 2 for loops are taking v1 and making sure to create a vector that is
		up to spec and capable of knowing when there is less than 3 per codon
		the inner for loop will be checking if the size is divisable by 3 and if not it
		will not be appending this new codon to the vecotr
	*/
	for (int i = 0; i < static_cast<int>(v1.size()); i++) {
		string working_sequence = v1.at(i);
		int start = 0;
		current = {};

		for (int j = 0; j < static_cast<int>(working_sequence.size()) / 3; ++j)	{
			current.push_back(working_sequence.substr(start, 3));
			start += 3;
		}
		return_vector.push_back(current);
	}
	return return_vector;
}

std::string Translate(const std::vector<std::string>& codon_sequence) {
	string result = "";
	vector<string> codons = { "GCU", "GCC", "GCA", "GCG", "CGU", "CGC", "CGA", "CGG",
		"AGA", "AGG",
"AAU", "AAC", "GAU", "GAC", "UGU", "UGC", "CAA", "CAG", "GAA", "GAG",
"GGU", "GGC", "GGA", "GGG", "CAU", "CAC", "AUU", "AUC", "AUA", "UUA",
"UUG", "CUU", "CUC", "CUA", "CUG", "AAA", "AAG", "AUG", "UUU", "UUC",
"CCU", "CCC", "CCA", "CCG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC",
"ACU", "ACC", "ACA", "ACG", "UGG", "UAU", "UAC", "GUU", "GUC", "GUA",
"GUG", "UAG", "UGA", "UAA" };
	vector<string> amino_acids = { "A", "A", "A", "A", "R", "R", "R", "R", "R", "R",
		"N", "N", "D", "D",
"C", "C", "Q", "Q", "E", "E", "G", "G", "G", "G", "H", "H", "I", "I",
"I", "L", "L", "L", "L", "L", "L", "K", "K", "M", "F", "F", "P", "P",
"P", "P", "S", "S", "S", "S", "S", "S", "T", "T", "T", "T", "W", "Y",
"Y", "V", "V", "V", "V", "*", "*", "*" };

	int possible_codon_size = static_cast<int>(codons.size()); 
	int codon_size = static_cast<int>(codon_sequence.size());

	/*
		these for loops will take a series of codons and convert them to amino acids
	*/
	for (int i = 0; i < codon_size; i++) {
		for (int j = 0; j < possible_codon_size; j++) {
			if (codon_sequence.at(i) == codons.at(j)) {
				result += amino_acids.at(j);
			}
		}
	}
	return result;
}

std::string GetLongestOpenReadingFrame(const std::string& DNA_sequence) {
	string temp_string = DNA_sequence;
	string working_string;
	string biggest_str = "";

	/*
		this for loop will using a boolean to find longest sequence this will loop
		thru the string and go true when an M is found and while true append to a
		new string then boolean will go false when a '*' is encountered and will
		check if the new string that has been appended is longer than the previous
	*/
	for (vector<string> frame : GetReadingFramesAsCodons(temp_string)) {
		working_string = Translate(frame);
		bool check_char = false;
		string current_str = "";
		for (auto c : working_string) {
			if (c == 'M') {
				check_char = true;
			}
			if (check_char == true) {
				current_str += c;
			}
			if (c == '*') {
				if (static_cast<int>(current_str.size()) > static_cast<int>(biggest_str.size())) {
					biggest_str = current_str;
				}
				current_str = "";
				check_char = false;
			}
		}
	}
	return biggest_str;
}