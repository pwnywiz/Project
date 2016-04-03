#ifndef SYLVESTER_H_
#define SYLVESTER_H_

#include "polys.h"

class Sylv
{
	private:
		int main1, main2;
		int hidden1, hidden2;
		int maxMain, maxHidden;
		UniPoly **func1;
		UniPoly **func2;
		UniPoly *zeros;
		UniPoly ***sylvester; // 2D Array of pointers to UniPoly object

	public:
		Sylv(BivPoly *clsA, BivPoly *clsB);
		~Sylv();

		void printSylv(); // Prints the Sylvester Matrix
		double specific(int row, int col, int part); // Ask for a specific Sylvester value
		int get_hidden() { return maxHidden; } // Returns max Hidden degree
		int get_main() { return maxMain; } // Returns max Main degree
		int get_total() { return main1+main2; } // Returns the Sylvester degree
};

#endif /* SYLVESTER_H_ */
