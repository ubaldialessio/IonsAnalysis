#include "definition.h"

int main(int argc, char *argv[]) {
	if (argc < 2) {
	    printf("Usage: \n");
	    printf("%s <input root file> \n", argv[0]);
	    return 1;
	}

	//check if is a valid input
    bool validInput = true;
	//check(); 

	//process
	if (validInput) {
		NAIA::NAIAChain chain;
		chain.Add(argv[1]);
		chain.SetupBranches();
	}
}
