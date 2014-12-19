#include "Sales_item.h"
using namespace std;


int main(int argc, char** args){
	long a = (float)sqrt(10.0f);

	// declare variables to hold running sum and data for the next record 
    Sales_item total, trans;

    // is there data to process?
    if (std::cin >> total) {
        // if so, read the transaction records 
        while (std::cin >> trans)
            if (total.same_isbn(trans)) 
                // match: update the running total 
                total = total + trans;
            else {   
                // no match: print & assign to total
                std::cout << total << std::endl;
                total = trans;
            }
        // remember to print last record
        std::clog << total << std::endl; 
    } else {
        // no input!, warn the user
        std::cerr << "No data?!" << std::endl;
        return -1;  // indicate failure
    }

    return 0;
}

//#include <vector>
//int main(){
//	vector<Sales_item> sis;
//	Sales_item si;	
//	if (cin >> si) {// is there data to process?
//		sis.push_back(si);
//		while (std::cin >> si){
//			 if (total.same_isbn(trans)) 
//		}
//	} else {// no input!, warn the user        
//        std::cout << "No data?!" << std::endl;
//        return -1;  // indicate failure
//    }
//	return 0;
//}

