#include <vector>
#include <fstream>
#include <iostream> 
#include <string>
#include <assert.h>
#include "string_util.h"

int main( int argc, char **argv )
{
	std::string seperator(" ");
	// write
	std::vector<float> arr;
	arr.push_back(1.0);	arr.push_back(2.0);	arr.push_back(3.0);	arr.push_back(4.0);
	std::ofstream ofs("output.txt");
	for (int i = 0; i < arr.size(); ++i)
	{
		ofs << arr[i] << seperator;
	}
	ofs << std::endl;
	ofs.close();
	
	// read
	std::vector<float> arr1;
	std::ifstream ifs("output.txt");
	std::string cur_line;
	while(std::getline(ifs,cur_line))
	{
		while (cur_line.size())
		{
			std::string tmp = jj::get_substr(cur_line, seperator);
			if (!tmp.compare(seperator))
				continue;
			arr1.push_back( jj::str2double(tmp) );
		}
	}
	ifs.close();

	// verify
	for (int i = 0; i < arr.size(); ++i)
	{
		assert(arr[i] == arr1[i]);
	}
}