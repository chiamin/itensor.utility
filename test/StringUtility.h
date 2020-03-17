#include <string>
using namespace std;

inline bool find (string str, string toFind)
{
    size_t pos = str.find (toFind);
    return !(pos == std::string::npos);
}

inline void eraseSubStr (string& mainStr, const string& toErase)
{
	// Search for the substring in string
	size_t pos = mainStr.find(toErase); 
	if (pos != string::npos)
	{
		// If found then erase it from string
		mainStr.erase(pos, toErase.length());
	}
}
