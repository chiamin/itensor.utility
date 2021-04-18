#ifndef __STRINGUTILITY_H_CMC__
#define __STRINGUTILITY_H_CMC__
#include <string>
using namespace std;

inline bool hasStr (string str, string key)
{
    size_t pos = str.find (key);
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

inline vector<string> splitStr (string str, const string& delimiter)
{
    vector<string> re;
	size_t pos = str.find (delimiter);
    while (pos != string::npos)
    {
        re.push_back (str.substr (0, pos));
        str = str.substr (pos + delimiter.size());
        pos = str.find (delimiter);
    }
    re.push_back (str);
    return re;
}
#endif
