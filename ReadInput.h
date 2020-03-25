#ifndef __READINPUT_H_CMC__
#define __READINPUT_H_CMC__
#include <algorithm>
#include <fstream>
#include <vector>
using namespace std;

template <typename T>
vector<T> split_str (const string& str, int skip=0)
{
    vector<T> re;
    std::istringstream iss (str);

    string drop;
    for(int i = 0; i < skip; i++)
        iss >> drop;

    T a;
    while (iss.good()) {
        iss >> a;
        re.push_back (a);
    }
    return re;
}

template <typename T>
vector<T> read_vector (ifstream& ifs, const string& key, int skip=2)
{
    // Go to the key-word line
    string line;
    while (getline(ifs, line))
    {
        vector<string> words = split_str<string> (line);
        if (std::find (words.begin(), words.end(), key) != words.end())
            break;
    }

    // Read data
    vector<T> re = split_str<T> (line, skip);
    if (re.size() == 0) {
        cout << "ReadInput.h: read_vector: " << key << " not found" << endl;
        cout << "             " << line << endl;
        throw;
    }
    return re;
}

template <typename T>
vector<T> read_vector (const string& fname, const string& key, int skip=2)
{
    ifstream ifs (fname);
    return read_vector<T> (ifs, key, skip);
}

bool fstream_goto (ifstream& ifs, const string& key, int skipline=0, string end="")
{
    // Go to the key-word line
    string line;
    while (getline(ifs, line))
    {
        vector<string> words = split_str<string> (line);
        if (std::find (words.begin(), words.end(), key) != words.end())
            break;
        if (end != "")
            if (std::find (words.begin(), words.end(), end) != words.end())
                return false;
        //auto n = line.find (key);
        //if (n != string::npos) break;
    }

    if (!ifs) return false;

    for(int i = 0; i < skipline; i++)
        getline(ifs, line);
    return true;
}

bool has_keyword (string fname, string key, string start="", string end="")
{
    ifstream ifs (fname);
    if (!ifs) {
        cout << "Error: cannot open file " << fname << endl;
        throw;
    }
    if (start != "")
        fstream_goto (ifs, start);
    return fstream_goto (ifs, key, 0, end);
}

vector<string> read_bracket (ifstream& ifs, string key, int skipline=0)
{
    vector<string> re;

    // Go to the key-word line
    bool good = fstream_goto (ifs, key);
    if (!good) return re;
    good = fstream_goto (ifs, "{", skipline);
    if (!good) return re;

    string line;
    while (getline(ifs, line))
    {
        auto n = line.find ("}");
        if (n != string::npos) break;

        re.push_back (line);
    }
    return re;
}

vector<string> read_bracket (string file, string key, int skipline=0)
{
    ifstream ifs (file);
    if (!ifs)
    {
        cout << "Error: read_bracket: cannot open file: " << file << endl;
        throw;
    }
    return read_bracket (ifs, key, skipline);
}
#endif
