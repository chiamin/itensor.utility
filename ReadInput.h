#ifndef __READINPUT_H_CMC__
#define __READINPUT_H_CMC__
#include <algorithm>
#include <fstream>
#include <vector>
#include <cctype>
#include <sstream>
using namespace std;

string raw_string (const string& str)
{
    // s is our escaped output string
    std::string s = "";
    // loop through all characters
    for(char c : str)
    {
        // check if a given character is printable
        // the cast is necessary to avoid undefined behaviour
        if(isprint((unsigned char)c))
            s += c;
        else
        {
            std::stringstream stream;
            // if the character is not printable
            // we'll convert it to a hex string using a stringstream
            // note that since char is signed we have to cast it to unsigned first
            stream << std::hex << (unsigned int)(unsigned char)(c);
            std::string code = stream.str();
            s += std::string("\\x")+(code.size()<2?"0":"")+code;
            // alternatively for URL encodings:
            //s += std::string("%")+(code.size()<2?"0":"")+code;
        }
    }
    return s;
}

// trim from left
inline void lstrip (std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
}

// trim from right
inline void rstrip (std::string& s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
}

// trim from left & right
inline void stripe (std::string& s, const char* t = " \t\n\r\f\v")
{
    rstrip(s, t);
    lstrip(s, t);
}

ifstream open_file (const string& file)
{
    ifstream ifs (file);
    if (!ifs)
    {
        cout << "Cannot open " << file << endl;
        throw;
    }
    return move(ifs);
}

template <typename T>
vector<T> split_str (const string& str, int skip=0)
{
    string str2 = str;
    rstrip (str2);
    vector<T> re;
    std::istringstream iss (str2);

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
