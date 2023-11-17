#ifndef CONFIG_PARSER
#define CONFIG_PARSER

#include <boost/multiprecision/cpp_int.hpp>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>

class ConfigParser
{
  private:
    std::string inputfileName;
    std::map<std::string, std::string> values;

    std::string findValue(const std::string &) const;

  public:
    ConfigParser() {}
    ConfigParser(const std::string &filename) { load(filename); }

    bool getBool(const std::string &) const;
    char getChar(const std::string &) const;
    double getDouble(const std::string &) const;
    float getFloat(const std::string &) const;
    int getInt(const std::string &) const;
    boost::multiprecision::cpp_int getMultiprecisionInt(const std::string &) const;
    short getShort(const std::string &) const;
    std::size_t getSizeT(const std::string &) const;
    std::string getString(const std::string &) const;

    void load(const std::string &);
};

#endif

