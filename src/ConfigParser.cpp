#include "ConfigParser.h"

//------------------------------------------------------------------------------
//
//   ConfigLine
//
//   Class used by ConfigParser to read lines of the config file
//
//------------------------------------------------------------------------------
class ConfigLine {
  private:
    std::string original;  // the original string given
    std::string formatted; // the string with comment removed and escape sequences handled
    std::string parameter; // first word of the line
    std::string value;     // the rest of the line excluding comments. this will
                           // only be 1 word unless quotes were used to give a string

    void format();
    std::string getNthWord(const std::string &, std::size_t) const;
    bool isComment(const std::string &) const;
    bool isWhitespace(const std::string &) const;
    std::size_t numWordsInString(const std::string &) const;
    std::string& ltrim(std::string &, const std::string & = "\t\n\v\f\r ");
    std::string& rtrim(std::string &, const std::string & = "\t\n\v\f\r ");
    std::string& trim(std::string &, const std::string & = "\t\n\v\f\r ");

  public:
    ConfigLine() {}
    ConfigLine(const std::string &str) : original(str) { format(); }

    std::string getFormatted() const { return formatted; }
    std::string getOriginal() const { return original; }
    std::string getParameter() const { return parameter; }
    std::string getValue() const { return value; }

    bool hasData() const { return !formatted.empty(); }
    void load(const std::string &);
};


//------------------------------------------------------------------------------
// Populates the formatted, parameter, and value strings
//------------------------------------------------------------------------------
void ConfigLine::format() {
  if (isWhitespace(original) || isComment(original)) {
    return;
  }

  std::ostringstream oss;
  bool escaped = false;
  bool quoteOpened = false;
  bool withinQuote = false;

  for (std::size_t i = 0; i < original.size(); ++i) {
    // *
    // * Character is backslash
    // *
    if (original.at(i) == '\\') {
      if (escaped) {
        oss << "\\";
        escaped = false;
      } else {
        escaped = true;
      }
    }

    // *
    // * Character is a quotation mark
    // *
    else if (original.at(i) == '"') {
      if (escaped) { // print a quote if escape sequence was started
        oss << '"';
        escaped = false;
      } else {
        oss << '"';

        if (withinQuote) { // close the quote
          withinQuote = false;
        } else {
          if (quoteOpened) { // line tried to open a quote twice
            std::string errorMessage = "ConfigParser: Invalid line: " + original;
            throw std::runtime_error(errorMessage);
          } else { // open the quote
            withinQuote = true;
            quoteOpened = true;
          }
        }
      }
    }

    // *
    // * Character is a hash
    // *
    else if (original.at(i) == '#') {
      if (escaped || withinQuote) { // print a hash if in escape sequence or a quote
        oss << '#';
        escaped = false;
      } else { // Comment has started
        break;
      }
    }

    // *
    // * Normal character
    // *
    else {
      if (escaped) {
        switch (original.at(i)) {
          case 'a':
            oss << '\a';
            break;
          case 'b':
            oss << '\b';
            break;
          case 'f':
            oss << '\f';
            break;
          case 'n':
            oss << '\n';
            break;
          case 'r':
            oss << '\r';
            break;
          case 't':
            oss << '\t';
            break;
          case 'v':
            oss << '\v';
            break;
          case '\'':
            oss << '\'';
            break;
          case '?':
            oss << '\?';
            break;
          default:
            oss << '\\';
            oss << original.at(i);
        }

        escaped = false;
      } else {
        oss << original.at(i);
      }
    }
  }

  formatted = oss.str();
  trim(formatted);
  if ((!quoteOpened && numWordsInString(formatted) > 2)           // quote wasn't opened for string argument
  ||  (quoteOpened && formatted.back() != '"')                    // quote wasn't the last character
  ||  (quoteOpened && withinQuote)                                // quote wasn't closed
  ||  (quoteOpened && getNthWord(formatted, 2).front() != '"')) { // quote wasn't opened immediately
    std::string errorMessage = "ConfigParser: Invalid line: " + original;
    throw std::runtime_error(errorMessage);
  }

  parameter = getNthWord(formatted, 1);

  if (numWordsInString(formatted) > 1) {
    if (!quoteOpened) {
      value = getNthWord(formatted, 2);
    } else {
      std::size_t openQuotePos = formatted.find_first_of('"', parameter.size());
      value = formatted.substr(openQuotePos + 1);
      value.pop_back();
    }
  }
}


//------------------------------------------------------------------------------
// Returns the nth word of the string (whitespace separated)
// Starts from 1, so the first word would be getNthWord(myString, 1)
// The last word would be getNthWord(myString, numWordsInString(myString))
// http://www.cplusplus.com/forum/beginner/99171/
//------------------------------------------------------------------------------
std::string ConfigLine::getNthWord(const std::string &str, std::size_t n) const {
  std::string word = str;
  std::istringstream iss(word);
  while (n-- > 0 && (iss >> word));
  return word;
}


//------------------------------------------------------------------------------
// Returns true if the entire line is a comment (i.e., the first non-whitespace
// character is #)
//------------------------------------------------------------------------------
bool ConfigLine::isComment(const std::string &str) const {
  std::string copy = str;
  std::string const chars = "\t\n\v\f\r ";
  copy.erase(0, copy.find_first_not_of(chars));
  std::size_t commentStart = copy.find_first_of('#');
  return (commentStart == 0);
}


//------------------------------------------------------------------------------
// Returns true if the string contains only whitespace
//------------------------------------------------------------------------------
bool ConfigLine::isWhitespace(const std::string &str) const {
  for (std::size_t i = 0; i < str.length(); i++) {
    if (str.at(i) != ' '
    &&  str.at(i) != '\n'
    &&  str.at(i) != '\t'
    &&  str.at(i) != '\r'
    &&  str.at(i) != '\f'
    &&  str.at(i) != '\v') {
      return false;
    }
  }

  return true;
}


//------------------------------------------------------------------------------
// Loads a string and formats it
//------------------------------------------------------------------------------
void ConfigLine::load(const std::string &str) {
  original = str;
  formatted.clear();
  parameter.clear();
  value.clear();
  format();
}

//------------------------------------------------------------------------------
// Returns the number of whitespace separated words in the string
//------------------------------------------------------------------------------
std::size_t ConfigLine::numWordsInString(const std::string &str) const {
  const char *whitespace = " \t\n\r\f\v";
  std::string copy = str;

  copy.erase(0, copy.find_first_not_of(whitespace));
  copy.erase(copy.find_last_not_of(whitespace) + 1);

  if (copy.length() == 0) {
    return 0;
  }

  std::size_t num = 0;
  int i = 0;

  while (copy.length() != 0) {
    if (i % 2 == 0) {
      copy.erase(0, copy.find_first_of(whitespace));
      num++;
    } else {
      copy.erase(0, copy.find_first_not_of(whitespace));
    }
    ++i;
  }

  return num;
}

//------------------------------------------------------------------------------
// ltrim - trim std::string from the left
// rtrim - from the right
// trim  - both sides
// http://www.martinbroadhurst.com/how-to-trim-a-stdstring.html
//------------------------------------------------------------------------------
std::string& ConfigLine::ltrim(std::string &str, const std::string &chars) {
  str.erase(0, str.find_first_not_of(chars));
  return str;
}
std::string& ConfigLine::rtrim(std::string &str, const std::string &chars) {
  str.erase(str.find_last_not_of(chars) + 1);
  return str;
}
std::string& ConfigLine::trim(std::string &str, const std::string &chars) {
  return ltrim(rtrim(str, chars), chars);
}

//------------------------------------------------------------------------------
//
//
//   ConfigParser member functions
//
//
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// Private helper function. Returns the value (as a string) mapped to the
// parameter name
//------------------------------------------------------------------------------
inline std::string ConfigParser::findValue(const std::string &parameterName) const {
  std::string value;
  auto iter = values.find(parameterName);
  if (iter != values.end()) {
    value = iter->second;
  } else {
    std::string errorMessage = "ConfigParser: '";
    errorMessage += parameterName;
    errorMessage += "' was not provided in ";
    errorMessage += inputfileName;
    throw std::runtime_error(errorMessage);
  }

  return value;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a bool
//------------------------------------------------------------------------------
bool ConfigParser::getBool(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  if (value == "true" || value == "1") {
    return true;
  } else if (value == "false" || value == "0") {
    return false;
  } else {
    std::istringstream stream(value);
    bool returnValue;
    stream >> returnValue;
    return returnValue;
  }
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a char
//------------------------------------------------------------------------------
char ConfigParser::getChar(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::istringstream stream(value);
  char returnValue;
  stream >> returnValue;
  return returnValue;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a double
//------------------------------------------------------------------------------
double ConfigParser::getDouble(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::istringstream stream(value);
  double returnValue;
  stream >> returnValue;
  return returnValue;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a float
//------------------------------------------------------------------------------
float ConfigParser::getFloat(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::istringstream stream(value);
  float returnValue;
  stream >> returnValue;
  return returnValue;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as an int
//------------------------------------------------------------------------------
int ConfigParser::getInt(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::istringstream stream(value);
  int returnValue;
  stream >> returnValue;
  return returnValue;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a boost::multiprecision::cpp_int
//------------------------------------------------------------------------------
// boost::multiprecision::cpp_int ConfigParser::getMultiprecisionInt(const std::string &parameterName) const {
//   std::string value = findValue(parameterName);
//   boost::multiprecision::cpp_int returnValue(value);
//   return returnValue;
// }

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a short
//------------------------------------------------------------------------------
short ConfigParser::getShort(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::istringstream stream(value);
  short returnValue;
  stream >> returnValue;
  return returnValue;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a size_t
//------------------------------------------------------------------------------
std::size_t ConfigParser::getSizeT(const std::string &parameterName) const {
  std::string value = findValue(parameterName);
  std::istringstream stream(value);
  std::size_t returnValue;
  stream >> returnValue;
  return returnValue;
}

//------------------------------------------------------------------------------
// Returns the value mapped to the parameter as a string
//------------------------------------------------------------------------------
std::string ConfigParser::getString(const std::string &parameterName) const {
  return findValue(parameterName);
}

//------------------------------------------------------------------------------
// Loads a config file and parses the information into the map
//------------------------------------------------------------------------------
void ConfigParser::load(const std::string &filename) {
  inputfileName = filename;
  std::ifstream inputfile;
  inputfile.open(inputfileName);
  if (!inputfile.is_open()) {
    std::string errorMessage = "ConfigParser: Could not open ";
    errorMessage += inputfileName;
    throw std::runtime_error(errorMessage);
  }

  while(inputfile.good()) {
    std::string str;
    std::getline(inputfile, str);

    ConfigLine line(str);
    if (line.hasData()) {
      auto status = values.insert(std::pair<std::string, std::string>(line.getParameter(), line.getValue()));
      if (status.second == false) {
        std::string errorMessage = "ConfigParser: '";
        errorMessage += line.getParameter();
        errorMessage += "' provided more than once in ";
        errorMessage += inputfileName;
        throw std::runtime_error(errorMessage);
      }
    }
  }

  inputfile.close();
}
