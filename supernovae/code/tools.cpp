#include "header.h"

void tokenize(std::string &str, char delim, std::vector<std::string> &out)
{
  out.clear();
  size_t start;
  size_t end = 0;
  
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
    {
      end = str.find(delim, start);
      out.push_back(str.substr(start, end - start));
    }
}
