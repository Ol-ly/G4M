#ifndef INT_TO_STR
#define INT_TO_STR

#include "intToStr.h"

std::string IntToStr(int i)
 {
  std::string res;
  bool minus = false;
  if (i == 0) {
    res = "0";
  } else {
    if (i < 0) {
      minus = true;
      i *= -1;
    }
    while (i!=0) {
      char t = i % 10;
      res += char(int('1') - 1 + t);
      i = (i - t)/10;
    }
  }
  char c;
  char len = res.length();
  for (int ii = 0; ii <= (len-ii-1); ii++)
   {
    c = res[ii];
    res[ii] = res[len-1-ii];
    res[len-1-ii] = c;
   }
  if (minus) res = std::string("-") + res;
  return res;
 }

void int2str(int i, char tmp[])
 {
  int in = i;
  int count=0;
  tmp[0]='\0';
    if (in == 0) {
      tmp[0] = '0';
      count++;
    }
  while (in!=0)
   {
    tmp[count]=char(int('1') - 1 + in%10);
    in = (in - in%10)/10;
    count++;
   }
  char c;
  for (int ii = 0; ii <= (count-ii-1); ii++)
   {
    c=tmp[ii];
    tmp[ii] = tmp[count-1-ii];
    tmp[count-1-ii] = c;
   } 
  tmp[count] = '\0';
 }

#endif