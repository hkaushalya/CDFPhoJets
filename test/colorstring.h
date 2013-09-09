#pragma once  //for precompiler optimization.
#ifndef COLORSTRING_HH
#define COLORSTRING_HH


/*************************************************************************
 * A simple color string. Has almost all std::string properties.
 *************************************************************************
 * Developer notes:
 * inheriting from std::string seems a bad idea as they are not designed 
 * to be inhereted(hence no virtual methods to be over written. this 
 * seems like a better and simpler way and still have all the 
 * functionality of std::string to use std::string methods.
 * Side effect of this prepending color characters is that
 * std::string.length() method count those too! need local method to
 * overcome this shortage.
 *************************************************************************
 * USAGE::
 * colorstring mystr("some text","red");
 * std::cout << mystr.str << std::endl;
 * mystr.str = "new text";
 * int length = mystr.str.length(); // this will include the hidden
 * 											// color characters too!
 *************************************************************************
 * Author: Samantha K. Hewamanage <samantha@fnal.gov> 12-03-2009
 *************************************************************************/

#include <string>

class colorstring
{
	public:
		std::string str;
		std::string strcolor;
		//std::string noattcolor;

		colorstring(std::string s, std::string color)
		{
			std::string pre;
			std::string noattcolor("\033[0m");

			if (color == "red")         { pre = "\033[31m"; } 
			else if (color == "yellow") { pre = "\033[33m"; } 
			else if (color == "green")  { pre = "\033[32m"; } 
			else if (color == "blue")   { pre = "\033[34m"; } 
			else if (color == "magenta"){ pre = "\033[35m"; } 
			else if (color == "cyan")   { pre = "\033[36m"; } 
			else                        { pre = "\033[0m"; } //deafult is black

			strcolor = pre;
			//std::cout << "strcolor length = " << strcolor.length() << std::endl;
			//std::cout << "strcolor = " << strcolor << std::endl;
			str = pre + s + noattcolor;	
		}

		//if you want only the original text length without the
		//color characters
		//unsigned int length() const
		//{
			//find locations of 'm'
		//	if (str.length()>0)
		//	{
				//std::cout << "str length = " << str.length() << std::endl;
				//std::cout << "str color = " << strcolor << std::endl;
				//const int s1 = strcolor.length();
				//const int s2 = noattcolor.length();
				//std::cout << "s1 length = " << s1 << std::endl;
				//std::cout << "s2 length = " << s2 << std::endl;
				//return (str.length() - s2 - s1);

		//	} else return 0;

		//};


};
#endif /* COLORSTRING_HH */

