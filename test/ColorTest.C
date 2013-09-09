//
//The ANSI Terminal specification gives programs (including BBS door games) the ability to change the text color or background color (as well as move the cursor, amoung other things)
//
//An ansi code begins with the ESC character* (ascii 27) and a left square bracket... followed by a number (or 2 or more separated by a semicolon) and a letter.
//
//In the case of colour codes, the trailing letter is "m"...
//
//So as an example, we have ESC[31m ... this will change the foreground colour to red.
//
//The codes are as follows:
//
//1m     -     Change text to hicolour (bold) mode
//4m     -        "    "   "  Underline (doesn't seem to work)
//5m     -        "    "   "  BLINK!!
//8m     -        "    "   "  Hidden (same colour as bg)
//30m    -        "    "   "  Black
//31m    -        "    "   "  Red
//32m    -        "    "   "  Green
//33m    -        "    "   "  Yellow
//34m    -        "    "   "  Blue
//35m    -        "    "   "  Magenta
//36m    -        "    "   "  Cyan
//37m    -        "    "   "  White
//40m    -     Change Background to Black
//41m    -        "       "      "  Red
//42m    -        "       "      "  Green
//43m    -        "       "      "  Yellow
//44m    -        "       "      "  Blue
//45m    -        "       "      "  Magenta
//46m    -        "       "      "  Cyan
//47m    -        "       "      "  White
//7m     -     Change to Black text on a White bg
//0m     -     Turn off all attributes.
//
//
//Now for example, say I wanted blinking, yellow text on a magenta background... I'd type ESC[45;33;5m
//
//* - The Escape character looks like an arrow pointing left, and can be produced in a hex editor, or with dos' EDIT command by typing CTRL+P and then ESC.



#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <ostream>
#include "colorstring.h"

inline std::ostream& red (std::ostream &s)
{
	s << "\033[31m";
	return s;
}
inline std::ostream& green (std::ostream &s)
{
	s << "\033[32m";
	return s;
}
inline std::ostream& yellow (std::ostream &s)
{
	s << "\033[33m";
	return s;
}

inline std::ostream& blue (std::ostream &s)
{
	s << "\033[34m";
	return s;
}



int main (void) {
	
	//std::string c_blk("black");
	//std::string c_mag("magenta");
	//printf("\033[0m%s\033[0m\n",c_blk.c_str());
	//printf("%s\033[0m\n",c_mag.c_str());
	//printf ("\033[%im%s, \033[1mbold %s\033[0m\n",
	//			0, "black", "black");
	//printf ("\033[%im%s, \033[1mbold %s\033[0m\n",
	//			16, "red", "red");

	//printf ("\033[31m%s  \033[1m%s  \n","red", "bold");
	//std::string str("\033[31mred  \033[1m bold");

	/*std::string str("some text");
	std::cout << red << str << std::endl;
	std::cout << green << str << std::endl;
	

	std::cout << "\033[0;" << 32 << "mHello!\033[0m" << std::endl;
	std::cout << "\033[0;" << 33 << "mHello World!\033[0m" << std::endl;
	std::cout << "\033[0;" << 34 << "mHello World Hello!\033[0m" << std::endl;
	std::cout << "\033[0;" << 35 << "mHello World Hello World!\033[0m" << std::endl;
*/

	std::string s("some new text");
	std::cout << "correct lenght = " << s.length() << std::endl;
	colorstring cs("some new text","red");
	colorstring cs2("some new text","yellow");
	std::cout << cs.str << std::endl;
	std::cout << cs2.str << std::endl;
//	std::cout << "length = " << cs2.str.length() << std::endl;
//	std::cout << "length = " << cs2.length() << std::endl;

	
	exit (0);
}

